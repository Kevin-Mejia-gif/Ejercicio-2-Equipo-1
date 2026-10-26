# Construcción de la base mensual de Colombia

setwd("~/Documents/Macroeconometría/Tareas/Ejercicio 2/Colombia")

library(readxl)
library(urca)

inicio <- as.Date("2005-01-01")
fin    <- as.Date("2026-07-01")

while (sink.number() > 0) sink()
sink("log/01_datos.log", split = TRUE)
cat("Construcción de la base:", format(Sys.time()), "\n\n")

leer_banrep <- function(archivo) {
  d <- read.csv(file.path("data/originales", archivo), header = FALSE, skip = 1,
                col.names = c("fecha", "valor"))
  d$fecha <- as.Date(d$fecha, format = "%Y/%m/%d")
  d
}

# Promedio mensual de una serie diaria
mensual <- function(d) {
  d$mes <- as.Date(format(d$fecha, "%Y-%m-01"))
  aggregate(valor ~ mes, data = d, FUN = mean)
}

leer_fred <- function(archivo) {
  d <- read.csv(file.path("data/originales", archivo), header = FALSE, skip = 1,
                col.names = c("fecha", "valor"))
  d$fecha <- as.Date(d$fecha)
  d
}

# Lleva la fecha al primer día del mes para poder unir las series
a_mes <- function(d) {
  data.frame(mes = as.Date(format(d$fecha, "%Y-%m-01")), valor = d$valor)
}

# Series diarias: TPM y TRM
tpm <- mensual(leer_banrep("TPM_BanRep_diaria.csv"))
trm <- mensual(leer_banrep("TRM_BanRep_diaria.csv"))

# Series mensuales: IPC e ITCR
ipc  <- a_mes(leer_banrep("IPC_DANE_mensual.csv"))
itcr <- a_mes(leer_banrep("ITCR_IPC_total_BanRep_mensual.csv"))

# Cartera bruta total en moneda legal, con ajuste por titularización (saldo a fin de mes)
cartera <- a_mes(leer_banrep("CARTERA_BRUTA_ML_BanRep_mensual.csv"))

# Variables externas (FRED): tasa de fondos federales efectiva y Brent, promedios mensuales
fed   <- a_mes(leer_fred("FEDFUNDS_FRED_mensual.csv"))
brent <- a_mes(leer_fred("BRENT_EIA_FRED_mensual.csv"))

# Expectativas de inflación a 12 meses: mediana de la Encuesta Mensual de Expectativas de Analistas (BanRep)
eme <- read.csv("data/originales/EME_expectativas_inflacion_BanRep_mensual.csv")
expect <- data.frame(mes = as.Date(paste0(eme$fecha, "-01")), valor = 100 * eme$exp12_mediana)

# ISE ajustado por estacionalidad y calendario (Cuadro 2 del anexo del DANE)
c2 <- read_excel("data/originales/anex-ISE-9actividades-jul2026.xlsx",
                 sheet = "Cuadro 2", col_names = FALSE, .name_repair = "minimal",
                 progress = FALSE)
filas <- grep("^Indicador de Seguimiento a la Econom.a$", trimws(c2[[1]]))
fila <- filas[!is.na(suppressWarnings(as.numeric(unlist(c2[filas, 2]))))][1]
valores <- as.numeric(unlist(c2[fila, -1]))
valores <- valores[!is.na(valores)]
ise <- data.frame(mes = seq(as.Date("2005-01-01"), by = "month", length.out = length(valores)),
                  valor = valores)

names(tpm)[2]  <- "tpm"
names(trm)[2]  <- "trm"
names(ipc)[2]  <- "ipc"
names(itcr)[2] <- "itcr"
names(cartera)[2] <- "cartera"
names(ise)[2]  <- "ise"
names(fed)[2]   <- "fed"
names(brent)[2] <- "brent"
names(expect)[2] <- "exp_inflacion_12m"

base <- Reduce(function(x, y) merge(x, y, by = "mes", all = TRUE),
               list(ise, ipc, tpm, trm, itcr, fed, brent, expect, cartera))
base <- base[order(base$mes), ]

# Transformaciones (antes de recortar, para no perder los primeros meses)
base$log_ise         <- 100 * log(base$ise)
base$inflacion_anual <- 100 * (log(base$ipc) - log(c(rep(NA, 12), head(base$ipc, -12))))
# Inflación mensual: la que entra al VAR. El IPC no está desestacionalizado, así que el VAR lleva dummies estacionales
# (ver las pruebas de raíz unitaria estacional al final y 03_comparacion.R)
base$inflacion_mensual <- 100 * (log(base$ipc) - log(c(NA, head(base$ipc, -1))))
base$log_itcr        <- 100 * log(base$itcr)
# Crédito real: cartera deflactada con el IPC, en variación anual
base$log_credito_real <- 100 * log(base$cartera / base$ipc)
base$credito_real_anual <- base$log_credito_real - c(rep(NA, 12), head(base$log_credito_real, -12))
base$d_log_itcr      <- c(NA, diff(base$log_itcr))
base$log_trm         <- 100 * log(base$trm)
base$log_brent       <- 100 * log(base$brent)
base$d_log_ise       <- c(NA, diff(base$log_ise))
base$d_log_brent     <- c(NA, diff(base$log_brent))
base$d_log_brent_l1  <- c(NA, head(base$d_log_brent, -1))
base$fed_l1          <- c(NA, head(base$fed, -1))

base <- base[base$mes >= inicio & base$mes <= fin, ]

# Dummies de impulso: marzo-diciembre de 2020 (pandemia) y mayo de 2021 (paro nacional)
meses_dummy <- c(seq(as.Date("2020-03-01"), as.Date("2020-12-01"), by = "month"), as.Date("2021-05-01"))
for (m in as.character(meses_dummy)) {
  base[[paste0("d", format(as.Date(m), "%Y%m"))]] <- as.integer(base$mes == as.Date(m))
}
# Dummies estacionales de febrero a diciembre (enero queda en la constante)
estacionales <- sprintf("est%02d", 2:12)
for (j in 2:12) base[[estacionales[j - 1]]] <- as.integer(as.integer(format(base$mes, "%m")) == j)

vars_var <- c("d_log_ise", "inflacion_mensual", "exp_inflacion_12m", "tpm", "d_log_itcr")
vars_ext <- c("d_log_brent", "d_log_brent_l1", "fed", "fed_l1")
dummies  <- paste0("d", format(meses_dummy, "%Y%m"))
cat("Muestra:", format(min(base$mes)), "a", format(max(base$mes)), "-", nrow(base), "meses\n")
cat("Faltantes en las variables del VAR y externas:\n")
print(colSums(is.na(base[, c(vars_var, vars_ext)])))

base$fecha <- format(base$mes, "%Y-%m")
write.csv(base[, c("fecha", "ise", "ipc", "tpm", "trm", "itcr", "fed", "brent", "cartera",
                   "log_ise", "d_log_ise", "inflacion_anual", "inflacion_mensual", "exp_inflacion_12m", "credito_real_anual",
                   "log_itcr", "d_log_itcr", "log_trm", "log_brent", "d_log_brent")],
          "data/colombia_mensual.csv", row.names = FALSE)
write.csv(base[, c("fecha", vars_var)], "data/colombia_var.csv", row.names = FALSE)
write.csv(base[, c("fecha", vars_ext, estacionales, dummies)], "data/colombia_exogenas.csv", row.names = FALSE)

# Gráficas
titulos <- c(log_ise = "Log ISE desestacionalizado (x100)",
             d_log_ise = "Variación mensual del log ISE (x100)",
             inflacion_mensual = "Inflación mensual (%), entra al VAR",
             inflacion_anual = "Inflación anual (%), referencia",
             exp_inflacion_12m = "Expectativa de inflación a 12 meses (%)",
             tpm = "Tasa de política monetaria (%)",
             d_log_itcr = "Variación mensual del log ITCR (x100)",
             log_itcr = "Log ITCR (x100)",
             log_trm = "Log TRM (x100)",
             log_brent = "Log Brent (x100)",
             d_log_brent = "Variación mensual del log Brent (x100)",
             fed = "Tasa de fondos federales (%)",
             credito_real_anual = "Crecimiento anual del crédito real (%)")
png("figuras/series_colombia.png", width = 1800, height = 3360, res = 170)
par(mfrow = c(7, 2), mar = c(3, 4, 2.5, 1))
for (v in names(titulos)) {
  plot(base$mes, base[[v]], type = "n", main = titulos[v], xlab = "", ylab = "", cex.main = 1)
  rect(as.Date("2020-03-01"), par("usr")[3], as.Date("2020-12-01"), par("usr")[4],
       col = adjustcolor("grey60", 0.25), border = NA)
  grid()
  lines(base$mes, base[[v]], lwd = 1.3, col = "#1f4e79")
}
invisible(dev.off())

# Pruebas de raíz unitaria: ADF (H0: raíz unitaria) y KPSS (H0: estacionariedad)
pruebas <- function(x, nombre, transformacion, tendencia) {
  x <- na.omit(x)
  adf  <- ur.df(x, type = if (tendencia) "trend" else "drift", lags = 12, selectlags = "AIC")
  kp   <- ur.kpss(x, type = if (tendencia) "tau" else "mu", lags = "short")
  data.frame(serie = nombre, transformacion = transformacion,
             deterministicos = if (tendencia) "constante y tendencia" else "constante",
             adf_estadistico = round(adf@teststat[1], 2),
             adf_critico_5 = adf@cval[1, "5pct"],
             adf_rechaza_raiz = adf@teststat[1] < adf@cval[1, "5pct"],
             kpss_estadistico = round(kp@teststat, 3),
             kpss_critico_5 = kp@cval[1, "5pct"],
             kpss_rechaza_estacionariedad = kp@teststat > kp@cval[1, "5pct"])
}

series <- list(log_ise = TRUE, inflacion_mensual = FALSE, inflacion_anual = FALSE, exp_inflacion_12m = FALSE, credito_real_anual = FALSE,
               tpm = FALSE, log_itcr = TRUE, log_brent = TRUE, fed = FALSE)
tabla <- do.call(rbind, lapply(names(series), function(v) rbind(
  pruebas(base[[v]], v, "nivel", series[[v]]),
  pruebas(diff(base[[v]]), v, "primera diferencia", FALSE))))
write.csv(tabla, "data/raices_unitarias.csv", row.names = FALSE)

cat("\nPruebas de raíz unitaria (valores críticos al 5%):\n")
print(tabla, row.names = FALSE)

# Pruebas de raíz unitaria estacional en el IPC. Deciden cómo tratar la estacionalidad (Kilian y Lütkepohl, sección 19.2):
# la diferencia de 12 meses (inflación anual) es la transformación adecuada si la estacionalidad tiene raíz unitaria;
# si es determinística, lo adecuado es la inflación mensual con dummies estacionales. Igual que ADF y KPSS, las dos
# pruebas tienen hipótesis nulas opuestas:
#   OCSB (Osborn, Chui, Smith y Birchenhall 1988), H0: hay raíz unitaria estacional (paquete forecast).
#   Canova y Hansen (1995), H0: la estacionalidad es determinística y estable. Se aplica a la inflación mensual con
#   las 11 frecuencias estacionales juntas; su valor crítico se simula de la distribución asintótica (puente browniano).
lipc <- ts(100 * log(base$ipc), start = c(2005, 1), frequency = 12)
canova_hansen <- function(y) {
  n <- length(y); tt <- 1:n
  Fm <- cbind(do.call(cbind, lapply(1:5, function(j) cbind(cos(2 * pi * j * tt / 12), sin(2 * pi * j * tt / 12)))), cos(pi * tt))
  X <- cbind(1, c(NA, y[-n]), Fm)                     # constante, un rezago de y y los 11 regresores trigonométricos
  ok <- complete.cases(X)
  e <- lm.fit(X[ok, ], y[ok])$residuals; nn <- length(e)
  G <- Fm[ok, ] * e; S <- apply(G, 2, cumsum)
  q <- round(4 * (nn / 100)^(2 / 9)); Om <- crossprod(G) / nn         # varianza de largo plazo (Newey-West)
  for (j in 1:q) { Gj <- crossprod(G[(j + 1):nn, ], G[1:(nn - j), ]) / nn; Om <- Om + (1 - j / (q + 1)) * (Gj + t(Gj)) }
  sum(diag(solve(Om, crossprod(S)))) / nn^2
}
set.seed(2026)
sim_ch <- replicate(5000, { W <- apply(matrix(rnorm(1000 * 11), 1000), 2, cumsum) / sqrt(1000)
  BB <- W - outer((1:1000) / 1000, W[1000, ]); mean(rowSums(BB^2)) })
vc_ch <- unname(quantile(sim_ch, 0.95))
estacional <- NULL
for (ini in c(2005, 2009)) {
  x <- if (ini == 2005) lipc else window(lipc, start = c(ini, 1))
  ch <- canova_hansen(as.numeric(diff(x)))
  oc <- if (requireNamespace("forecast", quietly = TRUE)) forecast::ocsb.test(x, lag.method = "AIC", maxlag = 12) else NULL
  estacional <- rbind(estacional, data.frame(
    muestra = paste0(ini, "-2026"),
    ocsb_estadistico = if (is.null(oc)) NA else round(oc$statistics, 2),
    ocsb_critico_5 = if (is.null(oc)) NA else round(oc$critical, 2),
    ocsb_rechaza_raiz_estacional = if (is.null(oc)) NA else oc$statistics < oc$critical,
    ch_estadistico = round(ch, 2), ch_critico_5 = round(vc_ch, 2), ch_rechaza_estacionalidad_estable = ch > vc_ch))
}
write.csv(estacional, "data/raices_estacionales.csv", row.names = FALSE)
cat("\nPruebas de raíz unitaria estacional (OCSB sobre 100 x log IPC; Canova-Hansen sobre la inflación mensual):\n")
print(estacional, row.names = FALSE)
if (!requireNamespace("forecast", quietly = TRUE)) cat("Para la prueba OCSB hay que instalar el paquete forecast.\n")

sink()
