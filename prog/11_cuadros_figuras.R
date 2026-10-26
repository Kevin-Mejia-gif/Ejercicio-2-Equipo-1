# Cuadros (.tex) y figuras de las secciones de datos y estimación y de sus anexos
# Se corre después de 01_datos.R y 02_var.R. Los cuadros quedan en cuadros/ y las figuras en figuras/.

setwd("~/Documents/Macroeconometría/Tareas/Ejercicio 2/Colombia")

dir.create("cuadros", showWarnings = FALSE)

base <- read.csv("data/colombia_mensual.csv")
base$mes <- as.Date(paste0(base$fecha, "-01"))
x_df <- read.csv("data/colombia_exogenas.csv")
v <- readRDS("data/var_reducido.rds")
m <- v$modelo
raices_ur <- read.csv("data/raices_unitarias.csv")
raices_est <- read.csv("data/raices_estacionales.csv")

inicio <- as.Date("2010-01-01")
fin    <- as.Date("2026-07-01")
cod_dummies <- grep("^d20", names(x_df), value = TRUE)
meses_dummy <- as.Date(paste0(substr(cod_dummies, 2, 5), "-", substr(cod_dummies, 6, 7), "-01"))
meses_es <- c("ene", "feb", "mar", "abr", "may", "jun", "jul", "ago", "sep", "oct", "nov", "dic")

etiquetas <- c(d_log_ise = "ISE", inflacion_mensual = "Inflación", exp_inflacion_12m = "Expectativas",
               tpm = "TPM", d_log_itcr = "ITCR")
K <- length(etiquetas)
num <- function(x, d = 3) ifelse(is.na(x), "", sprintf(paste0("$%.", d, "f$"), x))
fila <- function(...) paste0(paste(c(...), collapse = " & "), " \\\\")

escribir_cuadro <- function(archivo, titulo, etiqueta, columnas, encabezado, filas, nota, tamano = "\\small",
                            entorno = "tabular") {
  ancho <- if (entorno == "tabularx") "{\\textwidth}" else ""
  lineas <- c("\\begin{table}[H]", "\\centering",
              paste0("\\caption{", titulo, "}"), paste0("\\label{", etiqueta, "}"), tamano,
              paste0("\\begin{", entorno, "}", ancho, "{@{}", columnas, "@{}}"), "\\toprule", encabezado, "\\midrule",
              filas, "\\bottomrule", paste0("\\end{", entorno, "}"),
              paste0("\\notascuadro{\\textit{Nota.} ", nota, "}"), "\\end{table}")
  writeLines(lineas, file.path("cuadros", archivo), useBytes = TRUE)
}

# ---- Anexo A: datos ----------------------------------------------------------------------

fuentes <- rbind(
  c("ISE desestacionalizado y corregido por efecto calendario (2015 = 100)", "DANE", "Mensual", "ene 2005 -- jul 2026"),
  c("IPC total nacional (dic 2018 = 100)", "DANE, vía Banco de la República", "Mensual", "jul 1954 -- sep 2026"),
  c("Expectativa de inflación a 12 meses, mediana (Encuesta Mensual de Expectativas)", "Banco de la República", "Mensual", "sep 2003 -- sep 2026"),
  c("Tasa de política monetaria", "Banco de la República", "Diaria", "feb 1998 -- oct 2026"),
  c("ITCR según IPC, ponderaciones totales (2010 = 100)", "Banco de la República", "Mensual", "dic 1986 -- ago 2026"),
  c("Tasa de fondos federales efectiva", "Reserva Federal (H.15), vía FRED", "Mensual", "ene 2000 -- sep 2026"),
  c("Precio spot del Brent", "EIA, vía FRED", "Mensual", "ene 2000 -- sep 2026"))
escribir_cuadro("cuadro_fuentes.tex", "Series utilizadas y fuentes oficiales", "tab:fuentes",
                ">{\\raggedright\\arraybackslash}X >{\\raggedright\\arraybackslash}p{3.4cm} >{\\raggedright\\arraybackslash}p{1.9cm} >{\\raggedright\\arraybackslash}p{3.0cm}",
                fila("Serie", "Fuente", "Frecuencia original", "Cobertura descargada"),
                apply(fuentes, 1, function(r) fila(r)),
                "Las series originales se conservan sin modificar en \\texttt{data/originales/}. El ISE de 2024 es provisional y el de 2025--2026 preliminar. Las ponderaciones de los últimos meses del ITCR son estimaciones del Banco de la República.",
                entorno = "tabularx")

trans <- rbind(
  c("ISE", "Endógena 1", "$100 \\times (\\log ISE_t - \\log ISE_{t-1})$", "\\% mensual"),
  c("Inflación", "Endógena 2", "$100 \\times (\\log IPC_t - \\log IPC_{t-1})$", "\\% mensual"),
  c("Expectativas", "Endógena 3", "$100 \\times$ mediana de la expectativa a 12 meses", "\\%"),
  c("TPM", "Endógena 4", "Promedio mensual de la serie diaria", "\\%"),
  c("ITCR", "Endógena 5", "$100 \\times (\\log ITCR_t - \\log ITCR_{t-1})$", "\\% mensual"),
  c("Brent", "Exógena", "$100 \\times (\\log Brent_t - \\log Brent_{t-1})$, en $t$ y $t-1$", "\\% mensual"),
  c("Tasa Fed", "Exógena", "Nivel, en $t$ y $t-1$", "\\%"),
  c("Estacionales", "Exógena", "11 dummies, una por mes de febrero a diciembre", "---"),
  c("Puntuales", "Exógena", "Una dummy por mes: mar--dic 2020 y may 2021", "---"))
escribir_cuadro("cuadro_transformaciones.tex", "Variables del VAR y sus transformaciones", "tab:transformaciones",
                "llll", fila("Variable", "Papel en el VAR", "Transformación", "Unidad"),
                apply(trans, 1, function(r) fila(r)),
                "El orden de las endógenas es el del VAR. Un aumento del ITCR es una depreciación real del peso. Las respuestas al impulso del ISE y del ITCR se acumulan para leerse en nivel; la de la inflación se acumula 12 meses para leerse como inflación anual.")

muestra <- base[base$mes >= inicio & base$mes <= fin, ]
muestra$fed <- x_df$fed[match(muestra$fecha, x_df$fecha)]
desc_vars <- c(d_log_ise = "ISE (\\% mensual)", inflacion_mensual = "Inflación mensual (\\%)",
               inflacion_anual = "Inflación anual (\\%), referencia", exp_inflacion_12m = "Expectativas a 12 meses (\\%)",
               tpm = "TPM (\\%)", d_log_itcr = "ITCR (\\% mensual)", d_log_brent = "Brent (\\% mensual)", fed = "Tasa Fed (\\%)")
filas_desc <- sapply(names(desc_vars), function(k) {
  x <- muestra[[k]]
  fila(desc_vars[k], num(mean(x), 2), num(sd(x), 2), num(min(x), 2), num(max(x), 2))
})
escribir_cuadro("cuadro_descriptivos.tex", "Estadística descriptiva, 2010m1--2026m7", "tab:descriptivos",
                "lrrrr", fila("Variable", "Media", "Desv. est.", "Mínimo", "Máximo"), filas_desc,
                sprintf("Variables transformadas como entran al VAR; %d observaciones mensuales. La inflación anual no entra al VAR y se incluye solo como referencia. Fuente: DANE, Banco de la República y FRED.", nrow(muestra)))

nombres_ur <- c(log_ise = "Log ISE", inflacion_mensual = "Inflación mensual", exp_inflacion_12m = "Expectativas",
                tpm = "TPM", log_itcr = "Log ITCR", log_brent = "Log Brent", fed = "Tasa Fed")
ur <- raices_ur[raices_ur$serie %in% names(nombres_ur), ]
suppressMessages(library(urca))
infl_full <- base$inflacion_mensual[!is.na(base$inflacion_mensual)]
mes_full <- as.integer(substr(base$fecha[!is.na(base$inflacion_mensual)], 6, 7))
infl_sa <- resid(lm(infl_full ~ factor(mes_full)))
prueba_sa <- function(x, transformacion) {
  a <- ur.df(x, type = "drift", lags = 12, selectlags = "AIC"); k <- ur.kpss(x, type = "mu", lags = "short")
  data.frame(serie = "inflacion_sa", transformacion = transformacion, deterministicos = "constante",
             adf_estadistico = round(a@teststat[1], 2), adf_critico_5 = a@cval[1, "5pct"],
             adf_rechaza_raiz = a@teststat[1] < a@cval[1, "5pct"], kpss_estadistico = round(k@teststat, 3),
             kpss_critico_5 = k@cval[1, "5pct"], kpss_rechaza_estacionariedad = k@teststat > k@cval[1, "5pct"])
}
ur <- rbind(ur, prueba_sa(infl_sa, "nivel"), prueba_sa(diff(infl_sa), "primera diferencia"))
nombres_ur <- append(nombres_ur, c(inflacion_sa = "Inflación mensual desest."), after = 2)
ur <- ur[order(match(ur$serie, names(nombres_ur))), ]
si_no <- function(x) ifelse(x, "Sí", "No")
filas_ur <- sapply(seq_len(nrow(ur)), function(i) {
  r <- ur[i, ]
  fila(if (r$transformacion == "nivel") nombres_ur[r$serie] else "",
       if (r$transformacion == "nivel") "Nivel" else "Diferencia",
       num(r$adf_estadistico, 2), num(r$adf_critico_5, 2), si_no(r$adf_rechaza_raiz),
       num(r$kpss_estadistico, 3), num(r$kpss_critico_5, 3), si_no(r$kpss_rechaza_estacionariedad))
})
lineas <- c("\\begin{table}[H]", "\\centering", "\\caption{Pruebas de raíz unitaria}", "\\label{tab:raices}", "\\footnotesize",
            "\\begin{tabular}{@{}llrrcrrc@{}}", "\\toprule",
            " & & \\multicolumn{3}{c}{ADF} & \\multicolumn{3}{c}{KPSS} \\\\",
            "\\cmidrule(lr){3-5}\\cmidrule(l){6-8}",
            fila("Serie", "Transformación", "Estad.", "VC 5\\%", "Rechaza $H_0$", "Estad.", "VC 5\\%", "Rechaza $H_0$"),
            "\\midrule", filas_ur, "\\bottomrule", "\\end{tabular}",
            "\\notascuadro{\\textit{Nota.} ADF: $H_0$ = raíz unitaria (rezagos por AIC, máximo 12). KPSS: $H_0$ = estacionariedad. En nivel, log ISE, log ITCR y log Brent con constante y tendencia; las demás con constante. Las diferencias, con constante. La inflación mensual desestacionalizada (desest.) es el residuo de regresarla en 11 dummies estacionales. Muestra 2005--2026. VC = valor crítico.}",
            "\\end{table}")
writeLines(lineas, "cuadros/cuadro_raices.tex", useBytes = TRUE)

filas_est <- sapply(seq_len(nrow(raices_est)), function(i) {
  r <- raices_est[i, ]
  fila(gsub("-", "--", r$muestra), num(r$ocsb_estadistico, 2), num(r$ocsb_critico_5, 2), si_no(r$ocsb_rechaza_raiz_estacional),
       num(r$ch_estadistico, 2), num(r$ch_critico_5, 2), si_no(r$ch_rechaza_estacionalidad_estable))
})
lineas <- c("\\begin{table}[H]", "\\centering", "\\caption{Pruebas de raíz unitaria estacional en el IPC}", "\\label{tab:estacional}", "\\small",
            "\\begin{tabular}{@{}lrrcrrc@{}}", "\\toprule",
            " & \\multicolumn{3}{c}{OCSB} & \\multicolumn{3}{c}{Canova--Hansen} \\\\",
            "\\cmidrule(lr){2-4}\\cmidrule(l){5-7}",
            fila("Muestra", "Estad.", "VC 5\\%", "Rechaza $H_0$", "Estad.", "VC 5\\%", "Rechaza $H_0$"),
            "\\midrule", filas_est, "\\bottomrule", "\\end{tabular}",
            "\\notascuadro{\\textit{Nota.} OCSB (Osborn et al., 1988), sobre $100 \\times \\log IPC$: $H_0$ = hay raíz unitaria estacional. Canova y Hansen (1995), sobre la inflación mensual: $H_0$ = la estacionalidad es determinística y estable. Las dos apuntan a estacionalidad determinística, por lo que la inflación entra en variación mensual con dummies estacionales en lugar de en variación anual (Kilian y Lütkepohl, 2017, sección 19.2). VC = valor crítico.}",
            "\\end{table}")
writeLines(lineas, "cuadros/cuadro_estacional.tex", useBytes = TRUE)

# ---- Anexo A.6: especificaciones alternativas (resultados de 03_comparacion.R) -------------
comp <- read.csv("data/comparacion_especificaciones.csv")
comp_inf <- read.csv("data/comparacion_inflacion.csv")
d_final <- "mar-dic 2020 + may 2021"
nombres_esp <- c(base = "Base: sin expectativas", expectativas = "Base + expectativas (elegida)",
                 credito = "Base + crédito real (variación anual)", credito_dif = "Base + crédito real (cambio de la variación anual)",
                 exp_credito = "Base + expectativas + crédito real", exp_credito_dif = "Base + expectativas + cambio del crédito real")
nombres_dum <- c("mar-jun 2020" = "Mar--jun 2020", "mar-jun 2020 + may 2021" = "Mar--jun 2020 y may 2021",
                 "mar-dic 2020" = "Mar--dic 2020", "mar-dic 2020 + may 2021" = "Mar--dic 2020 y may 2021 (elegida)")
fila_comp <- function(r, etiqueta) fila(etiqueta, r$K, num(r$raiz_max, 3), sprintf("$%.0f$", r$cond_largo_plazo),
                                        num(r$portmanteau_h24, 3), r$lb_rechazos, r$arch_rechazos)
sel_v <- comp[comp$dummies == d_final & comp$p == v$p, ]
sel_d <- comp[comp$variables == "expectativas" & comp$p == v$p, ]
filas_comp <- c(paste0("\\multicolumn{7}{@{}l}{\\textit{A. Variables endógenas (dummies puntuales: mar--dic 2020 y may 2021)}} \\\\"),
                sapply(names(nombres_esp), function(k) fila_comp(sel_v[sel_v$variables == k, ], nombres_esp[k])),
                "\\addlinespace",
                paste0("\\multicolumn{7}{@{}l}{\\textit{B. Dummies puntuales (variables: base + expectativas)}} \\\\"),
                sapply(names(nombres_dum), function(k) fila_comp(sel_d[sel_d$dummies == k, ], nombres_dum[k])))
escribir_cuadro("cuadro_especificaciones.tex", "Especificaciones alternativas del VAR", "tab:especificaciones",
                ">{\\raggedright\\arraybackslash}X ccccccc",
                fila("Especificación", "$K$", "Raíz", "Cond. $I - A(1)$", "Portm. $h=24$ ($p$)", "LB", "ARCH"),
                filas_comp,
                sprintf("Todas las especificaciones con %d rezagos, inflación mensual, dummies estacionales, Brent y Fed en $t$ y $t-1$, y la muestra 2010m1--2026m7. Raíz: módulo de la raíz más grande del VAR. Cond.: número de condición de $I - A_1 - \\dots - A_p$; valores altos indican que los efectos de largo plazo, en los que se basa la identificación, se estiman con poca precisión. LB y ARCH: número de ecuaciones en que se rechaza, al 5\\%%, la ausencia de autocorrelación (Ljung-Box, 12 rezagos) y de heterocedasticidad condicional (ARCH-LM, 12 rezagos). Fuente: \\texttt{prog/03\\_comparacion.R}.", v$p),
                tamano = "\\footnotesize", entorno = "tabularx")
filas_inf <- sapply(sort(unique(comp_inf$p)), function(pp) {
  a <- comp_inf[comp_inf$inflacion == "anual" & comp_inf$p == pp, ]; mm <- comp_inf[comp_inf$inflacion == "mensual" & comp_inf$p == pp, ]
  fila(if (pp == v$p) paste0(pp, "$^{*}$") else pp, num(a$raiz_max, 3), sprintf("$%.0f$", a$cond_largo_plazo), num(a$portmanteau_h24, 3),
       num(mm$raiz_max, 3), sprintf("$%.0f$", mm$cond_largo_plazo), num(mm$portmanteau_h24, 3))
})
lineas <- c("\\begin{table}[H]", "\\centering", "\\caption{VAR con inflación anual frente a inflación mensual}", "\\label{tab:inflacion_anual_mensual}", "\\small",
            "\\begin{tabular}{@{}crrrrrr@{}}", "\\toprule",
            " & \\multicolumn{3}{c}{Inflación anual} & \\multicolumn{3}{c}{Inflación mensual} \\\\",
            "\\cmidrule(lr){2-4}\\cmidrule(l){5-7}",
            fila("Rezagos", "Raíz", "Cond.", "Portm. $h=24$", "Raíz", "Cond.", "Portm. $h=24$"),
            "\\midrule", filas_inf, "\\bottomrule", "\\end{tabular}",
            "\\notascuadro{\\textit{Nota.} Mismas variables, muestra, exógenas y dummies puntuales; solo cambia cómo entra la inflación (la mensual, con dummies estacionales). Raíz: módulo de la raíz más grande; un valor de 1 o mayor indica un VAR no estable. Cond.: número de condición de $I - A(1)$. Portm.: $p$-valor de la prueba Portmanteau ajustada a 24 meses. $^{*}$Especificación elegida. Fuente: \\texttt{prog/03\\_comparacion.R}.}",
            "\\end{table}")
writeLines(lineas, "cuadros/cuadro_inflacion_anual_mensual.tex", useBytes = TRUE)

# ---- Anexo B: estimación -------------------------------------------------------------------

cr <- v$criterios
minimos <- apply(cr[, c("AIC", "HQ", "BIC")], 2, which.min)
Y <- v$Y; X <- v$X
lr_p <- sapply(seq_len(nrow(cr)), function(i) {
  p <- cr[i, "p"]
  if (p == 1) return(NA)
  m0 <- v$estimar_var(Y, X, p - 1, inicio = v$t0); m1 <- v$estimar_var(Y, X, p, inicio = v$t0)
  LR <- (m1$n_obs - m1$n_par) * (log(det(m0$sigma_mv)) - log(det(m1$sigma_mv)))
  1 - pchisq(LR, K^2)
})
filas_cr <- sapply(seq_len(nrow(cr)), function(i) {
  celdas <- sapply(c("AIC", "HQ", "BIC"), function(k) {
    x <- sprintf("%.4f", cr[i, k]); if (minimos[k] == i) paste0("$\\mathbf{", x, "}$") else paste0("$", x, "$")
  })
  marca <- if (cr[i, "p"] == v$p) paste0(cr[i, "p"], "$^{*}$") else cr[i, "p"]
  fila(marca, celdas, num(lr_p[i], 3))
})
escribir_cuadro("cuadro_rezagos.tex", "Criterios de información y pruebas de razón de verosimilitud", "tab:rezagos",
                "crrrr", fila("Rezagos ($p$)", "AIC", "HQ", "BIC", "LR: $p$ vs. $p-1$ ($p$-valor)"), filas_cr,
                sprintf("Muestra común 2010m1--2026m7 ($T = %d$) para todos los rezagos. En negritas, el mínimo de cada criterio. LR: prueba de razón de verosimilitud de $p$ contra $p-1$ rezagos, con la corrección de Sims; $H_0$ = el rezago $p$ no hace falta. $^{*}$Especificación elegida.", m$n_obs))

dg <- v$diagnosticos
ps <- intersect(c("p4", "p5", "p6", "p8"), names(dg))
col_p <- function(f) sapply(ps, function(k) f(dg[[k]]))
rechazan <- function(tab) {
  r <- names(etiquetas)[tab[, "p_valor"] < 0.05]
  if (length(r) == 0) "Ninguna" else paste(etiquetas[r], collapse = ", ")
}
filas_dg <- c(
  fila("Raíz más grande (módulo)", col_p(function(d) num(d$raices[1], 3))),
  fila("Portmanteau, $h = 12$ ($p$)", col_p(function(d) num(d$portmanteau["h = 12", "p_valor"], 3))),
  fila("Portmanteau, $h = 24$ ($p$)", col_p(function(d) num(d$portmanteau["h = 24", "p_valor"], 3))),
  fila("LM, orden 1 ($p$)", col_p(function(d) num(d$lm["orden 1", "p_valor"], 3))),
  fila("LM, orden 12 ($p$)", col_p(function(d) num(d$lm["orden 12", "p_valor"], 3))),
  sapply(seq_len(K), function(j) fila(paste0("Ljung-Box (12), ", etiquetas[j], " ($p$)"),
                                       col_p(function(d) num(d$ljung_box[j, "p_valor"], 3)))),
  fila("Jarque-Bera: rechazan", col_p(function(d) rechazan(d$jarque_bera))),
  fila("ARCH-LM (12): rechazan", col_p(function(d) rechazan(d$arch))),
  fila("F estacionales: rechazan", col_p(function(d) rechazan(d$estacionales))))
encab_dg <- fila("Prueba", paste(sub("p", "", ps), "rezagos"))
escribir_cuadro("cuadro_diagnosticos.tex", "Diagnósticos del VAR según el número de rezagos", "tab:diagnosticos",
                paste0(">{\\raggedright\\arraybackslash}p{4.2cm}", paste(rep(">{\\centering\\arraybackslash}X", length(ps)), collapse = "")),
                encab_dg, filas_dg,
                sprintf("4 rezagos: HQ; 5: AIC; %d: especificación elegida; 8: robustez. El VAR es estable si la raíz más grande es menor que 1. Portmanteau ajustado, LM de autocorrelación (Breusch-Godfrey con la corrección de Edgerton y Shukur) y Ljung-Box: $H_0$ = sin autocorrelación. Jarque-Bera: $H_0$ = normalidad. ARCH-LM: $H_0$ = varianza constante. F: $H_0$ = las 11 dummies estacionales no hacen falta en la ecuación. Jarque-Bera y ARCH-LM excluyen los meses con dummy puntual. Nivel de significancia: 5\\%%.", v$p),
                tamano = "\\footnotesize", entorno = "tabularx")

filas_coef <- c(fila("Constante", num(m$constante)))
for (i in 1:v$p) {
  filas_coef <- c(filas_coef, "\\addlinespace")
  for (j in 1:K) filas_coef <- c(filas_coef, fila(paste0(etiquetas[j], "$_{t-", i, "}$"), num(m$A[, j, i])))
}
escribir_cuadro("cuadro_coeficientes.tex", "Coeficientes estimados del VAR: constante y rezagos", "tab:coeficientes",
                "lrrrrr", fila("Regresor", etiquetas), filas_coef,
                sprintf("Estimación por MCO ecuación por ecuación, 2010m1--2026m7 ($T = %d$), %d regresores por ecuación. Cada columna es una ecuación; el subíndice indica el rezago. Los coeficientes de las variables exógenas y las dummies están en el Cuadro~\\ref{tab:coeficientes_exogenas}.", m$n_obs, m$n_par),
                tamano = "\\footnotesize")

nom_x <- c(d_log_brent = "Brent$_t$", d_log_brent_l1 = "Brent$_{t-1}$", fed = "Fed$_t$", fed_l1 = "Fed$_{t-1}$")
filas_x <- sapply(names(nom_x), function(k) fila(nom_x[k], num(m$gamma[, k])))
filas_x <- c(filas_x, "\\addlinespace",
             sapply(sprintf("est%02d", 2:12), function(k) fila(paste("Estacional", meses_es[as.integer(substr(k, 4, 5))]), num(m$gamma[, k]))),
             "\\addlinespace",
             sapply(seq_along(cod_dummies), function(i) fila(paste("Puntual", meses_es[as.integer(format(meses_dummy[i], "%m"))], format(meses_dummy[i], "%Y")),
                                                            num(m$gamma[, cod_dummies[i]]))))
escribir_cuadro("cuadro_coeficientes_exogenas.tex", "Coeficientes estimados del VAR: exógenas y dummies", "tab:coeficientes_exogenas",
                "lrrrrr", fila("Regresor", etiquetas), filas_x,
                "Mismo VAR que el Cuadro~\\ref{tab:coeficientes}. Las dummies estacionales se miden respecto a enero, que queda en la constante.",
                tamano = "\\footnotesize")

S <- m$sigma; R <- cov2cor(S)
filas_s <- sapply(1:K, function(i) fila(etiquetas[i], sapply(1:K, function(j) {
  if (i == j) paste0("$\\mathbf{", sprintf("%.4f", S[i, j]), "}$") else if (j > i) num(S[i, j], 4) else num(R[i, j], 3)
})))
escribir_cuadro("cuadro_sigma.tex", "Matriz de covarianzas y correlaciones de los residuos", "tab:sigma",
                "lrrrrr", fila("", etiquetas), filas_s,
                sprintf("Diagonal (en negritas) y triángulo superior: covarianzas $\\hat\\Sigma_u$, con corrección por grados de libertad ($T - %d$). Triángulo inferior: correlaciones.", m$n_par))

# ---- Anexo B.4: de la forma reducida a la matriz de impacto ----------------------------------
P_chol <- t(chol(m$sigma))
A1 <- diag(K) - apply(m$A, c(1, 2), sum)
L0 <- solve(A1) %*% P_chol
theta <- atan2(-L0[5, 4], L0[5, 5])
Q <- diag(K); Q[4:5, 4:5] <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2)
H <- P_chol %*% Q
H <- H %*% diag(sign(diag(H)))
LP <- solve(A1) %*% H
etq_choque <- c("$\\varepsilon^{ISE}$", "$\\varepsilon^{\\pi}$", "$\\varepsilon^{\\pi^e}$", "$\\varepsilon^{TPM}$", "$\\varepsilon^{ITCR}$")
celda <- function(x) if (abs(x) < 5e-10) "$0$" else num(x, 3)
filas_H <- c("\\multicolumn{6}{@{}l}{\\textit{A. Impacto: $\\hat H$}} \\\\",
             sapply(1:K, function(i) fila(etiquetas[i], sapply(1:K, function(j) celda(H[i, j])))),
             "\\addlinespace",
             "\\multicolumn{6}{@{}l}{\\textit{B. Largo plazo: $(I - \\hat A(1))^{-1}\\hat H$}} \\\\",
             sapply(1:K, function(i) fila(etiquetas[i], sapply(1:K, function(j) celda(LP[i, j])))))
escribir_cuadro("cuadro_impacto.tex", "Matriz de impacto y efectos de largo plazo estimados", "tab:impacto",
                "lrrrrr", fila("Variable", etq_choque), filas_H,
                sprintf("Cada columna es un choque estructural de una desviación estándar y cada fila, la respuesta de una variable. Panel A: respuesta en el mes del choque; los ceros son las restricciones de corto plazo. Panel B: efecto acumulado de largo plazo; el cero en la fila del ITCR y la columna de la TPM es la restricción de largo plazo. El choque monetario (columna TPM) se obtiene rotando las dos últimas columnas de la descomposición de Cholesky un ángulo de %.1f grados. Un choque de una desviación estándar eleva la TPM en %.3f pp; las respuestas al impulso se reescalan a un alza de 25 pb.", theta * 180 / pi, H[4, 4]),
                tamano = "\\small")

# ---- Figuras ----------------------------------------------------------------------------
azul <- "#1f4e79"
gris_dummy <- "grey85"
sombra_dummies <- function() {
  for (d in meses_dummy) rect(d - 15, par("usr")[3], d + 15, par("usr")[4], col = gris_dummy, border = NA)
}
panel_serie <- function(fechas, x, titulo) {
  plot(fechas, x, type = "n", main = titulo, xlab = "", ylab = "", cex.main = 1, font.main = 1, las = 1)
  sombra_dummies()
  abline(h = 0, lty = 2, col = "grey40")
  lines(fechas, x, lwd = 1.3, col = azul)
}

titulos_series <- c(d_log_ise = "ISE (variación mensual, %)", inflacion_mensual = "Inflación mensual (%)",
                    exp_inflacion_12m = "Expectativa de inflación a 12 meses (%)", tpm = "TPM (%)",
                    d_log_itcr = "ITCR (variación mensual, %)", d_log_brent = "Brent (variación mensual, %)",
                    fed = "Tasa de fondos federales (%)")
png("figuras/fig_series_var.png", width = 2000, height = 2300, res = 220)
par(mfrow = c(4, 2), mar = c(2.6, 3.6, 2.4, 1), mgp = c(2, 0.6, 0))
for (k in names(titulos_series)) panel_serie(muestra$mes, muestra[[k]], titulos_series[k])
plot.new()
legend("center", c("Serie", "Meses con dummy puntual (mar-dic 2020, may 2021)"), lwd = c(1.3, 10),
       col = c(azul, gris_dummy), bty = "n", cex = 1.0)
invisible(dev.off())

# Versión horizontal con las endógenas, para las diapositivas
png("figuras/fig_series_presentacion.png", width = 2400, height = 1350, res = 220)
par(mfrow = c(2, 3), mar = c(2.6, 3.6, 2.4, 1), mgp = c(2, 0.6, 0))
for (k in names(titulos_series)[1:5]) panel_serie(muestra$mes, muestra[[k]], titulos_series[k])
plot.new()
legend("center", c("Serie", "Meses con dummy puntual"), lwd = c(1.3, 10), col = c(azul, gris_dummy), bty = "n", cex = 1.1)
invisible(dev.off())

# Estacionalidad de la inflación mensual: patrón por mes de cada año
infl <- muestra$inflacion_mensual
mes_num <- as.integer(format(muestra$mes, "%m"))
png("figuras/fig_estacionalidad.png", width = 2000, height = 1000, res = 220)
par(mar = c(3.2, 3.6, 1.5, 1), mgp = c(2, 0.6, 0))
plot(mes_num, infl, type = "n", xaxt = "n", xlab = "", ylab = "Inflación mensual (%)", las = 1)
axis(1, at = 1:12, labels = meses_es)
for (a in unique(format(muestra$mes, "%Y"))) {
  s <- format(muestra$mes, "%Y") == a
  lines(mes_num[s], infl[s], col = "grey75", lwd = 1)
}
promedio <- tapply(infl, mes_num, mean)
lines(1:12, promedio, col = azul, lwd = 2.5, type = "o", pch = 19)
abline(h = 0, lty = 2, col = "grey40")
legend("topright", c("Cada año, 2010-2026", "Promedio por mes"), lwd = c(1, 2.5), pch = c(NA, 19),
       col = c("grey75", azul), bty = "n", cex = 0.9)
invisible(dev.off())

U <- m$residuos
fechas_u <- as.Date(paste0(v$fechas[m$filas], "-01"))
png("figuras/fig_residuos.png", width = 2000, height = 1900, res = 220)
par(mfrow = c(3, 2), mar = c(2.6, 3.6, 2.4, 1), mgp = c(2, 0.6, 0))
for (j in 1:K) {
  panel_serie(fechas_u, U[, j], paste("Residuos:", etiquetas[j]))
  abline(h = c(-2, 2) * sqrt(S[j, j]), lty = 3, col = "grey30")
}
plot.new()
legend("center", c("Residuo", "± 2 desviaciones estándar", "Meses con dummy puntual"), lwd = c(1.3, 1, 10), lty = c(1, 3, 1),
       col = c(azul, "grey30", gris_dummy), bty = "n", cex = 1.0)
invisible(dev.off())

cat("Cuadros en cuadros/:", paste(list.files("cuadros"), collapse = ", "), "\n")
cat("Figuras: figuras/fig_series_var.png, figuras/fig_estacionalidad.png, figuras/fig_residuos.png\n")
