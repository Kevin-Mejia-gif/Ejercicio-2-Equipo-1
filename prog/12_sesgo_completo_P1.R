# Robustez de la inferencia de la P1: corrección de sesgo completa frente al atajo de 07_inferencia_P1.R
#
# Forma reducida: la de 02_var.R. Identificación: la P1 (tres ceros de impacto y neutralidad de largo plazo del ITCR).
# 07_inferencia_P1.R corrige el sesgo de muestra chica con el atajo de Kilian (1998), el "bootstrap después del
# bootstrap": estima el sesgo de los coeficientes una vez, en la muestra original, y se lo resta a todas las réplicas.
# Kilian y Lütkepohl (2017, pp. 364-365) describen el algoritmo completo, que estima el sesgo de cada réplica con un
# bootstrap anidado; el atajo servía para ahorrar cómputo. Este programa corre los dos sobre las mismas réplicas:
#   1. Etapa 1 (la de 07): sesgo estimado en la muestra original y A corregida, con el ajuste de estacionariedad.
#   2. Etapa 2: B réplicas por bloques de 12 meses, simuladas desde la A corregida con medias restadas (K&L p. 364).
#      A cada réplica se le resta
#        a) el sesgo de la etapa 1, igual para todas (atajo, como en 07), o
#        b) su propio sesgo: el promedio de B_anidado réplicas simuladas desde el VAR estimado en esa réplica, menos
#           su A (completo).
# Bandas de percentiles de Efron de 68%. Como referencia se usan las bandas por bloques sin corrección de 07.
# Advertencia (K&L p. 365): con heterocedasticidad condicional no hay resultados analíticos sobre el sesgo de MCO, y la
# justificación teórica del bootstrap con corrección de sesgo no está establecida. 02_var.R detecta ARCH en las
# ecuaciones del ISE y de la TPM, por eso estas bandas se reportan como robustez.
#
# Se corre desde la carpeta Colombia, después de 07_inferencia_P1.R. Usa varios núcleos con el paquete parallel (viene
# con R); con 7 núcleos tarda unos 17 minutos. Cada réplica tiene su propia semilla, así que el resultado no depende del
# número de núcleos.
# setwd("~/Documents/Macroeconometría/Tareas/Ejercicio 2/Colombia")

hmax      <- 48     # horizonte de las respuestas (meses)
B         <- 2000   # réplicas de la etapa 2 (como en 07)
B_anidado <- 500    # réplicas para estimar el sesgo de cada réplica
l_bloque  <- 12     # largo de los bloques (meses)
choque_pp <- 0.25   # alza de 25 pb en la TPM
nucleos   <- max(1, parallel::detectCores() - 1)

while (sink.number() > 0) sink()
sink("log/12_sesgo_completo_P1.log", split = TRUE)
cat("Corrección de sesgo completa para la Propuesta 1:", format(Sys.time()), "\n\n")

v  <- readRDS("data/var_reducido.rds")
x  <- readRDS("data/P1_inferencia.rds")             # de 07: sesgo de la etapa 1, A corregida y bandas
m  <- v$modelo
p  <- v$p
K  <- ncol(v$Y)
t0 <- v$t0
nombres     <- c("ISE", "Inflación", "Expectativas", "TPM", "ITCR")
estimar_var <- v$estimar_var
raices      <- v$raices
U0    <- m$residuos
filas <- m$filas
A_BC  <- x$A_BC                                     # A corregida de 07
sesgo <- x$sesgo                                    # sesgo de la etapa 1 de 07
# Medias restadas (K&L p. 364), como en 07: la A estimada no cambia, sólo la constante
Yd <- sweep(v$Y, 2, colMeans(v$Y[filas, ]))
Xd <- sweep(v$X, 2, colMeans(v$X[filas, ]))
md <- estimar_var(Yd, Xd, p, inicio = t0)

# ---- 1. Funciones -----------------------------------------------------------------------------------------------
# La P1 como función de la forma reducida, igual que en 07
anual <- function(z) { s <- cumsum(z); s - c(rep(0, 12), head(s, -12)) }
irf_P1 <- function(A, S) {
  P  <- t(chol(S))
  L  <- solve(diag(K) - apply(A, c(1, 2), sum)) %*% P
  th <- atan(-L[5, 4] / L[5, 5])
  r  <- matrix(0, hmax + 1, K)
  r[1, ] <- cos(th) * P[, 4] + sin(th) * P[, 5]
  for (h in 1:hmax) for (j in 1:min(h, p)) r[h + 1, ] <- r[h + 1, ] + A[, , j] %*% r[h + 1 - j, ]
  r <- r * choque_pp / r[1, 4]
  r[, 1] <- cumsum(r[, 1]); r[, 2] <- anual(r[, 2]); r[, 5] <- cumsum(r[, 5])
  list(r = r, angulo = th * 180 / pi)
}
raiz_max <- function(A) raices(A)[1]
# Corrección con ajuste de estacionariedad (K&L p. 364), igual que en 07
corregir <- function(A, s) {
  if (raiz_max(A) >= 1) { attr(A, "fraccion") <- 0; return(A) }
  k <- 0
  repeat { Ac <- A - s * 0.99^k; if (raiz_max(Ac) < 1) break; k <- k + 1 }
  attr(Ac, "fraccion") <- 0.99^k
  Ac
}
# Bloques de 12 meses con recentrado (K&L p. 356), igual que en 07
bloques <- function(U) {
  n <- nrow(U); s <- ceiling(n / l_bloque)
  ini <- sample(n - l_bloque + 1, s, replace = TRUE)
  idx <- as.vector(sapply(ini, function(a) a:(a + l_bloque - 1)))[1:n]
  media <- t(sapply(1:l_bloque, function(i) colMeans(U[i:(i + n - l_bloque), , drop = FALSE])))
  U[idx, ] - media[rep(1:l_bloque, s)[1:n], ]
}
# Simulación recursiva de 07, escrita para que corra más rápido: la parte determinista (constante, exógenas y dummies)
# se calcula una vez y los p rezagos entran en un solo producto de matrices. Abajo se comprueba que da lo mismo.
parte_det <- function(mod, X) sweep(X[filas, , drop = FALSE] %*% t(mod$gamma), 2, mod$constante, "+")
simular <- function(A, D, U, Y) {
  A_todas <- do.call(cbind, lapply(1:p, function(i) A[, , i]))     # [A1 ... Ap]
  for (j in seq_along(filas)) {
    t <- filas[j]
    Y[t, ] <- D[j, ] + A_todas %*% as.vector(t(Y[(t - 1):(t - p), ])) + U[j, ]
  }
  Y
}
estimar <- function(Y) estimar_var(Y, Xd, p, inicio = t0)
# Sesgo de una estimación: se simula desde el VAR estimado en esa muestra, con sus propios residuos por bloques
sesgo_de <- function(mod, Y) {
  D <- parte_det(mod, Xd)
  suma <- 0
  for (i in 1:B_anidado) suma <- suma + estimar(simular(mod$A, D, bloques(mod$residuos), Y))$A
  suma / B_anidado - mod$A
}
persistencia <- function(A) diag(apply(A, c(1, 2), sum))          # suma de coeficientes propios de cada ecuación
D_BC <- parte_det(md, Xd)

# Una réplica de la etapa 2: datos desde la A corregida; la estimación se corrige con el atajo y con su propio sesgo
replica <- function(semilla) {
  assign(".Random.seed", semilla, envir = globalenv())
  Yb <- simular(A_BC, D_BC, bloques(U0), Yd)
  mb <- estimar(Yb)
  if (raiz_max(mb$A) >= 1) return(NULL)             # sin estabilidad no existe el multiplicador de largo plazo
  s_r <- sesgo_de(mb, Yb)
  A_a <- corregir(mb$A, sesgo)                      # atajo
  A_c <- corregir(mb$A, s_r)                        # completo
  list(atajo = irf_P1(A_a, mb$sigma), completo = irf_P1(A_c, mb$sigma),
       raiz = c(sin = raiz_max(mb$A), atajo = raiz_max(A_a), completo = raiz_max(A_c)),
       fraccion = c(atajo = attr(A_a, "fraccion"), completo = attr(A_c, "fraccion")),
       sesgo_persistencia = persistencia(s_r))
}

# Comprobaciones: la simulación rápida reproduce la de 07 y la estimación corregida es la de 07
set.seed(1)
Ub <- bloques(U0)
lento <- v$Y
for (j in seq_along(filas)) {
  t <- filas[j]
  y <- m$constante + m$gamma %*% v$X[t, ] + Ub[j, ]
  for (i in 1:p) y <- y + m$A[, , i] %*% lento[t - i, ]
  lento[t, ] <- y
}
puntual    <- irf_P1(m$A, m$sigma)
puntual_BC <- irf_P1(A_BC, m$sigma)
cat(sprintf("Comprobaciones: simulación rápida frente a la de 07, diferencia máxima %.1e; P1 corregida frente a 07, %.1e\n\n",
            max(abs(simular(m$A, parte_det(m, v$X), Ub, v$Y) - lento)), max(abs(puntual_BC$r - x$puntual_BC$r))))

# ---- 2. Bootstrap -----------------------------------------------------------------------------------------------
# Una semilla por réplica (generador L'Ecuyer, el que usa parallel): el resultado no depende del número de núcleos.
# Al final se regresa al generador anterior, para que los demás scripts den lo mismo si se corren después en la sesión.
generador_previo <- RNGkind()
RNGkind("L'Ecuyer-CMRG")
set.seed(2026)
semillas <- vector("list", B)
semillas[[1]] <- .Random.seed
for (b in 2:B) semillas[[b]] <- parallel::nextRNGStream(semillas[[b - 1]])

t_ini <- Sys.time()
if (nucleos > 1) {
  cl <- parallel::makeCluster(nucleos)
  parallel::clusterExport(cl, c("hmax", "B_anidado", "l_bloque", "choque_pp", "p", "K", "t0", "estimar_var", "raices",
                                "U0", "filas", "A_BC", "sesgo", "Yd", "Xd", "D_BC", "anual", "irf_P1", "raiz_max",
                                "corregir", "bloques", "parte_det", "simular", "estimar", "sesgo_de", "persistencia"),
                          envir = environment())
  res <- parallel::parLapply(cl, semillas, replica)
  parallel::stopCluster(cl)
} else {
  res <- lapply(semillas, replica)
}
RNGkind(generador_previo[1], generador_previo[2], generador_previo[3])
validas <- !sapply(res, is.null)
res <- res[validas]
cat(sprintf("Etapa 2: %d réplicas x %d anidadas en %.1f minutos con %d núcleos. Réplicas con VAR estable: %d\n\n",
            B, B_anidado, as.numeric(difftime(Sys.time(), t_ini, units = "mins")), nucleos, sum(validas)))

irf <- list(atajo    = simplify2array(lapply(res, function(z) z$atajo$r)),      # (hmax + 1) x K x réplicas
            completo = simplify2array(lapply(res, function(z) z$completo$r)))
angulos <- list(atajo = sapply(res, function(z) z$atajo$angulo), completo = sapply(res, function(z) z$completo$angulo))
cuantil <- function(a, pr) apply(a, 1:2, quantile, pr)
bandas  <- lapply(irf, function(a) list(b05 = cuantil(a, .05), b16 = cuantil(a, .16), b50 = cuantil(a, .5),
                                        b84 = cuantil(a, .84), b95 = cuantil(a, .95)))

# ---- 3. Cuadros -------------------------------------------------------------------------------------------------
# Atajo (07) y Atajo (12) usan réplicas distintas: su diferencia es error de simulación.
# Atajo (12) y Completo (12) usan las mismas réplicas: su diferencia es sólo el sesgo que se resta.
todas <- list(bloques_07 = x$bandas$bloques, atajo_07 = x$bandas$sesgo, atajo_12 = bandas$atajo, completo_12 = bandas$completo)
celdas <- list(ITCR_0 = c(5, 0), ITCR_3 = c(5, 3), ITCR_12 = c(5, 12), Inflacion_12 = c(2, 12), Inflacion_24 = c(2, 24),
               Expectativas_6 = c(3, 6), ISE_12 = c(1, 12), ISE_24 = c(1, 24), ISE_36 = c(1, 36))
cuadro <- t(sapply(celdas, function(cc) {
  h <- cc[2] + 1; k <- cc[1]
  c(Estimacion = sprintf("%+.2f", puntual$r[h, k]), Estimacion_corregida = sprintf("%+.2f", puntual_BC$r[h, k]),
    sapply(todas, function(b) sprintf("[%+.2f, %+.2f]", b$b16[h, k], b$b84[h, k])))
}))
cat("Cuadro 1. Respuesta a un alza de 25 pb: estimación y bandas de 68% (variable_mes)\n")
cat("bloques_07: sin corrección | atajo_07: corrección de 07 | atajo_12 y completo_12: este programa, mismas réplicas\n")
print(noquote(cuadro))

ang_todos <- list(bloques_07 = x$angulos$bloques, atajo_07 = x$angulos$sesgo, atajo_12 = angulos$atajo,
                  completo_12 = angulos$completo)
imp_todos <- list(bloques_07 = x$impacto_ITCR$bloques, atajo_07 = x$impacto_ITCR$sesgo, atajo_12 = irf$atajo[1, 5, ],
                  completo_12 = irf$completo[1, 5, ])
cuadro2 <- t(sapply(names(ang_todos), function(k) {
  a <- ang_todos[[k]]
  c(angulo_p16 = unname(round(quantile(a, .16, na.rm = TRUE), 1)), angulo_p50 = round(median(a, na.rm = TRUE), 1),
    angulo_p84 = unname(round(quantile(a, .84, na.rm = TRUE), 1)), deprecian = round(mean(imp_todos[[k]] > 0, na.rm = TRUE), 2))
}))
cat(sprintf("\nCuadro 2. Ángulo de la P1 y réplicas que deprecian al impacto (estimación: %.1f grados; corregida: %.1f)\n",
            puntual$angulo, puntual_BC$angulo))
print(cuadro2)
cat(sprintf("ISE a 24 meses < 0: atajo_12 %.0f%% de las réplicas, completo_12 %.0f%%\n",
            100 * mean(irf$atajo[25, 1, ] < 0), 100 * mean(irf$completo[25, 1, ] < 0)))

# Qué tanto cambia el sesgo de una réplica a otra
sp <- t(sapply(res, `[[`, "sesgo_persistencia"))
s1 <- persistencia(sesgo)
cat("\nCuadro 3. Sesgo de la suma de coeficientes propios: etapa 1 (el que usa el atajo) frente al de cada réplica\n")
for (k in 1:K) cat(sprintf("%-13s etapa 1 %+.3f | réplicas: p16 %+.3f, mediana %+.3f, p84 %+.3f\n", nombres[k], s1[k],
                           quantile(sp[, k], .16), median(sp[, k]), quantile(sp[, k], .84)))
rz <- t(sapply(res, `[[`, "raiz"))
fr <- t(sapply(res, `[[`, "fraccion"))
cat(sprintf("\nRaíz más grande: MCO %.3f, corregida %.3f. En las réplicas, mediana (p84): sin corregir %.3f (%.3f),",
            raiz_max(m$A), raiz_max(A_BC), median(rz[, "sin"]), quantile(rz[, "sin"], .84)),
    sprintf("atajo %.3f (%.3f), completo %.3f (%.3f)\n", median(rz[, "atajo"]), quantile(rz[, "atajo"], .84),
            median(rz[, "completo"]), quantile(rz[, "completo"], .84)))
cat(sprintf("Réplicas en que el ajuste de estacionariedad encoge el sesgo: atajo %.0f%%, completo %.0f%%\n",
            100 * mean(fr[, "atajo"] < 1), 100 * mean(fr[, "completo"] < 1)))

saveRDS(list(puntual = puntual, puntual_BC = puntual_BC, bandas = bandas, angulos = angulos, raices = rz, fraccion = fr,
             sesgo_persistencia = sp, impacto_ITCR = list(atajo = irf$atajo[1, 5, ], completo = irf$completo[1, 5, ]),
             B = B, B_anidado = B_anidado), "data/P1_sesgo_completo.rds")
write.csv(data.frame(celda = names(celdas), cuadro, check.names = FALSE), "data/P1_bandas_sesgo_completo_25pb.csv",
          row.names = FALSE)

# ---- 4. Figura: bandas de 68% sin corrección (07), con el atajo y con la corrección completa ----------------------
titulos <- c("ISE (nivel, %)", "Inflación anual (pp)", "Expectativas de inflación (pp)", "TPM (pp)",
             "ITCR (nivel, %; + = depreciación)")
morado <- "#5b9bd5"
azul   <- "#1f4e79"
bq <- x$bandas$bloques; ba <- bandas$atajo; bc <- bandas$completo
png("figuras/irf_P1_sesgo_completo_25pb.png", width = 2000, height = 1900, res = 220)
par(mfrow = c(3, 2), mar = c(3.2, 3.6, 2.4, 1), mgp = c(2, 0.6, 0))
for (k in 1:K) {
  rango <- range(bq$b16[, k], bq$b84[, k], ba$b16[, k], ba$b84[, k], bc$b16[, k], bc$b84[, k], puntual$r[, k],
                 puntual_BC$r[, k], 0)
  plot(0:hmax, rep(0, hmax + 1), type = "n", ylim = rango, xlab = "Meses", ylab = "", main = titulos[k],
       cex.main = 1, font.main = 1, las = 1)
  polygon(c(0:hmax, hmax:0), c(ba$b16[, k], rev(ba$b84[, k])), col = adjustcolor(morado, 0.25), border = NA)
  abline(h = 0, lty = 2, col = "grey40")
  for (lim in list(bq$b16, bq$b84)) lines(0:hmax, lim[, k], col = "grey45", lwd = 1.5, lty = 2)
  for (lim in list(bc$b16, bc$b84)) lines(0:hmax, lim[, k], col = azul, lwd = 2)
  lines(0:hmax, puntual$r[, k], lwd = 2)
  lines(0:hmax, puntual_BC$r[, k], lwd = 2, lty = 3)
}
plot.new()
legend("center", c("Estimación (MCO)", "Estimación con corrección de sesgo", "Bloques sin corrección: banda 68%",
                   "Corrección con el atajo: banda 68%", "Corrección completa: banda 68%"),
       lwd = c(2, 2, 1.5, 8, 2), lty = c(1, 3, 2, 1, 1), col = c("black", "black", "grey45", adjustcolor(morado, 0.25), azul),
       bty = "n", cex = 1.0)
dev.off()

cat("\nListo: figuras/irf_P1_sesgo_completo_25pb.png; data/P1_sesgo_completo.rds y data/P1_bandas_sesgo_completo_25pb.csv\n")
sink()
