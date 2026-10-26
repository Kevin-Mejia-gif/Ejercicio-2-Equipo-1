# Propuestas 8 y 9 (combinación 2: corto plazo + signos) con la forma reducida de 02_var.R
#
# La forma reducida es la misma de 02_var.R: VAR(6), cinco variables (ISE, inflación mensual, expectativas, TPM, ITCR),
# 2010m1-2026m7, Brent y Fed exógenas, 11 dummies estacionales y 11 de impulso. Sólo cambia la identificación del
# choque monetario. La respuesta de la inflación se presenta como inflación anual (suma de 12 meses), igual que en 04.
#
#   Ceros (corto plazo): el choque monetario no mueve en el mes al ISE, a la inflación ni a las expectativas
#     (los mismos ceros del patrón de Rubén). Con esos ceros, la columna 4 de H es una mezcla de las columnas
#     4 y 5 de Cholesky con un ángulo theta:  H[, 4] = cos(theta) * P[, 4] + sin(theta) * P[, 5].
#   Signos: en los meses 0 a 5 la TPM no baja y el ITCR en nivel no sube (el peso se aprecia).
#     ISE, inflación y expectativas quedan libres: son lo que queremos aprender.
#   P8: sólo ceros y signos.
#   P9: además, una cota: en el impacto el peso no se aprecia más de 10% por cada pp que sube la TPM
#     (Kilian y Murphy, 2012; Kilian y Lütkepohl, sección 13.4.2).
#
# Inferencia bayesiana (Kilian y Lütkepohl, sección 13.6): se sortea la forma reducida de su posterior
# y, para cada sorteo, 40 ángulos al azar entre -90 y 90 grados; se guardan los que cumplen los signos.
#
# Se corre desde la carpeta Colombia, después de 02_var.R. Si usan RStudio, quiten el # y pongan su ruta:
# setwd("~/Documents/Macroeconometría/Tareas/Ejercicio 2/Colombia")

set.seed(2026)
hmax      <- 48     # horizonte de las respuestas (meses)
M         <- 2000   # sorteos de la forma reducida
n_ang     <- 40     # ángulos por sorteo
kappa     <- 10     # cota: apreciación máxima de impacto, % por cada pp de alza de la TPM
choque_pp <- 0.25   # para comparar con la Propuesta 1: alza de 25 pb en la TPM

while (sink.number() > 0) sink()
sink("log/05_propuestas_P8_P9.log", split = TRUE)
cat("Propuestas 8 y 9 con la forma reducida de 02_var.R:", format(Sys.time()), "\n\n")

v <- readRDS("data/var_reducido.rds")
Y <- v$Y
X <- v$X
p <- v$p
K <- ncol(Y)
nombres <- c("ISE", "Inflación", "Expectativas", "TPM", "ITCR")

# ---- 1. Forma reducida (la misma de 02_var.R), escrita como regresión -------------------------
filas <- v$t0:nrow(Y)
Z  <- cbind(1, do.call(cbind, lapply(1:p, function(i) Y[filas - i, ])), X[filas, ])
Yt <- Y[filas, ]
n  <- nrow(Yt)
m  <- ncol(Z)
ZtZi <- solve(crossprod(Z))
Bhat <- ZtZi %*% crossprod(Z, Yt)       # coeficientes por MCO (m x K)
U    <- Yt - Z %*% Bhat
S    <- crossprod(U)                     # suma de residuos al cuadrado (K x K)

# Matrices A1, ..., Ap a partir de B (en B, la fila 1 es la constante y luego vienen los rezagos)
matrices_A <- function(B) {
  A <- array(NA, c(K, K, p))
  for (i in 1:p) A[, , i] <- t(B[1 + (i - 1) * K + 1:K, ])
  A
}
cat(sprintf("Comprobación: diferencia máxima con los coeficientes de 02_var.R = %.1e\n\n",
            max(abs(matrices_A(Bhat) - v$modelo$A))))

# ---- 2. Respuestas al impulso ------------------------------------------------------------------
# Matrices Psi_h de la representación de medias móviles (igual que en 04_estructural.R)
psis <- function(A) {
  comp <- matrix(0, K * p, K * p)
  comp[1:K, ] <- do.call(cbind, lapply(1:p, function(i) A[, , i]))
  comp[(K + 1):(K * p), 1:(K * (p - 1))] <- diag(K * (p - 1))
  Psi <- list()
  Mc <- diag(K * p)
  for (h in 0:hmax) {
    Psi[[h + 1]] <- Mc[1:K, 1:K]
    Mc <- comp %*% Mc
  }
  Psi
}

# Respuesta a un choque con vector de impacto b. ISE e ITCR se acumulan para leerse en nivel; la inflación mensual se
# pasa a anual sumando las respuestas de los últimos 12 meses.
anual <- function(x) { s <- cumsum(x); s - c(rep(0, 12), head(s, -12)) }
respuesta <- function(Psi, b) {
  r <- t(sapply(Psi, function(Ps) Ps %*% b))   # (hmax + 1) x K
  r[, 1] <- cumsum(r[, 1])
  r[, 2] <- anual(r[, 2])
  r[, 5] <- cumsum(r[, 5])
  r
}

# ---- 3. Sorteos: forma reducida de la posterior y ángulos al azar --------------------------------
# Prior difuso: Sigma ~ Wishart inversa(S, n - m) y, dado Sigma, vec(B) ~ Normal(vec(Bhat), Sigma x (Z'Z)^-1)
LZ <- t(chol(ZtZi))
aceptadas <- list()   # respuestas que cumplen los signos (choque de 1 d.e.)
angulo    <- c()      # ángulo de cada respuesta aceptada, en grados
for (d in 1:M) {
  Sigma <- solve(rWishart(1, n - m, solve(S))[, , 1])
  B     <- Bhat + LZ %*% matrix(rnorm(m * K), m, K) %*% chol(Sigma)
  P     <- t(chol(Sigma))
  Psi   <- psis(matrices_A(B))
  r4 <- respuesta(Psi, P[, 4])    # respuesta al choque 4 de Cholesky
  r5 <- respuesta(Psi, P[, 5])    # respuesta al choque 5 de Cholesky
  for (j in 1:n_ang) {
    th <- runif(1, -pi / 2, pi / 2)
    r  <- cos(th) * r4 + sin(th) * r5          # respuesta al choque con ángulo th
    if (all(r[1:6, 4] >= 0) && all(r[1:6, 5] <= 0)) {
      aceptadas[[length(aceptadas) + 1]] <- r
      angulo <- c(angulo, th * 180 / pi)
    }
  }
}
irf_P8 <- simplify2array(aceptadas)                     # (hmax + 1) x K x aceptadas
elast <- irf_P8[1, 5, ] / irf_P8[1, 4, ]                 # % que se mueve el ITCR por cada pp de TPM, en el impacto
en_P9  <- elast >= -kappa
irf_P9 <- irf_P8[, , en_P9]

cat(sprintf("Sorteos: %d de la forma reducida x %d ángulos = %d candidatos\n", M, n_ang, M * n_ang))
cat(sprintf("P8 (signos):        %5d aceptados (%.0f%%); ángulo mediano %.0f grados [p16 %.0f, p84 %.0f]\n",
            dim(irf_P8)[3], 100 * dim(irf_P8)[3] / (M * n_ang), median(angulo), quantile(angulo, .16), quantile(angulo, .84)))
cat(sprintf("P9 (signos + cota): %5d aceptados (%.0f%%); ángulo mediano %.0f grados [p16 %.0f, p84 %.0f]\n",
            dim(irf_P9)[3], 100 * dim(irf_P9)[3] / (M * n_ang), median(angulo[en_P9]), quantile(angulo[en_P9], .16), quantile(angulo[en_P9], .84)))
cat(sprintf("Apreciación de impacto por cada pp de TPM: mediana %.1f%% en P8 (%.0f%% de los sorteos pasan de %d%%), %.1f%% en P9\n\n",
            -median(elast), 100 * mean(elast < -kappa), kappa, -median(elast[en_P9])))

# ---- 4. Cuadros ------------------------------------------------------------------------------------
cuantil <- function(a, pr) apply(a, c(1, 2), quantile, pr)
q_P8 <- list(p16 = cuantil(irf_P8, .16), p50 = cuantil(irf_P8, .5), p84 = cuantil(irf_P8, .84))
q_P9 <- list(p16 = cuantil(irf_P9, .16), p50 = cuantil(irf_P9, .5), p84 = cuantil(irf_P9, .84))

# Cuadro 1: choque de una desviación estándar (unidad natural cuando hay un conjunto de ángulos)
horizontes <- c(0, 12, 24)
cuadro <- function(q, a, etiqueta) {
  fila <- c(Version = etiqueta)
  for (k in c(1, 2, 3, 5)) for (h in horizontes) fila[paste0(nombres[k], "_", h)] <- sprintf("%.2f", q$p50[h + 1, k])
  fila["TPM_0"] <- sprintf("%.2f", q$p50[1, 4])
  fila["P(ISE_24<0)"] <- sprintf("%.2f", mean(a[25, 1, ] < 0))
  fila["P(Infl_24<0)"] <- sprintf("%.2f", mean(a[25, 2, ] < 0))
  fila
}
c1 <- as.data.frame(rbind(cuadro(q_P8, irf_P8, "P8: signos"), cuadro(q_P9, irf_P9, "P9: signos + cota")))
cat("Cuadro 1. Mediana posterior, choque de 1 d.e. (ISE e ITCR en nivel, %; inflación, expectativas y TPM en pp)\n")
print(c1, row.names = FALSE)

# Cuadro 2: P9 por cada 25 pb, contra la Propuesta 1 (04_estructural.R).
# Con la cota la TPM siempre sube claramente en el impacto, así que normalizar a 25 pb es estable.
# En la P8 no se hace: hay ángulos en los que la TPM casi no se mueve y, al dividir, las respuestas se disparan.
por_pp <- irf_P9
for (s in 1:dim(irf_P9)[3]) por_pp[, , s] <- irf_P9[, , s] / irf_P9[1, 4, s] * choque_pp
q_P9_25 <- list(p16 = cuantil(por_pp, .16), p50 = cuantil(por_pp, .5), p84 = cuantil(por_pp, .84))
p1 <- readRDS(paste0("data/identificacion_", round(100 * choque_pp), "pb.rds"))$irf$P1
c2 <- data.frame(Horizonte = rep(horizontes, each = 2), Esquema = rep(c("P9: signos + cota", "Propuesta 1"), 3))
for (k in c(1, 2, 3, 5)) c2[[nombres[k]]] <- sprintf("%.2f", c(rbind(q_P9_25$p50[horizontes + 1, k], p1[horizontes + 1, k])))
cat(sprintf("\nCuadro 2. Alza de %d pb en la TPM: P9 (mediana) contra la Propuesta 1 (estimación puntual)\n", round(100 * choque_pp)))
print(c2, row.names = FALSE)

sink()

write.csv(c1, "data/P8_P9_resumen_1de.csv", row.names = FALSE)
write.csv(c2, paste0("data/P9_vs_P1_", round(100 * choque_pp), "pb.csv"), row.names = FALSE)
saveRDS(list(P8 = q_P8, P9 = q_P9, P9_25pb = q_P9_25, angulo = angulo, elasticidad = elast, kappa = kappa),
        "data/P8_P9.rds")

# ---- 5. Figuras (mismo estilo que 04_estructural.R) ----------------------------------------------------
titulos <- c("ISE (nivel, %)", "Inflación anual (pp)", "Expectativas de inflación (pp)", "TPM (pp)",
             "ITCR (nivel, %; + = depreciación)")
sombra <- function(lo, hi, col) polygon(c(0:hmax, hmax:0), c(lo, rev(hi)), col = col, border = NA)
panel <- function(k, rangos) {
  plot(0:hmax, rep(0, hmax + 1), type = "n", ylim = range(rangos, 0), xlab = "Meses", ylab = "",
       main = titulos[k], cex.main = 1, font.main = 1, las = 1)
  abline(h = 0, lty = 2, col = "grey40")
}

# Figura 1: P8 contra P9, choque de 1 d.e.
png("figuras/irf_P8_vs_P9_1de.png", width = 2000, height = 1900, res = 220)
par(mfrow = c(3, 2), mar = c(3.2, 3.6, 2.4, 1), mgp = c(2, 0.6, 0))
for (k in 1:K) {
  panel(k, c(q_P8$p16[, k], q_P8$p84[, k], q_P9$p16[, k], q_P9$p84[, k]))
  sombra(q_P8$p16[, k], q_P8$p84[, k], adjustcolor("grey40", 0.25))
  sombra(q_P9$p16[, k], q_P9$p84[, k], adjustcolor("#9dc3e6", 0.6))
  abline(h = 0, lty = 2, col = "grey40")
  lines(0:hmax, q_P8$p50[, k], lwd = 2)
  lines(0:hmax, q_P9$p50[, k], lwd = 2, lty = 2, col = "#1f4e79")
}
plot.new()
legend("center", c("P8: signos (mediana)", "P9: signos + cota (mediana)", "P8: percentiles 16-84", "P9: percentiles 16-84"),
       lwd = c(2, 2, 8, 8), lty = c(1, 2, 1, 1),
       col = c("black", "#1f4e79", adjustcolor("grey40", 0.25), adjustcolor("#9dc3e6", 0.6)), bty = "n", cex = 1.05)
dev.off()

# Figura 2: P9 contra la Propuesta 1, alza de 25 pb
png(paste0("figuras/irf_P9_vs_P1_", round(100 * choque_pp), "pb.png"), width = 2000, height = 1900, res = 220)
par(mfrow = c(3, 2), mar = c(3.2, 3.6, 2.4, 1), mgp = c(2, 0.6, 0))
for (k in 1:K) {
  panel(k, c(q_P9_25$p16[, k], q_P9_25$p84[, k], p1[, k]))
  sombra(q_P9_25$p16[, k], q_P9_25$p84[, k], adjustcolor("#9dc3e6", 0.6))
  abline(h = 0, lty = 2, col = "grey40")
  lines(0:hmax, q_P9_25$p50[, k], lwd = 2, col = "#1f4e79")
  lines(0:hmax, p1[, k], lwd = 2, lty = 2, col = "grey25")
}
plot.new()
legend("center", c("P9: signos + cota (mediana)", "P9: percentiles 16-84", "Propuesta 1 (estimación puntual)"),
       lwd = c(2, 8, 2), lty = c(1, 1, 2), col = c("#1f4e79", adjustcolor("#9dc3e6", 0.6), "grey25"), bty = "n", cex = 1.05)
dev.off()

# Figura 3: sólo la P9, alza de 25 pb
png(paste0("figuras/irf_P9_", round(100 * choque_pp), "pb.png"), width = 2000, height = 1900, res = 220)
par(mfrow = c(3, 2), mar = c(3.2, 3.6, 2.4, 1), mgp = c(2, 0.6, 0))
for (k in 1:K) {
  panel(k, c(q_P9_25$p16[, k], q_P9_25$p84[, k]))
  sombra(q_P9_25$p16[, k], q_P9_25$p84[, k], adjustcolor("#9dc3e6", 0.6))
  abline(h = 0, lty = 2, col = "grey40")
  lines(0:hmax, q_P9_25$p50[, k], lwd = 2, col = "#1f4e79")
}
plot.new()
legend("center", c("P9: signos + cota (mediana)", "P9: percentiles 16-84"),
       lwd = c(2, 8), col = c("#1f4e79", adjustcolor("#9dc3e6", 0.6)), bty = "n", cex = 1.05)
dev.off()
