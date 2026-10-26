# Propuesta 10: la identificación de la P9 con nuestra forma reducida inicial (más las expectativas)
#
# Identificación: la misma de la P9 (05_propuestas_P8_P9.R).
#   Ceros: el choque monetario no mueve en el mes al ISE, al IPC ni a las expectativas.
#   Signos: en los meses 0 a 5 la TPM no baja y la TRM no sube (el peso se aprecia).
#   Cota: en el impacto el peso no se aprecia más de 10% por cada pp que sube la TPM.
#
# Forma reducida: la de nuestra propuesta inicial (guía de estudio, sección 2), con las expectativas añadidas
# igual que en la forma reducida principal (mediana a 12 meses de la encuesta del Banco, ordenada antes de la TPM).
#   Variables en niveles: 100 x log ISE, 100 x log IPC, expectativas, TPM, 100 x log TRM (nominal, promedio mensual).
#   12 rezagos; constante, 11 dummies mensuales (el IPC no está desestacionalizado) y dummies de impulso de
#   marzo a diciembre de 2020; Brent (100 x log) y tasa de la Fed, en t y t-1.
#   Estimación de 2006m1 a 2026m7; 2005 sólo sirve de rezagos iniciales.
# Todos los datos salen de data/colombia_mensual.csv (01_datos.R).
#
# Se corre desde la carpeta Colombia, después de 04_estructural.R y 05_propuestas_P8_P9.R
# (de ellos sólo se leen la P1 y la P9 para comparar). Si usan RStudio, quiten el # y pongan su ruta:
# setwd("~/Documents/Macroeconometría/Tareas/Ejercicio 2/Colombia")

set.seed(2026)
p         <- 12     # rezagos
hmax      <- 48     # horizonte de las figuras (meses)
hlargo    <- 240    # horizonte largo, para ver si las respuestas regresan a cero (20 años)
M         <- 2000   # sorteos de la forma reducida
n_ang     <- 40     # ángulos por sorteo
kappa     <- 10     # cota: apreciación máxima de impacto, % por cada pp de alza de la TPM (Inf = sin cota)
choque_pp <- 0.25   # alza de 25 pb en la TPM
inicio    <- "2006-01"
pandemia  <- sprintf("2020-%02d", 3:12)

while (sink.number() > 0) sink()
sink("log/06_propuesta_P10.log", split = TRUE)
cat("Propuesta 10 con nuestra forma reducida inicial:", format(Sys.time()), "\n\n")

# ---- 1. Datos -----------------------------------------------------------------------------------
b <- read.csv("data/colombia_mensual.csv")
Y <- cbind(ISE          = b$log_ise,              # ya viene como 100 x log
           IPC          = 100 * log(b$ipc),
           Expectativas = b$exp_inflacion_12m,
           TPM          = b$tpm,
           TRM          = b$log_trm)              # 100 x log de la TRM; si sube, el peso se deprecia
K <- ncol(Y)
nombres <- colnames(Y)
mes <- as.integer(substr(b$fecha, 6, 7))
estacionales <- sapply(2:12, function(j) as.numeric(mes == j))          # 11 dummies mensuales
impulsos     <- sapply(pandemia, function(f) as.numeric(b$fecha == f))  # 10 dummies de pandemia
rezago1 <- function(x) c(NA, head(x, -1))
oil <- 100 * log(b$brent)
X <- cbind(estacionales, impulsos, oil, rezago1(oil), b$fed, rezago1(b$fed))

# ---- 2. Forma reducida por MCO ---------------------------------------------------------------------
filas <- which(b$fecha == inicio):nrow(Y)
Z  <- cbind(1, do.call(cbind, lapply(1:p, function(i) Y[filas - i, ])), X[filas, ])
Yt <- Y[filas, ]
n  <- nrow(Yt)
m  <- ncol(Z)
ZtZi <- solve(crossprod(Z))
Bhat <- ZtZi %*% crossprod(Z, Yt)
U    <- Yt - Z %*% Bhat
S    <- crossprod(U)

matrices_A <- function(B) {        # en B, la fila 1 es la constante y luego vienen los rezagos
  A <- array(NA, c(K, K, p))
  for (i in 1:p) A[, , i] <- t(B[1 + (i - 1) * K + 1:K, ])
  A
}
raiz_max <- function(A) {          # raíz más grande de la matriz compañera (menor que 1 = estable)
  comp <- matrix(0, K * p, K * p)
  comp[1:K, ] <- do.call(cbind, lapply(1:p, function(i) A[, , i]))
  comp[(K + 1):(K * p), 1:(K * (p - 1))] <- diag(K * (p - 1))
  max(Mod(eigen(comp, only.values = TRUE)$values))
}
cat(sprintf("Muestra: %s a %s, n = %d meses, %d regresores por ecuación; raíz más grande %.3f\n\n",
            b$fecha[min(filas)], b$fecha[max(filas)], n, m, raiz_max(matrices_A(Bhat))))

# ---- 3. Respuestas al impulso -----------------------------------------------------------------------
# Recursión de la forma reducida: r_0 = b y r_h = A1 r_(h-1) + ... + Ap r_(h-p).
respuesta <- function(A, b, H) {
  r <- matrix(0, H + 1, K, dimnames = list(NULL, nombres))
  r[1, ] <- b
  for (h in 1:H) for (j in 1:min(h, p)) r[h + 1, ] <- r[h + 1, ] + A[, , j] %*% r[h + 1 - j, ]
  r
}
# Las variables ya están en niveles. Se añaden la inflación anual (para comparar con la P9) y la TRM real
# (TRM menos IPC; los precios de afuera no responden a la política colombiana).
completar <- function(r) {
  inflacion <- r[, 2] - c(rep(0, 12), head(r[, 2], -12))   # IPC hoy menos IPC hace 12 meses
  cbind(r, Inflacion = inflacion, TRM_real = r[, 5] - r[, 2])
}

# ---- 4. Sorteos: forma reducida de la posterior y ángulos al azar -----------------------------------------
LZ <- t(chol(ZtZi))
aceptadas <- list()     # respuestas aceptadas, choque de 1 d.e., meses 0 a hmax
largas    <- list()     # respuestas aceptadas por cada 25 pb, cada 12 meses hasta hlargo (sólo sorteos estables)
angulo    <- c()
estables  <- 0
for (d in 1:M) {
  Sigma <- solve(rWishart(1, n - m, solve(S))[, , 1])
  B     <- Bhat + LZ %*% matrix(rnorm(m * K), m, K) %*% chol(Sigma)
  A     <- matrices_A(B)
  P     <- t(chol(Sigma))
  estable <- raiz_max(A) < 1
  r4 <- completar(respuesta(A, P[, 4], hlargo))   # respuesta al choque 4 de Cholesky
  r5 <- completar(respuesta(A, P[, 5], hlargo))   # respuesta al choque 5 de Cholesky
  for (j in 1:n_ang) {
    th <- runif(1, -pi / 2, pi / 2)
    r  <- cos(th) * r4 + sin(th) * r5
    if (all(r[1:6, 4] >= 0) && all(r[1:6, 5] <= 0) && r[1, 5] >= -kappa * r[1, 4]) {
      aceptadas[[length(aceptadas) + 1]] <- r[1:(hmax + 1), ]
      angulo <- c(angulo, th * 180 / pi)
      if (estable) {
        largas[[length(largas) + 1]] <- r[seq(1, hlargo + 1, by = 12), ] * choque_pp / r[1, 4]
        estables <- estables + 1
      }
    }
  }
}
irf_P10 <- simplify2array(aceptadas)           # (hmax + 1) x 7 x aceptadas, choque de 1 d.e.
lar_P10 <- simplify2array(largas)              # 21 x 7 x aceptadas estables, por cada 25 pb
elast   <- irf_P10[1, 5, ] / irf_P10[1, 4, ]
cat(sprintf("Sorteos: %d de la forma reducida x %d ángulos = %d candidatos\n", M, n_ang, M * n_ang))
cat(sprintf("P10: %d aceptados (%.0f%%); ángulo mediano %.0f grados [p16 %.0f, p84 %.0f]; apreciación de impacto mediana %.1f%% por pp\n",
            dim(irf_P10)[3], 100 * dim(irf_P10)[3] / (M * n_ang), median(angulo), quantile(angulo, .16), quantile(angulo, .84), -median(elast)))
cat(sprintf("De los aceptados, %d (%.0f%%) vienen de un VAR estable; sólo esos se usan para el horizonte largo.\n\n",
            estables, 100 * estables / dim(irf_P10)[3]))

# ---- 5. Cuadros -----------------------------------------------------------------------------------------
cuantil <- function(a, pr) apply(a, c(1, 2), quantile, pr)
q_P10 <- list(p16 = cuantil(irf_P10, .16), p50 = cuantil(irf_P10, .5), p84 = cuantil(irf_P10, .84))
por_pp <- irf_P10
for (s in 1:dim(irf_P10)[3]) por_pp[, , s] <- irf_P10[, , s] / irf_P10[1, 4, s] * choque_pp
q_P10_25 <- list(p16 = cuantil(por_pp, .16), p50 = cuantil(por_pp, .5), p84 = cuantil(por_pp, .84))

# Cuadro 1: choque de una desviación estándar
horizontes <- c(0, 12, 24)
c1 <- c(Propuesta = "P10")
for (k in c("ISE", "Inflacion", "Expectativas", "TRM")) for (h in horizontes) c1[paste0(k, "_", h)] <- sprintf("%.2f", q_P10$p50[h + 1, k])
c1["TPM_0"] <- sprintf("%.2f", q_P10$p50[1, "TPM"])
c1["P(ISE_24<0)"] <- sprintf("%.2f", mean(irf_P10[25, "ISE", ] < 0))
c1["P(IPC_24<0)"] <- sprintf("%.2f", mean(irf_P10[25, "IPC", ] < 0))
c1 <- as.data.frame(t(c1))
cat("Cuadro 1. Mediana posterior, choque de 1 d.e. (ISE, IPC y TRM en nivel, %; inflación anual, expectativas y TPM en pp)\n")
print(c1, row.names = FALSE)

# Cuadro 2: alza de 25 pb, contra la P9 (forma reducida principal) y la P1. El tipo de cambio va en términos reales
# en los tres: TRM menos IPC en la P10, ITCR en la P9 y la P1.
P9 <- readRDS("data/P8_P9.rds")$P9_25pb$p50
P1 <- readRDS(paste0("data/identificacion_", round(100 * choque_pp), "pb.rds"))$irf$P1
c2 <- data.frame(Horizonte = rep(horizontes, each = 3), Propuesta = rep(c("P10", "P9", "P1"), 3))
c2$ISE          <- sprintf("%.2f", c(rbind(q_P10_25$p50[horizontes + 1, "ISE"], P9[horizontes + 1, 1], P1[horizontes + 1, 1])))
c2$Inflacion    <- sprintf("%.2f", c(rbind(q_P10_25$p50[horizontes + 1, "Inflacion"], P9[horizontes + 1, 2], P1[horizontes + 1, 2])))
c2$Expectativas <- sprintf("%.2f", c(rbind(q_P10_25$p50[horizontes + 1, "Expectativas"], P9[horizontes + 1, 3], P1[horizontes + 1, 3])))
c2$TC_real      <- sprintf("%.2f", c(rbind(q_P10_25$p50[horizontes + 1, "TRM_real"], P9[horizontes + 1, 5], P1[horizontes + 1, 5])))
cat(sprintf("\nCuadro 2. Alza de %d pb en la TPM (P10 y P9: mediana; P1: estimación puntual)\n", round(100 * choque_pp)))
print(c2, row.names = FALSE)
cat(sprintf("En la P10, la TRM nominal se mueve %.2f%% en el impacto, %.2f%% a los 12 meses y %.2f%% a los 24.\n",
            q_P10_25$p50[1, "TRM"], q_P10_25$p50[13, "TRM"], q_P10_25$p50[25, "TRM"]))

# Cuadro 3: ¿regresan las respuestas? Por cada 25 pb, mediana y percentiles 16-84, sólo sorteos estables.
anios <- c(2, 4, 5, 10, 20)
c3 <- NULL
for (k in c("ISE", "IPC", "TPM", "TRM", "TRM_real", "Expectativas")) {
  fila <- c(Variable = k)
  for (a in anios) fila[paste0(a, "_anios")] <- sprintf("%.2f [%.2f, %.2f]",
    median(lar_P10[a + 1, k, ]), quantile(lar_P10[a + 1, k, ], .16), quantile(lar_P10[a + 1, k, ], .84))
  c3 <- rbind(c3, fila)
}
c3 <- as.data.frame(c3)
cat(sprintf("\nCuadro 3. ¿Regresan a cero? Respuesta a un alza de %d pb a 2, 4, 5, 10 y 20 años: mediana [p16, p84], sólo sorteos con VAR estable\n", round(100 * choque_pp)))
print(c3, row.names = FALSE)

# Lo mismo con la estimación puntual (MCO), en el ángulo mediano de los aceptados
A_hat <- matrices_A(Bhat)
P_hat <- t(chol(S / (n - m)))
th_med <- median(angulo) * pi / 180
r_hat <- completar(respuesta(A_hat, cos(th_med) * P_hat[, 4] + sin(th_med) * P_hat[, 5], hlargo))
r_hat <- r_hat * choque_pp / r_hat[1, 4]
c3b <- data.frame(Variable = c("ISE", "IPC", "TPM", "TRM", "TRM_real", "Expectativas"))
for (a in anios) c3b[[paste0(a, "_anios")]] <- sprintf("%.2f", r_hat[12 * a + 1, c3b$Variable])
cat(sprintf("\nCuadro 3b. Lo mismo con la estimación puntual (raíz más grande %.3f) y el ángulo mediano (%.0f grados)\n", raiz_max(A_hat), median(angulo)))
print(c3b, row.names = FALSE)

sink()

write.csv(c1, "data/P10_resumen_1de.csv", row.names = FALSE)
write.csv(c2, paste0("data/P10_vs_P9_P1_", round(100 * choque_pp), "pb.csv"), row.names = FALSE)
write.csv(c3, paste0("data/P10_horizonte_largo_", round(100 * choque_pp), "pb.csv"), row.names = FALSE)
write.csv(c3b, paste0("data/P10_horizonte_largo_puntual_", round(100 * choque_pp), "pb.csv"), row.names = FALSE)
saveRDS(list(P10 = q_P10, P10_25pb = q_P10_25, largo_25pb = list(p16 = cuantil(lar_P10, .16), p50 = cuantil(lar_P10, .5), p84 = cuantil(lar_P10, .84)),
             angulo = angulo, elasticidad = elast, kappa = kappa), "data/P10.rds")

# ---- 6. Figuras (mismo estilo que 04 y 05) ------------------------------------------------------------------
sombra <- function(x, lo, hi, col) polygon(c(x, rev(x)), c(lo, rev(hi)), col = col, border = NA)
panel <- function(x, titulo, rangos) {
  plot(x, rep(0, length(x)), type = "n", ylim = range(rangos, 0), xlab = "Meses", ylab = "",
       main = titulo, cex.main = 1, font.main = 1, las = 1)
  abline(h = 0, lty = 2, col = "grey40")
}
meses <- 0:hmax
titulos <- c(ISE = "ISE (nivel, %)", Inflacion = "Inflación anual (pp)", Expectativas = "Expectativas de inflación (pp)",
             TPM = "TPM (pp)", TRM = "TRM (nivel, %; + = depreciación)")

# Figura 1: sólo la P10, alza de 25 pb
png(paste0("figuras/irf_P10_", round(100 * choque_pp), "pb.png"), width = 2000, height = 1900, res = 220)
par(mfrow = c(3, 2), mar = c(3.2, 3.6, 2.4, 1), mgp = c(2, 0.6, 0))
for (k in names(titulos)) {
  panel(meses, titulos[k], c(q_P10_25$p16[, k], q_P10_25$p84[, k]))
  sombra(meses, q_P10_25$p16[, k], q_P10_25$p84[, k], adjustcolor("#2a78d6", 0.30))
  abline(h = 0, lty = 2, col = "grey40")
  lines(meses, q_P10_25$p50[, k], lwd = 2, col = "#184f95")
}
plot.new()
legend("center", c("P10: signos + cota (mediana)", "P10: percentiles 16-84"),
       lwd = c(2, 8), col = c("#184f95", adjustcolor("#2a78d6", 0.30)), bty = "n", cex = 1.05)
dev.off()

# Figura 2: P10 contra P9 (misma identificación, distinta forma reducida), alza de 25 pb.
# El tipo de cambio va en términos reales en las dos: TRM menos IPC en la P10, ITCR en la P9.
P9q <- readRDS("data/P8_P9.rds")$P9_25pb
cols9 <- c(ISE = 1, Inflacion = 2, Expectativas = 3, TPM = 4, TRM_real = 5)
tit2 <- c(titulos[c("ISE", "Inflacion", "Expectativas", "TPM")], TRM_real = "Tipo de cambio real (nivel, %; + = depreciación)")
png(paste0("figuras/irf_P10_vs_P9_", round(100 * choque_pp), "pb.png"), width = 2000, height = 1900, res = 220)
par(mfrow = c(3, 2), mar = c(3.2, 3.6, 2.4, 1), mgp = c(2, 0.6, 0))
for (k in names(tit2)) {
  panel(meses, tit2[k], c(q_P10_25$p16[, k], q_P10_25$p84[, k], P9q$p16[, cols9[k]], P9q$p84[, cols9[k]]))
  sombra(meses, P9q$p16[, cols9[k]], P9q$p84[, cols9[k]], adjustcolor("#9dc3e6", 0.30))
  sombra(meses, q_P10_25$p16[, k], q_P10_25$p84[, k], adjustcolor("#2a78d6", 0.30))
  abline(h = 0, lty = 2, col = "grey40")
  lines(meses, P9q$p50[, cols9[k]], lwd = 2, lty = 2, col = "#1f4e79")
  lines(meses, q_P10_25$p50[, k], lwd = 2, col = "#184f95")
}
plot.new()
legend("center", c("P10: forma reducida en niveles (mediana)", "P9: forma reducida principal (mediana)", "P10: percentiles 16-84", "P9: percentiles 16-84"),
       lwd = c(2, 2, 8, 8), lty = c(1, 2, 1, 1),
       col = c("#184f95", "#1f4e79", adjustcolor("#2a78d6", 0.30), adjustcolor("#9dc3e6", 0.30)), bty = "n", cex = 1.0)
dev.off()

# Figura 3: horizonte largo de la P10 (20 años, un punto por año), alza de 25 pb, sólo sorteos estables
q_lar <- list(p16 = cuantil(lar_P10, .16), p50 = cuantil(lar_P10, .5), p84 = cuantil(lar_P10, .84))
x_lar <- seq(0, hlargo, by = 12)
tit3 <- c(ISE = "ISE (nivel, %)", IPC = "IPC (nivel, %)", Expectativas = "Expectativas de inflación (pp)",
          TPM = "TPM (pp)", TRM_real = "TRM real (nivel, %; + = depreciación)")
png(paste0("figuras/irf_P10_largo_", round(100 * choque_pp), "pb.png"), width = 2000, height = 1900, res = 220)
par(mfrow = c(3, 2), mar = c(3.2, 3.6, 2.4, 1), mgp = c(2, 0.6, 0))
for (k in names(tit3)) {
  panel(x_lar, tit3[k], c(q_lar$p16[, k], q_lar$p84[, k]))
  sombra(x_lar, q_lar$p16[, k], q_lar$p84[, k], adjustcolor("#2a78d6", 0.30))
  abline(h = 0, lty = 2, col = "grey40")
  lines(x_lar, q_lar$p50[, k], lwd = 2, col = "#184f95", type = "o", pch = 16, cex = 0.6)
}
plot.new()
legend("center", c("P10: mediana (un punto por año)", "P10: percentiles 16-84", "Sólo sorteos con VAR estable"),
       lwd = c(2, 8, NA), pch = c(16, NA, NA), col = c("#184f95", adjustcolor("#2a78d6", 0.30), NA), bty = "n", cex = 1.05)
dev.off()
