# Figuras didácticas para el documento de inferencia de la Propuesta 1 (Inferencia_P1.pdf)
#
# No estima nada nuevo: ilustra los conceptos del documento con los datos y los resultados de 07_inferencia_P1.R.
#   1. ARCH simulado contra ruido iid (mismo ruido de base, misma varianza incondicional).
#   2. Autocorrelaciones de los residuos y de sus cuadrados (ISE, TPM, ITCR).
#   3. Una réplica de los residuos de la TPM con cada método de remuestreo.
#   4. Sesgo de MCO en un AR(1) persistente y sensibilidad del multiplicador de largo plazo 1 / (1 - rho).
#   5. Por qué la corrección de sesgo mueve el ángulo de la P1 y casi no el conjunto de la P9.
#   6. Intervalos de Efron y de Hall sobre la distribución bootstrap de la respuesta de impacto del ITCR.
#   7. Bandas mes a mes contra la banda conjunta de Bonferroni (inflación y expectativas).
#
# Se corre desde la carpeta Colombia, después de 07_inferencia_P1.R (tarda unos segundos).
# setwd("~/Documents/Macroeconometría/Tareas/Ejercicio 2/Colombia")

set.seed(2026)
v  <- readRDS("data/var_reducido.rds")
x  <- readRDS("data/P1_inferencia.rds")
m  <- v$modelo
p  <- v$p
K  <- ncol(v$Y)
U0 <- m$residuos
Tn <- nrow(U0)
fechas  <- v$fechas[m$filas]
anio    <- as.numeric(substr(fechas, 1, 4)) + (as.numeric(substr(fechas, 6, 7)) - 1) / 12
activos <- rowSums(abs(U0)) > 1e-10            # meses sin dummy (residuo distinto de cero)
naranja <- "#1f4e79"; azul <- "grey25"; morado <- "#5b9bd5"
ac <- function(z, j) cor(z[-(1:j)], z[1:(length(z) - j)])

# ---- 1. ARCH simulado contra ruido iid --------------------------------------------------------------------
# Mismos epsilon para las dos series; ARCH(1) con omega = 0.5 y a = 0.5, de modo que la varianza incondicional es 1.
n_sim <- 199; eps <- rnorm(n_sim + 50)
u_arch <- numeric(n_sim + 50); s2 <- 1
for (t in 1:(n_sim + 50)) { u_arch[t] <- sqrt(s2) * eps[t]; s2 <- 0.5 + 0.5 * u_arch[t]^2 }
u_arch <- u_arch[-(1:50)]; u_iid <- eps[-(1:50)]
png("figuras/P1_ARCH_simulado.png", width = 2200, height = 1300, res = 220)
par(mfrow = c(2, 2), mar = c(3, 3.6, 2.4, 1), mgp = c(2, 0.6, 0))
for (s in list(list(u_iid, "Ruido iid"), list(u_arch, "ARCH(1), a = 0.5"))) {
  plot(s[[1]], type = "h", col = "grey45", las = 1, xlab = "", ylab = "", ylim = c(-6, 6),
       main = sprintf("%s: u_t   (autocorrelación %.2f)", s[[2]], ac(s[[1]], 1)), cex.main = 0.95, font.main = 1)
  abline(h = 0, col = "grey60")
}
for (s in list(list(u_iid, "Ruido iid"), list(u_arch, "ARCH(1), a = 0.5"))) {
  plot(s[[1]]^2, type = "h", col = naranja, las = 1, xlab = "Mes", ylab = "", ylim = c(0, 25),
       main = sprintf("%s: u_t al cuadrado   (autocorrelación %.2f)", s[[2]], ac(s[[1]]^2, 1)), cex.main = 0.95, font.main = 1)
}
dev.off()

# ---- 2. Autocorrelaciones de los residuos y de sus cuadrados ----------------------------------------------
png("figuras/P1_autocorrelaciones.png", width = 2200, height = 850, res = 220)
par(mfrow = c(1, 3), mar = c(3.4, 3.4, 2.4, 0.8), mgp = c(2, 0.6, 0))
lim <- 1.96 / sqrt(sum(activos))
for (k in c(1, 4, 5)) {
  u <- U0[activos, k]
  a <- rbind(sapply(1:12, function(j) ac(u, j)), sapply(1:12, function(j) ac(u^2, j)))
  barplot(a, beside = TRUE, names.arg = 1:12, col = c("grey60", naranja), border = NA, ylim = c(-0.25, 0.4), las = 1,
          xlab = "Rezago (meses)", main = c("ISE", "", "", "TPM", "ITCR")[k], cex.main = 1, font.main = 1, cex.names = 0.8)
  abline(h = c(-lim, lim), lty = 2, col = "grey30"); abline(h = 0)
  if (k == 1) legend("topright", c("residuos", "residuos al cuadrado", "±1.96 / raíz de n"), fill = c("grey60", naranja, NA),
                     border = NA, lty = c(NA, NA, 2), bty = "n", cex = 0.85)
}
dev.off()

# ---- 3. Una réplica de los residuos de la TPM con cada método -----------------------------------------------
# Las barras naranjas son los residuos que en los datos caen entre julio de 2021 y junio de 2023 (la racha de volatilidad
# de la TPM): se ve a dónde los manda cada método.
racha <- anio >= 2021.5 & anio < 2023.5
idx_iid <- sample(Tn, Tn, replace = TRUE)
signo   <- sample(c(-1, 1), Tn, replace = TRUE)
l <- 12; s <- ceiling(Tn / l); ini <- sample(Tn - l + 1, s, replace = TRUE)
idx_blq <- as.vector(sapply(ini, function(a) a:(a + l - 1)))[1:Tn]
media   <- sapply(1:l, function(i) mean(U0[i:(i + Tn - l), 4]))
replicas <- list(
  list(U0[, 4], 1:Tn, "Datos: residuos de la ecuación de la TPM"),
  list(U0[idx_iid, 4], idx_iid, "Réplica iid: cada mes se saca al azar"),
  list(U0[, 4] * signo, 1:Tn, "Réplica wild: cada residuo se queda en su mes y sólo cambia de signo"),
  list(U0[idx_blq, 4] - media[rep(1:l, s)[1:Tn]], idx_blq, "Réplica por bloques de 12 meses: se mueven bloques completos"))
png("figuras/P1_replicas_residuos.png", width = 2200, height = 1900, res = 220)
par(mfrow = c(4, 1), mar = c(2.4, 3.6, 2.2, 1), mgp = c(2, 0.6, 0))
for (r in replicas) {
  plot(anio, r[[1]], type = "h", lwd = 2, col = ifelse(racha[r[[2]]], naranja, "grey55"), las = 1, xlab = "", ylab = "",
       ylim = c(-1, 1) * max(abs(U0[, 4])) * 1.05, main = r[[3]], cex.main = 1, font.main = 1)
  abline(h = 0, col = "grey40")
}
dev.off()

# ---- 4. Sesgo de MCO en un AR(1) persistente y el multiplicador 1 / (1 - rho) ---------------------------------
rho <- round(v$raices(m$A)[1], 3); Tsim <- Tn; nsim <- 5000   # la raíz más grande del VAR estimado y su número de meses
rho_BC <- round(v$raices(x$A_BC)[1], 3)                          # la misma raíz con la corrección de sesgo (07)
est <- replicate(nsim, {
  y <- numeric(Tsim + 100); e <- rnorm(Tsim + 100)
  for (t in 2:(Tsim + 100)) y[t] <- rho * y[t - 1] + e[t]
  y <- y[-(1:100)]
  coef(lm(y[-1] ~ y[-Tsim]))[2]                 # MCO con constante, como en el VAR
})
aprox <- rho - (1 + 3 * rho) / Tsim               # aproximación de Marriott y Pope (1954)
png("figuras/P1_sesgo_AR1.png", width = 2200, height = 950, res = 220)
par(mfrow = c(1, 2), mar = c(3.6, 3.8, 2.4, 1), mgp = c(2.3, 0.6, 0))
hist(est, breaks = 60, col = "grey80", border = "white", freq = FALSE, xlim = c(0.8, 1), las = 1,
     xlab = "Estimación de MCO del coeficiente", ylab = "Densidad",
     main = sprintf("(a) AR(1) con rho = %.3f y T = %d", rho, Tsim), cex.main = 1, font.main = 1)
abline(v = rho, lwd = 2); abline(v = mean(est), lwd = 2, col = naranja, lty = 2)
legend("topleft", c(sprintf("Verdadero: %.3f", rho), sprintf("Promedio de MCO: %.3f", mean(est)),
                    sprintf("Marriott y Pope: %.3f", aprox)), lwd = c(2, 2, NA), lty = c(1, 2, NA),
       col = c("black", naranja, NA), bty = "n", cex = 0.85)
rr <- seq(0.85, 0.999, length.out = 400)
plot(rr, 1 / (1 - rr), type = "l", lwd = 2, log = "y", las = 1, xlab = "rho", ylab = "1 / (1 - rho), escala logarítmica",
     main = "(b) Efecto de largo plazo de un AR(1)", cex.main = 1, font.main = 1)
pts <- c(aprox, rho, rho_BC); cols_p <- c(naranja, "black", morado)
points(pts, 1 / (1 - pts), pch = 19, col = cols_p, cex = 1.2)
text(pts, 1 / (1 - pts) * c(1, 0.85, 1.35), sprintf("%.3f: %.0f", pts, 1 / (1 - pts)), pos = c(4, 2, 4), cex = 0.85, col = cols_p)
dev.off()

# ---- 5. El ángulo de la P1 contra el conjunto de la P9, con y sin corrección --------------------------------
P   <- t(chol(m$sigma))
irf <- function(A, b, H = 5) { r <- matrix(0, H + 1, K); r[1, ] <- b
  for (h in 1:H) for (j in 1:min(h, p)) r[h + 1, ] <- r[h + 1, ] + A[, , j] %*% r[h + 1 - j, ]
  r[, 5] <- cumsum(r[, 5]); r }
arco_P9 <- function(A) {
  g <- seq(-89.9, 89.9, by = 0.1); r4 <- irf(A, P[, 4]); r5 <- irf(A, P[, 5])
  ok <- sapply(g, function(th) { r <- cos(th * pi / 180) * r4 + sin(th * pi / 180) * r5
    all(r[, 4] >= 0) && all(r[, 5] <= 0) && r[1, 5] >= -10 * r[1, 4] })
  range(g[ok]) }
direccion <- function(A, th) { L <- solve(diag(K) - apply(A, c(1, 2), sum)) %*% P
  (cos(th) * L[5, 4] + sin(th) * L[5, 5]) / sqrt(L[5, 4]^2 + L[5, 5]^2) }
arco_ls <- arco_P9(m$A); arco_bc <- arco_P9(x$A_BC)
th_ls <- x$puntual$angulo; th_bc <- x$puntual_BC$angulo
g  <- seq(-60, 80, by = 0.25)
eta <- 0.25 * (P[5, 4] + P[5, 5] * tan(g * pi / 180)) / P[4, 4]      # ITCR al impacto por 25 pb (no depende de A)
eta_P1 <- 0.25 * (P[5, 4] + P[5, 5] * tan(c(th_ls, th_bc) * pi / 180)) / P[4, 4]
lim_a <- range(eta[g >= arco_ls[1] - 10 & g <= max(th_ls, th_bc) + 10], eta_P1, 0)
lim_a <- lim_a + c(-0.15, 0.25) * diff(lim_a)
png("figuras/P1_angulo_P1_P9.png", width = 2200, height = 950, res = 220)
par(mfrow = c(1, 2), mar = c(3.6, 3.8, 2.4, 1), mgp = c(2.3, 0.6, 0))
plot(g, eta, type = "n", ylim = lim_a, las = 1, xlab = "Ángulo (grados)", ylab = "ITCR al impacto (%, por 25 pb)",
     main = "(a) Respuesta de impacto según el ángulo", cex.main = 1, font.main = 1)
rect(arco_ls[1], lim_a[1] - 1, arco_ls[2], lim_a[2] + 1, col = adjustcolor("grey60", 0.3), border = NA)
abline(v = arco_bc[2], lty = 3, col = "grey30")
lines(g, eta, lwd = 2); abline(h = 0, col = "grey50")
points(c(th_ls, th_bc), eta_P1, pch = 19, col = c("black", morado), cex = 1.3)
text(c(th_ls, th_bc), eta_P1, c(sprintf("P1 con MCO\n%.1f°", th_ls), sprintf("P1 corregida\n%.1f°", th_bc)),
     pos = c(2, 4), cex = 0.8, col = c("black", morado))
text(mean(arco_ls), lim_a[2] - 0.08 * diff(lim_a), "Conjunto de la P9", cex = 0.85)
legend("bottomright", c("Conjunto de la P9 con A de MCO", "Su borde superior con A corregida"),
       fill = c(adjustcolor("grey60", 0.3), NA), border = NA, lty = c(NA, 3), col = c(NA, "grey30"), bty = "n", cex = 0.8)
plot(g, direccion(m$A, g * pi / 180), type = "l", lwd = 2, ylim = c(-1, 1), las = 1, xlab = "Ángulo (grados)",
     ylab = "Efecto de largo plazo sobre el ITCR (normalizado)", main = "(b) Dónde se anula el efecto de largo plazo",
     cex.main = 1, font.main = 1)
lines(g, direccion(x$A_BC, g * pi / 180), lwd = 2, lty = 2, col = morado)
abline(h = 0, col = "grey50"); abline(v = c(th_ls, th_bc), col = c("black", morado), lty = 3)
legend("topleft", c("Con A de MCO", "Con A corregida"), lwd = 2, lty = c(1, 2), col = c("black", morado), bty = "n", cex = 0.85)
dev.off()
cat(sprintf("Arco de la P9: MCO [%.1f, %.1f]; corregida [%.1f, %.1f]. Ángulo de la P1: %.1f -> %.1f\n",
            arco_ls[1], arco_ls[2], arco_bc[1], arco_bc[2], th_ls, th_bc))

# ---- 6. Efron contra Hall sobre la distribución bootstrap de la respuesta de impacto del ITCR -----------------
z <- x$impacto_ITCR$bloques; z <- z[!is.na(z)]; est0 <- x$puntual$r[1, 5]
q16 <- quantile(z, .16); q84 <- quantile(z, .84)
png("figuras/P1_Efron_Hall.png", width = 2200, height = 950, res = 220)
par(mar = c(3.6, 3.8, 1.5, 1), mgp = c(2.3, 0.6, 0))
d <- density(z[z > quantile(z, .005) & z < quantile(z, .995)])
plot(d, lwd = 2, col = naranja, main = "", las = 1, xlim = range(d$x, 2 * est0 - q84, 2 * est0 - q16, 0) + c(-1, 1.5),
     ylim = c(0, max(d$y) * 1.35),
     xlab = "ITCR al impacto (%, por 25 pb; + = depreciación)", ylab = "Densidad de las réplicas por bloques")
abline(v = est0, lwd = 2); abline(v = median(z), lty = 2, col = naranja); abline(v = 0, col = "grey60")
yE <- max(d$y) * 1.12; yH <- max(d$y) * 1.25
segments(q16, yE, q84, yE, lwd = 5, col = azul); segments(2 * est0 - q84, yH, 2 * est0 - q16, yH, lwd = 5, col = morado)
text(q16, yE, sprintf("Efron [%.2f, %.2f]", q16, q84), pos = 2, cex = 0.85, col = azul)
text(2 * est0 - q16, yH, sprintf("Hall [%.2f, %.2f]", 2 * est0 - q84, 2 * est0 - q16), pos = 4, cex = 0.85, col = morado)
legend("right", c(sprintf("Estimación: %.2f", est0), sprintf("Mediana de las réplicas: %.2f", median(z))),
       lwd = c(2, 1), lty = c(1, 2), col = c("black", naranja), bty = "n", cex = 0.85)
dev.off()

# ---- 7. Bandas mes a mes contra la banda conjunta de Bonferroni ----------------------------------------------
b <- x$bandas$bloques; hh <- 0:24
png("figuras/P1_Bonferroni.png", width = 2200, height = 950, res = 220)
par(mfrow = c(1, 2), mar = c(3.6, 3.8, 2.4, 1), mgp = c(2.3, 0.6, 0))
for (k in 2:3) {
  plot(hh, x$puntual$r[hh + 1, k], type = "n", ylim = range(b$b02[2:9, k], b$b98[2:9, k], b$b16[hh + 1, k], b$b84[hh + 1, k], 0),
       las = 1, xlab = "Meses", ylab = "pp, por 25 pb", main = c("", "Inflación anual", "Expectativas de inflación")[k],
       cex.main = 1, font.main = 1)
  polygon(c(1:8, 8:1), c(b$b02[2:9, k], rev(b$b98[2:9, k])), col = adjustcolor(morado, 0.18), border = NA)
  polygon(c(hh, rev(hh)), c(b$b16[hh + 1, k], rev(b$b84[hh + 1, k])), col = adjustcolor("grey40", 0.35), border = NA)
  lines(hh, x$puntual$r[hh + 1, k], lwd = 2); abline(h = 0, lty = 2, col = "grey30")
  if (k == 2) legend("topright", c("Estimación", "Banda de 68% mes a mes", "Banda conjunta de 68%, meses 1 a 8"),
                     lwd = c(2, 8, 8), col = c("black", adjustcolor("grey40", 0.35), adjustcolor(morado, 0.18)), bty = "n", cex = 0.8)
}
dev.off()
cat(sprintf("Bonferroni (bloques), meses 1-8: inflación percentil 2 mínimo %.3f, máximo %.3f; expectativas percentil 2 máximo %.3f\n",
            min(b$b02[2:9, 2]), max(b$b02[2:9, 2]), max(b$b02[2:9, 3])))
cat(sprintf("AR(1): promedio de MCO %.4f (Marriott y Pope %.4f); 1/(1-rho): %.1f, %.1f, %.1f\n",
            mean(est), aprox, 1 / (1 - aprox), 1 / (1 - rho), 1 / (1 - rho_BC)))
cat(sprintf("Simulación: autocorrelación de u y de u^2: iid %.2f y %.2f; ARCH %.2f y %.2f; curtosis ARCH %.1f\n",
            ac(u_iid, 1), ac(u_iid^2, 1), ac(u_arch, 1), ac(u_arch^2, 1), mean(u_arch^4) / mean(u_arch^2)^2))
cat(sprintf("Bloques, ITCR al impacto: Efron [%.2f, %.2f], Hall [%.2f, %.2f], mediana %.2f\n", q16, q84, 2 * est0 - q84, 2 * est0 - q16, median(z)))
