# ¿La corrección de sesgo de Kilian (1998) también cambia a la P9 y a la P10?
#
# El ángulo de la P1 sale de (I - A(1))^(-1), así que es sensible a la corrección (07_inferencia_P1.R). La P9 y la P10
# fijan su conjunto de ángulos con signos en los meses 0 a 5 y una cota de impacto. Aquí se repite su inferencia
# bayesiana (05_propuestas_P8_P9.R y 06_propuesta_P10.R) corrigiendo cada sorteo de la forma reducida con el sesgo
# estimado, como prueba de sensibilidad. La corrección de Kilian es una herramienta frecuentista: en la posterior se
# usa sólo para ver qué tanto dependen las respuestas de la persistencia.
#   P9: sesgo estimado en 07_inferencia_P1.R (bootstrap por bloques sobre la forma reducida de 02_var.R).
#   P10: sesgo estimado aquí, con 500 réplicas por bloques de 12 meses sobre la forma reducida en niveles.
# Se corre desde la carpeta Colombia, después de 07_inferencia_P1.R (tarda unos cinco minutos).
# setwd("~/Documents/Macroeconometría/Tareas/Ejercicio 2/Colombia")

set.seed(2026)
kappa <- 10; hmax <- 120; M <- 1000; n_ang <- 40

while (sink.number() > 0) sink()
sink("log/09_sesgo_P9_P10.log", split = TRUE)
cat("Corrección de sesgo en la P9 y la P10:", format(Sys.time()), "\n")

mA <- function(B, K, p) { A <- array(NA, c(K, K, p)); for (i in 1:p) A[, , i] <- t(B[1 + (i - 1) * K + 1:K, ]); A }
comp <- function(A) { K <- dim(A)[1]; p <- dim(A)[3]; C <- matrix(0, K * p, K * p)
  C[1:K, ] <- do.call(cbind, lapply(1:p, function(i) A[, , i]))
  if (p > 1) C[(K + 1):(K * p), 1:(K * (p - 1))] <- diag(K * (p - 1)); C }
raizmax <- function(A) max(Mod(eigen(comp(A), only.values = TRUE)$values))
corregir <- function(A, b) {            # ajuste de estacionariedad: el sesgo se encoge 1% a la vez (K&L p. 364)
  if (raizmax(A) >= 1) return(A); k <- 0
  repeat { Ac <- A - b * 0.99^k; if (raizmax(Ac) < 1) break; k <- k + 1 }
  attr(Ac, "fraccion") <- 0.99^k; Ac }
resp <- function(A, b, H = hmax) { K <- dim(A)[1]; p <- dim(A)[3]; r <- matrix(0, H + 1, K); r[1, ] <- b
  for (h in 1:H) for (j in 1:min(h, p)) r[h + 1, ] <- r[h + 1, ] + A[, , j] %*% r[h + 1 - j, ]; r }
bloques <- function(U, l = 12) { Tn <- nrow(U); s <- ceiling(Tn / l); ini <- sample(Tn - l + 1, s, replace = TRUE)
  idx <- as.vector(sapply(ini, function(a) a:(a + l - 1)))[1:Tn]
  media <- t(sapply(1:l, function(i) colMeans(U[i:(i + Tn - l), , drop = FALSE]))); U[idx, ] - media[rep(1:l, s)[1:Tn], ] }

# ---- Los dos modelos: datos, cómo leer niveles y qué variables se grafican ----
modelos <- list()
v <- readRDS("data/var_reducido.rds")
modelos$P9 <- list(Y = v$Y, X = v$X, p = v$p, filas = v$t0:nrow(v$Y),
  niveles = function(r) { r[, 1] <- cumsum(r[, 1]); r[, 5] <- cumsum(r[, 5]); r },        # ISE e ITCR acumulados
  series  = function(q) cbind(ISE = q[, 1], TC = q[, 5],                                   # inflación anual: suma de 12 meses
                              Infl = cumsum(q[, 2]) - c(rep(0, 12), head(cumsum(q[, 2]), -12))),
  sesgo   = readRDS("data/P1_inferencia.rds")$sesgo)
b <- read.csv("data/colombia_mensual.csv")
Y10 <- cbind(b$log_ise, 100 * log(b$ipc), b$exp_inflacion_12m, b$tpm, b$log_trm); mes <- as.integer(substr(b$fecha, 6, 7))
r1 <- function(x) c(NA, head(x, -1)); oil <- 100 * log(b$brent)
X10 <- cbind(sapply(2:12, function(j) as.numeric(mes == j)), sapply(sprintf("2020-%02d", 3:12), function(x) as.numeric(b$fecha == x)),
             oil, r1(oil), b$fed, r1(b$fed))
modelos$P10 <- list(Y = Y10, X = X10, p = 12, filas = which(b$fecha == "2006-01"):nrow(Y10),
  niveles = function(r) r,
  series  = function(q) cbind(ISE = q[, 1], TC = q[, 5] - q[, 2], Infl = q[, 2] - c(rep(0, 12), head(q[, 2], -12))),  # TC real
  sesgo   = NULL)

resultados <- list()
for (nm in names(modelos)) {
  mo <- modelos[[nm]]; Y <- mo$Y; X <- mo$X; p <- mo$p; f <- mo$filas; K <- ncol(Y)
  Z <- cbind(1, do.call(cbind, lapply(1:p, function(i) Y[f - i, ])), X[f, ]); n <- nrow(Z); m <- ncol(Z)
  ZtZi <- solve(crossprod(Z)); Bhat <- ZtZi %*% crossprod(Z, Y[f, ]); U <- Y[f, ] - Z %*% Bhat; S <- crossprod(U)
  A0 <- mA(Bhat, K, p)
  if (is.null(mo$sesgo)) {               # sesgo de la P10: bootstrap por bloques de 12 meses desde la A estimada
    Aest <- array(0, c(K, K, p)); R <- 500
    for (r in 1:R) { Ub <- bloques(U); Yb <- Y
      for (j in seq_along(f)) { t <- f[j]; y <- Bhat[1, ] + t(Bhat[(1 + K * p + 1):m, , drop = FALSE]) %*% X[t, ] + Ub[j, ]
        for (i in 1:p) y <- y + A0[, , i] %*% Yb[t - i, ]; Yb[t, ] <- y }
      Zb <- cbind(1, do.call(cbind, lapply(1:p, function(i) Yb[f - i, ])), X[f, ])
      Aest <- Aest + mA(solve(crossprod(Zb), crossprod(Zb, Yb[f, ])), K, p) / R }
    mo$sesgo <- Aest - A0 }
  ABC <- corregir(A0, mo$sesgo)
  frac <- if (is.null(attr(ABC, "fraccion"))) 1 else attr(ABC, "fraccion")
  cat(sprintf("\n========== %s ==========\nRaíz más grande: MCO %.3f -> corregida %.3f (se aplicó %.0f%% del sesgo)\n",
              nm, raizmax(A0), raizmax(ABC), 100 * frac))
  # Posterior con prior difuso, como en 05 y 06; cada sorteo se usa sin y con la corrección
  LZ <- t(chol(ZtZi)); out <- list(MCO = list(), corregida = list()); expl <- c(MCO = 0, corregida = 0)
  for (d in 1:M) {
    Sigma <- solve(rWishart(1, n - m, solve(S))[, , 1]); B <- Bhat + LZ %*% matrix(rnorm(m * K), m, K) %*% chol(Sigma)
    A <- mA(B, K, p); P <- t(chol(Sigma)); th <- runif(n_ang, -pi / 2, pi / 2)
    for (cual in c("MCO", "corregida")) {
      Au <- if (cual == "MCO") A else corregir(A, mo$sesgo)
      expl[cual] <- expl[cual] + (raizmax(Au) >= 1)
      r4 <- mo$niveles(resp(Au, P[, 4])); r5 <- mo$niveles(resp(Au, P[, 5]))
      for (t1 in th) { r <- cos(t1) * r4 + sin(t1) * r5
        if (all(r[1:6, 4] >= 0) && all(r[1:6, 5] <= 0) && r[1, 5] >= -kappa * r[1, 4])
          out[[cual]][[length(out[[cual]]) + 1]] <- mo$series(r * 0.25 / r[1, 4]) }
    }
  }
  for (cual in names(out)) {
    a <- simplify2array(out[[cual]])                           # (hmax + 1) x 3 x aceptados
    q <- lapply(c(p16 = .16, p50 = .5, p84 = .84), function(pr) apply(a, c(1, 2), quantile, pr))
    resultados[[nm]][[cual]] <- q
    cat(sprintf("%-9s aceptados %5d; VAR explosivo en %2.0f%% de los sorteos | mediana [p16, p84]:\n", cual, dim(a)[3], 100 * expl[cual] / M))
    for (h in c(0, 12, 24, 48, 120)) cat(sprintf("   mes %3d: ISE %+.2f [%+.2f, %+.2f] | TC %+.2f [%+.2f, %+.2f] | inflación %+.2f [%+.2f, %+.2f]\n", h,
      q$p50[h + 1, 1], q$p16[h + 1, 1], q$p84[h + 1, 1], q$p50[h + 1, 2], q$p16[h + 1, 2], q$p84[h + 1, 2],
      q$p50[h + 1, 3], q$p16[h + 1, 3], q$p84[h + 1, 3]))
  }
}

# ---- Figura: medianas y bandas de 68%, sin y con corrección, hasta 10 años ----
png("figuras/P9_P10_sesgo.png", width = 2200, height = 1800, res = 220)
par(mfrow = c(2, 2), mar = c(3.4, 3.8, 2.4, 1), mgp = c(2.2, 0.6, 0))
morado <- "#5b9bd5"; hh <- 0:hmax
titulos <- list(P9 = c("P9: ITCR (nivel, %)", "P9: ISE (nivel, %)"), P10 = c("P10: TC real (nivel, %)", "P10: ISE (nivel, %)"))
for (nm in c("P9", "P10")) for (j in 1:2) {
  k <- c(2, 1)[j]; q0 <- resultados[[nm]]$MCO; q1 <- resultados[[nm]]$corregida
  plot(hh, q0$p50[, k], type = "n", ylim = range(q0$p16[, k], q0$p84[, k], q1$p16[, k], q1$p84[, k], 0), las = 1,
       xlab = "Meses", ylab = "", main = titulos[[nm]][j], cex.main = 1, font.main = 1)
  polygon(c(hh, rev(hh)), c(q0$p16[, k], rev(q0$p84[, k])), col = adjustcolor("grey40", 0.35), border = NA)
  polygon(c(hh, rev(hh)), c(q1$p16[, k], rev(q1$p84[, k])), col = adjustcolor(morado, 0.22), border = NA)
  lines(hh, q0$p50[, k], lwd = 2); lines(hh, q1$p50[, k], lwd = 2, lty = 2, col = morado)
  abline(h = 0, lty = 2, col = "grey30")
  if (nm == "P9" && j == 1) legend("bottomleft", c("Sin corrección: mediana y banda 68%", "Con corrección: mediana y banda 68%"),
                                   lwd = 2, lty = c(1, 2), col = c("black", morado), bty = "n", cex = 0.85)
}
dev.off()
saveRDS(resultados, "data/P9_P10_sesgo.rds")
cat("\nListo: figuras/P9_P10_sesgo.png, data/P9_P10_sesgo.rds\n")
sink()
