# Identificación estructural: comparación de las siete propuestas y respuestas al impulso

setwd("~/Documents/Macroeconometría/Tareas/Ejercicio 2/Colombia")

set.seed(2026)
hmax  <- 48      # horizonte de las respuestas (meses)
B     <- 1000    # réplicas del bootstrap
n_rot <- 20000   # rotaciones sorteadas para las propuestas con signos
choque_pp <- 0.25  # tamaño del choque: alza sorpresa de la TPM en el impacto, en pp (0.25 = 25 pb)
sufijo <- paste0("_", round(100 * choque_pp), "pb")

while (sink.number() > 0) sink()
sink("log/04_estructural.log", split = TRUE)
cat("Identificación estructural:", format(Sys.time()), "\n\n")

v <- readRDS("data/var_reducido.rds")
m <- v$modelo
p <- v$p
K <- ncol(v$Y)
nombres <- c("ISE", "Inflación", "Expectativas", "TPM", "ITCR")

# Matrices Psi_h de la representación de medias móviles
psis <- function(A) {
  comp <- matrix(0, K * p, K * p)
  comp[1:K, ] <- do.call(cbind, lapply(1:p, function(i) A[, , i]))
  if (p > 1) comp[(K + 1):(K * p), 1:(K * (p - 1))] <- diag(K * (p - 1))
  M <- diag(K * p)
  lapply(0:hmax, function(h) { r <- M[1:K, 1:K]; M <<- comp %*% M; r })
}

# La inflación entra al VAR como variación mensual; su respuesta se presenta como inflación anual, es decir, la suma
# de las respuestas de los últimos 12 meses (antes del choque la respuesta es cero), para leerla igual que la meta.
anual <- function(x) { s <- cumsum(x); s - c(rep(0, 12), head(s, -12)) }

# Respuestas al choque monetario (columna 4), normalizadas a un alza de choque_pp en la TPM.
# ISE e ITCR se acumulan para leerse en nivel; la inflación mensual se pasa a anual.
respuestas <- function(Psi, H) {
  r <- t(sapply(Psi, function(Ps) (Ps %*% H)[, 4])) / H[4, 4] * choque_pp
  r[, 1] <- cumsum(r[, 1])
  r[, 2] <- anual(r[, 2])
  r[, 5] <- cumsum(r[, 5])
  r
}

# ---- Propuestas 1-4: identificación exacta ---------------------------------
# Cada columna j de Q (H = P Q, P = Cholesky de Sigma) debe ser ortogonal a sus
# restricciones y a las columnas ya resueltas; se resuelve columna por columna.
sr <- function(...) lapply(c(...), function(i) list(tipo = "sr", fila = i))
lr <- function(...) lapply(c(...), function(i) list(tipo = "lr", fila = i))

restricciones <- list(
  P1 = list(list(), sr(1), sr(1, 2), c(sr(1, 2, 3), lr(5)), sr(1, 2, 3)),
  P2 = list(list(), sr(1), sr(1, 2), c(sr(2, 3), lr(1, 5)), sr(1, 2, 3)),
  P3 = list(list(), sr(1), sr(1, 2), c(sr(1, 2, 3), lr(5)), c(sr(1, 2), lr(1))),
  P4 = list(list(), sr(1), sr(1, 2), c(sr(2, 3), lr(1, 5)), c(sr(1, 2), lr(1)))
)

identificar_exacta <- function(A, S, res) {
  P <- t(chol(S))
  L <- solve(diag(K) - apply(A, c(1, 2), sum)) %*% P
  Q <- matrix(0, K, K)
  hechas <- integer(0)
  for (j in order(-sapply(res, length))) {
    R <- do.call(rbind, lapply(res[[j]], function(z) if (z$tipo == "sr") P[z$fila, ] else L[z$fila, ]))
    R <- rbind(R, t(Q[, hechas, drop = FALSE]))
    if (nrow(R) != K - 1 || qr(R)$rank < K - 1) stop("restricciones no identifican la columna ", j)
    Q[, j] <- svd(R, nv = K)$v[, K]
    hechas <- c(hechas, j)
  }
  H <- P %*% Q
  H %*% diag(sign(diag(H)))   # orientación: diagonal positiva (TPM sube con el choque 4)
}

# ---- Propuestas 5-7: ceros y signos (identificación por conjuntos) ----------
ortogonal <- function(R, d) qr.Q(qr(cbind(t(R), diag(d))), complete = TRUE)[, -(1:nrow(R)), drop = FALSE]
unitario <- function(x) c(x) / sqrt(sum(x^2))

sortear_H <- function(P, propuesta) {
  R <- diag(K)
  if (propuesta == "P5") {            # 9 ceros: rotación libre de las columnas 4 y 5
    th <- runif(1, 0, 2 * pi)
    R[4:5, 4:5] <- matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2)
  } else if (propuesta == "P6") {     # 8 ceros: columnas 3-5 rotan con h34 = 0
    q4 <- unitario(ortogonal(matrix(P[3, 3:5], 1), 3) %*% rnorm(2))
    Z <- qr.Q(qr(cbind(q4, matrix(rnorm(9), 3))), complete = TRUE)
    R[3:5, 3:5] <- cbind(Z[, 2], q4, Z[, 3])
  } else {                            # P7, 7 ceros: columnas 2-5 rotan con h23 = h24 = h34 = 0
    q4 <- unitario(ortogonal(rbind(P[2, 2:5], P[3, 2:5]), 4) %*% rnorm(2))
    q3 <- unitario(ortogonal(rbind(P[2, 2:5], q4), 4) %*% rnorm(2))
    W  <- ortogonal(rbind(q4, q3), 4) %*% qr.Q(qr(matrix(rnorm(4), 2)))
    R[2:5, 2:5] <- cbind(W[, 1], q3, q4, W[, 2])
  }
  H <- P %*% R
  if (H[4, 4] < 0) H[, 4] <- -H[, 4]
  H
}

cumple_signos <- function(H, Psi, propuesta) {
  r <- t(sapply(Psi[1:13], function(Ps) (Ps %*% H)[, 4]))
  r[, 2] <- anual(r[, 2])                                 # inflación anual, como en las respuestas
  ok <- H[5, 4] < 0                                       # apreciación real al impacto
  if (propuesta %in% c("P5", "P7")) ok <- ok && all(r[2:4, 4] > 0) && all(r[c(7, 10, 13), 2] < 0)
  if (propuesta == "P6")            ok <- ok && all(r[c(4, 7, 10), 3] < 0)
  ok
}

# ---- Estimación puntual -----------------------------------------------------
Psi <- psis(m$A)
H_exactas <- lapply(restricciones, function(res) identificar_exacta(m$A, m$sigma, res))
irf_puntual <- lapply(H_exactas, function(H) respuestas(Psi, H))

cat("Verificación de las propuestas exactas:\n")
L_hat <- solve(diag(K) - apply(m$A, c(1, 2), sum))
for (nm in names(H_exactas)) {
  H <- H_exactas[[nm]]
  cat(sprintf("  %s: max|HH' - Sigma| = %.1e | largo plazo ITCR-choque 4 = %.1e | largo plazo ISE-choque 4 = %.1e\n",
              nm, max(abs(H %*% t(H) - m$sigma)), (L_hat %*% H)[5, 4], (L_hat %*% H)[1, 4]))
}
cat(sprintf("\nP3 vs P1, máxima diferencia en el choque monetario: %.1e\n", max(abs(H_exactas$P3[, 4] - H_exactas$P1[, 4]))))
cat(sprintf("P4 vs P2, máxima diferencia en el choque monetario: %.1e\n", max(abs(H_exactas$P4[, 4] - H_exactas$P2[, 4]))))
cat("\nMatriz de impacto H, Propuesta 1:\n")
print(round(H_exactas$P1, 3))
cat(sprintf("\nTamaño del choque de una desviación estándar: %.3f pp de TPM\n", H_exactas$P1[4, 4]))

# ---- Bootstrap wild para P1 y P2 ---------------------------------------------
boot_irf <- function(res) {
  out <- array(NA, c(B, hmax + 1, K))
  for (b in 1:B) {
    eta <- sample(c(-1, 1), nrow(m$residuos), replace = TRUE)
    U <- m$residuos * eta
    Yb <- v$Y
    for (j in seq_along(m$filas)) {
      t <- m$filas[j]
      y <- m$constante + m$gamma %*% v$X[t, ] + U[j, ]
      for (i in 1:p) y <- y + m$A[, , i] %*% Yb[t - i, ]
      Yb[t, ] <- y
    }
    mb <- v$estimar_var(Yb, v$X, p, inicio = v$t0)
    if (v$raices(mb$A)[1] >= 1) next
    out[b, , ] <- respuestas(psis(mb$A), identificar_exacta(mb$A, mb$sigma, res))
  }
  out
}
bandas <- function(bt, pr) apply(bt[!is.na(bt[, 1, 1]), , , drop = FALSE], c(2, 3), quantile, pr, na.rm = TRUE)

boot <- list(P1 = boot_irf(restricciones$P1), P2 = boot_irf(restricciones$P2))
for (nm in names(boot)) cat(sprintf("\nBootstrap %s: %d réplicas válidas de %d", nm, sum(!is.na(boot[[nm]][, 1, 1])), B))
cat("\n")
bd <- lapply(boot, function(bt) list(b05 = bandas(bt, .05), b16 = bandas(bt, .16), b84 = bandas(bt, .84), b95 = bandas(bt, .95)))

# ---- Conjuntos identificados con signos -----------------------------------------
P_hat <- t(chol(m$sigma))
conjuntos <- list()
for (prop in c("P5", "P6", "P7")) {
  acept <- list()
  for (i in 1:n_rot) {
    H <- sortear_H(P_hat, prop)
    if (cumple_signos(H, Psi, prop)) acept[[length(acept) + 1]] <- respuestas(Psi, H)
  }
  if (length(acept) > 0) conjuntos[[prop]] <- simplify2array(acept)   # (hmax+1) x K x aceptadas
  cat(sprintf("%s: %d de %d rotaciones cumplen los signos (%.0f%%)\n", prop, length(acept), n_rot, 100 * length(acept) / n_rot))
  if (length(acept) == 0) cat(sprintf("   %s queda vacía: con esta forma reducida ninguna rotación cumple sus signos; se omite en los cuadros\n", prop))
}
cuantil_conj <- function(a, pr) apply(a, c(1, 2), quantile, pr)

# ---- Cuadro resumen de respuestas -----------------------------------------------
horizontes <- c(0, 12, 24)
marca <- function(x, lo, hi) paste0(sprintf("%.2f", x), ifelse(lo > 0 | hi < 0, "*", ""))
filas <- list()
for (nm in c("P1", "P2")) {
  ir <- irf_puntual[[nm]]
  fila <- c(Propuesta = nm)
  for (k in c(1, 2, 3, 5)) for (h in horizontes)
    fila[paste0(nombres[k], "_", h)] <- marca(ir[h + 1, k], bd[[nm]]$b05[h + 1, k], bd[[nm]]$b95[h + 1, k])
  filas[[nm]] <- fila
}
for (nm in names(conjuntos)) {
  med <- cuantil_conj(conjuntos[[nm]], .5)
  fila <- c(Propuesta = nm)
  for (k in c(1, 2, 3, 5)) for (h in horizontes) fila[paste0(nombres[k], "_", h)] <- sprintf("%.2f", med[h + 1, k])
  filas[[nm]] <- fila
}
resumen <- as.data.frame(do.call(rbind, filas))
cat("\nRespuestas a un alza de", 100 * choque_pp, "pb en la TPM (P1-P2: estimación puntual, * = banda al 90% excluye cero; P5-P7: mediana del conjunto)\n")
print(resumen, row.names = FALSE)

sink()

write.csv(resumen, paste0("data/identificacion_resumen", sufijo, ".csv"), row.names = FALSE)
write.csv(round(H_exactas$P1, 6), "data/H_propuesta1.csv")
saveRDS(list(H = H_exactas, irf = irf_puntual, bandas = bd, conjuntos = lapply(conjuntos, function(a) list(p16 = cuantil_conj(a, .16), p50 = cuantil_conj(a, .5), p84 = cuantil_conj(a, .84), aceptadas = dim(a)[3])),
             identificar_exacta = identificar_exacta, restricciones = restricciones,
             respuestas = respuestas, psis = psis), paste0("data/identificacion", sufijo, ".rds"))

# ---- Figuras ---------------------------------------------------------------------
titulos <- c("ISE (nivel, %)", "Inflación anual (pp)", "Expectativas de inflación (pp)", "TPM (pp)",
             "ITCR (nivel, %; + = depreciación)")
panel <- function(k, rangos) {
  plot(0:hmax, rep(0, hmax + 1), type = "n", ylim = range(rangos, 0), xlab = "Meses", ylab = "",
       main = titulos[k], cex.main = 1, font.main = 1, las = 1)
  abline(h = 0, lty = 2, col = "grey40")
}
sombra <- function(lo, hi, col) polygon(c(0:hmax, hmax:0), c(lo, rev(hi)), col = col, border = NA)

png(paste0("figuras/irf_propuesta1", sufijo, ".png"), width = 2000, height = 1900, res = 220)
par(mfrow = c(3, 2), mar = c(3.2, 3.6, 2.4, 1), mgp = c(2, 0.6, 0))
for (k in 1:K) {
  panel(k, c(bd$P1$b05[, k], bd$P1$b95[, k]))
  sombra(bd$P1$b05[, k], bd$P1$b95[, k], "grey88")
  sombra(bd$P1$b16[, k], bd$P1$b84[, k], "grey70")
  abline(h = 0, lty = 2, col = "grey40")
  lines(0:hmax, irf_puntual$P1[, k], lwd = 2)
}
plot.new()
legend("center", c("Respuesta estimada", "Banda 68%", "Banda 90%"), lwd = c(2, 8, 8),
       col = c("black", "grey70", "grey88"), bty = "n", cex = 1.1)
dev.off()

png(paste0("figuras/irf_p1_vs_p2", sufijo, ".png"), width = 2000, height = 1900, res = 220)
par(mfrow = c(3, 2), mar = c(3.2, 3.6, 2.4, 1), mgp = c(2, 0.6, 0))
for (k in 1:K) {
  panel(k, c(bd$P1$b16[, k], bd$P1$b84[, k], bd$P2$b16[, k], bd$P2$b84[, k]))
  sombra(bd$P1$b16[, k], bd$P1$b84[, k], adjustcolor("grey40", 0.30))
  sombra(bd$P2$b16[, k], bd$P2$b84[, k], adjustcolor("#9dc3e6", 0.30))
  abline(h = 0, lty = 2, col = "grey40")
  lines(0:hmax, irf_puntual$P1[, k], lwd = 2)
  lines(0:hmax, irf_puntual$P2[, k], lwd = 2, lty = 2, col = "#1f4e79")
}
plot.new()
legend("center", c("Propuesta 1 (= P3)", "Propuesta 2 (= P4)", "Banda 68% P1", "Banda 68% P2"),
       lwd = c(2, 2, 8, 8), lty = c(1, 2, 1, 1),
       col = c("black", "#1f4e79", adjustcolor("grey40", 0.30), adjustcolor("#9dc3e6", 0.30)), bty = "n", cex = 1.05)
dev.off()

png(paste0("figuras/irf_propuesta5", sufijo, ".png"), width = 2000, height = 1900, res = 220)
par(mfrow = c(3, 2), mar = c(3.2, 3.6, 2.4, 1), mgp = c(2, 0.6, 0))
c16 <- cuantil_conj(conjuntos$P5, .16); c50 <- cuantil_conj(conjuntos$P5, .5); c84 <- cuantil_conj(conjuntos$P5, .84)
for (k in 1:K) {
  panel(k, c(c16[, k], c84[, k], irf_puntual$P1[, k]))
  sombra(c16[, k], c84[, k], "grey80")
  abline(h = 0, lty = 2, col = "grey40")
  lines(0:hmax, c50[, k], lwd = 2)
  lines(0:hmax, irf_puntual$P1[, k], lwd = 1.5, lty = 2, col = "#1f4e79")
}
plot.new()
legend("center", c("Propuesta 5: mediana del conjunto", "Propuesta 5: percentiles 16-84", "Propuesta 1 (referencia)"),
       lwd = c(2, 8, 1.5), lty = c(1, 1, 2), col = c("black", "grey80", "#1f4e79"), bty = "n", cex = 1.05)
dev.off()
