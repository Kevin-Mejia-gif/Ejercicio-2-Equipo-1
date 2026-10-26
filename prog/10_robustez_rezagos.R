# Robustez de la Propuesta 1 al número de rezagos: 4, 6 (especificación principal) y 8
#
# Misma forma reducida de 02_var.R (inflación mensual con dummies estacionales, Brent y Fed exógenas, dummies de impulso),
# estimada con 4, 6 y 8 rezagos. Para cada una:
#   - la P1 de 04_estructural.R (tres ceros de impacto y neutralidad de largo plazo del ITCR);
#   - bandas de 68% del bootstrap por bloques de 12 meses, como en 07_inferencia_P1.R;
#   - la estimación con la corrección de sesgo de Kilian (1998). Sirve para ver cuánto pesa el sesgo de muestra chica,
#     que crece con el número de parámetros: si la corrección mueve mucho las respuestas, sus magnitudes no son confiables.
# Para comparar, también se calcula el conjunto de ángulos de la P9 (signos en los meses 0 a 5 y la cota de 10% por pp)
# en la estimación puntual.
#
# Se corre desde la carpeta Colombia, después de 02_var.R (tarda unos minutos).
# setwd("~/Documents/Macroeconometría/Tareas/Ejercicio 2/Colombia")

set.seed(2026)
hmax      <- 48     # horizonte de las respuestas (meses)
B         <- 2000   # réplicas del bootstrap por cada número de rezagos
l_bloque  <- 12     # largo de los bloques (meses)
choque_pp <- 0.25   # alza de 25 pb en la TPM
kappa     <- 10     # cota de la P9

while (sink.number() > 0) sink()
sink("log/10_robustez_rezagos.log", split = TRUE)
cat("Robustez de la P1 al número de rezagos:", format(Sys.time()), "\n\n")

v <- readRDS("data/var_reducido.rds")
Y <- v$Y
X <- v$X
K <- ncol(Y)
rezagos <- sort(c(v$p, v$p_robustez))
nombres <- c("ISE", "Inflación anual", "Expectativas", "TPM", "ITCR", "IPC")

anual <- function(x) { s <- cumsum(x); s - c(rep(0, 12), head(s, -12)) }   # suma de los últimos 12 meses

# Respuesta a un vector de impacto b: recursión r_h = A1 r_(h-1) + ... + Ap r_(h-p)
recursion <- function(A, b, H) {
  p <- dim(A)[3]
  r <- matrix(0, H + 1, K)
  r[1, ] <- b
  for (h in 1:H) for (j in 1:min(h, p)) r[h + 1, ] <- r[h + 1, ] + A[, , j] %*% r[h + 1 - j, ]
  r
}

# P1 por cada 25 pb. Columnas: ISE (nivel), inflación anual, expectativas, TPM, ITCR (nivel) e IPC (nivel)
irf_P1 <- function(A, S) {
  P  <- t(chol(S))
  L  <- solve(diag(K) - apply(A, c(1, 2), sum)) %*% P
  th <- atan(-L[5, 4] / L[5, 5])
  r  <- recursion(A, cos(th) * P[, 4] + sin(th) * P[, 5], hmax)
  r  <- r * choque_pp / r[1, 4]
  list(r = cbind(cumsum(r[, 1]), anual(r[, 2]), r[, 3], r[, 4], cumsum(r[, 5]), cumsum(r[, 2])), angulo = th * 180 / pi)
}

# Conjunto de ángulos de la P9 en una estimación: TPM >= 0 e ITCR en nivel <= 0 en los meses 0 a 5, y la cota
arco_P9 <- function(A, S) {
  P <- t(chol(S)); r4 <- recursion(A, P[, 4], 5); r5 <- recursion(A, P[, 5], 5)
  g <- seq(-89.9, 89.9, by = 0.1)
  ok <- sapply(g, function(t) { r <- cos(t * pi / 180) * r4 + sin(t * pi / 180) * r5
    all(r[, 4] >= 0) && all(cumsum(r[, 5]) <= 0) && r[1, 5] >= -kappa * r[1, 4] })
  if (any(ok)) range(g[ok]) else c(NA, NA)
}

# Bootstrap por bloques de residuos, con recentrado (igual que 07_inferencia_P1.R)
bloques <- function(U, l) {
  n <- nrow(U); s <- ceiling(n / l)
  ini <- sample(n - l + 1, s, replace = TRUE)
  idx <- as.vector(sapply(ini, function(a) a:(a + l - 1)))[1:n]
  media <- t(sapply(1:l, function(i) colMeans(U[i:(i + n - l), , drop = FALSE])))
  U[idx, ] - media[rep(1:l, s)[1:n], ]
}
simular <- function(mod, U) {
  p <- mod$p; Yb <- Y
  for (j in seq_along(mod$filas)) {
    t <- mod$filas[j]
    y <- mod$constante + mod$gamma %*% X[t, ] + U[j, ]
    for (i in 1:p) y <- y + mod$A[, , i] %*% Yb[t - i, ]
    Yb[t, ] <- y
  }
  Yb
}
# Corrección de sesgo con ajuste de estacionariedad (K&L p. 364)
corregir <- function(A, sesgo) {
  if (v$raices(A)[1] >= 1) return(A)
  k <- 0
  repeat { Ac <- A - sesgo * 0.99^k; if (v$raices(Ac)[1] < 1) break; k <- k + 1 }
  attr(Ac, "fraccion") <- 0.99^k
  Ac
}

res <- list()
for (p in rezagos) {
  m  <- v$estimar_var(Y, X, p, inicio = v$t0)
  pt <- irf_P1(m$A, m$sigma)
  irf <- array(NA, c(B, hmax + 1, 6)); ang <- rep(NA, B); A_media <- array(0, dim(m$A))
  for (b in 1:B) {
    mb <- v$estimar_var(simular(m, bloques(m$residuos, l_bloque)), X, p, inicio = v$t0)
    A_media <- A_media + mb$A / B
    if (v$raices(mb$A)[1] >= 1) next
    s <- irf_P1(mb$A, mb$sigma)
    irf[b, , ] <- s$r
    ang[b] <- s$angulo
  }
  A_BC <- corregir(m$A, A_media - m$A)
  bc   <- irf_P1(A_BC, m$sigma)
  ok   <- !is.na(ang)
  res[[paste0("p", p)]] <- list(p = p, puntual = pt, corregida = bc, angulo = ang, irf = irf[ok, , ],
                                raiz = v$raices(m$A)[1], raiz_BC = v$raices(A_BC)[1],
                                fraccion = if (is.null(attr(A_BC, "fraccion"))) 1 else attr(A_BC, "fraccion"),
                                arco_P9 = arco_P9(m$A, m$sigma), par_ecuacion = m$n_par)
}

# ---- Cuadro ---------------------------------------------------------------------------------------------------
celdas <- list(c(5, 0), c(2, 6), c(2, 12), c(6, 12), c(6, 24), c(1, 24), c(5, 12))
etiq   <- c("ITCR al impacto", "Inflación anual, mes 6", "Inflación anual, mes 12", "IPC, mes 12", "IPC, mes 24",
            "ISE, mes 24", "ITCR, mes 12")
cuadro <- NULL
for (x in res) {
  q16 <- apply(x$irf, 2:3, quantile, .16); q84 <- apply(x$irf, 2:3, quantile, .84)
  fila <- data.frame(rezagos = x$p, regresores_por_ecuacion = x$par_ecuacion,
                     raiz_MCO = round(x$raiz, 3), raiz_corregida = round(x$raiz_BC, 3),
                     angulo = round(x$puntual$angulo, 1), angulo_p16 = round(unname(quantile(x$angulo, .16, na.rm = TRUE)), 1),
                     angulo_p84 = round(unname(quantile(x$angulo, .84, na.rm = TRUE)), 1), angulo_corregido = round(x$corregida$angulo, 1),
                     arco_P9 = sprintf("[%.1f, %.1f]", x$arco_P9[1], x$arco_P9[2]),
                     replicas_validas = dim(x$irf)[1], deprecian_al_impacto = round(mean(x$irf[, 1, 5] > 0), 2))
  for (i in seq_along(celdas)) {
    k <- celdas[[i]][1]; h <- celdas[[i]][2] + 1
    fila[[etiq[i]]] <- sprintf("%+.2f [%+.2f, %+.2f] -> %+.2f", x$puntual$r[h, k], q16[h, k], q84[h, k], x$corregida$r[h, k])
  }
  cuadro <- rbind(cuadro, fila)
}
rownames(cuadro) <- paste(cuadro$rezagos, "rezagos")
cat("Respuesta a un alza de 25 pb: estimación [banda de 68%, bloques de 12 meses] -> estimación con corrección de sesgo\n\n")
print(t(cuadro), quote = FALSE)
cat(sprintf("\nEn la especificación principal (%d rezagos) las bandas coinciden con las de 07_inferencia_P1.R salvo por el",
            v$p), "error de simulación: son réplicas distintas.\n")

# Meses en que la banda de 68% de la inflación anual excluye el cero, para cada número de rezagos
meses <- function(x) if (length(x)) paste(x, collapse = " ") else "ninguno"
cat("\nInflación anual: meses en que la banda de 68% excluye el cero\n")
for (x in res) {
  q16 <- apply(x$irf, 2:3, quantile, .16); q84 <- apply(x$irf, 2:3, quantile, .84)
  cat(sprintf("%d rezagos: sobre cero, %s; bajo cero, %s\n", x$p, meses(which(q16[, 2] > 0) - 1),
              meses(which(q84[, 2] < 0) - 1)))
}


sink()

write.csv(cuadro, "data/P1_robustez_rezagos_25pb.csv", row.names = FALSE)
saveRDS(lapply(res, function(x) x[c("p", "puntual", "corregida", "raiz", "raiz_BC", "fraccion", "arco_P9")]),
        "data/P1_robustez_rezagos.rds")

# ---- Figura: la P1 con 4, 6 y 8 rezagos; banda de 68% de la especificación principal ----------------------------
titulos <- c("ISE (nivel, %)", "Inflación anual (pp)", "Expectativas de inflación (pp)", "TPM (pp)",
             "ITCR (nivel, %; + = depreciación)", "IPC (nivel, %)")
estilo <- list(p4 = list(col = "grey35", lty = 2), p6 = list(col = "#1f4e79", lty = 1), p8 = list(col = "#5b9bd5", lty = 4))
principal <- res[[paste0("p", v$p)]]
b16 <- apply(principal$irf, 2:3, quantile, .16); b84 <- apply(principal$irf, 2:3, quantile, .84)
png("figuras/irf_P1_rezagos_25pb.png", width = 2000, height = 2000, res = 220)
par(mfrow = c(3, 2), mar = c(3.2, 3.6, 2.4, 1), mgp = c(2, 0.6, 0), oma = c(2.6, 0, 0, 0))
for (k in 1:6) {
  rango <- range(b16[, k], b84[, k], sapply(res, function(x) range(x$puntual$r[, k])), 0)
  plot(0:hmax, rep(0, hmax + 1), type = "n", ylim = rango, xlab = "Meses", ylab = "", main = titulos[k],
       cex.main = 1, font.main = 1, las = 1)
  polygon(c(0:hmax, hmax:0), c(b16[, k], rev(b84[, k])), col = "grey80", border = NA)
  abline(h = 0, lty = 2, col = "grey40")
  for (x in res) { e <- estilo[[paste0("p", x$p)]]; lines(0:hmax, x$puntual$r[, k], lwd = 2, col = e$col, lty = e$lty) }
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("bottom", c(paste(rezagos, "rezagos"), sprintf("Banda de 68%% con %d rezagos", v$p)), horiz = TRUE, bty = "n", cex = 0.9,
       lwd = c(rep(2, length(rezagos)), 8), lty = c(sapply(rezagos, function(p) estilo[[paste0("p", p)]]$lty), 1),
       col = c(sapply(rezagos, function(p) estilo[[paste0("p", p)]]$col), "grey80"))
dev.off()
cat("Listo: figuras/irf_P1_rezagos_25pb.png, data/P1_robustez_rezagos_25pb.csv y data/P1_robustez_rezagos.rds\n")
