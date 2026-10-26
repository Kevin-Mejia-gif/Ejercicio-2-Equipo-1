# VAR en forma reducida con variables exógenas

setwd("~/Documents/Macroeconometría/Tareas/Ejercicio 2/Colombia")

p_max <- 12
inicio_muestra <- "2010-01"   # primer mes estimado; los meses previos solo se usan como rezagos

while (sink.number() > 0) sink()
sink("log/02_var.log", split = TRUE)
cat("VAR reducido:", format(Sys.time()), "\n\n")

y_df <- read.csv("data/colombia_var.csv")
x_df <- read.csv("data/colombia_exogenas.csv")
datos <- merge(y_df, x_df, by = "fecha")
datos <- datos[complete.cases(datos), ]

# La inflación entra como variación mensual del IPC; como el IPC no está desestacionalizado, se agregan 11 dummies
# estacionales (est02 a est12) en todas las ecuaciones (Kilian y Lütkepohl, sección 19.2, pp. 663-668)
nombres_y <- c("d_log_ise", "inflacion_mensual", "exp_inflacion_12m", "tpm", "d_log_itcr")
nombres_x <- c("d_log_brent", "d_log_brent_l1", "fed", "fed_l1",
               grep("^est", names(x_df), value = TRUE),
               grep("^d20", names(x_df), value = TRUE))

Y <- as.matrix(datos[, nombres_y])
X <- as.matrix(datos[, nombres_x])
fechas <- datos$fecha
K <- ncol(Y)
t0 <- which(fechas == inicio_muestra)

# Estimación por MCO ecuación por ecuación.
# inicio: primera observación de la muestra efectiva (permite fijar la misma muestra para todos los p)
estimar_var <- function(Y, X, p, inicio = p + 1) {
  Tt <- nrow(Y)
  t_idx <- inicio:Tt
  rezagos <- do.call(cbind, lapply(1:p, function(i) Y[t_idx - i, , drop = FALSE]))
  Z <- cbind(1, rezagos, X[t_idx, , drop = FALSE])
  Yt <- Y[t_idx, , drop = FALSE]
  B <- solve(crossprod(Z), crossprod(Z, Yt))
  U <- Yt - Z %*% B
  n_obs <- nrow(Yt)
  n_par <- ncol(Z)
  K <- ncol(Y)
  A <- array(NA, c(K, K, p), dimnames = list(colnames(Y), colnames(Y), paste0("L", 1:p)))
  for (i in 1:p) A[, , i] <- t(B[1 + (i - 1) * K + 1:K, ])
  list(p = p, A = A,
       constante = B[1, ],
       gamma = t(B[(1 + K * p + 1):n_par, , drop = FALSE]),
       sigma = crossprod(U) / (n_obs - n_par),
       sigma_mv = crossprod(U) / n_obs,
       residuos = U, n_obs = n_obs, n_par = n_par, filas = t_idx)
}

# Matriz compañera y módulos de sus raíces
raices <- function(A) {
  K <- dim(A)[1]; p <- dim(A)[3]
  comp <- matrix(0, K * p, K * p)
  comp[1:K, ] <- do.call(cbind, lapply(1:p, function(i) A[, , i]))
  if (p > 1) comp[(K + 1):(K * p), 1:(K * (p - 1))] <- diag(K * (p - 1))
  sort(Mod(eigen(comp, only.values = TRUE)$values), decreasing = TRUE)
}

# Selección de rezagos con la misma muestra para todos los p
criterios <- t(sapply(1:p_max, function(p) {
  m <- estimar_var(Y, X, p, inicio = t0)
  ld <- log(det(m$sigma_mv))
  n <- m$n_obs
  c(p = p,
    AIC = ld + 2 * p * K^2 / n,
    HQ  = ld + 2 * log(log(n)) * p * K^2 / n,
    BIC = ld + log(n) * p * K^2 / n)
}))
cat("Criterios de información (muestra común desde", inicio_muestra, "):\n")
print(round(criterios, 4))
eleccion <- sapply(c("AIC", "HQ", "BIC"), function(k) criterios[which.min(criterios[, k]), "p"])
cat("\nRezagos elegidos:", paste(names(eleccion), eleccion, sep = " = ", collapse = ", "), "\n")

# Especificación principal: 6 rezagos. El AIC elige 5, que se toma como piso: en VAR mensuales el AIC da las respuestas
# con menor error (Kilian e Ivanov 2005), quedarse corto distorsiona más que pasarse y, en muestras chicas, hasta el AIC
# tiende a quedarse corto (Kilian y Lütkepohl, pp. 59-61); con restricciones de largo plazo el libro sugiere además subir
# el orden del VAR (sección 11.4.1). Se agrega un rezago y no más porque el sesgo de muestra chica crece con el número
# de parámetros. 4 y 8 rezagos quedan como robustez (10_robustez_rezagos.R).
p_elegido <- 6
p_robustez <- c(4, 8)

# Portmanteau ajustado (H0: sin autocorrelación hasta h)
portmanteau <- function(U, h, p) {
  n <- nrow(U)
  U <- scale(U, scale = FALSE)
  C0i <- solve(crossprod(U) / n)
  Q <- 0
  for (j in 1:h) {
    Cj <- crossprod(U[(j + 1):n, ], U[1:(n - j), ]) / n
    Q <- Q + sum(diag(t(Cj) %*% C0i %*% Cj %*% C0i)) / (n - j)
  }
  Q <- n^2 * Q
  gl <- ncol(U)^2 * (h - p)
  c(estadistico = Q, gl = gl, p_valor = 1 - pchisq(Q, gl))
}

# Jarque-Bera por ecuación (H0: normalidad)
jarque_bera <- function(u) {
  n <- length(u); u <- u - mean(u); s <- sqrt(mean(u^2))
  sk <- mean(u^3) / s^3; ku <- mean(u^4) / s^4
  jb <- n / 6 * (sk^2 + (ku - 3)^2 / 4)
  c(estadistico = jb, p_valor = 1 - pchisq(jb, 2))
}

# ARCH-LM por ecuación con q rezagos (H0: sin heterocedasticidad condicional)
# La constante va en la fórmula de lm: si va dentro de Z con "- 1", summary() reporta el R^2 sin centrar y la prueba
# rechaza casi siempre, aunque no haya ARCH.
arch_lm <- function(u, q = 12) {
  e2 <- u^2; n <- length(e2)
  Z <- sapply(1:q, function(i) e2[(q + 1 - i):(n - i)])
  y <- e2[(q + 1):n]
  r2 <- summary(lm(y ~ Z))$r.squared
  lm_est <- (n - q) * r2
  c(estadistico = lm_est, p_valor = 1 - pchisq(lm_est, q))
}

# Regresores de un modelo estimado (constante, rezagos y exógenas), para las pruebas que siguen
regresores <- function(m) cbind(1, do.call(cbind, lapply(1:m$p, function(i) Y[m$filas - i, , drop = FALSE])), X[m$filas, , drop = FALSE])

# LM de Breusch-Godfrey (H0: sin autocorrelación hasta el orden h), con la corrección F de Edgerton y Shukur
# (Lütkepohl 2005, sección 4.4.4): se regresan los residuos en los regresores del VAR y en h rezagos de los residuos
lm_autocorr <- function(m, h) {
  U <- m$residuos; n <- nrow(U); k <- ncol(U)
  Ul <- do.call(cbind, lapply(1:h, function(i) rbind(matrix(0, i, k), U[1:(n - i), , drop = FALSE])))
  W <- cbind(regresores(m), Ul)
  E <- U - W %*% solve(crossprod(W), crossprod(W, U))
  q <- k * h; r <- 0.5 * k * q - 1; N <- n - m$n_par - q - (k - q + 1) / 2
  s <- if (k^2 + q^2 - 5 > 0) sqrt((k^2 * q^2 - 4) / (k^2 + q^2 - 5)) else 1
  R <- (det(crossprod(U) / n) / det(crossprod(E) / n))^(1 / s) - 1
  Fs <- R * (N * s - r) / (k * q)
  c(estadistico = Fs, p_valor = 1 - pf(Fs, k * q, N * s - r))
}

# F por ecuación para las 11 dummies estacionales (H0: no hacen falta en esa ecuación)
f_estacionales <- function(m) {
  cols <- 1 + K * m$p + grep("^est", nombres_x)
  ZtZi <- solve(crossprod(regresores(m)))
  t(sapply(1:K, function(k) {
    b <- m$gamma[k, grep("^est", nombres_x)]
    Fs <- drop(t(b) %*% solve(m$sigma[k, k] * ZtZi[cols, cols], b)) / length(b)
    c(F = Fs, p_valor = 1 - pf(Fs, length(b), m$n_obs - m$n_par))
  }))
}

diagnosticos <- function(m, nombre) {
  U <- m$residuos
  r <- raices(m$A)
  pm <- rbind(portmanteau(U, 12, m$p), portmanteau(U, 24, m$p))
  rownames(pm) <- c("h = 12", "h = 24")
  # Las dummies dejan residuos exactamente en cero en sus meses; se excluyen de normalidad y ARCH
  activos <- rowSums(abs(U)) > 1e-10
  lb <- t(apply(U, 2, function(u) {
    b <- Box.test(u, lag = 12, type = "Ljung-Box")
    c(estadistico = unname(b$statistic), p_valor = b$p.value)
  }))
  jb <- t(apply(U[activos, ], 2, jarque_bera))
  ar <- t(apply(U[activos, ], 2, arch_lm))
  lmac <- rbind(lm_autocorr(m, 1), lm_autocorr(m, 12))
  rownames(lmac) <- c("orden 1", "orden 12")
  fe <- f_estacionales(m)
  cat("\n==========", nombre, "==========\n")
  cat("Observaciones efectivas:", m$n_obs, "| parámetros por ecuación:", m$n_par, "\n")
  cat("Muestra:", fechas[min(m$filas)], "a", fechas[max(m$filas)], "\n")
  cat("\nMódulo de la raíz más grande:", round(r[1], 4),
      if (r[1] < 1) "(estable)" else "(NO estable)", "\n")
  cat("\nPortmanteau ajustado:\n"); print(round(pm, 3))
  cat("\nLM de autocorrelación (Edgerton-Shukur):\n"); print(round(lmac, 3))
  cat("\nLjung-Box (12 rezagos) por ecuación:\n"); print(round(lb, 3))
  cat("\nJarque-Bera por ecuación:\n"); print(round(jb, 3))
  cat("\nARCH-LM (12 rezagos) por ecuación:\n"); print(round(ar, 3))
  cat("\nDummies estacionales, F por ecuación:\n"); print(round(fe, 3))
  invisible(list(raices = r, portmanteau = pm, lm = lmac, ljung_box = lb, jarque_bera = jb, arch = ar, estacionales = fe))
}

# Pruebas de razón de verosimilitud secuenciales, con la corrección de Sims (T - parámetros por ecuación)
cat("\nRazón de verosimilitud, p contra p + 1 (misma muestra):\n")
for (p in 4:7) {
  m0 <- estimar_var(Y, X, p, inicio = t0); m1 <- estimar_var(Y, X, p + 1, inicio = t0)
  LR <- (m1$n_obs - m1$n_par) * (log(det(m0$sigma_mv)) - log(det(m1$sigma_mv)))
  cat(sprintf("  VAR(%d) contra VAR(%d): LR = %.1f, gl = %d, p = %.3f\n", p, p + 1, LR, K^2, 1 - pchisq(LR, K^2)))
}

# Estimación con cada p sugerido por los criterios, con la principal y con las de robustez, con la muestra completa
modelos <- list()
diag_modelos <- list()
for (p in sort(unique(c(eleccion, p_elegido, p_robustez)))) {
  modelos[[paste0("p", p)]] <- estimar_var(Y, X, p, inicio = t0)
  diag_modelos[[paste0("p", p)]] <- diagnosticos(modelos[[paste0("p", p)]], paste("VAR con", p, "rezago(s)"))
}
var_final <- modelos[[paste0("p", p_elegido)]]
cat("\nEspecificación principal: VAR con", p_elegido, "rezagos\n")

sink()

saveRDS(list(Y = Y, X = X, fechas = fechas, nombres_y = nombres_y, nombres_x = nombres_x,
             inicio_muestra = inicio_muestra, t0 = t0, p = p_elegido, p_robustez = p_robustez, modelo = var_final, criterios = criterios,
             diagnosticos = diag_modelos, modelos_alternativos = modelos,
             estimar_var = estimar_var, raices = raices),
        "data/var_reducido.rds")

write.csv(round(criterios, 4), "data/var_criterios_rezagos.csv", row.names = FALSE)
write.csv(var_final$sigma, "data/var_sigma.csv")
write.csv(data.frame(fecha = fechas[var_final$filas], var_final$residuos),
          "data/var_residuos.csv", row.names = FALSE)
coefs <- rbind(constante = var_final$constante,
               do.call(rbind, lapply(1:p_elegido, function(i) {
                 a <- t(var_final$A[, , i]); rownames(a) <- paste0(rownames(a), "_L", i); a })),
               t(var_final$gamma))
write.csv(coefs, "data/var_coeficientes.csv")
