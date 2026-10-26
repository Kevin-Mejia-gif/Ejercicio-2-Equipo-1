# Comparación de especificaciones del VAR: variables, dummies, rezagos y forma de medir la inflación

setwd("~/Documents/Macroeconometría/Tareas/Ejercicio 2/Colombia")

p_max <- 12
inicio_muestra <- "2010-01"

while (sink.number() > 0) sink()
sink("log/03_comparacion.log", split = TRUE)
cat("Comparación de especificaciones:", format(Sys.time()), "\n\n")

base <- read.csv("data/colombia_mensual.csv")
x_df <- read.csv("data/colombia_exogenas.csv")
datos <- merge(base, x_df[, c("fecha", setdiff(names(x_df), names(base)))], by = "fecha")
datos <- datos[order(datos$fecha), ]
datos$d_credito_real_anual <- c(NA, diff(datos$credito_real_anual))

# Dummies de impulso para cualquier mes
dummy <- function(meses) {
  D <- sapply(meses, function(m) as.integer(datos$fecha == m))
  colnames(D) <- paste0("d", gsub("-", "", meses))
  D
}
meses <- function(desde, hasta) format(seq(as.Date(paste0(desde, "-01")), as.Date(paste0(hasta, "-01")), by = "month"), "%Y-%m")

conjuntos_y <- list(
  base           = c("d_log_ise", "inflacion_mensual", "tpm", "d_log_itcr"),
  expectativas   = c("d_log_ise", "inflacion_mensual", "exp_inflacion_12m", "tpm", "d_log_itcr"),
  credito        = c("d_log_ise", "inflacion_mensual", "credito_real_anual", "tpm", "d_log_itcr"),
  credito_dif    = c("d_log_ise", "inflacion_mensual", "d_credito_real_anual", "tpm", "d_log_itcr"),
  exp_credito    = c("d_log_ise", "inflacion_mensual", "exp_inflacion_12m", "credito_real_anual", "tpm", "d_log_itcr"),
  exp_credito_dif = c("d_log_ise", "inflacion_mensual", "exp_inflacion_12m", "d_credito_real_anual", "tpm", "d_log_itcr")
)

conjuntos_d <- list(
  "mar-jun 2020"         = meses("2020-03", "2020-06"),
  "mar-jun 2020 + may 2021" = c(meses("2020-03", "2020-06"), "2021-05"),
  "mar-dic 2020"         = meses("2020-03", "2020-12"),
  "mar-dic 2020 + may 2021" = c(meses("2020-03", "2020-12"), "2021-05")
)

exogenas <- c("d_log_brent", "d_log_brent_l1", "fed", "fed_l1")
estacionales <- grep("^est", names(x_df), value = TRUE)   # dummies estacionales, porque la inflación mensual no está desestacionalizada

estimar_var <- function(Y, X, p, inicio = p + 1) {
  Tt <- nrow(Y)
  t_idx <- inicio:Tt
  rezagos <- do.call(cbind, lapply(1:p, function(i) Y[t_idx - i, , drop = FALSE]))
  Z <- cbind(1, rezagos, X[t_idx, , drop = FALSE])
  Yt <- Y[t_idx, , drop = FALSE]
  B <- solve(crossprod(Z), crossprod(Z, Yt))
  U <- Yt - Z %*% B
  K <- ncol(Y)
  A <- array(NA, c(K, K, p))
  for (i in 1:p) A[, , i] <- t(B[1 + (i - 1) * K + 1:K, ])
  list(p = p, A = A, sigma = crossprod(U) / (nrow(Yt) - ncol(Z)), sigma_mv = crossprod(U) / nrow(Yt),
       residuos = U, n_obs = nrow(Yt), n_par = ncol(Z))
}

raices <- function(A) {
  K <- dim(A)[1]; p <- dim(A)[3]
  comp <- matrix(0, K * p, K * p)
  comp[1:K, ] <- do.call(cbind, lapply(1:p, function(i) A[, , i]))
  if (p > 1) comp[(K + 1):(K * p), 1:(K * (p - 1))] <- diag(K * (p - 1))
  max(Mod(eigen(comp, only.values = TRUE)$values))
}

portmanteau <- function(U, h, p) {
  n <- nrow(U); U <- scale(U, scale = FALSE)
  C0i <- solve(crossprod(U) / n)
  Q <- 0
  for (j in 1:h) {
    Cj <- crossprod(U[(j + 1):n, ], U[1:(n - j), ]) / n
    Q <- Q + sum(diag(t(Cj) %*% C0i %*% Cj %*% C0i)) / (n - j)
  }
  1 - pchisq(n^2 * Q, ncol(U)^2 * (h - p))
}

jb_p <- function(u) {
  n <- length(u); u <- u - mean(u); s <- sqrt(mean(u^2))
  jb <- n / 6 * ((mean(u^3) / s^3)^2 + (mean(u^4) / s^4 - 3)^2 / 4)
  1 - pchisq(jb, 2)
}

arch_p <- function(u, q = 12) {   # constante en la fórmula, para que summary() dé el R^2 centrado
  e2 <- u^2; n <- length(e2)
  Z <- sapply(1:q, function(i) e2[(q + 1 - i):(n - i)])
  r2 <- summary(lm(e2[(q + 1):n] ~ Z))$r.squared
  1 - pchisq((n - q) * r2, q)
}

resultados <- list()
for (ny in names(conjuntos_y)) for (nd in names(conjuntos_d)) {
  vars <- conjuntos_y[[ny]]
  X_all <- cbind(as.matrix(datos[, c(exogenas, estacionales)]), dummy(conjuntos_d[[nd]]))
  ok <- complete.cases(datos[, vars], X_all)
  Y <- as.matrix(datos[ok, vars]); X <- X_all[ok, , drop = FALSE]
  t0 <- which(datos$fecha[ok] == inicio_muestra)
  K <- length(vars)

  crit <- t(sapply(1:p_max, function(p) {
    m <- estimar_var(Y, X, p, inicio = t0)
    ld <- log(det(m$sigma_mv)); n <- m$n_obs
    c(AIC = ld + 2 * p * K^2 / n, HQ = ld + 2 * log(log(n)) * p * K^2 / n, BIC = ld + log(n) * p * K^2 / n)
  }))
  eleccion <- apply(crit, 2, which.min)

  for (p in sort(unique(c(eleccion, 6)))) {               # los que eligen los criterios y la especificación principal (6)
    m <- estimar_var(Y, X, p, inicio = t0)
    U <- m$residuos
    activos <- rowSums(abs(U)) > 1e-10
    lb <- apply(U, 2, function(u) Box.test(u, lag = 12, type = "Ljung-Box")$p.value)
    jb <- apply(U[activos, ], 2, jb_p)
    ar <- apply(U[activos, ], 2, arch_p)
    lr <- tryCatch(kappa(diag(K) - apply(m$A, c(1, 2), sum), exact = TRUE), error = function(e) Inf)
    resultados[[length(resultados) + 1]] <- data.frame(
      variables = ny, dummies = nd, K = K, p = p,
      elegido_por = if (any(eleccion == p)) paste(names(eleccion)[eleccion == p], collapse = "/") else "fijo",
      n_obs = m$n_obs, par_ecuacion = m$n_par,
      raiz_max = round(raices(m$A), 4),
      cond_largo_plazo = round(lr, 1),
      portmanteau_h12 = round(portmanteau(U, 12, p), 3),
      portmanteau_h24 = round(portmanteau(U, 24, p), 3),
      lb_rechazos = sum(lb < 0.05), lb_min = round(min(lb), 3),
      lb_inflacion = round(lb[vars == "inflacion_mensual"], 3),
      jb_rechazos = sum(jb < 0.05), arch_rechazos = sum(ar < 0.05),
      sd_res_ise = round(sqrt(m$sigma[1, 1]), 3),
      cor_tpm_itcr = round(cov2cor(m$sigma)[which(vars == "tpm"), which(vars == "d_log_itcr")], 3))
  }
}
tabla <- do.call(rbind, resultados)
rownames(tabla) <- NULL

# Una especificación es candidata si es estable con margen y no deja autocorrelación conjunta a 24 meses
tabla$candidata <- tabla$raiz_max < 0.98 & tabla$portmanteau_h24 > 0.05

cat("Muestra de estimación desde", inicio_muestra, "\n")
cat("cond_largo_plazo: número de condición de (I - A1 - ... - Ap); valores altos indican que C(1) es imprecisa\n\n")
print(tabla, row.names = FALSE)
cat("\nCandidatas (raíz < 0.98 y Portmanteau h24 > 0.05):\n")
print(tabla[tabla$candidata, c("variables", "dummies", "p", "elegido_por", "raiz_max", "cond_largo_plazo",
                               "portmanteau_h24", "lb_rechazos", "lb_inflacion", "par_ecuacion")], row.names = FALSE)

# ---- Inflación anual (diferencia de 12 meses) frente a mensual con dummies estacionales ---------------------------
# Mismas cinco variables y dummies de la especificación principal; sólo cambia cómo entra la inflación. La inflación
# anual es una suma móvil de 12 inflaciones mensuales: es muy persistente por construcción y un choque sale de ella
# exactamente 12 meses después. Las pruebas de raíz unitaria estacional de 01_datos.R indican estacionalidad
# determinística en el IPC, que se trata con dummies (Kilian y Lütkepohl, sección 19.2).
d_principal <- dummy(conjuntos_d[["mar-dic 2020 + may 2021"]])
formas <- list(
  anual   = list(vars = c("d_log_ise", "inflacion_anual", "exp_inflacion_12m", "tpm", "d_log_itcr"),
                 X = cbind(as.matrix(datos[, exogenas]), d_principal)),
  mensual = list(vars = c("d_log_ise", "inflacion_mensual", "exp_inflacion_12m", "tpm", "d_log_itcr"),
                 X = cbind(as.matrix(datos[, c(exogenas, estacionales)]), d_principal)))
comp_inflacion <- NULL
for (nm in names(formas)) {
  vars <- formas[[nm]]$vars
  ok <- complete.cases(datos[, vars], formas[[nm]]$X)
  Y <- as.matrix(datos[ok, vars]); X <- formas[[nm]]$X[ok, , drop = FALSE]
  t0 <- which(datos$fecha[ok] == inicio_muestra); K <- length(vars)
  for (p in c(4, 5, 6, 8, 10, 12)) {
    m <- estimar_var(Y, X, p, inicio = t0)
    U <- m$residuos
    activos <- rowSums(abs(U)) > 1e-10
    comp_inflacion <- rbind(comp_inflacion, data.frame(
      inflacion = nm, p = p, raiz_max = round(raices(m$A), 3),
      cond_largo_plazo = round(tryCatch(kappa(diag(K) - apply(m$A, c(1, 2), sum), exact = TRUE), error = function(e) Inf), 1),
      portmanteau_h24 = round(portmanteau(U, 24, p), 3),
      lb12_inflacion = round(Box.test(U[activos, 2], lag = 12, type = "Ljung-Box")$p.value, 3)))
  }
}
cat("\nInflación anual (sin dummies estacionales) frente a mensual (con dummies estacionales), resto igual a la especificación principal\n")
cat("raiz_max >= 1: VAR explosivo; lb12_inflacion: Ljung-Box a 12 rezagos de los residuos de la ecuación de inflación\n")
print(comp_inflacion, row.names = FALSE)

sink()

write.csv(tabla, "data/comparacion_especificaciones.csv", row.names = FALSE)
write.csv(comp_inflacion, "data/comparacion_inflacion.csv", row.names = FALSE)
