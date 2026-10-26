# Inferencia para la Propuesta 1: bootstrap por bloques de residuos y corrección de sesgo
#
# Forma reducida: la de 02_var.R. Identificación: la P1 de 04_estructural.R (tres ceros de impacto y neutralidad de
# largo plazo del ITCR). Este programa compara formas de ponerle bandas a la P1 (Kilian y Lütkepohl 2017, cap. 12):
#   1. Wild, el de 04_estructural.R: cada mes de residuos se multiplica por +1 o -1 al azar (sección 12.2.3).
#   2. Por bloques de residuos (Brüggemann, Jentsch y Trenkler 2016; sección 12.2.5): se remuestrean bloques de
#      12 meses seguidos de residuos. Es el que el libro considera válido para respuestas estructurales cuando hay
#      heterocedasticidad condicional (ARCH), y 02_var.R la detecta en las ecuaciones del ISE y de la TPM.
#   3. Por bloques con la corrección de sesgo de Kilian (1998) (sección 12.3), la que el libro recomienda con datos
#      persistentes.
# Bandas de percentiles de Efron: 16-84 (68%) y 5-95 (90%). Para comparar se reportan también bloques de 24 meses,
# el bootstrap iid y el intervalo de Hall.
#
# Se corre desde la carpeta Colombia, después de 04_estructural.R (tarda un par de minutos).
# Si usan RStudio, quiten el # y pongan su ruta:
# setwd("~/Documents/Macroeconometría/Tareas/Ejercicio 2/Colombia")

set.seed(2026)
hmax      <- 48     # horizonte de las respuestas (meses)
B_wild    <- 1000   # réplicas del wild: las mismas de 04_estructural.R, para reproducir sus bandas
B         <- 2000   # réplicas de los demás bootstraps
l_bloque  <- 12     # largo de los bloques (meses)
l_alterno <- 24     # otro largo, para ver la sensibilidad
choque_pp <- 0.25   # alza de 25 pb en la TPM

while (sink.number() > 0) sink()
sink("log/07_inferencia_P1.log", split = TRUE)
cat("Inferencia para la Propuesta 1:", format(Sys.time()), "\n\n")

v  <- readRDS("data/var_reducido.rds")
id <- readRDS("data/identificacion_25pb.rds")
m  <- v$modelo
p  <- v$p
K  <- ncol(v$Y)
nombres <- c("ISE", "Inflación", "Expectativas", "TPM", "ITCR")
U0 <- m$residuos                  # residuos de la forma reducida (199 meses x 5)
Tn <- nrow(U0)
fechas <- v$fechas[m$filas]
anio   <- as.numeric(substr(fechas, 1, 4)) + (as.numeric(substr(fechas, 6, 7)) - 1) / 12

# ---- 1. La P1 como función de la forma reducida ------------------------------------------------------------
# Con los tres ceros, la columna 4 de H es cos(th) P[, 4] + sin(th) P[, 5], con P = Cholesky de Sigma.
# La neutralidad de largo plazo del ITCR pide L[5, ] (cos(th), sin(th))' = 0, con L = (I - A1 - ... - Ap)^(-1) P,
# de modo que tan(th) = -L[5, 4] / L[5, 5]. Es la misma P1 de 04_estructural.R (se comprueba abajo).
anual <- function(x) { s <- cumsum(x); s - c(rep(0, 12), head(s, -12)) }
irf_P1 <- function(A, S) {
  P  <- t(chol(S))
  L  <- solve(diag(K) - apply(A, c(1, 2), sum)) %*% P
  th <- atan(-L[5, 4] / L[5, 5])
  r  <- matrix(0, hmax + 1, K, dimnames = list(NULL, nombres))
  r[1, ] <- cos(th) * P[, 4] + sin(th) * P[, 5]
  for (h in 1:hmax) for (j in 1:min(h, p)) r[h + 1, ] <- r[h + 1, ] + A[, , j] %*% r[h + 1 - j, ]
  r <- r * choque_pp / r[1, 4]            # alza de 25 pb en la TPM
  r[, 1] <- cumsum(r[, 1])                # ISE e ITCR en nivel
  r[, 2] <- anual(r[, 2])                 # inflación mensual pasada a anual (suma de 12 meses), como en 04
  r[, 5] <- cumsum(r[, 5])
  list(r = r, angulo = th * 180 / pi)
}
puntual <- irf_P1(m$A, m$sigma)
cat(sprintf("Comprobación con 04_estructural.R: diferencia máxima %.1e | ángulo de la P1 %.1f grados | raíz más grande %.3f\n\n",
            max(abs(puntual$r - id$irf$P1)), puntual$angulo, v$raices(m$A)[1]))

# ---- 2. Heterocedasticidad condicional (ARCH) en los residuos ---------------------------------------------------
# ARCH-LM (Engle 1982), igual que en 02_var.R: se regresa u_t^2 en 12 rezagos suyos; si la varianza es constante,
# (n - q) R^2 se distribuye chi-cuadrada con q grados de libertad. Los meses con dummy (residuo cero) se excluyen.
activos <- rowSums(abs(U0)) > 1e-10
arch_lm <- function(u, q = 12) {   # constante en la fórmula, para que summary() dé el R^2 centrado
  e2 <- u^2
  n  <- length(e2)
  Z  <- sapply(1:q, function(i) e2[(q + 1 - i):(n - i)])
  r2 <- summary(lm(e2[(q + 1):n] ~ Z))$r.squared
  c(estadistico = (n - q) * r2, p_valor = 1 - pchisq((n - q) * r2, q))
}
arch <- t(apply(U0[activos, ], 2, arch_lm))
ac1  <- function(x) cor(x[-1], x[-length(x)])
cuadro_arch <- data.frame(Ecuacion = nombres, ARCH_LM = round(arch[, 1], 1), p_valor = signif(arch[, 2], 2),
                          autocorr_u = round(apply(U0[activos, ], 2, ac1), 2),
                          autocorr_u2 = round(apply(U0[activos, ]^2, 2, ac1), 2))
cat("Cuadro 1. ARCH-LM (12 rezagos) y autocorrelación de orden 1 de los residuos y de sus cuadrados\n")
print(cuadro_arch, row.names = FALSE)

# Incertidumbre del estimador de la varianza del residuo del ITCR (sigma_55). Si los u_t^2 fueran independientes,
# su desviación estándar saldría de la varianza de u^2; con ARCH los cuadrados están autocorrelacionados y hace falta
# la varianza de largo plazo (Newey-West con 12 rezagos).
u2 <- U0[activos, 5]^2
na <- length(u2)
g  <- sapply(0:12, function(j) sum((u2[(j + 1):na] - mean(u2)) * (u2[1:(na - j)] - mean(u2))) / na)
de_indep <- sqrt(na * g[1]) / (Tn - m$n_par)
de_largo <- sqrt(na * (g[1] + 2 * sum((1 - (1:12) / 13) * g[-1]))) / (Tn - m$n_par)
cat(sprintf("\nsigma_55 estimada %.3f | desviación estándar con cuadrados independientes %.3f | con varianza de largo plazo %.3f\n\n",
            m$sigma[5, 5], de_indep, de_largo))

# ---- 3. Bootstrap -----------------------------------------------------------------------------------------------
# Diseño recursivo (como en 04_estructural.R): se arranca con los datos observados antes de la muestra, las
# exógenas quedan fijas y en cada mes se suma una innovación remuestreada.
simular <- function(A, U, mod = m, Y = v$Y, X = v$X) {
  Yb <- Y
  for (j in seq_along(mod$filas)) {
    t <- mod$filas[j]
    y <- mod$constante + mod$gamma %*% X[t, ] + U[j, ]
    for (i in 1:p) y <- y + A[, , i] %*% Yb[t - i, ]
    Yb[t, ] <- y
  }
  Yb
}
# Innovaciones de cada método
wild    <- function() U0 * sample(c(-1, 1), Tn, replace = TRUE)     # cada mes por +1 o -1
iid     <- function() U0[sample(Tn, Tn, replace = TRUE), ]           # meses sueltos, con reemplazo
bloques <- function(l) {                                             # bloques de l meses seguidos
  s   <- ceiling(Tn / l)                                             # bloques necesarios para llegar a Tn meses
  ini <- sample(Tn - l + 1, s, replace = TRUE)                       # dónde empieza cada bloque
  idx <- as.vector(sapply(ini, function(a) a:(a + l - 1)))[1:Tn]
  # Recentrado (K&L p. 356): la posición i de un bloque al azar tiene media igual al promedio de u_i, ..., u_(i+Tn-l),
  # que no es exactamente cero; se resta para que las innovaciones tengan media cero.
  media <- t(sapply(1:l, function(i) colMeans(U0[i:(i + Tn - l), , drop = FALSE])))
  U0[idx, ] - media[rep(1:l, s)[1:Tn], ]
}
# Corrección de sesgo con ajuste de estacionariedad (K&L p. 364): si al restar el sesgo la raíz más grande llega a 1,
# el sesgo se encoge 1% a la vez hasta que el VAR vuelva a ser estable.
corregir <- function(A, sesgo) {
  if (v$raices(A)[1] >= 1) return(A)
  k <- 0
  repeat {
    Ac <- A - sesgo * 0.99^k
    if (v$raices(Ac)[1] < 1) break
    k <- k + 1
  }
  attr(Ac, "fraccion") <- 0.99^k
  Ac
}
# Tamaño del multiplicador de largo plazo: valor singular mayor de (I - A(1))^(-1) = 1 / valor singular menor de I - A(1).
# Crece sin límite cuando I - A(1) se acerca a ser singular, es decir, cuando el VAR se acerca a una raíz unitaria.
multiplicador <- function(A) 1 / min(svd(diag(K) - apply(A, c(1, 2), sum))$d)
persistencia  <- function(A) diag(apply(A, c(1, 2), sum))[2:4]   # suma de coeficientes propios: inflación, expectativas, TPM

correr <- function(innov, R, A_dgp = m$A, sesgo = NULL, mod = m, Y = v$Y, X = v$X) {
  irf  <- array(NA, c(R, hmax + 1, K))
  Aest <- array(NA, c(R, K, K, p))
  angulo <- sigma55 <- mult <- rep(NA, R)
  pers <- matrix(NA, R, 3)
  for (b in 1:R) {
    mb <- v$estimar_var(simular(A_dgp, innov(), mod, Y, X), X, p, inicio = v$t0)
    Aest[b, , , ] <- mb$A
    sigma55[b] <- mb$sigma[5, 5]
    mult[b]    <- multiplicador(mb$A)
    pers[b, ]  <- persistencia(mb$A)
    A <- if (is.null(sesgo)) mb$A else corregir(mb$A, sesgo)
    if (v$raices(A)[1] >= 1) next                 # sin estabilidad no existe el multiplicador de largo plazo
    s <- irf_P1(A, mb$sigma)
    irf[b, , ] <- s$r
    angulo[b]  <- s$angulo
  }
  list(irf = irf, Aest = Aest, angulo = angulo, sigma55 = sigma55, mult = mult, pers = pers)
}

t0  <- Sys.time()
res <- list()
res$wild      <- correr(wild, B_wild)                       # primero: con la semilla de 04_estructural.R da sus bandas
res$iid       <- correr(iid, B)
res$bloques   <- correr(function() bloques(l_bloque), B)
res$bloques24 <- correr(function() bloques(l_alterno), B)

# Corrección de sesgo (Kilian 1998; K&L sección 12.3)
#  a) Sesgo de los coeficientes: promedio de A* en el bootstrap por bloques menos la A estimada.
sesgo <- apply(res$bloques$Aest, 2:4, mean) - m$A
#  b) A corregida, con el ajuste de estacionariedad.
A_BC <- corregir(m$A, sesgo)
#  c) Antes de simular se restan las medias (K&L p. 364). A, las Gamma y los residuos no cambian (sólo la constante),
#     pero así el VAR corregido conserva la media de los datos y las series simuladas no arrastran una transición.
Yd <- sweep(v$Y, 2, colMeans(v$Y[m$filas, ]))
Xd <- sweep(v$X, 2, colMeans(v$X[m$filas, ]))
md <- v$estimar_var(Yd, Xd, p, inicio = v$t0)
#  d) Bootstrap desde la A corregida; cada réplica se corrige con el mismo sesgo (bootstrap después del bootstrap).
res$sesgo  <- correr(function() bloques(l_bloque), B, A_dgp = A_BC, sesgo = sesgo, mod = md, Y = Yd, X = Xd)
puntual_BC <- irf_P1(A_BC, m$sigma)
cat(sprintf("Tiempo de los bootstraps: %.1f minutos\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat(sprintf("Con medias restadas, la A estimada no cambia: diferencia máxima %.1e\n\n", max(abs(md$A - m$A))))

# ---- 4. Cuadros -------------------------------------------------------------------------------------------------
cuantil <- function(a, pr) apply(a[!is.na(a[, 1, 1]), , , drop = FALSE], 2:3, quantile, pr)
bandas <- lapply(res, function(x) list(b02 = cuantil(x$irf, .02), b05 = cuantil(x$irf, .05), b16 = cuantil(x$irf, .16),
                                       b50 = cuantil(x$irf, .5), b84 = cuantil(x$irf, .84), b95 = cuantil(x$irf, .95),
                                       b98 = cuantil(x$irf, .98)))   # 2 y 98: banda conjunta de Bonferroni (8 meses)
# Intervalo de Hall (K&L p. 362): refleja la distribución del bootstrap alrededor de la estimación
hall <- list(b16 = 2 * puntual$r - bandas$bloques$b84, b84 = 2 * puntual$r - bandas$bloques$b16)
etiquetas <- c(wild = "Wild (04_estructural.R)", iid = "iid", bloques = "Bloques de 12 meses",
               bloques24 = "Bloques de 24 meses", sesgo = "Bloques de 12 + corrección de sesgo")

celdas <- list(c(5, 0), c(5, 3), c(5, 12), c(2, 3), c(2, 6), c(2, 12), c(2, 24), c(3, 3), c(3, 6), c(1, 12), c(1, 24), c(1, 36))
nombre_celda <- sapply(celdas, function(c) paste0(nombres[c[1]], "_", c[2]))
fila_banda <- function(lo, hi) sapply(celdas, function(c) sprintf("[%+.2f, %+.2f]", lo[c[2] + 1, c[1]], hi[c[2] + 1, c[1]]))
cuadro2 <- rbind(Estimacion = sapply(celdas, function(c) sprintf("%+.2f", puntual$r[c[2] + 1, c[1]])),
                 Estimacion_corregida = sapply(celdas, function(c) sprintf("%+.2f", puntual_BC$r[c[2] + 1, c[1]])),
                 t(sapply(names(bandas), function(k) fila_banda(bandas[[k]]$b16, bandas[[k]]$b84))),
                 Hall_bloques = fila_banda(hall$b16, hall$b84))
colnames(cuadro2) <- nombre_celda
cat("Cuadro 2. Respuesta a un alza de 25 pb: estimación y bandas de 68% (variable_mes)\n")
print(noquote(t(cuadro2)))

diag_metodo <- t(sapply(names(res), function(k) {
  x <- res[[k]]; ok <- !is.na(x$irf[, 1, 1])
  c(validas = sum(ok), deprecia = round(mean(x$irf[ok, 1, 5] > 0), 2),
    angulo_p16 = round(quantile(x$angulo, .16, na.rm = TRUE), 1), angulo_p50 = round(median(x$angulo, na.rm = TRUE), 1),
    angulo_p84 = round(quantile(x$angulo, .84, na.rm = TRUE), 1), multiplicador_p50 = round(median(x$mult), 1),
    de_sigma55 = round(sd(x$sigma55), 3), mediana_ITCR0 = round(bandas[[k]]$b50[1, 5], 2),
    mediana_Infl36 = round(bandas[[k]]$b50[37, 2], 2))
}))
cat(sprintf("\nCuadro 3. Diagnósticos por método (estimación: ángulo %.1f, multiplicador %.1f, sigma_55 %.3f, ITCR_0 %+.2f, Infl_36 %+.2f)\n",
            puntual$angulo, multiplicador(m$A), m$sigma[5, 5], puntual$r[1, 5], puntual$r[37, 2]))
print(diag_metodo)

cat("\nCuadro 4. Meses (0 a 12) en que la banda de 68% excluye el cero por arriba\n")
for (k in names(bandas)) {
  lo <- bandas[[k]]$b16
  lb <- cuantil(res[[k]]$irf, 0.32 / 8 / 2)            # Bonferroni: banda conjunta de 68% para los meses 1 a 8
  cat(sprintf("%-36s inflación: %-26s expectativas: %-26s ITCR: %-10s Bonferroni (inflación, meses 1-8): %s\n",
              etiquetas[k], paste(which(lo[2:13, 2] > 0), collapse = " "), paste(which(lo[2:13, 3] > 0), collapse = " "),
              paste(which(lo[1:13, 5] > 0) - 1, collapse = " "),
              if (isTRUE(any(lb[2:9, 2] > 0))) paste("excluye el cero en", paste(which(lb[2:9, 2] > 0), collapse = " ")) else "no excluye el cero"))
}

# Sesgo de muestra pequeña: en el mundo bootstrap la verdad es la A estimada; si las réplicas salen en promedio menos
# persistentes que ella, el MCO subestima la persistencia.
cat("\nSesgo de la persistencia (suma de coeficientes propios en A(1)); en el bootstrap por bloques la verdad es la A estimada\n")
cat(sprintf("%-13s MCO %.3f | promedio de las réplicas %.3f | corregida %.3f\n", c("Inflación", "Expectativas", "TPM"),
            persistencia(m$A), colMeans(res$bloques$pers), persistencia(A_BC)), sep = "")
cat(sprintf("Multiplicador de largo plazo 1/(valor singular menor de I - A(1)): MCO %.1f | mediana de las réplicas %.1f | corregido %.1f\n",
            multiplicador(m$A), median(res$bloques$mult), multiplicador(A_BC)))
cat(sprintf("Corrección: se aplicó %.0f%% del sesgo estimado; raíz más grande %.3f -> %.3f; ángulo de la P1 %.1f -> %.1f grados\n",
            100 * attr(A_BC, "fraccion"), v$raices(m$A)[1], v$raices(A_BC)[1], puntual$angulo, puntual_BC$angulo))

# ---- 5. Figuras -------------------------------------------------------------------------------------------------
titulos <- c("ISE (nivel, %)", "Inflación anual (pp)", "Expectativas de inflación (pp)", "TPM (pp)",
             "ITCR (nivel, %; + = depreciación)")
sombra <- function(lo, hi, col) polygon(c(0:hmax, hmax:0), c(lo, rev(hi)), col = col, border = NA)
panel  <- function(k, rangos) {
  plot(0:hmax, rep(0, hmax + 1), type = "n", ylim = range(rangos, 0), xlab = "Meses", ylab = "",
       main = titulos[k], cex.main = 1, font.main = 1, las = 1)
}
naranja <- "#1f4e79"

# Figura 1: residuos al cuadrado, en múltiplos de su promedio, con su promedio móvil de 12 meses.
# Los meses con dummy no tienen residuo y quedan fuera del promedio móvil. Las barras se cortan en 10.
tope <- 10
png("figuras/P1_residuos_ARCH.png", width = 2000, height = 1900, res = 220)
par(mfrow = c(3, 2), mar = c(3.2, 3.6, 2.4, 1), mgp = c(2, 0.6, 0))
for (k in 1:K) {
  e2 <- U0[, k]^2 / mean(U0[activos, k]^2)
  e2[!activos] <- NA
  ec <- pmin(e2, tope)                              # barras cortadas en el tope
  pm <- sapply(seq_len(Tn), function(t) {
    w <- ec[max(1, t - 5):min(Tn, t + 6)]
    if (sum(!is.na(w)) >= 6) mean(w, na.rm = TRUE) else NA
  })
  plot(anio, ec, type = "h", col = "grey55", ylim = c(0, tope), las = 1, xlab = "", ylab = "",
       main = sprintf("%s: ARCH-LM p = %s", nombres[k], format(signif(arch[k, 2], 2), scientific = TRUE)), cex.main = 1, font.main = 1)
  arriba <- which(e2 > tope)
  if (length(arriba)) {
    points(anio[arriba], rep(tope, length(arriba)), pch = 17, col = "grey25", cex = 0.8)
    grupo <- cumsum(c(1, diff(arriba) > 6))         # barras cortadas a menos de 6 meses se rotulan juntas
    for (g in unique(grupo)) {
      i <- arriba[grupo == g]
      text(anio[i[1]], tope - 0.4, paste(sprintf("%.0f", e2[i]), collapse = ", "), pos = 2, cex = 0.75, col = "grey25")
    }
  }
  lines(anio, pm, col = naranja, lwd = 2)
  abline(h = 1, lty = 2, col = "grey30")
}
plot.new()
legend("center", c("Residuo al cuadrado / su promedio", "Promedio móvil de 12 meses", "Varianza promedio (= 1)",
                   "Barra cortada en 10 (el número es su valor)"),
       lwd = c(1.5, 2, 1, NA), lty = c(1, 1, 2, NA), pch = c(NA, NA, NA, 17), col = c("grey55", naranja, "grey30", "grey25"),
       bty = "n", cex = 1.0)
dev.off()

# Figura 2: lo que ve cada bootstrap. (a) varianza del residuo del ITCR, centrada en su promedio de cada bootstrap
# (los residuos de MCO tienen varianza promedio Sigma (T - k) / T, un factor común que no cambia la P1), con las
# dos curvas teóricas: cuadrados independientes y varianza de largo plazo; (b) respuesta de impacto del ITCR.
png("figuras/P1_bootstrap_distribuciones.png", width = 2200, height = 1000, res = 220)
par(mfrow = c(1, 2), mar = c(3.6, 3.6, 2.4, 1), mgp = c(2.2, 0.6, 0))
cols <- c(iid = "grey65", wild = "grey25", bloques = naranja, sesgo = "#5b9bd5")
d_a <- lapply(c("iid", "wild", "bloques"), function(k) density(res[[k]]$sigma55 - mean(res[[k]]$sigma55)))
xs  <- seq(-3, 3, length.out = 400)
plot(NA, xlim = c(-3, 3), ylim = c(0, max(sapply(d_a, function(d) max(d$y)))), las = 1,
     xlab = "Varianza del residuo del ITCR, menos su promedio", ylab = "Densidad",
     main = "(a) Incertidumbre sobre Sigma", font.main = 1, cex.main = 1)
for (i in 1:3) lines(d_a[[i]], col = cols[c("iid", "wild", "bloques")[i]], lwd = 2)
lines(xs, dnorm(xs, 0, de_indep), lty = 3, lwd = 1.5)
lines(xs, dnorm(xs, 0, de_largo), lty = 2, lwd = 1.5)
legend("topleft", c("iid", "Wild", "Bloques (12)", "Teoría: cuadrados indep.", "Teoría: largo plazo"),
       col = c(cols[1:3], "black", "black"), lwd = c(2, 2, 2, 1.5, 1.5), lty = c(1, 1, 1, 3, 2), bty = "n", cex = 0.8)
d_b <- lapply(c("wild", "bloques", "sesgo"), function(k) {
  x <- res[[k]]$irf[, 1, 5]; x <- x[!is.na(x)]
  density(x[x > quantile(x, .01) & x < quantile(x, .99)])            # sin el 1% de cada cola, para la densidad
})
x_b <- range(sapply(c("wild", "bloques", "sesgo"), function(k) quantile(res[[k]]$irf[, 1, 5], c(.02, .98), na.rm = TRUE)), 0)
plot(NA, xlim = x_b, ylim = c(0, max(sapply(d_b, function(d) max(d$y)))), las = 1,
     xlab = "ITCR al impacto (%, por 25 pb; + = depreciación)", ylab = "Densidad", main = "(b) Respuesta de impacto del ITCR",
     font.main = 1, cex.main = 1)
for (i in 1:3) lines(d_b[[i]], col = cols[c("wild", "bloques", "sesgo")[i]], lwd = 2)
abline(v = puntual$r[1, 5], lty = 2); abline(v = puntual_BC$r[1, 5], lty = 3, col = cols["sesgo"], lwd = 1.5)
abline(v = 0, col = "grey60")
legend("topright", c("Wild", "Bloques (12)", "Bloques + sesgo", "Estimación", "Estimación corregida"),
       col = c(cols[c("wild", "bloques", "sesgo")], "black", cols["sesgo"]), lwd = c(2, 2, 2, 1, 1.5), lty = c(1, 1, 1, 2, 3),
       bty = "n", cex = 0.9)
dev.off()

# Figura 3: bandas por bloques (68% y 90%) y, en línea punteada, la banda de 68% del wild
png("figuras/irf_P1_bloques_25pb.png", width = 2000, height = 1900, res = 220)
par(mfrow = c(3, 2), mar = c(3.2, 3.6, 2.4, 1), mgp = c(2, 0.6, 0))
bq <- bandas$bloques; bw <- bandas$wild
for (k in 1:K) {
  panel(k, c(bq$b05[, k], bq$b95[, k], bw$b16[, k], bw$b84[, k]))
  sombra(bq$b05[, k], bq$b95[, k], "grey88")
  sombra(bq$b16[, k], bq$b84[, k], "grey70")
  lines(0:hmax, bw$b16[, k], lty = 2, col = "grey30", lwd = 1.5)
  lines(0:hmax, bw$b84[, k], lty = 2, col = "grey30", lwd = 1.5)
  abline(h = 0, lty = 2, col = "grey40")
  lines(0:hmax, puntual$r[, k], lwd = 2, col = naranja)
}
plot.new()
legend("center", c("Respuesta estimada (P1)", "Bloques de 12 meses: banda 68%", "Bloques de 12 meses: banda 90%",
                   "Wild (04_estructural.R): banda 68%"),
       lwd = c(2, 8, 8, 1.5), lty = c(1, 1, 1, 2), col = c(naranja, "grey70", "grey88", "grey30"), bty = "n", cex = 1.0)
dev.off()

# Figura 4: bloques sin y con corrección de sesgo (bandas de 68%)
png("figuras/irf_P1_sesgo_25pb.png", width = 2000, height = 1900, res = 220)
par(mfrow = c(3, 2), mar = c(3.2, 3.6, 2.4, 1), mgp = c(2, 0.6, 0))
bs <- bandas$sesgo
for (k in 1:K) {
  panel(k, c(bq$b16[, k], bq$b84[, k], bs$b16[, k], bs$b84[, k], puntual_BC$r[, k]))
  sombra(bq$b16[, k], bq$b84[, k], "grey70")
  sombra(bs$b16[, k], bs$b84[, k], adjustcolor("#5b9bd5", 0.25))
  abline(h = 0, lty = 2, col = "grey40")
  lines(0:hmax, puntual$r[, k], lwd = 2, col = naranja)
  lines(0:hmax, puntual_BC$r[, k], lwd = 2, lty = 2, col = "#5b9bd5")
}
plot.new()
legend("center", c("Estimación (MCO)", "Estimación con corrección de sesgo", "Bloques: banda 68%",
                   "Bloques + corrección de sesgo: banda 68%"),
       lwd = c(2, 2, 8, 8), lty = c(1, 2, 1, 1), col = c(naranja, "#5b9bd5", "grey70", adjustcolor("#5b9bd5", 0.25)),
       bty = "n", cex = 1.0)
dev.off()

# ---- 6. Se guardan los resultados -------------------------------------------------------------------------------
saveRDS(list(puntual = puntual, puntual_BC = puntual_BC, bandas = bandas, hall = hall, A_BC = A_BC, sesgo = sesgo,
             diagnosticos = diag_metodo, arch = cuadro_arch, angulos = lapply(res, `[[`, "angulo"),
             sigma55 = lapply(res, `[[`, "sigma55"), multiplicador = lapply(res, `[[`, "mult"),
             impacto_ITCR = lapply(res, function(x) x$irf[, 1, 5])), "data/P1_inferencia.rds")
write.csv(data.frame(celda = nombre_celda, t(cuadro2), check.names = FALSE), "data/P1_bandas_25pb.csv", row.names = FALSE)
cat("\nListo: figuras/P1_residuos_ARCH.png, figuras/P1_bootstrap_distribuciones.png, figuras/irf_P1_bloques_25pb.png,",
    "figuras/irf_P1_sesgo_25pb.png; data/P1_inferencia.rds y data/P1_bandas_25pb.csv\n")
sink()
