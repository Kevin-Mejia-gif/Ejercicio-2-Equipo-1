# Cuadros y figuras del Anexo D (inferencia, resultados y robustez)
#
# Arma los cuadros del Anexo D (cuadros/cuadro_D_*.tex) y dibuja dos figuras de respuestas con bandas de la P1
# (figuras/irf_P1_bloques_solo_25pb.png e irf_P1_metodos_25pb.png), con las salidas de 02, 07 y 12.
# Se corre desde la carpeta Colombia/, después de 01 a 12.

ruta_colombia <- "."

# ---- Datos ------------------------------------------------------------------------------------------------------
v <- readRDS(file.path(ruta_colombia, "data/var_reducido.rds"))       # 02_var.R
x <- readRDS(file.path(ruta_colombia, "data/P1_inferencia.rds"))      # 07_inferencia_P1.R
s <- readRDS(file.path(ruta_colombia, "data/P1_sesgo_completo.rds"))  # 12_sesgo_completo_P1.R
A <- v$modelo$A
K <- dim(A)[1]

# ---- Formato ----------------------------------------------------------------------------------------------------
num   <- function(z, d = 2) { r <- sprintf(paste0("%.", d, "f"), z); sub(paste0("^-(0\\.0{", d, "})$"), "\\1", r) }
mat   <- function(z, d = 2) paste0("$", num(z, d), "$")
banda <- function(lo, hi) if (abs(hi - lo) < 1e-10) "---" else paste0("$[", num(lo), ",\\ ", num(hi), "]$")
valor_p <- function(p) if (p < 0.01) { e <- floor(log10(p)); sprintf("$%.1f\\times10^{%d}$", p / 10^e, e) } else mat(p)
escribir <- function(archivo, lineas) {
  con <- file(file.path("cuadros", archivo), "wb")
  writeLines(enc2utf8(lineas), con, useBytes = TRUE)
  close(con)
}
fuente <- "Fuente: Elaboración propia con datos de DANE y Banco de la República."
unidades <- "Alza de 25 pb en la TPM; ITCR e ISE en nivel (\\%), inflación anual en pp."
tabla <- function(titulo, etiqueta, letra, columnas, cuerpo, nota) {
  c("\\begin{table}[H]", "\\centering", paste0("\\caption{", titulo, "}"), paste0("\\label{", etiqueta, "}"), letra,
    paste0("\\begin{tabular}{", columnas, "}"), "\\toprule", cuerpo, "\\bottomrule", "\\end{tabular}",
    paste0("\\notascuadro{", nota), "", paste0(fuente, "}"), "\\end{table}")
}

# Cuadros con las cuatro celdas del Cuadro 2 de la guía: ITCR mes 0, inflación meses 6 y 24, ISE mes 24
celdas <- list(c(1, 5), c(7, 2), c(25, 2), c(25, 1))                  # (mes + 1, variable)
encabezado <- c(" & ITCR, mes 0 & Inflación, mes 6 & Inflación, mes 24 & ISE, mes 24 \\\\", "\\midrule")
fila_estimacion <- function(nombre, r)
  paste(nombre, "&", paste(sapply(celdas, function(c) mat(r[c[1], c[2]])), collapse = " & "), "\\\\")
fila_banda <- function(nombre, b)
  paste(nombre, "&", paste(sapply(celdas, function(c) banda(b$b16[c[1], c[2]], b$b84[c[1], c[2]])), collapse = " & "), "\\\\")

# ---- D.1.2: prueba ARCH-LM --------------------------------------------------------------------------------------
arch <- x$arch
escribir("cuadro_D_arch.tex", tabla(
  "Prueba ARCH-LM (12 rezagos) y autocorrelación de orden 1 de los residuos y de sus cuadrados", "tab:D_arch", "\\small",
  "@{}lcccc@{}",
  c("Ecuación & ARCH-LM & Valor $p$ & Autocorrelación de $\\hat u_t$ & Autocorrelación de $\\hat u_t^2$ \\\\", "\\midrule",
    sapply(seq_len(nrow(arch)), function(i) sprintf("%s & %s & %s & %s & %s \\\\", arch$Ecuacion[i], num(arch$ARCH_LM[i], 1),
                                                     valor_p(arch$p_valor[i]), mat(arch$autocorr_u[i]), mat(arch$autocorr_u2[i])))),
  "\\textit{Nota.} Sin los 11 meses con \\textit{dummy} puntual, cuyos residuos son cero."))

# ---- D.1.3: estimación y bandas por bloques (Cuadro 2 de la guía: MCO, bloques de 12 y de 24 meses) -------------
escribir("cuadro_D_bandas_bloques.tex", tabla(
  "Estimación y bandas de 68\\% del bootstrap por bloques", "tab:D_bandas_bloques", "\\small", "@{}lcccc@{}",
  c(encabezado,
    fila_estimacion("Estimación (MCO)", x$puntual$r),
    fila_banda("Bloques de 12 meses", x$bandas$bloques),
    fila_banda("Bloques de 24 meses", x$bandas$bloques24)),
  paste("\\textit{Nota.}", unidades, "Percentiles 16 y 84 de 2,000 réplicas.")))

# ---- D.2: resultados generales ----------------------------------------------------------------------------------
horizontes <- c(0, 6, 12, 24, 36)
variables  <- c("ISE (nivel, \\%)", "Inflación anual (pp)", "Expectativas (pp)", "TPM (pp)", "ITCR (nivel, \\%)")
cuerpo <- c(paste0("Variable & ", paste0("Mes ", horizontes, collapse = " & "), " \\\\"), "\\midrule")
for (k in 1:K) {
  cuerpo <- c(cuerpo,
              paste(variables[k], "&", paste(sapply(horizontes, function(h) mat(x$puntual$r[h + 1, k])), collapse = " & "), "\\\\"),
              paste(" &", paste(sapply(horizontes, function(h) banda(x$bandas$bloques$b16[h + 1, k], x$bandas$bloques$b84[h + 1, k])),
                                collapse = " & "), "\\\\"),
              if (k < K) "\\addlinespace")
}
escribir("cuadro_D_resultados_generales.tex", tabla(
  "Respuestas a un alza de 25 pb en la TPM: estimación y banda de 68\\%", "tab:D_resultados", "\\footnotesize",
  "@{}lccccc@{}", cuerpo,
  paste("\\textit{Nota.} Primera fila de cada variable: estimación de MCO; segunda: percentiles 16 y 84 del bootstrap",
        "por bloques de 12 meses (2,000 réplicas). ---: fijado por construcción (ceros de impacto y alza de 25 pb).")))

# ---- D.3.2: persistencia y multiplicador (Cuadro 3 de la guía) --------------------------------------------------
persistencia  <- function(A) diag(apply(A, c(1, 2), sum))                     # suma de coeficientes propios
multiplicador <- function(A) 1 / min(svd(diag(K) - apply(A, c(1, 2), sum))$d)
p_mco <- persistencia(A)[2:4]
p_rep <- persistencia(A + x$sesgo)[2:4]       # promedio de las réplicas: el sesgo es su promedio menos la A de MCO
p_bc  <- persistencia(x$A_BC)[2:4]
nombres_p <- c("Persistencia de la inflación mensual", "Persistencia de las expectativas", "Persistencia de la TPM")
escribir("cuadro_D_persistencia.tex", tabla(
  "Persistencia y multiplicador de largo plazo: MCO, réplicas por bloques y corrección", "tab:D_persistencia", "\\small",
  "@{}lccc@{}",
  c(" & MCO & Réplicas por bloques & Corregido \\\\", "\\midrule",
    sprintf("%s & %s & %s (promedio) & %s \\\\", nombres_p, num(p_mco, 3), num(p_rep, 3), num(p_bc, 3)),
    sprintf("Multiplicador de largo plazo & %s & %s (mediana) & %s \\\\", num(multiplicador(A), 1),
            num(median(x$multiplicador$bloques), 1), num(multiplicador(x$A_BC), 1)),
    sprintf("Ángulo de la P1 & $%s^\\circ$ & $%s^\\circ$ (mediana) & $%s^\\circ$ \\\\", num(x$puntual$angulo, 1),
            num(median(x$angulos$bloques, na.rm = TRUE), 1), num(x$puntual_BC$angulo, 1))),
  paste("\\textit{Nota.} Persistencia: suma de los coeficientes propios en $A(1) = A_1 + \\dots + A_6$.",
        "Multiplicador: $1/\\sigma_{\\min}(I - A(1))$. 2,000 réplicas por bloques de 12 meses.")))

# ---- D.3.2: corrección de sesgo con el atajo y con el algoritmo completo ----------------------------------------
escribir("cuadro_D_sesgo_completo.tex", tabla(
  "Bandas de 68\\% sin corrección de sesgo, con el atajo y con la corrección completa", "tab:D_sesgo_completo",
  "\\small", "@{}lcccc@{}",
  c(encabezado,
    fila_estimacion("Estimación (MCO)", x$puntual$r),
    fila_estimacion("Estimación corregida", x$puntual_BC$r),
    "\\midrule",
    fila_banda("Sin corrección", x$bandas$bloques),
    fila_banda("Con el atajo", s$bandas$atajo),
    fila_banda("Con la corrección completa", s$bandas$completo)),
  paste("\\textit{Nota.}", unidades, "2,000 réplicas por bloques de 12 meses; en la corrección completa, 500 réplicas",
        "anidadas por réplica.")))

# ---- D.3.3: otros métodos de inferencia -------------------------------------------------------------------------
escribir("cuadro_D_otros_metodos.tex", tabla(
  "Bandas de 68\\% con distintos métodos de inferencia", "tab:D_otros_metodos", "\\small", "@{}lcccc@{}",
  c(encabezado,
    fila_estimacion("Estimación (MCO)", x$puntual$r),
    fila_banda("Estándar (iid)", x$bandas$iid),
    fila_banda("\\textit{Wild}", x$bandas$wild),
    fila_banda("Bloques de 12 meses", x$bandas$bloques),
    fila_banda("Hall (bloques de 12 meses)", x$hall)),
  paste("\\textit{Nota.}", unidades, "\\textit{Wild}: 1,000 réplicas; iid y bloques: 2,000.")))

# ---- Figuras de la P1 con las bandas de 07 (mismo estilo que 07) --------------------------------------------------
titulos <- c("ISE (nivel, %)", "Inflación anual (pp)", "Expectativas de inflación (pp)", "TPM (pp)",
             "ITCR (nivel, %; + = depreciación)")
r  <- x$puntual$r
bq <- x$bandas$bloques; bw <- x$bandas$wild; bh <- x$hall
meses <- 0:(nrow(r) - 1)
naranja <- "grey30"; azul <- "#1f4e79"
sombra <- function(lo, hi, col) polygon(c(meses, rev(meses)), c(lo, rev(hi)), col = col, border = NA)
panel  <- function(k, rangos) {
  plot(meses, rep(0, length(meses)), type = "n", ylim = range(rangos, 0), xlab = "Meses", ylab = "",
       main = titulos[k], cex.main = 1, font.main = 1, las = 1)
}

# D.1.3: sólo el bootstrap por bloques, bandas de 68% y 90%
png("figuras/irf_P1_bloques_solo_25pb.png", width = 2000, height = 1900, res = 220)
par(mfrow = c(3, 2), mar = c(3.2, 3.6, 2.4, 1), mgp = c(2, 0.6, 0))
for (k in 1:K) {
  panel(k, c(bq$b05[, k], bq$b95[, k], r[, k]))
  sombra(bq$b05[, k], bq$b95[, k], "grey88")
  sombra(bq$b16[, k], bq$b84[, k], "grey70")
  abline(h = 0, lty = 2, col = "grey40")
  lines(meses, r[, k], lwd = 2)
}
plot.new()
legend("center", c("Respuesta estimada (P1)", "Bloques de 12 meses: banda 68%", "Bloques de 12 meses: banda 90%"),
       lwd = c(2, 8, 8), col = c("black", "grey70", "grey88"), bty = "n", cex = 1.05)
dev.off()

# D.3.3: banda de 68% por bloques contra las de wild y Hall
png("figuras/irf_P1_metodos_25pb.png", width = 2000, height = 1900, res = 220)
par(mfrow = c(3, 2), mar = c(3.2, 3.6, 2.4, 1), mgp = c(2, 0.6, 0))
for (k in 1:K) {
  panel(k, c(bq$b16[, k], bq$b84[, k], bw$b16[, k], bw$b84[, k], bh$b16[, k], bh$b84[, k], r[, k]))
  sombra(bq$b16[, k], bq$b84[, k], "grey80")
  for (lim in list(bw$b16, bw$b84)) lines(meses, lim[, k], lty = 2, col = naranja, lwd = 1.5)
  for (lim in list(bh$b16, bh$b84)) lines(meses, lim[, k], col = azul, lwd = 1.5)
  abline(h = 0, lty = 2, col = "grey40")
  lines(meses, r[, k], lwd = 2)
}
plot.new()
legend("center", c("Respuesta estimada (P1)", "Bloques de 12 meses: banda 68%", "Wild: banda 68%",
                   "Hall (bloques de 12 meses): banda 68%"),
       lwd = c(2, 8, 1.5, 1.5), lty = c(1, 1, 2, 1), col = c("black", "grey80", naranja, azul), bty = "n", cex = 1.0)
dev.off()

cat("Listo: 2 figuras dibujadas y 6 cuadros en cuadros/\n")
