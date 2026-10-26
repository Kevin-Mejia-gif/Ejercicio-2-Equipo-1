# Efectos de un choque de política monetaria en Colombia

Ejercicio Empírico 2 de Macroeconometría (El Colegio de México, otoño 2026).
Equipo 1: Oscar Daniel Hernández Romero, Kevin Mejía Rodríguez y Rubén Perea Cervantes.

Se identifica un choque de política monetaria en Colombia con un VAR estructural de cinco
variables (ISE, inflación, expectativas de inflación, TPM e ITCR), datos mensuales de
2010 a 2026 y restricciones de corto y largo plazo, siguiendo a Bjørnland (2009).

## Estructura
- `prog/`: código en R (scripts 01 a 13).
- `data/`: datos originales (`data/originales/`) y generados.
- `cuadros/`, `figuras/` y `log/`: salidas de los scripts.
- `trabajo_escrito/`: proyecto en LaTeX del reporte (`main.tex`, `cuadros/`, `figuras/`) y su PDF.
- `presentacion/`: presentación en Beamer (tema `colmex`) y su PDF.

## Cómo replicar
Desde la raíz del repositorio, correr en orden los scripts de `prog/` (01 a 13).
En los scripts que traen `setwd(...)`, cambiar la ruta por la de esta carpeta.
Los pasos y las salidas de cada uno están en `data/LEEME.txt`.
Paquetes de R: `readxl`, `urca` y `forecast`; `12_sesgo_completo_P1.R` usa `parallel`
(viene con R) y tarda unos 17 minutos con 7 núcleos.
