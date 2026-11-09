# ============================================================
# 00_main.R  -  Script principal: UNO o VARIOS genomas, todo en una sola corrida
# ------------------------------------------------------------
# USO (RStudio):
#   1. Abre el proyecto (git) cuya carpeta contiene TODOS estos .R
#   2. Edita SOLO 01_parametros.R: lista `genomas`, `genomas_run` (cuales correr) y umbrales
#   3. Corre este archivo (Source).
#
# QUE HACE POR GENOMA (genomas_run):
#   modo "completo" (el genoma tiene GFF): modulos 02 a 15 + SAT350 + microsatelites -> 1 Excel + figuras
#   modo "ligero"   (sin GFF)            : cromosomas + SAT350 + microsatelites (para comparar)
#   Cada genoma deja sus resultados en resultados_<id>_<ventana>/ y sus objetos en objetos/objetos_<id>.rds
# DESPUES DEL BUCLE: si hay >= 2 genomas (corridos ahora o con .rds previos, ver `genomas_comparar`),
#   16_comparacion.R arma las tablas comparativas (+ comparacion_genomas.xlsx) y 15_figuras.R dibuja la
#   comparacion. Con UN solo genoma la comparacion se omite y todo lo demas corre igual.
#
# ORDEN DE MODULOS (completo):
#   01 parametros (incluye el genoma activo) -> 02 datos_base -> 03 intergenico -> 04 ventanas -> 05 telomeros
#   -> 06 subtelomeros -> 06b sat350 (GenBank/BLASTN/TRF) -> 07 satelites -> 07b microsatelites
#   -> 08 intergenico_metricas -> 09 filtro_banderas -> 10 repeticiones -> 11 extremos -> 12 estadistica
#   -> 13 zonas_criticas -> 14 excel -> 90 guardar_genoma -> 15 figuras (UNICO script de graficas)
#
# TRABAJAR UN SOLO MODULO: pon `hasta <- "NN_modulo.R"` (p. ej. "13_zonas_criticas.R"), corre el main y luego
# edita y corre solo ese archivo con Source (los objetos del ultimo genoma siguen en memoria).
# ============================================================
rm(list = ls())

# -------------------------------------------------------------------------


paquetes <- c("GenomicFeatures", "GenomicRanges", "GenomeInfoDb", "Biostrings", "rtracklayer",
              "txdbmaker", "ggplot2", "patchwork", "dplyr", "openxlsx", "BSgenome", "scales")
faltan <- paquetes[!vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltan)) stop("Faltan paquetes: ", paste(faltan, collapse = ", "),
                         "\n  CRAN: install.packages(...) | Bioconductor: BiocManager::install(...)")
suppressPackageStartupMessages({
  library(GenomicFeatures); library(GenomicRanges); library(GenomeInfoDb); library(Biostrings)
  library(rtracklayer);     library(txdbmaker);     library(ggplot2);      library(patchwork)
  library(grid);            library(dplyr);         library(openxlsx);     library(BSgenome)
})

# Carpeta de los scripts (absoluta: 01_parametros.R cambia el directorio de trabajo)
dir_scripts <- tryCatch(dirname(rstudioapi::getSourceEditorContext()$path), error = function(e) "")
if (!nzchar(dir_scripts) || dir_scripts == ".") dir_scripts <- getwd()
dir_scripts <- normalizePath(dir_scripts)

modulos_completo <- c(
  "01_parametros.R","01b_genoma_activo.R", "02_datos_base.R", "03_intergenico.R", "04_ventanas.R",
  "05_telomeros.R", "06_subtelomeros.R", "06b_sat350.R", "07_satelites.R", "07b_microsatelites.R",
  "08_intergenico_metricas.R", "09_filtro_banderas.R", "10_repeticiones.R", "11_extremos.R",
  "12_estadistica.R", "13_zonas_criticas.R", "14_excel.R", "90_guardar_genoma.R", "15_figuras.R","17_preparar_proteinas.R,
  18_graficas_pca_proteinas.R")
modulos_ligero <- c(
  "01_parametros.R","01b_genoma_activo.R", "02l_base_ligera.R", "06b_sat350.R",
  "07b_microsatelites.R","15F08_idiograma_ligero.R", "90_guardar_genoma.R,17_preparar_proteinas.R,
  18_graficas_pca_proteinas.R")

# -------------------------------------------------------------------------


hasta  <- NULL            # nombre del modulo donde parar (p. ej. "13_zonas_criticas.R"); NULL = todos
saltar <- character(0)    # p. ej. c("07_satelites.R") -- ojo: otros modulos dependen de lo que crea

correr <- function(f) {
  cat("\n>>> ", f, "\n", sep = "")
  t0 <- Sys.time()
  source(file.path(dir_scripts, f), encoding = "UTF-8")      # se evalua en el entorno global
  cat("    listo en ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s\n", sep = "")
}

# Se lee 01 una vez para conocer `genomas` y `genomas_run`
source(file.path(dir_scripts, "01_parametros.R"), encoding = "UTF-8")
genomas_run_main <- if (is.null(genomas_run)) names(genomas) else genomas_run
stopifnot(all(genomas_run_main %in% names(genomas)))
modos_gen <- vapply(genomas, function(x) {
  if (!is.null(x[["modo"]])) x[["modo"]]
  else if (!is.null(x[["gff"]]) && file.exists(x[["gff"]])) "completo" else "ligero"
}, character(1))

.keep <- c("paquetes", "faltan", "dir_scripts", "modulos_completo", "modulos_ligero", "hasta", "saltar",
           "correr", "genomas", "genomas_run_main", "modos_gen", "fig_modo", ".keep", ".gen",
           "genoma_activo", "resumen_corridas")
resumen_corridas <- list()

# ---- Bucle por genoma: cada uno parte de un entorno limpio ----
for (.gen in genomas_run_main) {
  rm(list = setdiff(ls(all.names = TRUE), .keep))
  genoma_activo <- .gen
  fig_modo <- "genoma"
  mods <- if (modos_gen[[.gen]] == "completo") modulos_completo else modulos_ligero
  if (!is.null(hasta)) mods <- mods[seq_len(match(hasta, mods, nomatch = length(mods)))]
  cat("\n################ GENOMA: ", .gen, " (modo ", modos_gen[[.gen]], ") ################\n", sep = "")
  ok <- tryCatch({ for (f in mods) if (!f %in% saltar) correr(f); TRUE },
                 error = function(e) { message("\n*** El genoma '", .gen, "' fallo: ", conditionMessage(e)); FALSE })
  resumen_corridas[[.gen]] <- c(modo = modos_gen[[.gen]], ok = ok,
                                salida = get0("carpeta_salida", ifnotfound = NA_character_))
}

# ---- Comparacion entre genomas (solo si hay >= 2) ----
fig_modo <- "comparacion"
cmp_ok <- FALSE
correr("16_comparacion.R")
if (isTRUE(cmp_ok)) correr("15_figuras.R")
fig_modo <- "genoma"      # deja el modo por defecto para correr 15_figuras.R a mano

cat("\n=== Resumen de la corrida ===\n")
for (n in names(resumen_corridas)) with(as.list(resumen_corridas[[n]]),
                                        cat(sprintf("%-18s modo %-8s %s  -> %s\n", n, modo, ifelse(ok == "TRUE", "OK    ", "FALLO "), salida)))
if (isTRUE(cmp_ok)) cat("Comparacion en:", carpeta_cmp, "\n")