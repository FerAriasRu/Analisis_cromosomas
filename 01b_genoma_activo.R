# ============================================================
# 01b_genoma_activo.R  -  Aplica al entorno los datos del genoma activo (de la lista `genomas`)
# ------------------------------------------------------------
# REQUIERE   : 01 (genomas, genomas_run, dir_salida, bin); `genoma_activo` lo fija el main
#              (si se corre a mano y no existe, toma el primero de genomas_run)
# CREA       : genoma_id, fna_file, gff_file, fai_file, rep_file, trf_file, eq_manual,
#              carpeta_externos, ventana_label, carpeta_salida
# ============================================================
if (!exists("genoma_activo"))
  genoma_activo <- if (!is.null(genomas_run)) genomas_run[1] else names(genomas)[1]
stopifnot(genoma_activo %in% names(genomas))
.G <- genomas[[genoma_activo]]                       # [[ ]] = coincidencia exacta de nombres

genoma_id <- genoma_activo
fna_file  <- .G[["fna"]]; gff_file <- .G[["gff"]]; fai_file <- .G[["fai"]]
rep_file  <- .G[["rep"]]; trf_file <- .G[["trf"]]; eq_manual <- .G[["eq_manual"]]
rm(.G)

ventana_label <- if (bin >= 1e6) {
  paste0(bin / 1e6, "mb")
} else {
  paste0(bin / 1e3, "kb")
}

carpeta_salida <- paste0(
  file.path(dir_salida, paste0("resultados_", genoma_id, "_", ventana_label)),
  "/"
)

carpeta_externos    <- file.path(carpeta_salida, paste0("externos_", genoma_id))
carpeta_referencias <- file.path(carpeta_salida, "externos_referencias")
carpeta_objetos     <- file.path(carpeta_salida, "objetos")

dir.create(carpeta_salida, showWarnings = FALSE, recursive = TRUE)
dir.create(carpeta_externos, showWarnings = FALSE, recursive = TRUE)
dir.create(carpeta_referencias, showWarnings = FALSE, recursive = TRUE)
dir.create(carpeta_objetos, showWarnings = FALSE, recursive = TRUE)