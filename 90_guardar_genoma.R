# ============================================================
# 90_guardar_genoma.R  -  Guarda los objetos del genoma para la comparacion entre genomas
# ------------------------------------------------------------
# REQUIERE   : eq, sl (02/02L) y, si existen, resumen_cromosomas, carac_extremos, sat_bloques,
#              sat350_regiones, microsat_loci, ext_vs_int
# CREA       : <dir_salida>/objetos/objetos_<genoma_id>.rds    (lo lee 16_comparacion.R)
# ============================================================
g0 <- function(n) get0(n, envir = globalenv(), inherits = FALSE)
dir_obj <- carpeta_objetos
dir.create(dir_obj, showWarnings = FALSE, recursive = TRUE)
zonas_comp <- g0("zonas_df")
if (!is.null(zonas_comp) && nrow(zonas_comp) > 0) {
  zonas_comp <- as.data.frame(zonas_comp)
  zonas_comp <- zonas_comp[
    as.character(zonas_comp$tipo) == "Bajo GC + baja densidad génica",
    c("chrom", "zona_id", "start_pb", "end_pb", "longitud_kb")
  ]
}

genes_comp <- NULL
genes_gr_obj <- g0("genes_gr")
if (!is.null(genes_gr_obj) && length(genes_gr_obj) > 0) {
  genes_comp <- data.frame(
    chrom = nombre_chr(as.character(seqnames(genes_gr_obj))),
    start = start(genes_gr_obj),
    end = end(genes_gr_obj),
    strand = as.character(strand(genes_gr_obj))
  )
}

repeticiones_comp <- NULL
rep_red_obj <- g0("rep_red")
if (!is.null(rep_red_obj) && length(rep_red_obj) > 0) {
  repeticiones_comp <- data.frame(
    chrom = nombre_chr(as.character(seqnames(rep_red_obj))),
    start = start(rep_red_obj),
    end = end(rep_red_obj)
  )
}

obj_genoma <- list(
  zonas = zonas_comp,
  genes = genes_comp,
  repeticiones = repeticiones_comp,
  genoma = genoma_id, modo = if (!is.null(g0("tabla"))) "completo" else "ligero",
  ventana = ventana_label, fecha = Sys.time(),
  longitudes = data.frame(chrom = eq$nombre, accession = eq$acc, long_pb = as.numeric(sl), stringsAsFactors = FALSE),
  resumen_cromosomas = g0("resumen_cromosomas"), carac_extremos = g0("carac_extremos"),
  sat_bloques = g0("sat_bloques"), sat350_regiones = g0("sat350_regiones"),
  microsat_loci = g0("microsat_loci"), ext_vs_int = g0("ext_vs_int"),
  parametros_st = c(w_st_max = w_st_max, gap_st = gap_st, k_gaps_st = k_gaps_st, min_st = min_st, dens_st_max = dens_st_max))
archivo_objeto <- file.path(
  dir_obj,
  paste0("objetos_", genoma_id, ".rds")
)

saveRDS(obj_genoma, archivo_objeto)
cat("Objetos guardados en", archivo_objeto, "\n")