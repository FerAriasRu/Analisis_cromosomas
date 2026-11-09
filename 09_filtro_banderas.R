# ============================================================
# 09_filtro_banderas.R  -  Filtro de ventanas, outliers, banderas y conjuntos d / db
# ------------------------------------------------------------
# REQUIERE   : 01 (bin, max_N_pct, min_ancho_frac, k_mad, mad_vars, excluir_*, umbral_pct),
#              02 (genes_gr, tiene_rep), 04 (tabla, bins, wn), 05 (tabla$telomere), 06 (tabla$subtelomero)
# CREA       : fuera_mad(x, k)
#              tabla$valida, outlier_MAD, base, analisis
#              tabla$bajo_GC, baja_dens_genica, alto_repeat, bajo_repeat, alto_intron, alto_intergenico
#              tabla$region_interes
#              umbrales   (percentiles por cromosoma)
#              d          = tabla[analisis]   conjunto de ANÁLISIS (interior filtrado)
#              db         = tabla[base]       conjunto BASE (incluye extremos; para análisis de extremos)
#              ventanas_revision   ventanas excluidas/marcadas con motivo y genes
# USADO POR  : 10 (banderas), 11 (db), 12 (d), 13 (banderas), 14, 15, 16, 17
# ============================================================

# ============================================================
# 10. FILTRO DE VENTANAS Y BANDERAS (una sola vez)
# ============================================================

tabla$valida <- !is.na(tabla$GC_percent) & tabla$N_pct <= max_N_pct &
  tabla$width >= min_ancho_frac * bin & wn > 0

fuera_mad <- function(x, k) {
  s <- mad(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(FALSE, length(x)))
  abs(x - median(x, na.rm = TRUE)) / s > k
}
tabla$outlier_MAD <- FALSE
for (ch in levels(tabla$chrom)) {
  i <- which(tabla$chrom == ch & tabla$valida)
  if (length(i) < 10) next
  o <- rep(FALSE, length(i))
  for (v in mad_vars) o <- o | fuera_mad(tabla[[v]][i], k_mad)
  tabla$outlier_MAD[i] <- o
}
tabla$base     <- tabla$valida & (!excluir_outliers_MAD | !tabla$outlier_MAD)
tabla$analisis <- tabla$base & (!excluir_extremos_analisis | !tabla$ext_tecnico)

for (v in c("bajo_GC", "baja_dens_genica", "alto_repeat", "bajo_repeat", "alto_intron", "alto_intergenico"))
  tabla[[v]] <- FALSE
umbrales <- list()
for (ch in levels(tabla$chrom)) {
  i <- which(tabla$chrom == ch & tabla$analisis)
  if (length(i) < 10) next
  qlo <- function(v) unname(quantile(tabla[[v]][i], umbral_pct, na.rm = TRUE))
  qhi <- function(v) unname(quantile(tabla[[v]][i], 1 - umbral_pct, na.rm = TRUE))
  u <- data.frame(chrom = ch, n_ventanas = length(i),
                  umbral_GC_bajo = qlo("GC_percent"), umbral_genes_Mb_bajo = qlo("genes_per_Mb"),
                  umbral_repeat_alto = qhi("repeat_pct"), umbral_repeat_bajo = qlo("repeat_pct"),
                  umbral_intron_alto = qhi("intron_pct"), umbral_intergenico_alto = qhi("intergenic_pct"))
  umbrales[[ch]] <- u
  tabla$bajo_GC[i]          <- tabla$GC_percent[i]     <= u$umbral_GC_bajo
  tabla$baja_dens_genica[i] <- tabla$genes_per_Mb[i]   <= u$umbral_genes_Mb_bajo   # incluye ventanas con 0 genes
  tabla$alto_repeat[i]      <- tabla$repeat_pct[i]     >= u$umbral_repeat_alto
  tabla$bajo_repeat[i]      <- tabla$repeat_pct[i]     <= u$umbral_repeat_bajo
  tabla$alto_intron[i]      <- tabla$intron_pct[i]     >= u$umbral_intron_alto
  tabla$alto_intergenico[i] <- tabla$intergenic_pct[i] >= u$umbral_intergenico_alto
}
umbrales <- do.call(rbind, umbrales); if (!is.null(umbrales)) rownames(umbrales) <- NULL

# Sin información de repeticiones las banderas de repetición no tienen sentido
if (!tiene_rep) { tabla$alto_repeat <- FALSE; tabla$bajo_repeat <- FALSE }

tabla$region_interes <- case_when(
  tabla$telomere == "Si"                                  ~ "Región telomérica",
  tabla$subtelomero == "Si"                               ~ "Región subtelomérica",
  tabla$bajo_GC & tabla$baja_dens_genica                  ~ "Bajo GC + baja densidad génica",
  tabla$alto_repeat & tabla$baja_dens_genica              ~ "Alta repetición + baja densidad génica",
  tabla$alto_repeat & tabla$bajo_GC                       ~ "Alta repetición + bajo GC",
  tabla$alto_intergenico & tabla$baja_dens_genica         ~ "Intergénica con pocos genes",
  TRUE                                                    ~ "Sin característica extrema")

d <- tabla[tabla$analisis, ]          # conjunto de análisis (interior filtrado)
d$grupo <- factor(ifelse(d$bajo_GC, "Bajo GC", "Resto"), levels = c("Bajo GC", "Resto"))
db <- tabla[tabla$base, ]             # conjunto base (incluye extremos) para el análisis de extremos

cat("\n=== Filtro de ventanas ===\n")
cat("Ventanas totales:", nrow(tabla), "\n")
cat("  No válidas (N > ", max_N_pct, "%, ventana corta o sin ACGT): ", sum(!tabla$valida), "\n", sep = "")
cat("  Outliers MAD en ", paste(mad_vars, collapse = ","), " (entre válidas): ", sum(tabla$outlier_MAD),
    if (excluir_outliers_MAD) " -> EXCLUIDAS" else " -> solo marcadas", "\n", sep = "")
cat("  Ventanas terminales técnicas (primeros/últimos ", w_excluir_ext, " pb): ", sum(tabla$valida & tabla$ext_tecnico),
    if (excluir_extremos_analisis) " -> EXCLUIDAS del análisis" else " -> incluidas", "\n", sep = "")
cat("  Conjunto base:", nrow(db), "| conjunto de análisis:", nrow(d), "\n")
cat("  Bajo GC:", sum(d$bajo_GC), "| baja dens. génica:", sum(d$baja_dens_genica),
    "| alta repetición:", sum(d$alto_repeat), "\n")

# Ventanas excluidas/marcadas con los genes que contienen
rev_i <- which(!tabla$analisis | tabla$outlier_MAD)
mm <- cbind(ifelse(!is.na(tabla$N_pct[rev_i]) & tabla$N_pct[rev_i] > max_N_pct, "N > umbral", ""),
            ifelse(tabla$width[rev_i] < min_ancho_frac * bin, "ventana terminal corta", ""),
            ifelse(is.na(tabla$GC_percent[rev_i]), "sin ACGT", ""),
            ifelse(tabla$outlier_MAD[rev_i], "Outlier MAD", ""),
            ifelse(tabla$valida[rev_i] & tabla$ext_tecnico[rev_i] & excluir_extremos_analisis, "Ventana terminal técnica", ""))
ov_x <- findOverlaps(bins[rev_i], genes_gr, ignore.strand = TRUE)
gl_x <- tapply(genes_gr$gene_id[subjectHits(ov_x)], factor(queryHits(ov_x), levels = seq_along(rev_i)),
               paste, collapse = ";")
ventanas_revision <- data.frame(
  chrom = as.character(tabla$chrom[rev_i]), start = tabla$start[rev_i], end = tabla$end[rev_i],
  width = tabla$width[rev_i], GC_percent = round(tabla$GC_percent[rev_i], 2),
  N_pct = round(tabla$N_pct[rev_i], 2), repeat_pct = round(tabla$repeat_pct[rev_i], 2),
  genes_per_Mb = round(tabla$genes_per_Mb[rev_i], 2),
  motivo = apply(mm, 1, function(x) paste(x[x != ""], collapse = " + ")),
  excluida_del_analisis = ifelse(tabla$analisis[rev_i], "No (solo marcada)", "Si"),
  n_genes = tabulate(queryHits(ov_x), nbins = length(rev_i)),
  genes = ifelse(is.na(gl_x), "", gl_x), stringsAsFactors = FALSE)
