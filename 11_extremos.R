# ============================================================
# 11_extremos.R  -  ¿Repeticiones concentradas en los extremos? (extremo vs interior)
# ------------------------------------------------------------
# REQUIERE   : 01 (ext_umbrales_kb), 09 (db), 10 (ext_ref_kb)
# CREA       : vars_ext, defs_ext
#              ext_vs_int       Wilcoxon extremo vs interior (varias definiciones de extremo)
#              ext_chr          lo mismo por cromosoma, con FDR y veredicto del patrón
#              cor_dist         Spearman de cada variable vs distancia al extremo
#              perfil_dist      mediana por clase de distancia
# USADO POR  : 14_excel (hojas Ext_*), 15_figuras (usa db y mid_d directamente)
# ============================================================

# ============================================================
# 12. ¿REPETICIONES CONCENTRADAS EN LOS EXTREMOS?
#     Conjunto `db` (base): incluye los extremos pero sin N ni outliers.
# ============================================================

vars_ext <- c(repeat_pct = "% repeticiones", repeat_sinTelo_pct = "% repeticiones (sin bloques teloméricos)",
              genes_per_Mb = "Genes/Mb", gene_pct = "% cobertura génica", intergenic_pct = "% intergénico",
              intron_pct = "% intrones", GC_percent = "GC (%)")

defs_ext <- c(setNames(lapply(ext_umbrales_kb, function(u) db$dist_extremo_pb <= u * 1000),
                       paste0("<= ", ext_umbrales_kb, " kb del borde")),
              list("10% terminal de cada cromosoma" = db$region == "Extremos (10%)"))

ext_vs_int <- do.call(rbind, lapply(names(defs_ext), function(nd) {
  g <- defs_ext[[nd]]
  do.call(rbind, lapply(names(vars_ext), function(v) {
    x <- db[[v]][g]; y <- db[[v]][!g]; x <- x[!is.na(x)]; y <- y[!is.na(y)]
    if (length(x) < 3 || length(y) < 3) return(NULL)
    w <- wilcox.test(x, y, exact = FALSE)
    data.frame(definicion_extremo = nd, variable = vars_ext[v], n_extremo = length(x), n_interior = length(y),
               mediana_extremo = round(median(x), 2), mediana_interior = round(median(y), 2),
               razon_ext_int = ifelse(median(y) > 0, round(median(x) / median(y), 2), NA),
               cliff_delta = round(2 * unname(w$statistic) / (length(x) * length(y)) - 1, 3),
               p_value = signif(w$p.value, 3), stringsAsFactors = FALSE)
  }))
}))
ext_vs_int$p_FDR <- signif(p.adjust(ext_vs_int$p_value, "BH"), 3)
rownames(ext_vs_int) <- NULL

ext_chr <- do.call(rbind, lapply(levels(db$chrom), function(ch) {
  s <- db[db$chrom == ch, ]; g <- s$dist_extremo_pb <= ext_ref_kb * 1000
  if (sum(g) < 3 || sum(!g) < 3) return(NULL)
  pr <- wilcox.test(s$repeat_pct[g], s$repeat_pct[!g], exact = FALSE)$p.value
  pg <- wilcox.test(s$genes_per_Mb[g], s$genes_per_Mb[!g], exact = FALSE)$p.value
  data.frame(chrom = ch, n_ext = sum(g), n_int = sum(!g),
             repeat_med_ext = round(median(s$repeat_pct[g], na.rm = TRUE), 2),
             repeat_med_int = round(median(s$repeat_pct[!g], na.rm = TRUE), 2), p_repeat = signif(pr, 3),
             genes_Mb_med_ext = round(median(s$genes_per_Mb[g], na.rm = TRUE), 2),
             genes_Mb_med_int = round(median(s$genes_per_Mb[!g], na.rm = TRUE), 2), p_genes = signif(pg, 3),
             stringsAsFactors = FALSE)
}))
ext_chr$p_repeat_FDR <- signif(p.adjust(ext_chr$p_repeat, "BH"), 3)
ext_chr$p_genes_FDR  <- signif(p.adjust(ext_chr$p_genes, "BH"), 3)
ext_chr$patron <- ifelse(ext_chr$repeat_med_ext > ext_chr$repeat_med_int & ext_chr$genes_Mb_med_ext < ext_chr$genes_Mb_med_int,
                         ifelse(ext_chr$p_repeat_FDR < 0.05 & ext_chr$p_genes_FDR < 0.05,
                                "Más repetición y menos genes en extremos (significativo)",
                                "Más repetición y menos genes en extremos (no significativo)"),
                         "No cumple el patrón")
rownames(ext_chr) <- NULL

cor_dist <- do.call(rbind, lapply(names(vars_ext), function(v) {
  ok <- !is.na(db[[v]])
  if (length(unique(db[[v]][ok])) < 2) return(NULL)   # variable constante: no hay correlación
  ct <- cor.test(db$dist_extremo_pb[ok], db[[v]][ok], method = "spearman", exact = FALSE)
  por_chr <- sapply(levels(db$chrom), function(ch) {
    s <- db[db$chrom == ch & ok, ]
    if (nrow(s) < 10 || length(unique(s[[v]])) < 2) return(NA_real_)
    unname(cor(s$dist_extremo_pb, s[[v]], method = "spearman"))
  })
  data.frame(variable = vars_ext[v], rho_global_vs_distancia = round(unname(ct$estimate), 3),
             p_value = signif(ct$p.value, 3), rho_mediana_por_chr = round(median(por_chr, na.rm = TRUE), 3),
             chr_rho_neg = sum(por_chr < 0, na.rm = TRUE), chr_rho_pos = sum(por_chr > 0, na.rm = TRUE),
             interpretacion = ifelse(unname(ct$estimate) < 0, "Aumenta hacia los extremos", "Disminuye hacia los extremos"),
             stringsAsFactors = FALSE)
}))
rownames(cor_dist) <- NULL

perfil_dist <- db %>% group_by(dist_clase) %>%
  summarise(n_ventanas = n(),
            repeat_mediana = round(median(repeat_pct, na.rm = TRUE), 2),
            repeat_Q1 = round(quantile(repeat_pct, 0.25, na.rm = TRUE), 2),
            repeat_Q3 = round(quantile(repeat_pct, 0.75, na.rm = TRUE), 2),
            repeat_sinTelo_mediana = round(median(repeat_sinTelo_pct, na.rm = TRUE), 2),
            genes_Mb_mediana = round(median(genes_per_Mb, na.rm = TRUE), 2),
            genes_Mb_Q1 = round(quantile(genes_per_Mb, 0.25, na.rm = TRUE), 2),
            genes_Mb_Q3 = round(quantile(genes_per_Mb, 0.75, na.rm = TRUE), 2),
            gene_pct_mediana = round(median(gene_pct, na.rm = TRUE), 2),
            intron_pct_mediana = round(median(intron_pct, na.rm = TRUE), 2),
            GC_mediana = round(median(GC_percent, na.rm = TRUE), 2), .groups = "drop") %>% as.data.frame()
perfil_dist$dist_clase <- as.character(perfil_dist$dist_clase)
cat("\n=== Extremos vs interior ===\n"); print(ext_vs_int[ext_vs_int$variable %in% c("% repeticiones", "Genes/Mb"), ])
