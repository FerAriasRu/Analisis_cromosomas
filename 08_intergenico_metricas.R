# ============================================================
# 08_intergenico_metricas.R  -  Métricas del intergénico por región y por orientación
# ------------------------------------------------------------
# REQUIERE   : 01 (gap_st), 02 (fa, exon_red, intrones_red, prom_red, rep_red),
#              03 (inter_all, inter_gr, ok_i), 04 (gc_de, cobertura_bp), 06 (st_gr)
# CREA       : n_ig
#              inter_all  (REEMPLAZADO: añade long_pb, GC, N_pct, repeat_pct, promoter_pct, en_subtelomero)
#              ib         solo tramos entre_genes con dist_pb > 0
#              res_orient, res_dist, kw (Kruskal-Wallis), pw (pairwise Wilcoxon)
#              elem, resumen_elem   GC y longitud de exón / intrón / intergénico
# USADO POR  : 14_excel (inter_all, res_*, resumen_elem, pw), 12_estadistica (nada directo),
#              15_figuras (ib)
# ============================================================

# ============================================================
# 9. INTERGÉNICO: MÉTRICAS POR REGIÓN Y POR ORIENTACIÓN
# ============================================================

n_ig <- sum(ok_i)
inter_all$long_pb <- ifelse(ok_i, inter_all$dist_pb, NA)
for (v in c("GC", "N_pct", "repeat_pct", "promoter_pct")) inter_all[[v]] <- NA_real_
inter_all$en_subtelomero <- NA_character_
if (n_ig > 0) {
  s_ig <- getSeq(fa, inter_gr)
  inter_all$GC[ok_i]           <- round(gc_de(inter_gr), 2)
  inter_all$N_pct[ok_i]        <- round(100 * letterFrequency(s_ig, "N")[, 1] / width(inter_gr), 2)
  inter_all$repeat_pct[ok_i]   <- round(100 * cobertura_bp(inter_gr, rep_red) / width(inter_gr), 2)
  inter_all$promoter_pct[ok_i] <- round(100 * cobertura_bp(inter_gr, prom_red) / width(inter_gr), 2)
  inter_all$en_subtelomero[ok_i] <- ifelse(overlapsAny(inter_gr, st_gr, ignore.strand = TRUE), "Si", "No")
}
inter_all <- inter_all[, c("chrom", "accession", "tipo", "start", "end", "dist_pb", "long_pb", "orientacion",
                           "gen_izq", "gen_der", "GC", "N_pct", "repeat_pct", "promoter_pct", "en_subtelomero")]
ib <- inter_all[inter_all$tipo == "entre_genes" & inter_all$dist_pb > 0, ]

res_orient <- inter_all %>% group_by(orientacion) %>%
  summarise(n_pares = n(), n_con_region = sum(dist_pb > 0),
            long_mediana_pb = median(dist_pb[dist_pb > 0]),
            long_Q1_pb = quantile(dist_pb[dist_pb > 0], 0.25, names = FALSE),
            long_Q3_pb = quantile(dist_pb[dist_pb > 0], 0.75, names = FALSE),
            long_total_Mb = round(sum(dist_pb[dist_pb > 0]) / 1e6, 3),
            GC_mediana = round(median(GC, na.rm = TRUE), 2),
            repeat_mediana_pct = round(median(repeat_pct, na.rm = TRUE), 2),
            .groups = "drop") %>% as.data.frame()

res_dist <- data.frame(
  estadistico = c("n tramos entre genes", "percentil 5 (pb)", "percentil 25 (pb)", "mediana (pb)",
                  "percentil 75 (pb)", "percentil 95 (pb)",
                  paste0("% tramos > ", gap_st, " pb (umbral subtelómero)")),
  valor = c(nrow(ib), quantile(ib$dist_pb, c(.05, .25, .5, .75, .95), names = FALSE),
            round(100 * mean(ib$dist_pb > gap_st), 2)))

kw <- tryCatch(kruskal.test(dist_pb ~ orientacion, data = ib), error = function(e) NULL)
pw <- tryCatch(pairwise.wilcox.test(ib$dist_pb, ib$orientacion, p.adjust.method = "BH", exact = FALSE),
               error = function(e) NULL)
cat("\n=== Intergénico por orientación de los genes flanqueantes ===\n"); print(res_orient)
if (!is.null(kw)) cat("Kruskal-Wallis (longitud ~ orientación): p =", signif(kw$p.value, 3), "\n")

elem <- rbind(
  data.frame(tipo = "Exón",        GC = gc_de(exon_red),     long = width(exon_red)),
  data.frame(tipo = "Intrón",      GC = gc_de(intrones_red), long = width(intrones_red)),
  data.frame(tipo = "Intergénico", GC = inter_all$GC[ok_i],  long = width(inter_gr)))
elem <- elem[is.finite(elem$GC), ]
elem$tipo <- factor(elem$tipo, levels = c("Exón", "Intrón", "Intergénico"))
resumen_elem <- elem %>% group_by(tipo) %>%
  summarise(n = n(), GC_mediana = round(median(GC), 2), long_mediana_pb = median(long),
            total_Mb = round(sum(long) / 1e6, 2), .groups = "drop") %>% as.data.frame()
p_intr_exon <- tryCatch(wilcox.test(GC ~ tipo, data = droplevels(elem[elem$tipo != "Intergénico", ]),
                                    exact = FALSE)$p.value, error = function(e) NA_real_)
resumen_elem$p_intron_vs_exon_wilcoxon <- ifelse(resumen_elem$tipo == "Intrón", signif(p_intr_exon, 3), NA)
