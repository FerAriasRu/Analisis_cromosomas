# ============================================================
# 12_estadistica.R  -  Resumen por cromosoma y estadística sobre el conjunto d
# ------------------------------------------------------------
# REQUIERE   : 02 (fa, sl, eq, genes_gr), 04 (tabla), 05 (tel_df), 06 (carac_extremos, estado),
#              07 (sat_bloques), 01 (patron_sat), 09 (d, tabla$chrom), 10 (L_sl, rep_bp_chr, gene_bp_chr)
# CREA       : gc_chr, n_gen_chr, resumen_cromosomas
#              cor_chr (GC vs genes/Mb por cromosoma), mh_res (Mantel-Haenszel)
#              vars                 variables comparadas (nombre -> etiqueta)
#              comparar(g1, g2, e1, e2)   función: Wilcoxon + Cliff's delta + FDR
#              cmp_bajoGC, cmp_bajoDen, cmp_intron, cmp_rep_AB, cmp_rep_alto, cmp_repeticion
#              perfil_rep, perfil, cor_vs(), cor_tab, het_gc
#              mod, or_tab, bg, res_perfil, mod_gc, coef_tab
#              rho_ig, tab_ig, fi_ig, jacc, concordancia
# MODIFICA   : d$grupo_rep (nueva columna)
# ROBUSTEZ   : los modelos (mod, mod_gc) usan solo variables numéricas con variación (helpers
#              preds_validos / formula_dinamica); cor_vs() omite variables vacías o constantes.
# USADO POR  : 14_excel (hojas Est_*, Repeticion_*, Inter_concordancia, Resumen_cromosomas),
#              15_figuras (rho_ig, jacc, vars)
# ============================================================

# ============================================================
# 13. RESUMEN POR CROMOSOMA Y ESTADÍSTICA (conjunto `d`)
# ============================================================

gc_chr <- letterFrequency(fa, "GC")[, 1] / rowSums(letterFrequency(fa, c("A", "C", "G", "T"))) * 100
n_gen_chr <- as.integer(table(factor(as.character(seqnames(genes_gr)), levels = names(sl))))

resumen_cromosomas <- do.call(rbind, lapply(seq_along(sl), function(k) {
  a  <- names(sl)[k]; ch <- eq$nombre[k]
  i  <- carac_extremos[carac_extremos$accession == a & carac_extremos$extremo == "inicio", ]
  f  <- carac_extremos[carac_extremos$accession == a & carac_extremos$extremo == "final", ]
  fila <- data.frame(
    chrom = ch, accession = a, long_pb = L_sl[k], length_Mb = round(L_sl[k] / 1e6, 3),
    GC_percent = round(as.numeric(gc_chr[k]), 2),
    n_genes = n_gen_chr[k], genes_per_Mb = round(n_gen_chr[k] / (L_sl[k] / 1e6), 2),
    gene_cobertura_pct = round(100 * gene_bp_chr[k] / L_sl[k], 2),
    intergenic_pct = round(100 - 100 * gene_bp_chr[k] / L_sl[k], 2),
    repeat_pct = round(100 * rep_bp_chr[k] / L_sl[k], 2),
    telomeros = estado(i$tiene_telomero, f$tiene_telomero),
    telo_inicio_pb = i$telo_long_pb, telo_final_pb = f$telo_long_pb,
    bloques_telo_internos = sum(tel_df$seqnames == a & tel_df$extremo == "interno" & tel_df$es_telomero),
    subtelomeros = estado(i$ST_detectada, f$ST_detectada),
    ST_inicio_pb = i$ST_long_pb, ST_final_pb = f$ST_long_pb,
    ST_inicio_dens_genes_kb = i$densidad_ST_genes_kb, ST_final_dens_genes_kb = f$densidad_ST_genes_kb,
    stringsAsFactors = FALSE)
  for (sn in names(patron_sat)) {
    b <- sat_bloques[sat_bloques$chrom == ch & sat_bloques$satelite == sn, ]
    fila[[paste0(sn, "_bloques")]] <- nrow(b)
    fila[[paste0(sn, "_pb")]]      <- sum(b$longitud_pb)
  }
  fila
}))
cat("\n=== Resumen por cromosoma ===\n"); print(resumen_cromosomas)

# ---- 13.1 Correlación GC vs genes/Mb por cromosoma + FDR ----
cor_chr <- do.call(rbind, lapply(levels(tabla$chrom), function(ch) {
  s <- d[d$chrom == ch, ]
  ct <- if (nrow(s) >= 3 && length(unique(s$GC_percent)) > 1 && length(unique(s$genes_per_Mb)) > 1)
    tryCatch(cor.test(s$GC_percent, s$genes_per_Mb, method = "spearman", exact = FALSE), error = function(e) NULL) else NULL
  data.frame(chrom = ch, n_ventanas_usadas = nrow(s), n_excluidas = sum(tabla$chrom == ch) - nrow(s),
             rho = if (is.null(ct)) NA else round(unname(ct$estimate), 3),
             p_value = if (is.null(ct)) NA else signif(ct$p.value, 3))
}))
cor_chr$p_FDR <- NA
ix <- !is.na(cor_chr$p_value)
cor_chr$p_FDR[ix] <- signif(p.adjust(cor_chr$p_value[ix], "BH"), 3)

# ---- 13.2 Mantel-Haenszel ----
tabs <- lapply(levels(tabla$chrom), function(ch) {
  s <- d[d$chrom == ch, ]
  table(bajo_GC = factor(s$bajo_GC, levels = c(TRUE, FALSE)),
        baja_dens = factor(s$baja_dens_genica, levels = c(TRUE, FALSE)))
})
mh <- tryCatch(mantelhaen.test(array(unlist(tabs), dim = c(2, 2, length(tabs))), alternative = "greater"),
               error = function(e) NULL)
mh_res <- if (!is.null(mh)) data.frame(OR_MH = round(unname(mh$estimate), 2), p_value = signif(mh$p.value, 3),
                                       hipotesis = "bajo GC enriquecido en baja densidad génica") else
                                         data.frame(nota = "Mantel-Haenszel no calculable")

# ---- 13.3 Comparaciones de grupos ----
vars <- c(GC_percent = "GC (%)", repeat_pct = "% repeticiones", genes_per_Mb = "Genes/Mb",
          gene_pct = "% cobertura génica", exon_pct = "% exones", intron_pct = "% intrones",
          intron_in_gene = "% intrón dentro del gen", introns_por_gen = "Intrones por gen",
          intergenic_pct = "% intergénico", inter_rep_pct = "% intergénico repetitivo",
          promoter_pct = "% promotores", GC_exon = "GC de exones (%)", GC_intron = "GC de intrones (%)",
          GC_inter = "GC intergénico (%)")

comparar <- function(g1, g2, e1, e2, dat = d, vv = vars) {
  res <- do.call(rbind, lapply(names(vv), function(v) {
    x <- dat[[v]][g1]; y <- dat[[v]][g2]; x <- x[!is.na(x)]; y <- y[!is.na(y)]
    if (length(x) < 3 || length(y) < 3) return(NULL)
    w <- wilcox.test(x, y, exact = FALSE)
    data.frame(grupo_1 = e1, grupo_2 = e2, variable = vv[v], n_1 = length(x), n_2 = length(y),
               mediana_1 = round(median(x), 2), mediana_2 = round(median(y), 2),
               cliff_delta = round(2 * unname(w$statistic) / (length(x) * length(y)) - 1, 3),
               p_value = signif(w$p.value, 3), stringsAsFactors = FALSE)
  }))
  if (!is.null(res)) { res$p_FDR <- signif(p.adjust(res$p_value, "BH"), 3); rownames(res) <- NULL }
  res
}
cmp_bajoGC   <- comparar(d$bajo_GC, !d$bajo_GC, "Bajo GC", "Resto")
cmp_bajoDen  <- comparar(d$baja_dens_genica, !d$baja_dens_genica, "Baja densidad génica", "Resto")
cmp_intron   <- comparar(d$alto_intron, !d$alto_intron, "Alto % intrones", "Resto")
cmp_rep_AB   <- comparar(d$alto_repeat, d$bajo_repeat, "ALTA repetición", "BAJA repetición")   # pregunta: % intrones en alta vs baja repetición
cmp_rep_alto <- comparar(d$alto_repeat, !d$alto_repeat, "ALTA repetición", "Resto")
cmp_repeticion <- rbind(cmp_rep_AB, cmp_rep_alto)
cat("\nALTA vs BAJA repetición:\n"); print(cmp_rep_AB)

d$grupo_rep <- factor(ifelse(d$alto_repeat, "Alta repetición", ifelse(d$bajo_repeat, "Baja repetición", "Intermedia")),
                      levels = c("Alta repetición", "Intermedia", "Baja repetición"))
perfil_rep <- d %>% group_by(grupo_rep) %>%
  summarise(n_ventanas = n(), across(all_of(names(vars)), ~ round(median(.x, na.rm = TRUE), 2)), .groups = "drop") %>%
  as.data.frame()
perfil <- d %>% group_by(grupo) %>%
  summarise(across(all_of(names(vars)), ~ round(median(.x, na.rm = TRUE), 2)), .groups = "drop") %>% as.data.frame()

# ---- 13.4 Correlaciones ----
cor_vs <- function(ref, etiqueta) {
  do.call(rbind, lapply(setdiff(names(vars), ref), function(v) {
    ok <- !is.na(d[[v]]) & !is.na(d[[ref]])
    if (!is.numeric(d[[v]]) || !is.numeric(d[[ref]]) || sum(ok) < 10 ||
        length(unique(d[[v]][ok])) < 2 || length(unique(d[[ref]][ok])) < 2) {
      message("cor_vs: se omite '", v, "' vs '", ref, "' (sin datos numéricos o sin variación)")
      return(NULL)
    }
    glob <- cor.test(d[[ref]][ok], d[[v]][ok], method = "spearman", exact = FALSE)
    por_chr <- sapply(levels(d$chrom), function(ch) {
      s <- d[d$chrom == ch & ok, ]
      if (nrow(s) < 10 || length(unique(s[[v]])) < 2 || length(unique(s[[ref]])) < 2) return(NA_real_)
      unname(cor(s[[ref]], s[[v]], method = "spearman"))
    })
    data.frame(referencia = etiqueta, variable = vars[v], rho_global = round(unname(glob$estimate), 3),
               p_global = signif(glob$p.value, 3), rho_mediana_por_chr = round(median(por_chr, na.rm = TRUE), 3),
               chr_rho_neg = sum(por_chr < 0, na.rm = TRUE), chr_rho_pos = sum(por_chr > 0, na.rm = TRUE))
  }))
}
cor_tab <- rbind(cor_vs("GC_percent", "GC (%)"), cor_vs("genes_per_Mb", "Genes/Mb"), cor_vs("repeat_pct", "% repeticiones"))
cor_tab$p_FDR <- signif(p.adjust(cor_tab$p_global, "BH"), 3)
rownames(cor_tab) <- NULL

het_gc <- d %>% group_by(chrom) %>%
  summarise(n = n(), GC_mediana = median(GC_percent), GC_sd = sd(GC_percent),
            GC_min = min(GC_percent), GC_max = max(GC_percent), rango_IQR = IQR(GC_percent),
            .groups = "drop") %>% as.data.frame()

# ---- 13.5 Modelos ----
# Variables numéricas con variación (si no hay repeticiones/intrones se excluyen solas)
preds_validos <- function(dat, preds)
  preds[vapply(preds, function(v) is.numeric(dat[[v]]) && length(unique(dat[[v]][!is.na(dat[[v]])])) > 1, logical(1))]
formula_dinamica <- function(resp, preds, dat) {
  t <- paste0("scale(", preds, ")")
  if (nlevels(droplevels(dat$chrom)) > 1) t <- c(t, "chrom")
  as.formula(paste(resp, "~", paste(t, collapse = " + ")))
}
preds_or <- preds_validos(d, c("repeat_pct", "genes_per_Mb", "intron_pct", "intergenic_pct"))
mod <- glm(formula_dinamica("bajo_GC", preds_or, d), data = d, family = binomial)
ci_or <- exp(confint.default(mod))
or_tab <- data.frame(termino = names(coef(mod)), OR = exp(coef(mod)), IC_inf = ci_or[, 1], IC_sup = ci_or[, 2],
                     p_value = summary(mod)$coefficients[, 4])
or_tab <- or_tab[grepl("scale", or_tab$termino), ]
or_tab$termino <- sub("scale\\((.*)\\)", "\\1", or_tab$termino)
or_tab[, 2:4] <- round(or_tab[, 2:4], 2); or_tab$p_value <- signif(or_tab$p_value, 3)
rownames(or_tab) <- NULL

bg <- d[d$bajo_GC, ]
bg$perfil <- case_when(
  bg$alto_repeat & bg$baja_dens_genica ~ "Repetitiva y pobre en genes",
  bg$alto_repeat                       ~ "Repetitiva con genes",
  bg$baja_dens_genica                  ~ "Pobre en genes, no repetitiva",
  TRUE                                 ~ "Génica, no repetitiva")
res_perfil <- as.data.frame(table(bg$perfil)); colnames(res_perfil) <- c("perfil_bajoGC", "n_ventanas")
res_perfil$pct <- round(100 * res_perfil$n_ventanas / sum(res_perfil$n_ventanas), 1)

preds_gc <- preds_validos(d, c("intron_in_gene", "gene_pct", "genes_per_Mb", "repeat_pct"))
dm <- d[complete.cases(d[, c("GC_percent", preds_gc)]), ]
preds_gc <- preds_validos(dm, preds_gc)
mod_gc <- lm(formula_dinamica("GC_percent", preds_gc, dm), data = dm)
co <- summary(mod_gc)$coefficients
coef_tab <- data.frame(termino = sub("scale\\((.*)\\)", "\\1", rownames(co)),
                       efecto_GC_por_SD = round(co[, 1], 3), error = round(co[, 2], 3), p_value = signif(co[, 4], 3))
coef_tab <- coef_tab[!grepl("^chrom", coef_tab$termino) & coef_tab$termino != "(Intercept)", ]
coef_tab$R2_modelo <- round(summary(mod_gc)$r.squared, 3)
rownames(coef_tab) <- NULL

# ---- 13.6 Intergénico vs baja densidad génica ----
rho_ig <- cor.test(d$genes_per_Mb, d$intergenic_pct, method = "spearman", exact = FALSE)
tab_ig <- table(baja_dens_genica = d$baja_dens_genica, alto_intergenico = d$alto_intergenico)
fi_ig  <- fisher.test(tab_ig)
jacc   <- sum(d$baja_dens_genica & d$alto_intergenico) / sum(d$baja_dens_genica | d$alto_intergenico)
concordancia <- data.frame(rho = unname(rho_ig$estimate), p_rho = rho_ig$p.value,
                           OR_fisher = unname(fi_ig$estimate), p_fisher = fi_ig$p.value, jaccard = jacc)
