# ============================================================
# 04_ventanas.R  -  Ventanas: composición, coberturas (sobre pb no-N) y conteos
# ------------------------------------------------------------
# REQUIERE   : 01 (bin, w_excluir_ext, ext_*, min_bp_*), 02 (sl, eq, fa, nombre_chr, restar,
#              N_red, rep_red, genes_red, genes_gr, exon_red, intrones_red, prom, prom_red),
#              03 (inter_gr)
# CREA       : chrs, bins            GRanges de cromosomas y de ventanas
#              tabla                 TABLA MAESTRA (una fila por ventana). Columnas:
#                chrom_acc, chrom, start, end, width, *_Mb, pos_rel, region, dist_extremo_pb,
#                lado_extremo, ext_tecnico, dist_clase, GC_percent, AT_percent, N_pct, width_noN,
#                repeat_pct, gene_pct, intergenic_pct, exon_pct, intron_pct, promoter_pct,
#                inter_rep_pct, has_intron, has_promoter, n_genes, genes_per_Mb, n_introns,
#                introns_por_gen, n_promoters, intron_in_gene, GC_intron, GC_exon, GC_inter
#              wn                    pb no-N por ventana (vector)
#              L_win, centro, cortes, etq_d, mid_d   (distancia al extremo y clases de distancia)
#              seqs, genes_mid, win_gene, gene_bp, intron_bp
# CREA (funciones):
#   cobertura_bp(bins_, gr_reducido)  pb de cada ventana cubiertos por un GRanges
#   cob_bp(gr) / cob_pct(gr)          cobertura en pb / % sobre pb no-N (resta N_red)
#   agg_win(valor)                    suma un valor por gen dentro de cada ventana
#   gc_de(gr)                         %GC de cada rango
#   gc_ventana(gr, min_bp)            %GC por ventana restringido a un GRanges
# USADO POR  : 05 al 17 (es el módulo central; casi todo cuelga de `tabla`).
# ============================================================

# ============================================================
# 5. VENTANAS: COMPOSICIÓN, COBERTURAS (sobre pb NO-N) Y CONTEOS
# ============================================================

chrs <- GRanges(names(sl), IRanges(1, sl), seqlengths = sl)
bins <- unlist(tile(chrs, width = bin))
seqlevels(bins) <- names(sl)
strand(bins) <- "*"
names(bins) <- NULL

tabla <- data.frame(
  chrom_acc = as.character(seqnames(bins)),
  chrom     = factor(nombre_chr(as.character(seqnames(bins))), levels = eq$nombre),
  start = start(bins), end = end(bins), width = width(bins),
  stringsAsFactors = FALSE)
tabla$start_Mb    <- tabla$start / 1e6
tabla$end_Mb      <- tabla$end / 1e6
tabla$position_Mb <- (tabla$start + tabla$end) / 2 / 1e6
tabla$pos_rel     <- tabla$position_Mb * 1e6 / as.numeric(sl[tabla$chrom_acc])
tabla$region      <- ifelse(tabla$pos_rel < 0.10 | tabla$pos_rel > 0.90, "Extremos (10%)", "Interior")

# Distancia de cada ventana (centro) al extremo cromosómico más cercano
L_win  <- as.numeric(sl[tabla$chrom_acc])
centro <- (tabla$start + tabla$end) / 2
tabla$dist_extremo_pb <- pmin(centro, L_win - centro)
tabla$lado_extremo    <- ifelse(centro <= L_win / 2, "inicio", "final")
tabla$ext_tecnico     <- tabla$start <= w_excluir_ext | tabla$end >= L_win - w_excluir_ext + 1
cortes <- c(seq(0, ext_max_kb, by = ext_paso_kb), Inf)
etq_d  <- c(paste0(head(cortes, -2), "-", cortes[-c(1, length(cortes))], " kb"), paste0("≥", ext_max_kb, " kb"))
mid_d  <- c(head(cortes, -2) + ext_paso_kb / 2, ext_max_kb + ext_paso_kb / 2)
tabla$dist_clase <- cut(tabla$dist_extremo_pb / 1000, breaks = cortes, labels = etq_d, right = FALSE)

# ---- 5.1 %GC, %AT (sobre A+C+G+T) y N ----
seqs    <- getSeq(fa, bins)
lf_base <- letterFrequency(seqs, c("A", "C", "G", "T"))
lf_gc   <- letterFrequency(seqs, "GC")[, 1]
n_base  <- rowSums(lf_base)
tabla$GC_percent <- ifelse(n_base > 0, 100 * lf_gc / n_base, NA_real_)
tabla$AT_percent <- ifelse(n_base > 0, 100 * (lf_base[, "A"] + lf_base[, "T"]) / n_base, NA_real_)

# ---- 5.2 Cobertura ----
cobertura_bp <- function(bins_, gr_reducido) {
  ov   <- findOverlaps(bins_, gr_reducido, ignore.strand = TRUE)
  out  <- numeric(length(bins_))
  if (length(ov) > 0) {
    ints <- pintersect(bins_[queryHits(ov)], gr_reducido[subjectHits(ov)], ignore.strand = TRUE)
    s <- tapply(width(ints), queryHits(ov), sum)
    out[as.integer(names(s))] <- s
  }
  out
}
N_bp <- cobertura_bp(bins, N_red)
tabla$N_pct     <- 100 * N_bp / tabla$width
tabla$width_noN <- tabla$width - N_bp
wn <- tabla$width_noN
cob_bp  <- function(gr) cobertura_bp(bins, restar(gr, N_red))
cob_pct <- function(gr) ifelse(wn > 0, 100 * cob_bp(gr) / wn, NA_real_)   # % de pb NO-N

tabla$repeat_pct     <- cob_pct(rep_red)
tabla$gene_pct       <- cob_pct(genes_red)
tabla$intergenic_pct <- cob_pct(inter_gr)                 # gene_pct + intergenic_pct = 100 (sobre pb no-N)
tabla$exon_pct       <- cob_pct(exon_red)
tabla$intron_pct     <- cob_pct(intrones_red)
tabla$promoter_pct   <- cob_pct(prom_red)
tabla$inter_rep_pct  <- cob_pct(intersect(inter_gr, rep_red, ignore.strand = TRUE))
tabla$has_intron     <- ifelse(tabla$intron_pct > 0, "Si", "No")
tabla$has_promoter   <- ifelse(tabla$promoter_pct > 0, "Si", "No")

# ---- 5.3 Conteos por punto medio ----
genes_mid <- resize(genes_gr, width = 1, fix = "center")
win_gene  <- findOverlaps(genes_mid, bins, select = "first", ignore.strand = TRUE)
ov_ig <- findOverlaps(intrones_red, genes_gr, type = "within", ignore.strand = TRUE)
genes_gr$n_introns <- tabulate(subjectHits(ov_ig), nbins = length(genes_gr))

agg_win <- function(valor) {
  out <- numeric(length(bins))
  s <- tapply(valor, win_gene, sum)
  out[as.integer(names(s))] <- s
  out
}
tabla$n_genes         <- agg_win(rep(1, length(genes_gr)))
tabla$genes_per_Mb    <- ifelse(wn > 0, tabla$n_genes / (wn / 1e6), NA_real_)
tabla$n_introns       <- agg_win(genes_gr$n_introns)
tabla$introns_por_gen <- ifelse(tabla$n_genes > 0, tabla$n_introns / tabla$n_genes, NA)
tabla$n_promoters     <- countOverlaps(bins, resize(prom, 1, fix = "center"), ignore.strand = TRUE)

gene_bp   <- cob_bp(genes_red)
intron_bp <- cob_bp(intrones_red)
tabla$intron_in_gene <- ifelse(gene_bp >= min_bp_gen_ventana, 100 * intron_bp / gene_bp, NA)

# ---- 5.4 GC por componente ----
gc_de <- function(gr) {
  s <- getSeq(fa, gr)
  100 * letterFrequency(s, "GC")[, 1] / rowSums(letterFrequency(s, c("A", "C", "G", "T")))
}
gc_ventana <- function(gr, min_bp = min_bp_gc_comp) {
  out <- rep(NA_real_, length(bins))
  ov  <- findOverlaps(bins, gr, ignore.strand = TRUE)
  if (length(ov) == 0) return(out)
  ints <- pintersect(bins[queryHits(ov)], gr[subjectHits(ov)], ignore.strand = TRUE)
  s    <- getSeq(fa, ints)
  gc   <- rowsum(letterFrequency(s, "GC")[, 1], queryHits(ov))
  at   <- rowsum(rowSums(letterFrequency(s, c("A", "C", "G", "T"))), queryHits(ov))
  out[as.integer(rownames(gc))] <- ifelse(at[, 1] >= min_bp, 100 * gc[, 1] / at[, 1], NA)
  out
}
tabla$GC_intron <- gc_ventana(intrones_red)
tabla$GC_exon   <- gc_ventana(exon_red)
tabla$GC_inter  <- gc_ventana(inter_gr)
