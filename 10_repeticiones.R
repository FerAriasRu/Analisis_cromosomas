# ============================================================
# 10_repeticiones.R  -  Familias, clases y ubicación (exón/intrón/intergénico) de las repeticiones
# ------------------------------------------------------------
# REQUIERE   : 01 (top_n_rep, max_seq_chars, ext_umbrales_kb), 02 (fa, sl, eq, tiene_rep, rep_gr, rep_red,
#              genes_red, exon_red, intrones_red, N_red, nombre_chr, restar),
#              03 (inter_gr), 04 (chrs, bins, cobertura_bp), 09 (banderas alto/bajo_repeat, base, analisis)
# CREA       : L_sl, rep_bp_chr, gene_bp_chr, rep_total_bp
#              por_nombre, bp_nombre, fam_nombre, div_nombre, nchr_nombre, top_nm
#              rep_top, rep_familias
#              clases, red_cl, rep_clase_bp, rep_clases_pb, rep_clases_pct, rep_clases_resumen
#              comps (Exón/Intrón/Génico otro/Intergénico), tabla_comp(), reg_vent(), ext_ref_kb, subconj
#              rep_comp, rep_comp_chr, rep_comp_clase, rep_comp_top
#              rep_hits, rep_hits_res
# USADO POR  : 11 (ext_ref_kb), 12 (L_sl, rep_bp_chr, gene_bp_chr), 14_excel (todas las rep_*),
#              15_figuras (clases, rep_clase_bp, L_sl)
# SIN REPETICIONES (tiene_rep = FALSE): crea L_sl, rep_bp_chr, gene_bp_chr, ext_ref_kb y deja las
#              tablas rep_* en NULL (14_excel las omite, 15_figuras omite F y G).
# ============================================================

# ============================================================
# 11. REPETICIONES: FAMILIAS, CLASES Y UBICACIÓN (exón / intrón / intergénico)
# ============================================================

L_sl <- as.numeric(sl)
rep_bp_chr  <- vapply(names(sl), function(a) sum(width(rep_red[as.character(seqnames(rep_red)) == a])), numeric(1))
gene_bp_chr <- vapply(names(sl), function(a) sum(width(genes_red[as.character(seqnames(genes_red)) == a])), numeric(1))
rep_total_bp <- sum(width(rep_red))

ext_ref_kb <- max(ext_umbrales_kb)

if (tiene_rep) {

# ---- 11.1 Top de repeticiones (por nombre) con la secuencia de su copia más larga ----
por_nombre <- split(rep_gr, rep_gr$nombre)
bp_nombre  <- sum(width(reduce(por_nombre, ignore.strand = TRUE)))
bp_nombre  <- setNames(as.numeric(bp_nombre), names(por_nombre))
fam_nombre <- tapply(rep_gr$familia, rep_gr$nombre, function(x) names(sort(table(x), decreasing = TRUE))[1])
fam_nombre <- fam_nombre[names(por_nombre)]
div_nombre <- tapply(rep_gr$div, rep_gr$nombre, median, na.rm = TRUE)[names(por_nombre)]
nchr_nombre <- vapply(por_nombre, function(x) length(unique(as.character(seqnames(x)))), integer(1))

top_nm <- names(sort(bp_nombre, decreasing = TRUE))[seq_len(min(top_n_rep, length(bp_nombre)))]
rep_top <- do.call(rbind, lapply(seq_along(top_nm), function(k) {
  nm <- top_nm[k]; g <- por_nombre[[nm]]
  gi <- g[which.max(width(g))]
  s  <- subseq(fa[[as.character(seqnames(gi))]], start(gi), end(gi))
  if (as.character(strand(gi)) == "-") s <- reverseComplement(s)
  s <- as.character(s)
  data.frame(rank = k, repeat_name = nm, clase_familia = unname(fam_nombre[nm]),
             clase = sub("\\?$", "", sub("/.*", "", unname(fam_nombre[nm]))),
             n_copias = length(g), n_cromosomas = unname(nchr_nombre[nm]),
             bp_total_fusionado = bp_nombre[nm], pct_genoma = round(100 * bp_nombre[nm] / sum(L_sl), 3),
             pct_de_las_repeticiones = round(100 * bp_nombre[nm] / rep_total_bp, 2),
             divergencia_mediana_pct = round(unname(div_nombre[nm]), 2),
             copia_mas_larga_pb = width(gi),
             copia_mas_larga_coord = paste0(nombre_chr(as.character(seqnames(gi))), ":", start(gi), "-", end(gi),
                                            " (", ifelse(as.character(strand(gi)) == "-", "-", "+"), ")"),
             secuencia_truncada = ifelse(nchar(s) > max_seq_chars, "Si", "No"),
             secuencia = substr(s, 1, max_seq_chars),
             stringsAsFactors = FALSE)
}))
rownames(rep_top) <- NULL
cat("\n=== Top de repeticiones (por bp fusionados) ===\n")
print(rep_top[, c("rank", "repeat_name", "clase_familia", "n_copias", "bp_total_fusionado", "pct_genoma")])

# ---- 11.2 Resumen por clase/familia (todas) ----
por_fam <- split(rep_gr, rep_gr$familia)
bp_fam  <- setNames(as.numeric(sum(width(reduce(por_fam, ignore.strand = TRUE)))), names(por_fam))
rep_familias <- data.frame(
  clase_familia = names(por_fam), clase = sub("\\?$", "", sub("/.*", "", names(por_fam))),
  n_copias = as.integer(lengths(por_fam)), bp_fusionado = bp_fam,
  pct_genoma = round(100 * bp_fam / sum(L_sl), 3),
  pct_de_las_repeticiones = round(100 * bp_fam / rep_total_bp, 2),
  n_nombres = vapply(names(por_fam), function(f) sum(fam_nombre == f, na.rm = TRUE), integer(1)),
  top3_repeticiones = vapply(names(por_fam), function(f) {
    nm <- names(fam_nombre)[fam_nombre == f & !is.na(fam_nombre)]
    paste(head(nm[order(-bp_nombre[nm])], 3), collapse = "; ")
  }, character(1)), stringsAsFactors = FALSE)
rep_familias <- rep_familias[order(-rep_familias$bp_fusionado), ]; rownames(rep_familias) <- NULL

# ---- 11.3 % por clase (LTR, LINE, SINE, DNA, Satellite, Simple_repeat...) por cromosoma ----
clases <- sort(unique(rep_gr$clase))
red_cl <- lapply(setNames(clases, clases), function(cl) reduce(rep_gr[rep_gr$clase == cl], ignore.strand = TRUE))
rep_clase_bp <- sapply(red_cl, function(r) {
  s <- tapply(width(r), factor(as.character(seqnames(r)), levels = names(sl)), sum)
  s[is.na(s)] <- 0
  as.numeric(s)
})
if (is.null(dim(rep_clase_bp))) rep_clase_bp <- matrix(rep_clase_bp, ncol = length(clases), dimnames = list(NULL, clases))
orden_cl <- order(-colSums(rep_clase_bp)); rep_clase_bp <- rep_clase_bp[, orden_cl, drop = FALSE]
clases <- colnames(rep_clase_bp)

rep_clases_pb <- data.frame(chrom = c(eq$nombre, "Genoma"), long_pb = c(L_sl, sum(L_sl)),
                            rbind(rep_clase_bp, colSums(rep_clase_bp)),
                            Total_repeticion_pb = c(rep_bp_chr, sum(rep_bp_chr)),
                            check.names = FALSE, stringsAsFactors = FALSE)
rep_clases_pct <- rep_clases_pb
rep_clases_pct[, c(clases, "Total_repeticion_pb")] <-
  round(100 * rep_clases_pb[, c(clases, "Total_repeticion_pb")] / rep_clases_pb$long_pb, 2)
names(rep_clases_pct)[names(rep_clases_pct) == "Total_repeticion_pb"] <- "Total_repeticion_pct"
rep_clases_pct$long_pb <- NULL

rep_clases_resumen <- data.frame(
  clase = clases, n_copias = as.integer(table(factor(rep_gr$clase, levels = clases))),
  bp_fusionado = colSums(rep_clase_bp), pct_genoma = round(100 * colSums(rep_clase_bp) / sum(L_sl), 3),
  pct_de_las_repeticiones = round(100 * colSums(rep_clase_bp) / rep_total_bp, 2),
  stringsAsFactors = FALSE)
rownames(rep_clases_resumen) <- NULL
cat("\n% de cada cromosoma ocupado por cada clase de repetición:\n"); print(rep_clases_pct)

# ---- 11.4 ¿En qué compartimento caen las repeticiones? (exclusivos: Exón > Intrón > Génico otro / Intergénico) ----
comp_exon   <- restar(exon_red, N_red)
comp_intron <- restar(restar(intrones_red, exon_red), N_red)
ex_int      <- reduce(c(comp_exon, comp_intron), ignore.strand = TRUE)
comp_otro   <- restar(restar(genes_red, ex_int), N_red)
comp_inter  <- restar(restar(inter_gr, ex_int), N_red)
comps <- list("Exón" = comp_exon, "Intrón" = comp_intron,
              "Génico otro (UTR/no clasificado)" = comp_otro, "Intergénico" = comp_inter)

tabla_comp <- function(region, etiqueta, rep_obj = rep_red) {
  out <- do.call(rbind, lapply(names(comps), function(cn) {
    cr <- intersect(comps[[cn]], region, ignore.strand = TRUE)
    data.frame(subconjunto = etiqueta, compartimento = cn, bp_compartimento = sum(width(cr)),
               bp_repeticion = sum(width(intersect(cr, rep_obj, ignore.strand = TRUE))),
               stringsAsFactors = FALSE)
  }))
  out$pct_del_compartimento_que_es_repeticion <- round(100 * out$bp_repeticion / out$bp_compartimento, 2)
  out$pct_de_las_repeticiones <- round(100 * out$bp_repeticion / sum(out$bp_repeticion), 2)
  out$pct_del_subconjunto_compartimento <- round(100 * out$bp_compartimento / sum(out$bp_compartimento), 2)
  out$enriquecimiento <- round(out$pct_de_las_repeticiones / out$pct_del_subconjunto_compartimento, 2)
  out
}
reg_vent <- function(flag) if (any(flag)) reduce(bins[flag], ignore.strand = TRUE) else NULL
ext_ref_kb <- max(ext_umbrales_kb)
subconj <- list()
subconj[["Genoma completo"]]                 <- reduce(chrs, ignore.strand = TRUE)
subconj[["Ventanas ALTO % repetición"]]      <- reg_vent(tabla$alto_repeat)
subconj[["Ventanas BAJO % repetición"]]      <- reg_vent(tabla$bajo_repeat)
subconj[["Ventanas intermedias (análisis)"]] <- reg_vent(tabla$analisis & !tabla$alto_repeat & !tabla$bajo_repeat)
subconj[[paste0("Extremos (<= ", ext_ref_kb, " kb)")]] <- reg_vent(tabla$base & tabla$dist_extremo_pb <= ext_ref_kb * 1000)
subconj[[paste0("Interior (> ", ext_ref_kb, " kb)")]]  <- reg_vent(tabla$base & tabla$dist_extremo_pb > ext_ref_kb * 1000)
subconj <- Filter(Negate(is.null), subconj)
stopifnot(length(subconj) > 0, all(nchar(names(subconj)) > 0))
rep_comp <- do.call(rbind, lapply(names(subconj), function(n) tabla_comp(subconj[[n]], n)))
rownames(rep_comp) <- NULL

rep_comp_chr <- do.call(rbind, lapply(seq_along(sl), function(k)
  tabla_comp(reduce(chrs[k], ignore.strand = TRUE), eq$nombre[k])))
rownames(rep_comp_chr) <- NULL

rep_comp_clase <- do.call(rbind, lapply(clases, function(cl) {
  r <- red_cl[[cl]]
  bp <- vapply(comps, function(cm) sum(width(intersect(r, cm, ignore.strand = TRUE))), numeric(1))
  data.frame(clase = cl, compartimento = names(comps), bp = bp,
             pct_de_la_clase = round(100 * bp / sum(bp), 2), stringsAsFactors = FALSE)
}))
rownames(rep_comp_clase) <- NULL

rep_comp_top <- do.call(rbind, lapply(top_nm, function(nm) {
  r <- reduce(por_nombre[[nm]], ignore.strand = TRUE)
  bp <- vapply(comps, function(cm) sum(width(intersect(r, cm, ignore.strand = TRUE))), numeric(1))
  data.frame(repeat_name = nm, clase_familia = unname(fam_nombre[nm]), compartimento = names(comps), bp = bp,
             pct_de_la_repeticion = round(100 * bp / sum(bp), 2), stringsAsFactors = FALSE)
}))
rownames(rep_comp_top) <- NULL
cat("\nUbicación de las repeticiones (genoma completo):\n")
print(rep_comp[rep_comp$subconjunto == "Genoma completo", ])

# ---- 11.5 Cada hit de RepeatMasker clasificado por coordenadas ----
bpm <- sapply(comps, function(cm) cobertura_bp(rep_gr, cm))
if (is.null(dim(bpm))) bpm <- matrix(bpm, nrow = length(rep_gr), dimnames = list(NULL, names(comps)))
rep_hits <- data.frame(chrom = nombre_chr(as.character(seqnames(rep_gr))), start = start(rep_gr), end = end(rep_gr),
                       strand = as.character(strand(rep_gr)), repeat_name = rep_gr$nombre,
                       clase_familia = rep_gr$familia, clase = rep_gr$clase, divergencia_pct = rep_gr$div,
                       pb_exon = bpm[, 1], pb_intron = bpm[, 2], pb_genico_otro = bpm[, 3], pb_intergenico = bpm[, 4],
                       compartimento_dominante = ifelse(rowSums(bpm) == 0, "Sin clasificar",
                                                        colnames(bpm)[max.col(bpm, ties.method = "first")]),
                       stringsAsFactors = FALSE)
rep_hits_res <- as.data.frame.matrix(table(factor(rep_hits$clase, levels = clases), rep_hits$compartimento_dominante))
rep_hits_res <- cbind(clase = rownames(rep_hits_res), rep_hits_res, stringsAsFactors = FALSE)
rownames(rep_hits_res) <- NULL


} else {
  message("10_repeticiones: sin repeticiones, se omite el análisis de familias/clases/ubicación.")
  rep_top <- rep_familias <- rep_clases_resumen <- rep_clases_pct <- rep_clases_pb <- NULL
  rep_comp <- rep_comp_chr <- rep_comp_clase <- rep_comp_top <- rep_hits <- rep_hits_res <- NULL
  clases <- character(0); rep_clase_bp <- NULL
}
