# ============================================================
# 06_subtelomeros.R  -  Subtélomeros por extremo + tabla de telómeros por cromosoma
# ------------------------------------------------------------
# REQUIERE   : 01 (w_st_max, gap_st, k_gaps_st, min_st, dens_st_max), 02 (sl, eq, genes_gr,
#              rep_red, armonizar, nombre_chr), 03, 04 (bins, tabla, gc_de, cobertura_bp),
#              05 (tel_df, res_telo)
# CREA       : genes_df                      genes como data.frame simple
#              subtelomero(acc, ext)         función: ST de un extremo (longitud + densidad génica)
#              carac_extremos                UNA FILA POR EXTREMO: telómero (res_telo) + ST + GC_ST,
#                                            repeat_pct_ST
#              ok_len, st_all, st_gr         st_all = todos los tramos; st_gr = solo ST válidos
#              ST_genes                      genes dentro de cada tramo ST
#              estado(a, b)                  "Ambos extremos"/"Solo inicio"/"Solo final"/"Ninguno"
#              ci, cf, telo_por_cromosoma    telómero inicio/final por cromosoma
#              tabla$subtelomero
# USADO POR  : 07_satelites (st_gr), 08 (st_gr), 12 (carac_extremos, estado), 13 (st_gr),
#              14_excel (carac_extremos, ST_genes, telo_por_cromosoma), 15_figuras (st_gr)
# ============================================================

# ============================================================
# 7. SUBTÉLOMEROS
#    Tramo desde el borde hasta que aparecen `k_gaps_st` intervalos consecutivos <= gap_st.
#    ST válido ("Si") si: longitud >= min_st  Y  densidad génica < dens_st_max (genes/kb).
# ============================================================

genes_df <- data.frame(acc = as.character(seqnames(genes_gr)), s = start(genes_gr),
                       e = end(genes_gr), id = genes_gr$gene_id, stringsAsFactors = FALSE)

subtelomero <- function(acc, ext) {
  L <- as.numeric(sl[acc])
  g <- genes_df[genes_df$acc == acc, ]
  if (ext == "inicio") { s <- g$s; e <- g$e } else { s <- L - g$e + 1; e <- L - g$s + 1 }
  id <- g$id
  
  d200  <- sum(s <= 200000) / 200
  d300  <- sum(s <= 300000) / 300
  d_chr <- nrow(g) / (L / 1000)
  
  keep <- s <= w_st_max
  s <- s[keep]; e <- e[keep]; id <- id[keep]
  o <- order(s); s <- s[o]; e <- e[o]; id <- id[o]
  n <- length(s)
  
  idx <- NA_integer_; gaps <- numeric(0)
  if (n >= 2) {
    gaps  <- s[-1] - cummax(e)[-n] - 1
    denso <- gaps <= gap_st
    if (length(denso) >= k_gaps_st)
      for (i in seq_len(length(denso) - k_gaps_st + 1))
        if (all(denso[i:(i + k_gaps_st - 1)])) { idx <- i; break }
  }
  
  if (is.na(idx)) { st_len <- w_st_max; n_st <- n; gaps_in <- gaps }
  else            { st_len <- s[idx] - 1; n_st <- idx - 1; gaps_in <- gaps[seq_len(idx - 1)] }
  
  if (!is.na(idx) && idx == 1 && incluir_primer_gen) {
    st_len  <- e[1]
    n_st    <- 1
    gaps_in <- numeric(0)
  }
  
  st_len <- max(0, min(st_len, L))
  
  
  
  if (st_len > 0) {
    st_start <- if (ext == "inicio") 1 else L - st_len + 1
    st_end   <- if (ext == "inicio") st_len else L
  } else { st_start <- NA; st_end <- NA }
  
  dens_st   <- if (st_len > 0) n_st / (st_len / 1000) else NA_real_
  crit_long  <- st_len >= min_st
  dens_aplica <- st_len >= dens_min_len
  crit_dens  <- if (!dens_aplica) TRUE else (!is.na(dens_st) && dens_st < dens_st_max)
  motivo_no <- if (crit_long && crit_dens) "" else
    paste(c(if (!crit_long) paste0("longitud < ", min_st, " pb"),
            if (!crit_dens) paste0("densidad >= ", dens_st_max, " genes/kb")), collapse = " + ")
  
  data.frame(
    ST_detectada = ifelse(crit_long && crit_dens, "Si", "No"),
    ST_crit_longitud = ifelse(crit_long, "Si", "No"),
    ST_crit_densidad = ifelse(!dens_aplica, "No aplica", ifelse(crit_dens, "Si", "No")),    ST_motivo_no = motivo_no,
    ST_zona_densa_alcanzada = ifelse(is.na(idx), "No (truncado a w_st_max)", "Si"),
    ST_start = as.integer(st_start), ST_end = as.integer(st_end), ST_long_pb = as.integer(st_len),
    n_genes_ST = as.integer(n_st),
    densidad_ST_genes_kb = round(dens_st, 4),
    dist_borde_primer_gen_pb = ifelse(n > 0, s[1] - 1, NA),
    gap_max_pb   = if (length(gaps_in)) max(gaps_in) else NA,
    gap_medio_pb = if (length(gaps_in)) round(mean(gaps_in)) else NA,
    n_gaps_mayores_umbral = sum(gaps_in > gap_st),
    dens_200kb_genes_kb = round(d200, 4),
    dens_300kb_genes_kb = round(d300, 4),
    dens_cromosoma_genes_kb = round(d_chr, 4),
    genes_ST = if (n_st > 0) paste(id[seq_len(n_st)], collapse = "; ") else "",
    stringsAsFactors = FALSE)
}






carac_extremos <- do.call(rbind, lapply(names(sl), function(acc)
  do.call(rbind, lapply(c("inicio", "final"), function(ext)
    cbind(data.frame(chrom = nombre_chr(acc), accession = acc, extremo = ext,
                     long_cromosoma_pb = as.numeric(sl[acc]), stringsAsFactors = FALSE),
          res_telo(acc, ext), subtelomero(acc, ext))))))

carac_extremos$GC_ST <- NA_real_
carac_extremos$repeat_pct_ST <- NA_real_
ok_len <- !is.na(carac_extremos$ST_start) & carac_extremos$ST_long_pb > 0
st_all <- armonizar(GRanges(carac_extremos$accession[ok_len],
                            IRanges(carac_extremos$ST_start[ok_len], carac_extremos$ST_end[ok_len])))
if (length(st_all) > 0) {
  carac_extremos$GC_ST[ok_len]         <- round(gc_de(st_all), 2)
  carac_extremos$repeat_pct_ST[ok_len] <- round(100 * cobertura_bp(st_all, rep_red) / width(st_all), 2)
}
st_gr <- st_all[carac_extremos$ST_detectada[ok_len] == "Si"]   # solo ST que cumplen longitud Y densidad
tabla$subtelomero <- "No"
if (length(st_gr) > 0) tabla$subtelomero[unique(queryHits(findOverlaps(bins, st_gr)))] <- "Si"

# Genes dentro de cada tramo subtelomérico (con descripción: familias génicas de interés)
ov_stg <- findOverlaps(st_all, genes_gr, ignore.strand = TRUE)
ce_ok  <- carac_extremos[ok_len, ]
ST_genes <- data.frame(
  chrom = ce_ok$chrom[queryHits(ov_stg)], extremo = ce_ok$extremo[queryHits(ov_stg)],
  ST_detectada = ce_ok$ST_detectada[queryHits(ov_stg)],
  gene_id = genes_gr$gene_id[subjectHits(ov_stg)], gene_nombre = genes_gr$gene_nombre[subjectHits(ov_stg)],
  descripcion = genes_gr$descripcion[subjectHits(ov_stg)], biotipo = genes_gr$biotipo[subjectHits(ov_stg)],
  start = start(genes_gr)[subjectHits(ov_stg)], end = end(genes_gr)[subjectHits(ov_stg)],
  strand = as.character(strand(genes_gr))[subjectHits(ov_stg)], stringsAsFactors = FALSE)

cat("\n=== Telómeros y subtélomeros por extremo ===\n")
print(carac_extremos[, c("chrom", "extremo", "tiene_telomero", "telo_long_pb", "ST_detectada", "ST_motivo_no",
                         "ST_long_pb", "n_genes_ST", "densidad_ST_genes_kb")])


# ---- Telómeros por cromosoma: ¿en qué extremo hay? ----
estado <- function(a, b)
  ifelse(a == "Si" & b == "Si", "Ambos extremos",
         ifelse(a == "Si", "Solo inicio", ifelse(b == "Si", "Solo final", "Ninguno")))

ci <- carac_extremos[carac_extremos$extremo == "inicio", ]
cf <- carac_extremos[carac_extremos$extremo == "final",  ]
ci <- ci[match(eq$acc, ci$accession), ]
cf <- cf[match(eq$acc, cf$accession), ]

telo_por_cromosoma <- data.frame(
  chrom             = eq$nombre,
  accession         = eq$acc,
  long_cromosoma_pb = as.numeric(sl),
  telo_inicio       = ci$tiene_telomero,
  inicio_start      = ci$telo_start,
  inicio_end        = ci$telo_end,
  inicio_long_pb    = ci$telo_long_pb,
  inicio_copias     = ci$telo_copias,
  telo_final        = cf$tiene_telomero,
  final_start       = cf$telo_start,
  final_end         = cf$telo_end,
  final_long_pb     = cf$telo_long_pb,
  final_copias      = cf$telo_copias,
  estado            = estado(ci$tiene_telomero, cf$tiene_telomero),
  stringsAsFactors  = FALSE)

cat("\n=== Telómeros por cromosoma ===\n")
print(telo_por_cromosoma[, c("chrom", "telo_inicio", "telo_final", "estado",
                             "inicio_long_pb", "final_long_pb")])
cat("Con telómero en inicio:", sum(telo_por_cromosoma$telo_inicio == "Si"),
    "| en final:", sum(telo_por_cromosoma$telo_final == "Si"),
    "de", nrow(telo_por_cromosoma), "cromosomas\n")

### Validacion de la longitud del ST contra el articulo (solo T. gondii; crea `paper` y `grid` para F15)
if (isTRUE(validar_st_articulo) && any(carac_extremos$chrom %in% c("Ia", "Ib", "II"))) {
###Prueba para saber variables de ST 
paper <- data.frame(
  chrom   = c("Ia","Ia","Ib","Ib","II","II","III","III","IV","IV","V","V","VI","VI",
              "VIIa","VIIa","VIIb","VIII","IX","IX","X","X","XI","XI","XII","XII"),
  extremo = c("inicio","final","inicio","final","inicio","final","inicio","final","inicio","final",
              "inicio","final","inicio","final","inicio","final","inicio","final","inicio","final",
              "inicio","final","inicio","final","inicio","final"),
  len_art = c(15.6,10.5,18.4,40,22.2,109.4,35.7,232.4,77,40.2,112.6,177.9,8.1,53.7,
              58.4,76.3,82.9,53.1,179.1,49.6,139.8,27.5,88.2,13.2,47.6,215) * 1000)
pm <- merge(paper, carac_extremos[, c("chrom","extremo","accession")], by = c("chrom","extremo"))
if (nrow(pm) > 0) {

g0 <- gap_st; k0 <- k_gaps_st            # guardar originales
grid <- expand.grid(k = 2:5, gap = c(2000, 2500, 3000, 3500, 4500, 6000))
grid$n_ok <- NA; grid$err_med <- NA
for (i in seq_len(nrow(grid))) {
  k_gaps_st <- grid$k[i]; gap_st <- grid$gap[i]
  mine <- mapply(function(a, e) subtelomero(a, e)$ST_long_pb, pm$accession, pm$extremo)
  r <- mine / pm$len_art
  grid$n_ok[i]    <- sum(r > 0.75 & r < 1.25)
  grid$err_med[i] <- round(median(abs(log(r))), 3)
}
gap_st <- g0; k_gaps_st <- k0            # restaurar
}
}
