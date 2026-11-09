# ============================================================
# 07_satelites.R  -  sat240 / sat350: RepeatMasker, TRF, de novo y secuencias GenBank
# ------------------------------------------------------------
# REQUIERE   : 01 (patron_sat, periodos_sat, gap_sat, trf_file, trf_min_copias, buscar_denovo, dn_*,
#              gb_*, usar_blast, w_st_max), 02 (fa, sl, eq, rm_df, fna_file, nombre_chr, armonizar),
#              04 (bins, tabla), 05 (motivos, tel_gr), 06 (st_gr)
#              Opcional: blastn/makeblastdb instalados; si no, usa k-mers en R.
# CREA (funciones): leer_trf(path), buscar_sat_denovo()   (GenBank/BLASTN de SAT350: ahora en 06b_sat350.R)
# CREA (datos): hits_rm, trf, trf_sat, hits_trf, hits_dn
#               (gb_hits, gb_mat, gb_presencia, metodo_gb, hits_gb los crea 06b_sat350.R)
#               hits                 todos los hits (acc, start, end, satelite, fuente)
#               sat_bloques          bloques fusionados por satélite, con ubicación y fuentes
#               trf_consensos        consensos TRF de sat240/sat350
#               tabla$sat240, tabla$sat350
# OPCIONAL   : sin trf_file se omiten los satélites por TRF; sin gb_accesiones se omite GenBank;
#              sin patron_sat/periodos_sat el módulo no busca satélites y crea tablas vacías.
#              (gb_mat y gb_presencia quedan NULL si no hay secuencias GenBank).
# USADO POR  : 12_estadistica (sat_bloques), 14_excel (sat_bloques, trf_consensos, gb_*),
#              15_figuras (gb_hits, gb_mat, metodo_gb)
# ============================================================

# ============================================================
# 8. SATÉLITES
# ============================================================


hits_rm <- do.call(rbind, lapply(names(patron_sat), function(sn) {
  x <- rm_df[rm_df$V5 %in% names(sl) & grepl(patron_sat[[sn]], rm_df$V10, ignore.case = TRUE), ]
  if (nrow(x) == 0) return(NULL)
  data.frame(acc = x$V5, start = x$V6, end = x$V7, satelite = sn, fuente = "RepeatMasker",
             stringsAsFactors = FALSE)
}))

leer_trf <- function(path) {
  if (is.null(path) || is.na(path) || !file.exists(path)) {
    message("TRF no disponible: se omiten los satélites por TRF."); return(NULL) }
  dat <- readLines(path); acc <- NA; filas <- list()
  for (l in dat) {
    if (startsWith(l, "Sequence:")) acc <- sub("Sequence: (\\S+).*", "\\1", l)
    else if (grepl("^[0-9]+ [0-9]+ ", l)) {
      f <- strsplit(l, " ")[[1]]
      filas[[length(filas) + 1]] <- data.frame(
        seqnames = acc, start = as.integer(f[1]), end = as.integer(f[2]),
        periodo = as.integer(f[3]), copias = as.numeric(f[4]),
        pct_match = as.numeric(f[6]), consenso = f[14], stringsAsFactors = FALSE)
    }
  }
  if (length(filas) == 0) return(NULL)
  do.call(rbind, filas)
}

trf <- leer_trf(trf_file); hits_trf <- NULL; trf_sat <- NULL
if (!is.null(trf)) {
  k_telo <- nchar(unidad_telo[1])
  trf$es_telo <- vapply(trf$consenso, function(s)
    !is.na(s) && nchar(s) %% k_telo == 0 && substr(s, 1, k_telo) %in% motivos &&
      s == strrep(substr(s, 1, k_telo), nchar(s) / k_telo), logical(1))
  trf$satelite <- NA_character_
  for (sn in names(periodos_sat)) trf$satelite[trf$periodo %in% periodos_sat[[sn]]] <- sn
  trf_sat <- trf[!is.na(trf$satelite) & !trf$es_telo & trf$copias >= trf_min_copias &
                   trf$seqnames %in% names(sl), ]
  if (nrow(trf_sat) > 0)
    hits_trf <- data.frame(acc = trf_sat$seqnames, start = trf_sat$start, end = trf_sat$end,
                           satelite = trf_sat$satelite, fuente = "TRF", stringsAsFactors = FALSE)
}

buscar_sat_denovo <- function(acc, periodos_sat, win = 60, min_ident = 0.60,
                              min_len = 500, gap_fusion = 200, max_micro = 0.80) {
  x <- charToRaw(toupper(as.character(fa[[acc]]))); n <- length(x)
  valido <- x != charToRaw("N")
  lags <- sort(unique(unlist(periodos_sat)))
  best <- numeric(n); bestp <- integer(n)
  
  # identidad local (ventana `win`) a distancia p; se queda el mejor p en cada posición
  for (p in lags) {
    m  <- (x[1:(n - p)] == x[(p + 1):n]) & valido[1:(n - p)] & valido[(p + 1):n]
    cs <- c(0L, cumsum(m))
    nw <- n - p - win + 1
    id <- (cs[(win + 1):(win + nw)] - cs[1:nw]) / win
    idx <- which(id > best[1:nw])
    best[idx] <- id[idx]; bestp[idx] <- p
  }
  
  ok <- best >= min_ident
  r  <- rle(ok); fin <- cumsum(r$lengths); ini <- fin - r$lengths + 1
  if (!any(r$values)) return(NULL)
  ir <- reduce(IRanges(ini[r$values], fin[r$values]), min.gapwidth = gap_fusion + 1)
  
  res <- lapply(seq_along(ir), function(j) {
    s0 <- start(ir)[j]; e0 <- end(ir)[j]
    bp <- bestp[s0:e0][ok[s0:e0]]
    p_mod <- as.integer(names(which.max(table(bp))))          # periodo dominante del bloque
    s <- s0; e <- min(n, e0 + win - 1 + p_mod)                # la ventana cubre win-1+p bases más allá
    if (e - s + 1 < min_len) return(NULL)
    sub <- x[s:e]; L <- length(sub)
    micro <- max(vapply(1:6, function(q) mean(sub[1:(L - q)] == sub[(q + 1):L]), numeric(1)))
    if (micro > max_micro) return(NULL)                       # microsatélite / baja complejidad
    sn <- names(periodos_sat)[vapply(periodos_sat, function(v) p_mod %in% v, logical(1))]
    if (!length(sn)) return(NULL)
    data.frame(acc = acc, start = s, end = e, satelite = sn[1], periodo = p_mod,
               ident = round(mean(best[s0:e0]), 3), fuente = "de novo", stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}

hits_dn <- NULL
if (buscar_denovo) {
  cat("\nBúsqueda de novo de satélites (lenta)...\n")
  dn <- do.call(rbind, lapply(names(sl), buscar_sat_denovo, periodos_sat = periodos_sat,
                              win = dn_win, min_ident = dn_min_ident, min_len = dn_min_len))
  if (!is.null(dn) && nrow(dn) > 0) {
    gr_dn <- armonizar(GRanges(dn$acc, IRanges(dn$start, dn$end)))
    tel_red <- reduce(tel_gr, ignore.strand = TRUE)
    ov <- findOverlaps(gr_dn, tel_red, ignore.strand = TRUE)
    solape <- pintersect(ranges(gr_dn)[queryHits(ov)], ranges(tel_red)[subjectHits(ov)])
    cob <- tapply(width(solape), queryHits(ov), sum)
    frac_telo <- numeric(length(gr_dn))
    frac_telo[as.integer(names(cob))] <- cob / width(gr_dn)[as.integer(names(cob))]
    keep <- frac_telo < 0.5                                   # solo descarta bloques que son mayormente telómero
    hits_dn <- dn[keep, c("acc", "start", "end", "satelite", "fuente")]
  }
}

# (SAT350 por GenBank/BLASTN/TRF: ver 06b_sat350.R; aporta `hits_gb`)
#hits <- do.call(rbind, list(hits_rm, hits_trf, hits_dn, hits_gb))
hits <- do.call(rbind, list(if (usar_rm_trf_sat) hits_rm, if (usar_rm_trf_sat) hits_trf, hits_dn, hits_gb))

if (is.null(hits))
  hits <- data.frame(acc = character(), start = integer(), end = integer(),
                     satelite = character(), fuente = character(), stringsAsFactors = FALSE)

vacio_sat <- data.frame(satelite = character(), chrom = character(), accession = character(),
                        start = integer(), end = integer(), longitud_pb = integer(),
                        n_fragmentos = integer(), fuentes = character(), ubicacion = character(),
                        en_subtelomero = character(), stringsAsFactors = FALSE)

sat_bloques <- do.call(rbind, c(list(vacio_sat), lapply(names(patron_sat), function(sn) {
  h <- hits[hits$satelite == sn & hits$acc %in% names(sl), ]
  if (nrow(h) == 0) return(NULL)
  g   <- armonizar(GRanges(h$acc, IRanges(h$start, h$end)))
  red <- reduce(g, min.gapwidth = gap_sat + 1, ignore.strand = TRUE)
  ov  <- findOverlaps(red, g)
  fu  <- tapply(h$fuente[subjectHits(ov)], queryHits(ov), function(x) paste(sort(unique(x)), collapse = " + "))
  fuentes <- character(length(red)); fuentes[as.integer(names(fu))] <- as.character(fu)
  a <- as.character(seqnames(red)); L <- as.numeric(sl[a])
  data.frame(satelite = sn, chrom = nombre_chr(a), accession = a, start = start(red), end = end(red),
             longitud_pb = width(red), n_fragmentos = countOverlaps(red, g, ignore.strand = TRUE),
             fuentes = fuentes,
             ubicacion = ifelse(start(red) <= w_st_max, "Extremo inicio",
                                ifelse(end(red) >= L - w_st_max + 1, "Extremo final", "Interno")),
             en_subtelomero = ifelse(overlapsAny(red, st_gr, ignore.strand = TRUE), "Si", "No"),
             stringsAsFactors = FALSE)
})))
sat_bloques <- sat_bloques[order(factor(sat_bloques$chrom, levels = eq$nombre),
                                 sat_bloques$satelite, sat_bloques$start), ]
rownames(sat_bloques) <- NULL

for (sn in names(patron_sat)) {
  tabla[[sn]] <- "No"
  b <- sat_bloques[sat_bloques$satelite == sn, ]
  if (nrow(b) > 0) {
    gb <- armonizar(GRanges(b$accession, IRanges(b$start, b$end)))
    tabla[[sn]][unique(queryHits(findOverlaps(bins, gb)))] <- "Si"
  }
}

trf_consensos <- NULL
if (!is.null(trf_sat) && nrow(trf_sat) > 0) {
  trf_consensos <- data.frame(
    satelite = trf_sat$satelite, chrom = nombre_chr(trf_sat$seqnames),
    start = trf_sat$start, end = trf_sat$end, longitud_pb = trf_sat$end - trf_sat$start + 1,
    periodo = trf_sat$periodo, copias = trf_sat$copias, pct_match = trf_sat$pct_match,
    consenso = trf_sat$consenso, stringsAsFactors = FALSE)
  trf_consensos <- trf_consensos[order(trf_consensos$satelite, -trf_consensos$longitud_pb), ]
  rownames(trf_consensos) <- NULL
}
if (nrow(sat_bloques) == 0) cat("No se encontraron satélites (", paste(names(patron_sat), collapse = "/"),
                                "). Revisa rep_file/trf_file o activa `buscar_denovo`.\n")



