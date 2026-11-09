# ============================================================
# 06b_sat350.R  -  Familia SAT350: GenBank -> BLASTN -> regiones -> TRF   (antes: Sat350.R)
# ------------------------------------------------------------
# DONDE CORRE : despues de 06_subtelomeros.R y ANTES de 07_satelites.R (07 suma `hits_gb` a `hits`)
#               Tambien corre en modo "ligero" (genomas sin GFF): solo necesita fa, sl, eq.
# REQUIERE   : 01 (sat350_refs, gb_accesiones, gb_satelite, gb_fasta_local, usar_blast, gb_evalue,
#              gb_min_ident, gb_min_len, gap_sat, sat350_flanco, sat350_periodo, trf_params_sat,
#              n_hilos, rehacer_externos, carpeta_externos, dir_salida),
#              02/02L (fa, sl, eq, nombre_chr), entradas_genoma.R (correr_trf_fresh, parse_trf_full,
#              asegurar_fasta_chr)
#              Programas: blastn + makeblastdb (obligatorios para este modulo), trf (opcional)
# CREA       : gb_hits, gb_mat, gb_presencia, metodo_gb, hits_gb     (mismos nombres que usaba 07/14/15)
#              sat350_refs_tab      las secuencias de referencia + longitud GenBank vs articulo
#              sat350_resumen_q     hits por secuencia de referencia
#              sat350_trf_ref       TRF sobre las referencias
#              sat350_regiones      regiones fusionadas (gap_sat) + resumen TRF (periodo, copias, monomero)
#              sat350_trf_regiones  todas las repeticiones TRF dentro de las regiones
# SIN BLAST   : el modulo se omite (todos los objetos quedan vacios) y el pipeline sigue.
# CACHE       : BLAST, base de datos y TRF se reutilizan si ya existen (rehacer_externos = TRUE los rehace).
# USADO POR  : 07 (hits_gb), 07b (sat350_regiones), 14_excel, 15_figuras (F12, F20), 90_guardar_genoma
# ============================================================

gb_empty <- data.frame(query = character(), acc = character(), start = integer(), end = integer(),
                       strand = character(), longitud = numeric(), identidad_pct = numeric(),
                       score = numeric(), cobertura_query_pct = numeric(), metodo = character(),
                       chrom = character(), stringsAsFactors = FALSE)
sat350_refs_tab <- NULL; sat350_resumen_q <- NULL; sat350_trf_ref <- NULL
sat350_regiones <- NULL; sat350_trf_regiones <- NULL
gb_hits <- gb_empty; gb_mat <- NULL; gb_presencia <- NULL; metodo_gb <- NA_character_; hits_gb <- NULL

usa_gb    <- length(gb_accesiones) > 0
hay_blast <- isTRUE(usar_blast) && nzchar(Sys.which("blastn")) && nzchar(Sys.which("makeblastdb"))
hay_trf   <- nzchar(Sys.which("trf"))
if (usa_gb && !hay_blast)
  message("06b_sat350: blastn/makeblastdb no disponibles (o usar_blast = FALSE): se omite SAT350 por GenBank.")

# ---- Descarga (o lectura local) de las secuencias de referencia ----
leer_genbank <- function(accs, fasta_local) {
  if (file.exists(fasta_local)) return(readDNAStringSet(fasta_local))
  url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=",
                paste(accs, collapse = ","), "&rettype=fasta&retmode=text")
  tmp <- tempfile(fileext = ".fasta")
  ok <- tryCatch({ suppressWarnings(download.file(url, tmp, quiet = TRUE, mode = "wb"))
    file.exists(tmp) && file.size(tmp) > 0 }, error = function(e) FALSE)
  x <- if (ok) tryCatch(readDNAStringSet(tmp), error = function(e) NULL) else NULL
  if (is.null(x) || length(x) == 0) {
    message("No se pudieron descargar las secuencias de GenBank. Guardalas en '", fasta_local,
            "' (FASTA) y vuelve a correr.")
    return(NULL)
  }
  writeXStringSet(x, fasta_local)
  x
}

gb_q <- NULL
if (usa_gb && hay_blast) {
  dir_ref <- carpeta_referencias      # las referencias no dependen del genoma
  dir.create(dir_ref, showWarnings = FALSE, recursive = TRUE)
  gb_q <- leer_genbank(gb_accesiones, file.path(dir_ref, basename(gb_fasta_local)))
}

if (!is.null(gb_q)) {
  gb_desc <- sub("^[^ ]+ ?", "", names(gb_q))
  names(gb_q) <- sub("\\.[0-9]+$", "", sub(" .*", "", names(gb_q)))
  keep <- names(gb_q) %in% gb_accesiones
  gb_q <- gb_q[keep]; gb_desc <- gb_desc[keep]
  faltan <- setdiff(gb_accesiones, names(gb_q))
  if (length(faltan) > 0) message("No se encontraron en el FASTA de GenBank: ", paste(faltan, collapse = ", "))
  cat("\nSecuencias GenBank cargadas:", paste0(names(gb_q), " (", width(gb_q), " pb)", collapse = ", "), "\n")

  # Tabla de referencias: longitud descargada vs la del articulo
  ref <- sat350_refs
  ref$longitud_genbank   <- as.integer(width(gb_q)[match(ref$acceso, names(gb_q))])
  ref$dif_longitud       <- ref$longitud_genbank - ref$longitud_paper
  ref$definicion_genbank <- gb_desc[match(ref$acceso, names(gb_q))]
  sat350_refs_tab <- ref

  dir.create(carpeta_externos, showWarnings = FALSE, recursive = TRUE)
  qf <- file.path(carpeta_externos, "sat350_queries.fasta")
  writeXStringSet(gb_q, qf)

  # ---- TRF sobre las referencias (solo como referencia) ----
  if (hay_trf) {
    p <- correr_trf_fresh(qf, carpeta_externos, trf_params_sat)
    if (!is.null(p)) sat350_trf_ref <- parse_trf_full(p)
  }

  # ---- BLASTN contra los cromosomas analizados (una sola vez, con cache) ----
  metodo_gb <- "BLASTN"
  fasta_chr <- asegurar_fasta_chr(fa, carpeta_externos)
  dir.create(file.path(carpeta_externos, "blast_db"), showWarnings = FALSE, recursive = TRUE)
  db_pref <- file.path(carpeta_externos, "blast_db", genoma_id)
  if (isTRUE(rehacer_externos)) unlink(Sys.glob(paste0(db_pref, "*")))
  if (length(Sys.glob(paste0(db_pref, "*.nsq"))) == 0)
    system2("makeblastdb", c("-in", shQuote(fasta_chr), "-dbtype", "nucl", "-parse_seqids",
                             "-out", shQuote(db_pref)), stdout = FALSE, stderr = FALSE)
  of <- file.path(carpeta_externos, "sat350_vs_genoma.tsv")
  cols_b <- c("qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send",
              "evalue", "bitscore", "qlen")
  if (isTRUE(rehacer_externos) || !file.exists(of)) {
    message("Ejecutando BLASTN (SAT350)...")
    system2("blastn", c("-query", shQuote(qf), "-db", shQuote(db_pref), "-task", "blastn",
                        "-evalue", format(gb_evalue, scientific = TRUE),
                        "-perc_identity", as.character(gb_min_ident),
                        "-dust", "no", "-soft_masking", "false",
                        "-max_target_seqs", "100000",          # como texto: numerico se escribe "1e+05"
                        "-num_threads", as.character(n_hilos),
                        "-outfmt", shQuote(paste("6", paste(cols_b, collapse = " "))),
                        "-out", shQuote(of)))
  } else message("BLASTN: se reutiliza ", of)

  gb_raw <- NULL
  if (file.exists(of) && file.size(of) > 0) {
    b <- read.delim(of, header = FALSE, col.names = cols_b, quote = "", stringsAsFactors = FALSE)
    ids <- unique(b$sseqid)
    s2  <- sub("^(lcl|gnl\\|[^|]+)\\|", "", ids)
    map <- vapply(seq_along(ids), function(k) {
      if (s2[k] %in% names(sl)) return(s2[k])
      h <- names(sl)[vapply(names(sl), grepl, logical(1), x = ids[k], fixed = TRUE)]
      if (length(h)) h[which.max(nchar(h))] else NA_character_
    }, character(1))
    b$acc <- unname(map[match(b$sseqid, ids)])
    b <- b[!is.na(b$acc) & b$pident >= gb_min_ident & b$length >= gb_min_len, ]
    if (nrow(b) > 0)
      gb_raw <- data.frame(query = b$qseqid, acc = b$acc, start = pmin(b$sstart, b$send), end = pmax(b$sstart, b$send),
                           strand = ifelse(b$send >= b$sstart, "+", "-"), longitud = b$length,
                           identidad_pct = b$pident, score = b$bitscore,
                           cobertura_query_pct = round(100 * b$length / b$qlen, 1),
                           metodo = "BLASTN", stringsAsFactors = FALSE)
  }
  gb_hits <- if (is.null(gb_raw) || nrow(gb_raw) == 0) gb_empty else { gb_raw$chrom <- nombre_chr(gb_raw$acc); gb_raw }
  gb_hits <- gb_hits[order(factor(gb_hits$query, levels = gb_accesiones),
                           factor(gb_hits$chrom, levels = eq$nombre), gb_hits$start), ]
  rownames(gb_hits) <- NULL
  cat("SAT350 (BLASTN):", nrow(gb_hits), "hits\n")

  # ---- Presencia por secuencia y cromosoma (pb cubiertos, fusionados) ----
  gb_mat <- do.call(rbind, lapply(gb_accesiones, function(q)
    vapply(names(sl), function(a) {
      h <- gb_hits[gb_hits$query == q & gb_hits$acc == a, ]
      if (nrow(h) == 0) 0 else sum(width(reduce(IRanges(h$start, h$end))))
    }, numeric(1))))
  dimnames(gb_mat) <- list(gb_accesiones, eq$nombre)
  gb_presencia <- data.frame(
    accesion = gb_accesiones,
    longitud_referencia_pb = as.integer(width(gb_q)[match(gb_accesiones, names(gb_q))]),
    descripcion = gb_desc[match(gb_accesiones, names(gb_q))],
    gb_mat, total_pb = rowSums(gb_mat), n_cromosomas = rowSums(gb_mat > 0),
    n_hits = as.integer(table(factor(gb_hits$query, levels = gb_accesiones))),
    presente = ifelse(rowSums(gb_mat) > 0, "Si", "No"), metodo = metodo_gb,
    check.names = FALSE, stringsAsFactors = FALSE)
  rownames(gb_presencia) <- NULL
  cat("\n=== Secuencias GenBank de", gb_satelite, "en el genoma ===\n")
  print(gb_presencia[, c("accesion", "longitud_referencia_pb", "total_pb", "n_cromosomas", "n_hits", "presente")])

  # ---- Resumen por secuencia de referencia ----
  rq <- do.call(rbind, lapply(gb_accesiones, function(q) {
    h <- gb_hits[gb_hits$query == q, ]
    data.frame(acceso = q, n_hits = nrow(h), n_cromosomas = length(unique(h$acc)),
               pident_media = if (nrow(h)) round(mean(h$identidad_pct), 1) else NA_real_,
               mejor_bitscore = if (nrow(h)) max(h$score) else NA_real_, stringsAsFactors = FALSE)
  }))
  sat350_resumen_q <- merge(sat350_refs_tab[, c("nombre", "acceso", "longitud_paper", "longitud_genbank", "dif_longitud")],
                            rq, by = "acceso", all.y = TRUE)
  sat350_resumen_q <- sat350_resumen_q[, c("nombre", "acceso", "longitud_paper", "longitud_genbank",
                                           "dif_longitud", "n_hits", "n_cromosomas", "pident_media", "mejor_bitscore")]

  # ---- Fusion de hits en regiones (gap_sat) ----
  if (nrow(gb_hits) > 0) {
    red  <- reduce(GRanges(gb_hits$acc, IRanges(gb_hits$start, gb_hits$end)), min.gapwidth = gap_sat + 1,
                   ignore.strand = TRUE, with.revmap = TRUE)
    rmap <- as.list(mcols(red)$revmap)
    ac   <- as.character(seqnames(red))
    sat350_regiones <- data.frame(
      region_id = sprintf("%s_%d_%d", ac, start(red), end(red)),
      chrom = nombre_chr(ac), accession = ac, start = start(red), end = end(red), ancho = width(red),
      n_hits = lengths(rmap),
      mejor_pident   = vapply(rmap, function(ix) max(gb_hits$identidad_pct[ix]), numeric(1)),
      mejor_bitscore = vapply(rmap, function(ix) max(gb_hits$score[ix]), numeric(1)),
      mejor_query    = vapply(rmap, function(ix) gb_hits$query[ix][which.max(gb_hits$score[ix])], character(1)),
      n_rep_trf = NA_integer_, periodo = NA_integer_, copias = NA_real_, pct_matches = NA_integer_,
      score_trf = NA_integer_, consenso = NA_character_, monomero_sat350 = NA,
      stringsAsFactors = FALSE)
    sat350_regiones <- sat350_regiones[order(match(sat350_regiones$chrom, eq$nombre), sat350_regiones$start), ]
    rownames(sat350_regiones) <- NULL

    # ---- Extraer las regiones del genoma y correr TRF ----
    if (hay_trf) {
      L_ac <- as.numeric(sl[sat350_regiones$accession])
      a <- as.integer(pmax(1, sat350_regiones$start - sat350_flanco))
      b <- as.integer(pmin(L_ac, sat350_regiones$end + sat350_flanco))
      seqs <- subseq(fa[sat350_regiones$accession], start = a, end = b)
      names(seqs) <- sat350_regiones$region_id
      fr <- file.path(carpeta_externos, "sat350_regiones.fasta")
      writeXStringSet(seqs, fr)
      cat("Extraidas", length(seqs), "regiones SAT350 (", sum(width(seqs)), "pb ) -> TRF\n")
      p <- correr_trf_fresh(fr, carpeta_externos, trf_params_sat)
      trf_r <- if (!is.null(p)) parse_trf_full(p) else parse_trf_full(NULL)
      if (nrow(trf_r) > 0) {
        names(trf_r)[1] <- "region_id"
        o   <- order(-trf_r$score)
        dom <- trf_r[o[!duplicated(trf_r$region_id[o])], ]
        i   <- match(sat350_regiones$region_id, dom$region_id)
        nr  <- as.integer(table(trf_r$region_id)[sat350_regiones$region_id]); nr[is.na(nr)] <- 0L
        sat350_regiones$n_rep_trf  <- nr
        sat350_regiones$periodo    <- dom$periodo[i]
        sat350_regiones$copias     <- dom$copias[i]
        sat350_regiones$pct_matches <- dom$pct_matches[i]
        sat350_regiones$score_trf  <- dom$score[i]
        sat350_regiones$consenso   <- dom$consenso[i]
        sat350_regiones$monomero_sat350 <- ifelse(is.na(sat350_regiones$periodo), NA,
                                                  sat350_regiones$periodo >= sat350_periodo[1] &
                                                    sat350_regiones$periodo <= sat350_periodo[2])
        sat350_trf_regiones <- trf_r
        cat("Regiones con monomero de periodo ", sat350_periodo[1], "-", sat350_periodo[2], " pb: ",
            sum(sat350_regiones$monomero_sat350, na.rm = TRUE), " de ", nrow(sat350_regiones), "\n", sep = "")
      } else sat350_regiones$n_rep_trf <- 0L
    }
  }
  hits_gb <- if (nrow(gb_hits) > 0)
    data.frame(acc = gb_hits$acc, start = gb_hits$start, end = gb_hits$end, satelite = gb_satelite,
               fuente = paste0("GenBank:", gb_hits$query), stringsAsFactors = FALSE) else NULL
}
