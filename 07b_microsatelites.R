# ============================================================
# 07b_microsatelites.R  -  Microsatelites: primers -> seqkit amplicon -> TRF   (antes: microsatelites.R)
# ------------------------------------------------------------
# DONDE CORRE : despues de 07_satelites.R (usa st_gr y sat350_regiones para dar contexto a cada locus).
#               Tambien corre en modo "ligero" (sin GFF): en ese caso no hay st_gr ni tabla.
# REQUIERE   : 01 (primers_micro, micro_mismatches, micro_max_amplicon, micro_trf_params, n_hilos,
#              rehacer_externos, carpeta_externos), 02/02L (fa, sl, eq, nombre_chr, armonizar),
#              entradas_genoma.R (rc, norm_lab, correr_trf_fresh, parse_trf_full, asegurar_fasta_chr)
#              Opcional: 06 (st_gr), 04 (tabla), 06b (sat350_regiones)
#              Programas: seqkit (obligatorio para este modulo), trf (opcional: sin TRF no hay repeticiones)
# CREA       : micro_amps       amplicones encontrados (coordenadas, hebra, secuencia)
#              micro_trf        todas las repeticiones TRF dentro de los amplicones
#              micro_resumen    una fila por locus: TRF principal, motivo, comparacion con el articulo
#              microsat_loci    una fila por amplicon para graficar: posicion, distancia al borde,
#                               en_subtelomero, en_sat350
#              tabla$microsat   "Si"/"No" por ventana (si existe `tabla`)
# SIN seqkit o SIN primers_micro: el modulo se omite y los objetos quedan en NULL.
# CACHE      : los amplicones se guardan en carpeta_externos/microsat_amplicones.rds
# USADO POR  : 14_excel, 15_figuras (F18, F19), 90_guardar_genoma
# ============================================================

micro_amps <- NULL; micro_trf <- NULL; micro_resumen <- NULL; microsat_loci <- NULL
hay_seqkit <- nzchar(Sys.which("seqkit")); hay_trf <- nzchar(Sys.which("trf"))

if (is.null(primers_micro) || nrow(primers_micro) == 0) {
  message("07b_microsatelites: sin primers_micro: se omite.")
} else if (!hay_seqkit) {
  message("07b_microsatelites: seqkit no esta en el PATH: se omite (conda install -c bioconda seqkit).")
} else {

res_m <- local({
  primers <- primers_micro
  primers$cromosoma_paper_lab <- norm_lab(sub("^.*chr", "", primers$cromosoma_paper))
  dir.create(carpeta_externos, showWarnings = FALSE, recursive = TRUE)
  fasta_chr <- asegurar_fasta_chr(fa, carpeta_externos)

  hamming <- function(a, b) sum(strsplit(a, "")[[1]] != strsplit(b, "")[[1]])
  # Motivo canonico: minimo alfabetico entre rotaciones y complemento inverso (TG = GT = CA = AC)
  canon <- function(m) {
    if (is.na(m) || m == "") return(NA_character_)
    m <- toupper(m); n <- nchar(m)
    cand <- c(m, rc(m))
    rots <- unlist(lapply(cand, function(s)
      vapply(0:(n - 1), function(i) paste0(substring(s, i + 1, n), substring(s, 1, i)), "")))
    sort(rots)[1]
  }
  canon_esp <- function(motivo) {
    if (is.na(motivo)) return(NA_character_)
    unique(vapply(strsplit(motivo, "[;/]")[[1]], canon, ""))
  }

  # ---- seqkit amplicon (con cache) ----
  buscar <- function(i) {
    p <- primers[i, ]
    message(sprintf("[%d/%d] %s ...", i, nrow(primers), p$marcador))
    run <- function(m) {
      l <- suppressWarnings(system2("seqkit", c("amplicon", "-F", p$fwd, "-R", p$rev, "-m", as.character(m),
                                                "-j", as.character(n_hilos), "--bed", shQuote(fasta_chr)),
                                    stdout = TRUE, stderr = FALSE))
      l[nzchar(l)]
    }
    lin <- run(0)                                              # 1) exacta (rapida)
    if (length(lin) == 0 && micro_mismatches > 0) lin <- run(micro_mismatches)   # 2) con mismatches
    if (length(lin) == 0) return(NULL)
    partes <- strsplit(lin, "\t", fixed = TRUE)
    acc <- vapply(partes, `[`, "", 1)
    ini <- as.integer(vapply(partes, `[`, "", 2)); fin <- as.integer(vapply(partes, `[`, "", 3))
    ok <- acc %in% names(sl) & (fin - ini) <= micro_max_amplicon
    if (any(!ok)) message("   descartados ", sum(!ok), " amplicon(es) fuera de cromosomas o > ", micro_max_amplicon, " pb")
    acc <- acc[ok]; ini <- ini[ok]; fin <- fin[ok]
    if (length(acc) == 0) return(NULL)
    do.call(rbind, lapply(seq_along(acc), function(h) {
      s_plus <- as.character(subseq(fa[[acc[h]]], ini[h] + 1, fin[h]))
      strand <- if (hamming(substr(s_plus, 1, nchar(p$fwd)), p$fwd) <= micro_mismatches) "+" else "-"
      s_fin  <- if (strand == "+") s_plus else rc(s_plus)
      data.frame(marcador = p$marcador, hit = h, id_amplicon = paste0(p$marcador, "_h", h),
                 accession = acc[h], amp_ini = ini[h] + 1L, amp_fin = fin[h], hebra = strand,
                 longitud_amplicon = nchar(s_fin), secuencia = s_fin, stringsAsFactors = FALSE)
    }))
  }

  cache <- file.path(carpeta_externos, "microsat_amplicones.rds")
  clave <- list(primers = primers[, c("marcador", "fwd", "rev")], mm = micro_mismatches, max = micro_max_amplicon)
  amps <- NULL
  if (!isTRUE(rehacer_externos) && file.exists(cache)) {
    cc <- readRDS(cache)
    if (identical(cc$clave, clave)) { amps <- cc$amps; message("Microsatelites: se reutilizan los amplicones en cache") }
  }
  if (is.null(amps)) {
    message("Buscando amplicones con seqkit...")
    amps <- do.call(rbind, lapply(seq_len(nrow(primers)), buscar))
    saveRDS(list(clave = clave, amps = amps), cache)
  }
  if (is.null(amps) || nrow(amps) == 0) {
    warning("Microsatelites: ningun primer dio amplicon en este genoma.")
    return(NULL)
  }
  amps$chrom <- nombre_chr(amps$accession)
  sin_hits <- setdiff(primers$marcador, amps$marcador)
  if (length(sin_hits)) message("Sin amplicon: ", paste(sin_hits, collapse = ", "))

  # ---- TRF sobre los amplicones ----
  trf_m <- parse_trf_full(NULL)
  if (hay_trf) {
    fl <- file.path(carpeta_externos, "microsat_loci.fasta")
    writeXStringSet(DNAStringSet(setNames(amps$secuencia, amps$id_amplicon)), fl)
    p <- correr_trf_fresh(fl, carpeta_externos, micro_trf_params)
    if (!is.null(p)) trf_m <- parse_trf_full(p)
  } else message("Microsatelites: sin TRF en el PATH, no se analizan las repeticiones.")
  names(trf_m)[1] <- "id_amplicon"

  # ---- Integrar TRF + amplicones + primers ----
  trf_full <- trf_m
  if (nrow(trf_full) > 0) {
    mk <- match(trf_full$id_amplicon, amps$id_amplicon)
    trf_full$marcador <- amps$marcador[mk]; trf_full$chrom <- amps$chrom[mk]
    trf_full$amp_ini  <- amps$amp_ini[mk];  trf_full$amp_fin <- amps$amp_fin[mk]; trf_full$hebra <- amps$hebra[mk]
    trf_full$motivo_esperado  <- primers$motivo[match(trf_full$marcador, primers$marcador)]
    trf_full$motivo_canonico  <- vapply(trf_full$consenso, canon, "")
    esp <- setNames(lapply(primers$motivo, canon_esp), primers$marcador)
    trf_full$coincide_motivo <- vapply(seq_len(nrow(trf_full)), function(k) {
      e <- esp[[trf_full$marcador[k]]]
      if (is.null(e) || anyNA(e)) NA else trf_full$motivo_canonico[k] %in% e
    }, logical(1))
    plus <- trf_full$hebra == "+"
    trf_full$gen_ini <- ifelse(plus, trf_full$amp_ini + trf_full$rep_ini - 1L, trf_full$amp_fin - trf_full$rep_fin + 1L)
    trf_full$gen_fin <- ifelse(plus, trf_full$amp_ini + trf_full$rep_fin - 1L, trf_full$amp_fin - trf_full$rep_ini + 1L)
    o <- order(-trf_full$score)
    trf_full$principal <- seq_len(nrow(trf_full)) %in% o[!duplicated(trf_full$id_amplicon[o])]
  }
  pr <- if (nrow(trf_full) > 0) trf_full[trf_full$principal, ] else trf_full
  ir <- match(amps$id_amplicon, pr$id_amplicon)
  ip <- match(amps$marcador, primers$marcador)
  g  <- function(v) if (nrow(pr) > 0) pr[[v]][ir] else rep(NA, nrow(amps))
  n_rep <- as.integer(table(factor(trf_full$id_amplicon, levels = amps$id_amplicon)))
  resumen <- data.frame(
    marcador = amps$marcador, id_amplicon = amps$id_amplicon, tipo = primers$tipo[ip],
    cromosoma = amps$chrom, accession = amps$accession, amp_ini = amps$amp_ini, amp_fin = amps$amp_fin,
    hebra = amps$hebra, longitud_amplicon = amps$longitud_amplicon, n_repeticiones_TRF = n_rep,
    rep_ini = g("rep_ini"), rep_fin = g("rep_fin"), gen_ini = g("gen_ini"), gen_fin = g("gen_fin"),
    periodo = g("periodo"), copias = g("copias"), consenso = g("consenso"),
    motivo_esperado = primers$motivo[ip], coincide_motivo = g("coincide_motivo"),
    pct_matches = g("pct_matches"), pct_indels = g("pct_indels"), score = g("score"),
    entropia = g("entropia"), secuencia_repetida = g("secuencia_repetida"),
    cromosoma_paper = primers$cromosoma_paper[ip], pos_ini_paper = primers$pos_ini_paper[ip],
    pos_fin_paper = primers$pos_fin_paper[ip], stringsAsFactors = FALSE)
  resumen$longitud_paper <- resumen$pos_fin_paper - resumen$pos_ini_paper + 1
  resumen$dif_longitud   <- resumen$longitud_amplicon - resumen$longitud_paper
  resumen$mismo_cromosoma <- norm_lab(resumen$cromosoma) == primers$cromosoma_paper_lab[ip]
  resumen$desplazamiento_vs_paper <- resumen$amp_ini - resumen$pos_ini_paper
  resumen$observacion <- ifelse(!hay_trf, "TRF no disponible",
                         ifelse(resumen$n_repeticiones_TRF == 0, "TRF no detecto repeticion (baja el score minimo o revisa)",
                         ifelse(!is.na(resumen$coincide_motivo) & !resumen$coincide_motivo, "Motivo detectado distinto al esperado", "OK")))
  if (length(sin_hits)) {
    fila <- resumen[rep(NA_integer_, length(sin_hits)), ]
    fila$marcador <- sin_hits; fila$observacion <- "Sin amplicon en el genoma"
    resumen <- rbind(resumen, fila)
  }
  rownames(resumen) <- NULL
  if (nrow(trf_full) > 0) trf_full <- trf_full[, c("marcador", "id_amplicon", "principal", "rep_ini", "rep_fin", "periodo",
                                                   "copias", "tam_consenso", "pct_matches", "pct_indels", "score", "pct_A",
                                                   "pct_C", "pct_G", "pct_T", "entropia", "consenso", "secuencia_repetida",
                                                   "motivo_esperado", "coincide_motivo", "chrom", "gen_ini", "gen_fin")]
  list(amps = amps[, setdiff(names(amps), "secuencia")], trf = trf_full, resumen = resumen)
})

if (!is.null(res_m)) {
  micro_amps <- res_m$amps; micro_trf <- res_m$trf; micro_resumen <- res_m$resumen
  m  <- micro_amps
  microsat_loci <- data.frame(
    chrom = m$chrom, accession = m$accession, marcador = m$marcador,
    nombre = gsub("_", "-", sub("_(multiplex|intron)$", "", m$marcador)),
    tipo = primers_micro$tipo[match(m$marcador, primers_micro$marcador)],
    ini = m$amp_ini, fin = m$amp_fin, pos_mb = (m$amp_ini + m$amp_fin) / 2 / 1e6, stringsAsFactors = FALSE)
  L  <- as.numeric(sl[microsat_loci$accession])
  microsat_loci$dist_borde_pb <- pmin(microsat_loci$ini - 1, L - microsat_loci$fin)
  gm <- armonizar(GRanges(microsat_loci$accession, IRanges(microsat_loci$ini, microsat_loci$fin)))
  microsat_loci$en_subtelomero <- if (exists("st_gr") && length(st_gr) > 0)
    ifelse(overlapsAny(gm, st_gr, ignore.strand = TRUE), "Si", "No") else "No"
  microsat_loci$en_sat350 <- if (!is.null(sat350_regiones) && nrow(sat350_regiones) > 0)
    ifelse(overlapsAny(gm, armonizar(GRanges(sat350_regiones$accession, IRanges(sat350_regiones$start, sat350_regiones$end))),
                       ignore.strand = TRUE), "Si", "No") else "No"
  if (exists("tabla") && exists("bins")) {
    tabla$microsat <- "No"
    tabla$microsat[unique(queryHits(findOverlaps(bins, gm)))] <- "Si"
  }
  cat("\nMicrosatelites:", nrow(microsat_loci), "loci |", sum(microsat_loci$en_subtelomero == "Si"), "en subtelomero |",
      sum(microsat_loci$en_sat350 == "Si"), "dentro de SAT350 |",
      sum(micro_resumen$observacion == "OK", na.rm = TRUE), "con TRF 'OK'\n")
}
}

if (!exists("sat_bloques")) {
  sat_bloques <- data.frame(
    chrom = character(),
    satelite = character(),
    start = integer(),
    end = integer(),
    longitud_pb = numeric(),
    fuentes = character(),
    stringsAsFactors = FALSE
  )
}


###################Resumen satelites encontrados 

# Una fila por cromosoma y familia, incluso si no hubo detecciones
tabla_satelites_chr <- expand.grid(
  chrom = eq$nombre,
  satelite = names(patron_sat),
  stringsAsFactors = FALSE
)

tabla_satelites_chr$n_bloques <- 0L
tabla_satelites_chr$longitud_total_kb <- 0
tabla_satelites_chr$ubicaciones <- "No detectado"

for (i in seq_len(nrow(tabla_satelites_chr))) {
  chr <- tabla_satelites_chr$chrom[i]
  fam <- tabla_satelites_chr$satelite[i]
  
  bloques <- sat_bloques[
    sat_bloques$chrom == chr & sat_bloques$satelite == fam,
  ]
  
  if (nrow(bloques) > 0) {
    tabla_satelites_chr$n_bloques[i] <- nrow(bloques)
    tabla_satelites_chr$longitud_total_kb[i] <-
      sum(bloques$longitud_pb, na.rm = TRUE) / 1000
    
    tabla_satelites_chr$ubicaciones[i] <- paste(
      paste0(
        bloques$start, "-", bloques$end, " pb",
        " [", bloques$fuentes, "]"
      ),
      collapse = "; "
    )
  }
}

View(tabla_satelites_chr)

################# satelites encontrados con primers 

if (is.null(microsat_loci) || nrow(microsat_loci) == 0) {
  message("No hay loci detectados con los primers.")
} else {
  tabla_primer_hits <- microsat_loci[
    order(
      match(microsat_loci$chrom, eq$nombre),
      microsat_loci$ini
    ),
    c(
      "chrom", "marcador", "nombre", "tipo",
      "accession", "ini", "fin", "pos_mb",
      "en_subtelomero", "en_sat350"
    )
  ]
  
  View(tabla_primer_hits)
}
