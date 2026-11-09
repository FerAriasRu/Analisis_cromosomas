# ============================================================
# entradas_genoma.R  -  FUNCIONES para trabajar con cualquier genoma / GFF
# (no ejecuta nada al cargarse)
# ------------------------------------------------------------
# LLAMADO POR: 02_datos_base.R (source)
# FUNCIONES  :
#   col_txt(gr, nm)                          columna de mcols como texto (NA si falta o está vacía)
#   col_primera(gr, nombres)                 primera columna no vacía entre varias candidatas
#   restar(a, b)                             setdiff sin hebra
#   nombre_corto(acc)                        "TGME49_chrIa" -> "Ia"; sin "chr" deja el accession
#   seleccionar_cromosomas(...)              decide qué secuencias se analizan y cómo se llaman
#   tipos_gen_detectar(gff_all, ...)         detecta qué tipos de feature son genes
#   construir_genes(gff_all, ...)            GRanges de genes (plan B: desde transcritos)
#   construir_exones_intrones(...)           exones/intrones (TxDb; plan B directo del GFF)
#   leer_repeatmasker(f), rm_vacio(), rm_a_granges(df)
#   repeticiones_desde_softmask(fna, accs)   repeticiones = tramos en minúscula del FASTA
#   correr_trf(...), correr_repeatmasker(...) ejecutan el programa si está instalado
# ============================================================

col_txt <- function(gr, nm) {
  if (!nm %in% names(mcols(gr))) return(rep(NA_character_, length(gr)))
  x <- mcols(gr)[[nm]]
  out <- unname(vapply(x, function(v)
    if (length(v) == 0 || all(is.na(v))) NA_character_ else paste(v, collapse = ","), character(1)))
  out[!is.na(out) & out == ""] <- NA_character_
  out
}

col_primera <- function(gr, nombres) {
  out <- rep(NA_character_, length(gr))
  for (nm in nombres) out <- ifelse(is.na(out), col_txt(gr, nm), out)
  out
}

restar <- function(a, b) setdiff(a, b, ignore.strand = TRUE)

nombre_corto <- function(acc) {
  x <- sub("^.*[Cc]hr(omosome)?[_ ]?", "", acc)
  ifelse(x == acc | x == "", acc, x)
}

seleccionar_cromosomas <- function(sl, min_len = 1e6, patron_incluir = NULL, patron_excluir = NULL,
                                   eq_manual = NULL, orden = "archivo") {
  if (!is.null(eq_manual)) {
    stopifnot(all(c("acc", "nombre") %in% names(eq_manual)))
    faltan <- setdiff(eq_manual$acc, names(sl))
    if (length(faltan))
      stop("eq_manual: estos accessions no están en el FASTA/.fai: ", paste(faltan, collapse = ", "),
           "\n  Nombres disponibles (primeros 10): ", paste(head(names(sl), 10), collapse = ", "))
    return(data.frame(acc = eq_manual$acc, nombre = eq_manual$nombre, stringsAsFactors = FALSE))
  }
  sel <- sl >= min_len
  if (!is.null(patron_incluir)) sel <- sel & grepl(patron_incluir, names(sl), ignore.case = TRUE)
  if (!is.null(patron_excluir)) sel <- sel & !grepl(patron_excluir, names(sl), ignore.case = TRUE)
  if (!any(sel)) {
    warning("Ninguna secuencia cumple min_len_cromosoma/patrones; se usan TODAS las secuencias.")
    sel <- rep(TRUE, length(sl))
  }
  s <- sl[sel]
  if (orden == "longitud") s <- sort(s, decreasing = TRUE)
  data.frame(acc = names(s), nombre = make.unique(nombre_corto(names(s))), stringsAsFactors = FALSE)
}

# ---------- Genes ----------
tipos_gen_detectar <- function(gff_all, incluir_pseudogenes = TRUE) {
  tp <- unique(as.character(gff_all$type))
  g  <- grep("(^|_)gene$|^pseudogene$", tp, value = TRUE)
  if (!incluir_pseudogenes) g <- g[!grepl("pseudogene", g)]
  g
}

construir_genes <- function(gff_all, tipos_gen = NULL, incluir_pseudogenes = TRUE) {
  if (is.null(tipos_gen)) {
    tipos  <- tipos_gen_detectar(gff_all, incluir_pseudogenes)
    metodo <- "autodetectado"
  } else {
    tipos  <- if (incluir_pseudogenes) tipos_gen else setdiff(tipos_gen, "pseudogene")
    metodo <- "definido en 01_parametros.R"
  }
  g <- gff_all[gff_all$type %in% tipos]
  origen <- paste0("features ", paste(tipos, collapse = " + "), " (", metodo, ")")
  if (length(g) == 0) {                      # plan B: el GFF no tiene features de gen
    tx_tipos <- c("mRNA", "transcript", "ncRNA", "lnc_RNA", "tRNA", "rRNA", "miRNA", "snRNA",
                  "snoRNA", "primary_transcript")
    tx <- gff_all[gff_all$type %in% tx_tipos]
    if (length(tx) > 0) {
      padre <- col_primera(tx, c("Parent", "gene_id", "ID"))
      gid   <- ifelse(is.na(padre), paste0("tx_", seq_along(tx)), padre)
      g <- unlist(range(split(tx, gid)), use.names = TRUE)
      g$ID <- names(g); names(g) <- NULL
      tipos  <- tx_tipos
      origen <- "PLAN B: transcritos agrupados por Parent (el GFF no tiene features de gen)"
    }
  }
  if (length(g) == 0) {
    tt <- table(as.character(gff_all$type))
    stop("No se encontraron genes en el GFF. Tipos de feature disponibles: ",
         paste0(names(tt), " (", as.integer(tt), ")", collapse = ", "),
         "\n  Define `tipos_gen` en 01_parametros.R con los tipos que correspondan a genes.")
  }
  list(genes = g, tipos = tipos, origen = origen)
}

# ---------- Exones e intrones ----------
construir_exones_intrones <- function(gff_file, gff_all, armonizar, gr_vacio, usar_txdb = TRUE) {
  tx_tipos <- c("mRNA", "transcript", "ncRNA", "lnc_RNA", "tRNA", "rRNA", "miRNA", "snRNA",
                "snoRNA", "primary_transcript")
  ex <- NULL; intr <- NULL; metodo <- "TxDb"
  if (usar_txdb) {
    txdb <- tryCatch(suppressWarnings(suppressMessages(makeTxDbFromGFF(gff_file, format = "auto"))),
                     error = function(e) NULL)
    if (!is.null(txdb)) {
      ex   <- tryCatch(armonizar(exons(txdb)), error = function(e) NULL)
      intr <- tryCatch(armonizar(unlist(intronsByTranscript(txdb))), error = function(e) NULL)
    }
  }
  if (is.null(ex) || length(ex) == 0) {
    metodo <- "directo del GFF"
    ex <- armonizar(gff_all[gff_all$type == "exon"])
    if (length(ex) == 0)
      ex <- armonizar(gff_all[gff_all$type %in% c("CDS", "five_prime_UTR", "three_prime_UTR")])
  }
  if (is.null(intr) || length(intr) == 0) {
    tx <- armonizar(gff_all[gff_all$type %in% tx_tipos])
    intr <- if (length(tx) > 0 && length(ex) > 0)
      armonizar(setdiff(reduce(tx, ignore.strand = TRUE), reduce(ex, ignore.strand = TRUE))) else gr_vacio()
  }
  if (length(ex) == 0) ex <- gr_vacio()
  list(exon_red = reduce(ex, ignore.strand = TRUE), intrones = intr, metodo = metodo)
}

# ---------- Repeticiones ----------
leer_repeatmasker <- function(rep_file) {
  rm_df <- read.table(rep_file, skip = 3, header = FALSE, sep = "", fill = TRUE,
                      stringsAsFactors = FALSE, comment.char = "", quote = "",
                      row.names = NULL, check.names = FALSE)
  rm_df$V6 <- suppressWarnings(as.integer(rm_df$V6))
  rm_df$V7 <- suppressWarnings(as.integer(rm_df$V7))
  rm_df    <- rm_df[!is.na(rm_df$V6) & !is.na(rm_df$V7), ]
  rm_df$clase <- sub("\\?$", "", sub("/.*", "", rm_df$V11))
  rm_df
}

rm_vacio <- function()
  data.frame(V2 = numeric(0), V5 = character(0), V6 = integer(0), V7 = integer(0), V9 = character(0),
             V10 = character(0), V11 = character(0), clase = character(0), stringsAsFactors = FALSE)

rm_a_granges <- function(rm_df)
  GRanges(seqnames = rm_df$V5, IRanges(start = rm_df$V6, end = rm_df$V7),
          strand = ifelse(rm_df$V9 == "C", "-", "+"),
          nombre = rm_df$V10, familia = rm_df$V11, clase = rm_df$clase,
          div = suppressWarnings(as.numeric(rm_df$V2)))

# Genomas "soft-masked" (NCBI, Ensembl...): las repeticiones están en minúscula en el FASTA
repeticiones_desde_softmask <- function(fna_file, accs, min_len = 20) {
  bs <- readBStringSet(fna_file)                 # BStringSet conserva mayúsculas/minúsculas
  names(bs) <- sub(" .*", "", names(bs))
  bs <- bs[names(bs) %in% accs]
  res <- lapply(names(bs), function(a) {
    r  <- charToRaw(as.character(bs[[a]]))
    minus <- r >= as.raw(0x61) & r != as.raw(0x6e)        # a-z excepto 'n'
    ir <- as(Rle(minus), "IRanges")
    ir <- ir[width(ir) >= min_len]
    if (length(ir)) data.frame(seqnames = a, start = start(ir), end = end(ir), stringsAsFactors = FALSE)
  })
  do.call(rbind, res)                              # NULL si no hay minúsculas
}

# ---------- Programas externos (opcionales) ----------
correr_trf <- function(fna_file, dir_out, params = c("2", "7", "7", "80", "10", "50", "500")) {
  exe <- Sys.which("trf")
  if (!nzchar(exe)) return(NULL)
  dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
  dir_out <- normalizePath(dir_out)
  out <- file.path(dir_out, paste0(basename(fna_file), ".", paste(params, collapse = "."), ".dat"))
  if (file.exists(out)) return(out)
  fna_abs <- normalizePath(fna_file)
  old <- setwd(dir_out); on.exit(setwd(old), add = TRUE)
  message("Ejecutando TRF (puede tardar)...")
  system2(exe, c(shQuote(fna_abs), params, "-d", "-h"), stdout = FALSE, stderr = FALSE)
  if (file.exists(out)) out else NULL
}

correr_repeatmasker <- function(fna_file, dir_out, species, threads = 4) {
  exe <- Sys.which("RepeatMasker")
  if (!nzchar(exe) || is.null(species)) return(NULL)
  dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
  dir_out <- normalizePath(dir_out)
  out <- file.path(dir_out, paste0(basename(fna_file), ".out"))
  if (file.exists(out)) return(out)
  message("Ejecutando RepeatMasker (puede tardar horas)...")
  system2(exe, c("-pa", as.character(threads), "-species", shQuote(species), "-dir", shQuote(dir_out),
                 shQuote(normalizePath(fna_file))), stdout = FALSE, stderr = FALSE)
  if (file.exists(out)) out else NULL
}


# ---------- Utilidades compartidas por 02L, 06b, 07b (SAT350 / microsatelites) ----------
# "CM065576.1 ... chromosome Ia, ..." -> "chrIa" | "TGME49_chrIX" -> "chrIX" | scaffold -> su accesion
canon_chr <- function(h) {
  h   <- as.character(h)
  tok <- sub("[[:space:]|].*$", "", h)
  g <- function(x, pat) {
    m <- regmatches(x, regexec(pat, x, ignore.case = TRUE, perl = TRUE))
    vapply(m, function(v) if (length(v) > 1) v[2] else NA_character_, character(1))
  }
  m1  <- g(tok, "(?:^|_)chr(?:omosome)?_?([A-Za-z0-9]+)$")
  m2  <- g(h,   "\\bchromosome[[:space:]_]+([A-Za-z0-9]+)")
  lab <- ifelse(!is.na(m1), m1, m2)
  ifelse(!is.na(lab), paste0("chr", lab), tok)
}
norm_lab <- function(x) sub("/.*$", "", as.character(x))        # "VIIb/VIII" (RH, fusionado) -> "VIIb"
rc <- function(s) as.character(reverseComplement(DNAString(s)))

# FASTA solo con los cromosomas analizados (nombres = accesion); lo comparten BLAST y seqkit
asegurar_fasta_chr <- function(fa, dir_out) {
  dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(dir_out, "cromosomas.fasta")
  if (!file.exists(f) || isTRUE(get0("rehacer_externos", ifnotfound = FALSE))) writeXStringSet(fa, f)
  f
}

# TRF sobre un FASTA SIN reutilizar una salida vieja (correr_trf devuelve la existente si la hay)
correr_trf_fresh <- function(fasta, dir_out, params) {
  dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
  dat <- file.path(normalizePath(dir_out), paste0(basename(fasta), ".", paste(params, collapse = "."), ".dat"))
  unlink(dat)
  correr_trf(fasta, dir_out, as.character(params))
}

# Lee la salida .dat de TRF completa (una fila por repeticion; la 1a columna es el id de la secuencia)
parse_trf_full <- function(dat) {
  molde <- data.frame(secuencia = character(), rep_ini = integer(), rep_fin = integer(), periodo = integer(),
                      copias = numeric(), tam_consenso = integer(), pct_matches = integer(), pct_indels = integer(),
                      score = integer(), pct_A = integer(), pct_C = integer(), pct_G = integer(), pct_T = integer(),
                      entropia = numeric(), consenso = character(), secuencia_repetida = character(),
                      stringsAsFactors = FALSE)
  if (is.null(dat) || length(dat) != 1 || is.na(dat) || !file.exists(dat)) return(molde)
  id <- NA_character_; ids <- character(); filas <- list()
  for (l in readLines(dat)) {
    if (startsWith(l, "Sequence:")) id <- trimws(sub("^Sequence:\\s*", "", l))
    else if (grepl("^[0-9]+ [0-9]+ ", l)) {
      f <- strsplit(trimws(l), "\\s+")[[1]]
      if (length(f) >= 15) { filas[[length(filas) + 1]] <- f[1:15]; ids[length(ids) + 1] <- id }
    }
  }
  if (length(filas) == 0) return(molde)
  M <- do.call(rbind, filas)
  data.frame(secuencia = ids, rep_ini = as.integer(M[, 1]), rep_fin = as.integer(M[, 2]), periodo = as.integer(M[, 3]),
             copias = as.numeric(M[, 4]), tam_consenso = as.integer(M[, 5]), pct_matches = as.integer(M[, 6]),
             pct_indels = as.integer(M[, 7]), score = as.integer(M[, 8]), pct_A = as.integer(M[, 9]),
             pct_C = as.integer(M[, 10]), pct_G = as.integer(M[, 11]), pct_T = as.integer(M[, 12]),
             entropia = as.numeric(M[, 13]), consenso = M[, 14], secuencia_repetida = M[, 15], stringsAsFactors = FALSE)
}
