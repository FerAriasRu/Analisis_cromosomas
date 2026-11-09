# ============================================================
# 02_datos_base.R  -  Genoma, cromosomas, anotación (GFF), N y repeticiones
# ------------------------------------------------------------
# REQUIERE   : 01_parametros.R (fna_file, gff_file, fai_file, rep_file, trf_file, selección de
#              cromosomas, tipos_gen, attr_*, mapa_gff_a_fasta, correr_*, usar_softmask, prom_up...)
#              + entradas_genoma.R (funciones auxiliares, misma carpeta)
# CREA (datos):
#   eq, sl, n_filas, max_length_Mb, lab      cromosomas analizados (acc, nombre) y longitudes
#   fa            DNAStringSet de los cromosomas analizados
#   gff_all, genes_gr (gene_id, gene_nombre, descripcion, biotipo, gene_lab), genes_red, tipos_g, origen_genes
#   intrones_red, exon_red, metodo_exones, prom, prom_red
#   N_df, N_red   regiones de N (huecos del ensamblaje)
#   rm_df, rep_gr (nombre, familia, clase, div), rep_red, tiene_rep, fuente_rep   <- repeticiones
#   fuente_fai, fuente_trf, trf_file (queda NA si no hay TRF)
# CREA (funciones): nombre_chr(acc), armonizar(gr), gr_vacio()
#   (col_txt, col_primera y restar vienen de entradas_genoma.R)
# FUENTE DE REPETICIONES (en este orden): RepeatMasker (archivo o ejecutado) > minúsculas del FASTA
#   (softmask) > ninguna. Con "ninguna" los módulos de repeticiones se omiten sin romper el pipeline.
# USADO POR  : prácticamente todos los módulos.
# ============================================================

source(file.path(dir_scripts, "entradas_genoma.R"), encoding = "UTF-8")

# ============================================================
# 2.1 GENOMA Y CROMOSOMAS (el .fai es opcional)
# ============================================================
stopifnot(file.exists(fna_file), file.exists(gff_file))
fa_all <- readDNAStringSet(fna_file)
names(fa_all) <- sub(" .*", "", names(fa_all))
if (anyDuplicated(names(fa_all))) stop("El FASTA tiene nombres de secuencia duplicados.")

if (!is.null(fai_file) && file.exists(fai_file)) {
  fai    <- read.table(fai_file, sep = "\t", quote = "", stringsAsFactors = FALSE)[, 1:2]
  sl_all <- setNames(as.numeric(fai$V2), fai$V1)
  stopifnot(all(names(sl_all) %in% names(fa_all)))
  if (any(sl_all != width(fa_all)[match(names(sl_all), names(fa_all))]))
    stop("El .fai no coincide con el FASTA (¿es de otra versión del genoma?).")
  fuente_fai <- "archivo .fai"
} else {
  if (!is.null(fai_file)) message("No se encontró fai_file: las longitudes se calculan desde el FASTA.")
  sl_all <- setNames(as.numeric(width(fa_all)), names(fa_all))
  fuente_fai <- "calculado desde el FASTA"
}

eq <- seleccionar_cromosomas(sl_all, min_len_cromosoma, patron_incluir, patron_excluir,
                             eq_manual, orden_cromosomas)
sl <- sl_all[eq$acc]
fa <- fa_all[eq$acc]
rm(fa_all)
n_filas <- length(sl)
max_length_Mb <- max(sl) / 1e6
cat("Cromosomas analizados:", n_filas, "| longitud total:", round(sum(sl) / 1e6, 2), "Mb | .fai:", fuente_fai, "\n")
print(data.frame(eq, long_Mb = round(sl / 1e6, 3), row.names = NULL))
if (n_filas > 40) message("Aviso: ", n_filas, " secuencias analizadas; sube `min_len_cromosoma` si son scaffolds.")

lab <- setNames(paste0(eq$nombre, "\n", format(round(sl / 1e6, 2), nsmall = 2), " Mb"), eq$acc)
nombre_chr <- function(acc) eq$nombre[match(acc, eq$acc)]

# ============================================================
# 2.2 FUNCIONES QUE DEPENDEN DE `sl`
# ============================================================
armonizar <- function(gr) {
  gr <- gr[as.character(seqnames(gr)) %in% names(sl)]
  seqlevels(gr, pruning.mode = "coarse") <- names(sl)
  tryCatch(seqlengths(gr) <- sl, error = function(e) NULL)
  gr
}
gr_vacio <- function() {
  g <- GRanges(factor(character(0), levels = names(sl)), IRanges(integer(0), integer(0)))
  seqlengths(g) <- sl
  g
}

# ============================================================
# 2.3 GENES (autodetecta el tipo de feature y los atributos del GFF)
# ============================================================
gff_all <- import(gff_file)
if (!is.null(mapa_gff_a_fasta)) {
  m <- mapa_gff_a_fasta[names(mapa_gff_a_fasta) %in% seqlevels(gff_all)]
  gff_all <- renameSeqlevels(gff_all, m)
}
comunes <- intersect(as.character(seqlevelsInUse(gff_all)), names(sl))
if (length(comunes) == 0)
  stop("Ningún nombre de secuencia del GFF coincide con el FASTA.\n  GFF   : ",
       paste(head(seqlevels(gff_all), 8), collapse = ", "), "\n  FASTA : ",
       paste(head(names(sl), 8), collapse = ", "),
       "\n  Usa `mapa_gff_a_fasta` en 01_parametros.R para traducir los nombres.")

res_g <- construir_genes(gff_all, tipos_gen, incluir_pseudogenes)
tipos_g      <- res_g$tipos
origen_genes <- res_g$origen
genes_gr <- armonizar(res_g$genes)
if (length(genes_gr) == 0) stop("Hay genes en el GFF pero ninguno cae en los cromosomas analizados.")
genes_gr <- sort(genes_gr)
genes_gr$gene_id <- col_primera(genes_gr, attr_id)
sin_id <- is.na(genes_gr$gene_id)
if (any(sin_id)) genes_gr$gene_id[sin_id] <- paste0("gen_", which(sin_id))
if (anyDuplicated(genes_gr$gene_id)) {
  message("Hay gene_id repetidos; se hacen únicos con make.unique().")
  genes_gr$gene_id <- make.unique(genes_gr$gene_id)
}
genes_gr$gene_nombre <- col_primera(genes_gr, attr_nombre)
genes_gr$descripcion <- col_primera(genes_gr, attr_desc)
bt <- col_primera(genes_gr, c("gene_biotype", "biotype"))
genes_gr$biotipo     <- ifelse(is.na(bt), col_txt(genes_gr, "type"), bt)
genes_gr$gene_lab    <- ifelse(!is.na(genes_gr$gene_nombre) & genes_gr$gene_nombre != "",
                               genes_gr$gene_nombre, genes_gr$gene_id)
cat("Genes usados:", length(genes_gr), "|", origen_genes, "\n")

# ---- Intrones y exones ----
ei <- construir_exones_intrones(gff_file, gff_all, armonizar, gr_vacio, usar_txdb = is.null(mapa_gff_a_fasta))
exon_red      <- ei$exon_red
intrones_red  <- reduce(ei$intrones, ignore.strand = TRUE)
metodo_exones <- ei$metodo
cat("Exones:", length(exon_red), "| intrones:", length(intrones_red), "| método:", metodo_exones, "\n")

# ---- Promotores ----
prom     <- armonizar(suppressWarnings(trim(promoters(genes_gr, upstream = prom_up, downstream = 0))))
prom_red <- reduce(prom, ignore.strand = TRUE)

genes_red <- reduce(genes_gr, ignore.strand = TRUE)

# ---- Regiones de N (huecos del ensamblaje) ----
N_df <- do.call(rbind, lapply(names(sl), function(acc) {
  m <- gregexpr("N+", as.character(fa[[acc]]))[[1]]
  if (m[1] == -1) return(NULL)
  data.frame(seqnames = acc, start = as.integer(m), end = as.integer(m) + attr(m, "match.length") - 1L,
             stringsAsFactors = FALSE)
}))
N_red <- if (is.null(N_df)) gr_vacio() else
  reduce(armonizar(GRanges(N_df$seqnames, IRanges(N_df$start, N_df$end))), ignore.strand = TRUE)
cat("Regiones de N:", length(N_red), "| total N (pb):", sum(width(N_red)), "\n")

# ============================================================
# 2.4 REPETICIONES: RepeatMasker > softmask del FASTA > ninguna
# ============================================================
dir.create(carpeta_externos, showWarnings = FALSE, recursive = TRUE)
rm_df <- rm_vacio(); fuente_rep <- "ninguna"

if (!is.null(rep_file) && file.exists(rep_file)) {
  rm_df <- leer_repeatmasker(rep_file); fuente_rep <- "RepeatMasker (archivo)"
} else {
  if (!is.null(rep_file)) message("No se encontró rep_file: ", rep_file)
  if (ejecutar_repeatmasker) {
    f <- correr_repeatmasker(fna_file, carpeta_externos, rm_species, rm_threads)
    if (!is.null(f)) { rep_file <- f; rm_df <- leer_repeatmasker(f); fuente_rep <- "RepeatMasker (ejecutado)" }
    else message("RepeatMasker no está disponible o falló (¿definiste rm_species?).")
  }
}
if (nrow(rm_df) > 0) {
  rep_gr <- armonizar(rm_a_granges(rm_df))
} else {
  rep_gr <- gr_vacio()
  if (usar_softmask) {
    sm <- repeticiones_desde_softmask(fna_file, names(sl))
    if (!is.null(sm) && nrow(sm) > 0) {
      rep_gr <- armonizar(GRanges(sm$seqnames, IRanges(sm$start, sm$end), strand = "*",
                                  nombre = "softmask", familia = "Unknown", clase = "Unknown", div = NA_real_))
      fuente_rep <- "minúsculas (softmask) del FASTA"
    }
  }
}
if (length(rep_gr) == 0) {
  rep_gr$nombre <- character(0); rep_gr$familia <- character(0)
  rep_gr$clase  <- character(0); rep_gr$div <- numeric(0)
  message("Sin información de repeticiones: se omiten los análisis de repeticiones.")
}
rep_red   <- reduce(rep_gr, ignore.strand = TRUE)
tiene_rep <- length(rep_gr) > 0
cat("Repeticiones:", length(rep_gr), "hits | fuente:", fuente_rep, "\n")

# ---- TRF (solo se necesita para satélites) ----
fuente_trf <- "no disponible"
trf_ok <- !is.null(trf_file) && !is.na(trf_file) && file.exists(trf_file)
if (trf_ok) {
  fuente_trf <- "TRF (archivo)"
} else if (length(periodos_sat) > 0 && ejecutar_trf) {
  f <- correr_trf(fna_file, carpeta_externos)
  if (!is.null(f)) { trf_file <- f; trf_ok <- TRUE; fuente_trf <- "TRF (ejecutado)" }
}
if (!trf_ok) {
  trf_file <- NA_character_
  if (length(periodos_sat) > 0) message("Sin salida de TRF (ni programa instalado): satélites solo vía RepeatMasker/GenBank/de novo.")
}
