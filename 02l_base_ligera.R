# ============================================================
# 02L_base_ligera.R  -  Base minima para genomas SIN GFF (modo "ligero")
# ------------------------------------------------------------
# REQUIERE   : 01 (fna_file, fai_file, min_len_cromosoma, patron_incluir, patron_excluir, eq_manual,
#              orden_cromosomas), entradas_genoma.R
# CREA       : eq, sl, fa, n_filas, nombre_chr(), armonizar()   (mismas definiciones que 02_datos_base.R)
# NOMBRES    : se deducen del encabezado completo del FASTA ("... chromosome Ia ..." o "TGME49_chrIa" -> "Ia");
#              lo que no sea cromosoma conserva su accesion. Con eq_manual se usan los nombres manuales.
# ============================================================
source(file.path(dir_scripts, "entradas_genoma.R"), encoding = "UTF-8")
stopifnot(file.exists(fna_file))
fa_all <- readDNAStringSet(fna_file)
hdr <- names(fa_all)
names(fa_all) <- sub(" .*", "", hdr)
if (anyDuplicated(names(fa_all))) stop("El FASTA tiene nombres de secuencia duplicados.")
sl_all <- setNames(as.numeric(width(fa_all)), names(fa_all))

eq <- seleccionar_cromosomas(sl_all, min_len_cromosoma, patron_incluir, patron_excluir, eq_manual, orden_cromosomas)
if (is.null(eq_manual)) {
  lab <- canon_chr(hdr)[match(eq$acc, names(fa_all))]
  eq$nombre <- make.unique(ifelse(grepl("^chr", lab), sub("^chr", "", lab), eq$acc))
}
sl <- sl_all[eq$acc]; fa <- fa_all[eq$acc]; rm(fa_all)
n_filas <- length(sl)
cat("Cromosomas analizados (modo ligero):", n_filas, "| longitud total:", round(sum(sl) / 1e6, 2), "Mb\n")
print(data.frame(eq, long_Mb = round(sl / 1e6, 3), row.names = NULL))

nombre_chr <- function(acc) eq$nombre[match(acc, eq$acc)]
armonizar <- function(gr) {
  gr <- gr[as.character(seqnames(gr)) %in% names(sl)]
  seqlevels(gr, pruning.mode = "coarse") <- names(sl)
  tryCatch(seqlengths(gr) <- sl, error = function(e) NULL)
  gr
}
