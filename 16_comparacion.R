# ============================================================
# 16_comparacion.R  -  Tablas comparativas entre genomas (SIN graficas: las dibuja 15_figuras.R)
# ------------------------------------------------------------
# REQUIERE   : 01 (genomas_comparar, dir_salida, w_st_max), main (genomas_run_main),
#              objetos/objetos_<id>.rds de cada genoma (los crea 90_guardar_genoma.R)
# GENOMAS    : genomas_comparar (NULL = los corridos ahora; o c("RH","RH88") = los que ya tengan .rds)
# CREA       : cmp_ok (TRUE si hay >= 2 genomas), cmp_nombres, carpeta_cmp, niveles_chr, col_gen
#              cmp_long, cmp_sat350, cmp_micro, cmp_extremos, cmp_chr, cmp_sat_bloques, cmp_resumen
#              <dir_salida>/comparacion_genomas/comparacion_genomas.xlsx
# NOTA       : "VIIb/VIII" (RH, cromosoma fusionado) se normaliza a "VIIb" para alinearlo con los demas.
# ============================================================
cmp_ok <- FALSE
dir_obj <- paste0(dir_salida, "objetos/")
cmp_nombres <- unique(if (is.null(genomas_comparar)) get0("genomas_run_main", ifnotfound = character(0)) else genomas_comparar)
ventana_cmp <- if (bin >= 1e6) {
  paste0(bin / 1e6, "mb")
} else {
  paste0(bin / 1e3, "kb")
}

archivos <- file.path(
  dir_salida,
  paste0("resultados_", cmp_nombres, "_", ventana_cmp),
  "objetos",
  paste0("objetos_", cmp_nombres, ".rds")
)

falta <- !file.exists(archivos)
if (any(falta)) message("Comparacion: sin objetos para ", paste(cmp_nombres[falta], collapse = ", "),
                        " (corre primero ese genoma).")
cmp_nombres <- cmp_nombres[!falta]; archivos <- archivos[!falta]

if (length(cmp_nombres) < 2) {
  message("Comparacion omitida: se necesitan >= 2 genomas con objetos guardados (hay ", length(cmp_nombres), ").")
} else {
  cmp_ok <- TRUE
  carpeta_cmp <- paste0(dir_salida, "comparacion_genomas/")
  dir.create(carpeta_cmp, showWarnings = FALSE, recursive = TRUE)
  objs <- setNames(lapply(archivos, readRDS), cmp_nombres)
  apilar <- function(f) { x <- do.call(rbind, lapply(names(objs), function(n) f(objs[[n]], n))); if (!is.null(x)) rownames(x) <- NULL; x }

  cmp_zonas <- apilar(function(o, n) {
    z <- o$zonas
    if (is.null(z) || nrow(z) == 0) return(NULL)
    z$chrom <- norm_lab(z$chrom)
    cbind(genoma = n, z)
  })
  
  cmp_genes <- apilar(function(o, n) {
    g <- o$genes
    if (is.null(g) || nrow(g) == 0) return(NULL)
    g$chrom <- norm_lab(g$chrom)
    cbind(genoma = n, g)
  })
  
  cmp_repeticiones <- apilar(function(o, n) {
    r <- o$repeticiones
    if (is.null(r) || nrow(r) == 0) return(NULL)
    r$chrom <- norm_lab(r$chrom)
    cbind(genoma = n, r)
  })
  
  cmp_long <- apilar(function(o, n) data.frame(genoma = n, chrom = norm_lab(o$longitudes$chrom),
                                               accession = o$longitudes$accession, long_pb = o$longitudes$long_pb,
                                               stringsAsFactors = FALSE))
  cmp_sat350 <- apilar(function(o, n) { r <- o$sat350_regiones; if (is.null(r) || nrow(r) == 0) return(NULL)
    data.frame(genoma = n, chrom = norm_lab(r$chrom), accession = r$accession, ini = r$start, fin = r$end,
               n_hits = r$n_hits, mejor_pident = r$mejor_pident, periodo = r$periodo,
               monomero_sat350 = r$monomero_sat350, stringsAsFactors = FALSE) })
  cmp_micro <- apilar(function(o, n) { m <- o$microsat_loci; if (is.null(m) || nrow(m) == 0) return(NULL)
    cbind(genoma = n, chrom_n = norm_lab(m$chrom), m[, setdiff(names(m), "chrom")], stringsAsFactors = FALSE) })
  if (!is.null(cmp_micro)) names(cmp_micro)[names(cmp_micro) == "chrom_n"] <- "chrom"
  cmp_extremos <- apilar(function(o, n) { e <- o$carac_extremos; if (is.null(e)) return(NULL)
   
  data.frame(
    genoma = n,
    chrom = norm_lab(e$chrom),
    extremo = e$extremo,
    tiene_telomero = e$tiene_telomero,
    telo_long_pb = e$telo_long_pb,
    ST_detectada = e$ST_detectada,
    ST_start = e$ST_start,
    ST_end = e$ST_end,
    ST_long_pb = e$ST_long_pb,
    n_genes_ST = e$n_genes_ST,
    densidad_ST_genes_kb = e$densidad_ST_genes_kb,
    stringsAsFactors = FALSE
  ) })
  cmp_chr <- apilar(function(o, n) { r <- o$resumen_cromosomas; if (is.null(r)) return(NULL)
    r$chrom <- norm_lab(r$chrom); cbind(genoma = n, r, stringsAsFactors = FALSE) })
  cmp_sat_bloques <- apilar(function(o, n) { b <- o$sat_bloques; if (is.null(b) || nrow(b) == 0) return(NULL)
    b$chrom <- norm_lab(b$chrom); cbind(genoma = n, b, stringsAsFactors = FALSE) })

  # ---- Resumen por genoma ----
  cmp_resumen <- do.call(rbind, lapply(names(objs), function(n) {
    o <- objs[[n]]; r <- o$resumen_cromosomas; ce <- o$carac_extremos
    lon <- cmp_long[cmp_long$genoma == n, ]; tot <- sum(lon$long_pb)
    s <- cmp_sat350[cmp_sat350$genoma == n, ]; kb_s <- sum(s$fin - s$ini + 1) / 1e3
    Ls <- lon$long_pb[match(s$chrom, lon$chrom)]
    w <- function(v) if (!is.null(r)) round(sum(r[[v]] * r$long_pb, na.rm = TRUE) / sum(r$long_pb), 2) else NA
    data.frame(genoma = n, modo = o$modo, n_cromosomas = nrow(lon), tamano_Mb = round(tot / 1e6, 2),
               GC_pct = w("GC_percent"), repeat_pct = w("repeat_pct"), n_genes = if (!is.null(r)) sum(r$n_genes) else NA,
               extremos_con_telomero = if (!is.null(ce)) sum(ce$tiene_telomero == "Si") else NA,
               extremos_con_ST_valido = if (!is.null(ce)) sum(ce$ST_detectada == "Si") else NA,
               n_regiones_SAT350 = nrow(s), kb_SAT350 = round(kb_s, 1), pct_genoma_SAT350 = round(kb_s * 1e3 / tot * 100, 3),
               SAT350_en_extremos = sum(s$ini <= w_st_max | s$fin >= Ls - w_st_max + 1, na.rm = TRUE),
               n_microsat = sum(cmp_micro$genoma == n), stringsAsFactors = FALSE)
  }))
  cat("\n=== Comparacion entre genomas ===\n"); print(cmp_resumen)

  orden_toxo  <- c("Ia", "Ib", "II", "III", "IV", "V", "VI", "VIIa", "VIIb", "VIII", "IX", "X", "XI", "XII")
  pres        <- unique(cmp_long$chrom)
  niveles_chr <- c(orden_toxo[orden_toxo %in% pres], sort(setdiff(pres, orden_toxo)))
  col_gen     <- setNames(scales::hue_pal()(length(cmp_nombres)), cmp_nombres)

  # ---- Excel de la comparacion ----
  wb_c <- createWorkbook(); est_c <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "bottom")
  hoja_c <- function(nm, df) if (!is.null(df) && nrow(df) > 0) {
    addWorksheet(wb_c, nm); writeData(wb_c, nm, df, headerStyle = est_c, withFilter = TRUE)
    freezePane(wb_c, nm, firstRow = TRUE); setColWidths(wb_c, nm, cols = seq_len(ncol(df)), widths = "auto") }
  hoja_c("Resumen_genomas", cmp_resumen);  hoja_c("Longitudes", cmp_long)
  hoja_c("SAT350_regiones", cmp_sat350);   hoja_c("Microsatelites", cmp_micro)
  hoja_c("Extremos_ST_telomero", cmp_extremos); hoja_c("Resumen_cromosomas", cmp_chr)
  hoja_c("Satelites_bloques", cmp_sat_bloques)
  saveWorkbook(wb_c, paste0(carpeta_cmp, "comparacion_genomas.xlsx"), overwrite = TRUE)
  cat("Excel de comparacion:", paste0(carpeta_cmp, "comparacion_genomas.xlsx"), "\n")
}
