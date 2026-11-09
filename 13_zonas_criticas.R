# ============================================================
# 13_zonas_criticas.R  -  Zonas contiguas con bajo GC y/o baja densidad génica
# ------------------------------------------------------------
# REQUIERE   : 01 (umbral_pct), 02 (sl, eq, genes_gr, nombre_chr), 04 (bins, tabla),
#              06 (st_gr), 09 (banderas bajo_GC, baja_dens_genica, alto_repeat)
# CREA       : niveles_zona, colores_zona
#              crear_zonas(flag, tipo), resumir_zonas(zonas, prefijo)   funciones
#              lista_z, zonas_df, genes_zonas
#              zr_gr, zonas_rep_df, genes_zonas_rep      (baja densidad + alta repetición)
# USADO POR  : 14_excel (hojas Zonas_*), 15_figuras (zonas_df, niveles_zona, colores_zona)
# ============================================================

# ============================================================
# 14. ZONAS CRÍTICAS (ventanas contiguas con la misma combinación de banderas)
# ============================================================

niveles_zona <- c("Bajo GC + baja densidad génica", "Bajo GC", "Baja densidad génica")
colores_zona <- c("Bajo GC + baja densidad génica" = "#D55E00", "Bajo GC" = "#56B4E9",
                  "Baja densidad génica" = "#009E73")

crear_zonas <- function(flag, tipo) {
  if (!any(flag)) return(NULL)
  z <- reduce(bins[flag]); mcols(z) <- NULL; z$tipo <- tipo; z
}

resumir_zonas <- function(zonas, prefijo) {
  zch <- nombre_chr(as.character(seqnames(zonas)))
  zonas$zona_id <- paste0(zch, prefijo, ave(seq_along(zonas), zch, FUN = seq_along))
  ov  <- findOverlaps(zonas, bins)
  idx <- split(subjectHits(ov), factor(queryHits(ov), levels = seq_along(zonas)))
  est <- function(col, f) unname(vapply(idx, function(i) f(tabla[[col]][i], na.rm = TRUE), numeric(1)))
  ovg <- findOverlaps(zonas, genes_gr, ignore.strand = TRUE)
  qz <- queryHits(ovg); sg <- subjectHits(ovg)
  genes_z <- data.frame(
    chrom = zch[qz], zona_id = zonas$zona_id[qz], tipo = zonas$tipo[qz],
    zona_start = start(zonas)[qz], zona_end = end(zonas)[qz],
    gene_id = genes_gr$gene_id[sg], gene_nombre = genes_gr$gene_nombre[sg],
    descripcion = genes_gr$descripcion[sg], biotipo = genes_gr$biotipo[sg],
    gene_start = start(genes_gr)[sg], gene_end = end(genes_gr)[sg],
    strand = as.character(strand(genes_gr))[sg], stringsAsFactors = FALSE)
  genes_z <- genes_z[order(factor(genes_z$chrom, levels = eq$nombre), genes_z$zona_start, genes_z$gene_start), ]
  lista_genes <- tapply(genes_gr$gene_lab[sg], factor(qz, levels = seq_along(zonas)), paste, collapse = "; ")
  Lz <- as.numeric(sl[as.character(seqnames(zonas))])
  zdf <- data.frame(
    chrom = zch, zona_id = zonas$zona_id, tipo = zonas$tipo,
    region = paste0(zch, ":", format(start(zonas), big.mark = ","), "-", format(end(zonas), big.mark = ",")),
    start_pb = start(zonas), end_pb = end(zonas),
    start_Mb = round(start(zonas) / 1e6, 3), end_Mb = round(end(zonas) / 1e6, 3),
    longitud_kb = width(zonas) / 1e3, n_ventanas = unname(lengths(idx)),
    dist_al_extremo_kb = round(pmin(start(zonas) - 1, Lz - end(zonas)) / 1e3, 1),
    GC_medio = round(est("GC_percent", mean), 2), AT_medio = round(est("AT_percent", mean), 2),
    genes_per_Mb_medio = round(est("genes_per_Mb", mean), 1), genes_per_Mb_max = round(est("genes_per_Mb", max), 1),
    repeat_pct_medio = round(est("repeat_pct", mean), 1), intron_pct_medio = round(est("intron_pct", mean), 1),
    intergenic_pct_medio = round(est("intergenic_pct", mean), 1),
    solapa_subtelomero = ifelse(overlapsAny(zonas, st_gr, ignore.strand = TRUE), "Si", "No"),
    n_genes = tabulate(qz, nbins = length(zonas)), genes = unname(lista_genes), stringsAsFactors = FALSE)
  zdf$genes[is.na(zdf$genes)] <- ""
  zdf$chrom <- factor(zdf$chrom, levels = eq$nombre)
  zdf <- zdf[order(zdf$chrom, zdf$start_pb), ]; rownames(zdf) <- NULL
  list(zonas = zdf, genes = genes_z)
}

lista_z <- Filter(Negate(is.null), list(
  crear_zonas(tabla$bajo_GC &  tabla$baja_dens_genica, niveles_zona[1]),
  crear_zonas(tabla$bajo_GC & !tabla$baja_dens_genica, niveles_zona[2]),
  crear_zonas(!tabla$bajo_GC & tabla$baja_dens_genica, niveles_zona[3])))
if (length(lista_z) == 0) stop("No hay zonas críticas con estos umbrales: sube `umbral_pct`.")
rz <- resumir_zonas(sort(do.call(c, lista_z)), "_z")
zonas_df <- rz$zonas; genes_zonas <- rz$genes
zonas_df$tipo <- factor(zonas_df$tipo, levels = niveles_zona)
cat("\nZonas críticas por cromosoma y tipo:\n"); print(table(zonas_df$chrom, zonas_df$tipo))

# Zonas de baja densidad génica + alta repetición
zr_gr <- crear_zonas(tabla$baja_dens_genica & tabla$alto_repeat, "Baja densidad génica + alta repetición")
if (!is.null(zr_gr)) {
  rzr <- resumir_zonas(sort(zr_gr), "_rep")
  zonas_rep_df <- rzr$zonas; genes_zonas_rep <- rzr$genes
} else {
  zonas_rep_df <- data.frame(nota = "Sin ventanas con baja densidad génica y alta repetición")
  genes_zonas_rep <- data.frame(nota = "Sin genes")
}


#########Tabla de las zonas criticas segun variables graficadas ########
#GC (%) por cromosoma y densidad de genes.
#Porcentaje de repeticiones y densidad de genes.
#Porcentaje intergénico y densidad de genes.
#Porcentaje de intrones y densidad de genes.
########################################################################


tabla_zonas <- zonas_df[
  as.character(zonas_df$tipo) == "Bajo GC + baja densidad génica",
  c("chrom", "zona_id", "region", "start_pb", "end_pb", "longitud_kb",
    "n_ventanas", "GC_medio", "genes_per_Mb_medio",
    "repeat_pct_medio", "intergenic_pct_medio", "intron_pct_medio")
]

View(tabla_zonas)
####################### Genes en zonas criticas
View(genes_zonas[genes_zonas$zona_id == "Ia_z6", ])
genes_zonas[
  genes_zonas$gene_id == "TgRH_000018100",
  c("chrom", "zona_id", "gene_id", "gene_nombre", "descripcion", "biotipo")
]
