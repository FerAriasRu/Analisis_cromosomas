# ============================================================
# 14_excel.R  -  Un solo Excel con todas las tablas (hoja "Indice" al inicio)
# ------------------------------------------------------------
# REQUIERE   : prácticamente TODAS las tablas de los módulos 01 al 13. Si falla con
#              "objeto no encontrado", corre antes el módulo que lo crea (ver cabeceras).
#              Mínimos: tabla, ventanas_revision, umbrales, resumen_cromosomas, carac_extremos,
#              telo_por_cromosoma, tel_df, ST_genes, sat_bloques, trf_consensos, gb_presencia, gb_hits,
#              rep_*, ext_*, cor_dist, perfil_dist, res_orient, res_dist, resumen_elem, concordancia,
#              tab_ig, pw, inter_all, cor_chr, mh_res, cmp_*, perfil*, cor_tab, het_gc, or_tab,
#              res_perfil, coef_tab, zonas_df, genes_zonas, zonas_rep_df, genes_zonas_rep
# CREA       : tabla_out, parametros, telo_bloques, telo_term, mejor, telo_secuencia, extremos, st_detalle
#              wb, indice, estilo_enc, agregar_hoja(nombre, df, descripcion, ancho_max)
#              archivo_xlsx
# PARA AGREGAR UNA HOJA: llama agregar_hoja("Nombre", df, "descripción") antes de la hoja "Indice".
# EFECTO     : escribe carpeta_salida/analisis_genomico_<ventana>.xlsx
# ============================================================

# ============================================================
# 15. UN SOLO EXCEL CON TODO
# ============================================================

tabla_out <- tabla
num_cols <- c("GC_percent", "AT_percent", "N_pct", "width_noN", "repeat_pct", "repeat_sinTelo_pct", "gene_pct",
              "intergenic_pct", "exon_pct", "intron_pct", "promoter_pct", "inter_rep_pct", "genes_per_Mb",
              "introns_por_gen", "intron_in_gene", "GC_intron", "GC_exon", "GC_inter", "start_Mb", "end_Mb",
              "position_Mb", "pos_rel", "dist_extremo_pb")
tabla_out[num_cols] <- lapply(tabla_out[num_cols], round, 3)
tabla_out$dist_clase <- as.character(tabla_out$dist_clase)

parametros <- data.frame(
  parametro = c("bin", "incluir_pseudogenes", "prom_up", "max_N_pct", "min_ancho_frac", "k_mad", "mad_vars",
                "excluir_outliers_MAD", "excluir_extremos_analisis", "w_excluir_ext", "umbral_pct",
                "min_bp_gen_ventana", "min_bp_gc_comp", "w_telo", "gap_fusion", "min_copias", "unidad_telo",
                "w_st_max", "gap_st", "k_gaps_st", "min_st", "dens_st_max (genes/kb)", "gap_sat", "trf_file",
                "trf_min_copias", "buscar_denovo", "top_n_rep", "max_seq_chars", "ext_umbrales_kb",
                "ext_paso_kb", "ext_max_kb"),
  valor = as.character(c(bin, incluir_pseudogenes, prom_up, max_N_pct, min_ancho_frac, k_mad,
                         paste(mad_vars, collapse = ","), excluir_outliers_MAD, excluir_extremos_analisis,
                         w_excluir_ext, umbral_pct, min_bp_gen_ventana, min_bp_gc_comp, w_telo, gap_fusion,
                         min_copias, paste(unidad_telo, collapse = ","), w_st_max, gap_st, k_gaps_st, min_st, dens_st_max, gap_sat,
                         trf_file, trf_min_copias, buscar_denovo, top_n_rep, max_seq_chars,
                         paste(ext_umbrales_kb, collapse = ","), ext_paso_kb, ext_max_kb)),
  stringsAsFactors = FALSE)

parametros <- rbind(parametros, data.frame(
  parametro = c("genoma_id", "n_cromosomas", "fuente_fai", "fuente_repeticiones", "fuente_trf",
                "tipos_de_gen_usados", "origen_genes", "metodo_exones_intrones", "unidad_telo_usada"),
  valor = c(genoma_id, nrow(eq), fuente_fai, fuente_rep, fuente_trf, paste(tipos_g, collapse = ","),
            origen_genes, metodo_exones, unidad_telo[1]), stringsAsFactors = FALSE))

parametros <- rbind(parametros, data.frame(
  parametro = c("dens_min_len (pb)", "incluir_primer_gen"),
  valor     = as.character(c(dens_min_len, incluir_primer_gen)),
  stringsAsFactors = FALSE))

parametros <- rbind(parametros, data.frame(
  parametro = c("genomas_run", "n_hilos", "rehacer_externos", "gb_evalue", "gb_min_ident (%)", "gb_min_len (pb)",
                "sat350_periodo", "sat350_flanco", "trf_params_sat", "micro_mismatches", "micro_max_amplicon",
                "micro_trf_params", "n_secuencias_sat350"),
  valor = as.character(c(paste(get0("genomas_run_main", ifnotfound = genoma_id), collapse = ","), n_hilos, rehacer_externos, gb_evalue, gb_min_ident,
                         gb_min_len, paste(sat350_periodo, collapse = "-"), sat350_flanco,
                         paste(trf_params_sat, collapse = " "), micro_mismatches, micro_max_amplicon,
                         paste(micro_trf_params, collapse = " "), length(gb_accesiones))),
  stringsAsFactors = FALSE))

telo_bloques <- data.frame(
  chrom = nombre_chr(tel_df$seqnames), accession = tel_df$seqnames, start = tel_df$start, end = tel_df$end,
  longitud_pb = tel_df$longitud, copias = tel_df$copias, ubicacion = tel_df$extremo,
  dist_borde_pb = tel_df$dist_borde, es_telomero = ifelse(tel_df$es_telomero, "Si", "No"),
  secuencia = tel_df$secuencia, hebra = tel_df$hebra, n_motivos = tel_df$n_motivos,
  orientacion_esperada_ok = ifelse(tel_df$orientacion_ok, "Si", "No"), calidad = tel_df$calidad,
  stringsAsFactors = FALSE)




telo_bloques <- telo_bloques[order(factor(telo_bloques$chrom, levels = eq$nombre), telo_bloques$start), ]

telo_term <- telo_bloques[telo_bloques$ubicacion != "interno" & telo_bloques$es_telomero == "Si", ]
telo_term <- telo_term[order(telo_term$chrom, telo_term$ubicacion, -telo_term$n_motivos), ]
mejor <- telo_term[!duplicated(telo_term[, c("chrom", "ubicacion")]), ]

telo_secuencia <- data.frame(chrom = carac_extremos$chrom, extremo = carac_extremos$extremo,
                             tiene_telomero = carac_extremos$tiene_telomero, stringsAsFactors = FALSE)
idx <- match(paste(telo_secuencia$chrom, telo_secuencia$extremo), paste(mejor$chrom, mejor$ubicacion))

telo_secuencia$secuencia_mejor_bloque <- mejor$secuencia[idx]
telo_secuencia$hebra_mejor_bloque <- mejor$hebra[idx]
tel_f <- unidad_telo[1]; tel_r <- rc(tel_f)      # motivo y su complemento inverso
telo_secuencia$lectura <- ifelse(is.na(idx), NA,
                                 ifelse(mejor$hebra[idx] == "+",
                                        paste0(tel_f, " en la hebra + del FASTA"),
                                        paste0(tel_r, " en la hebra + del FASTA (= ", tel_f, " en la hebra complementaria)")))
telo_secuencia$secuencia_esperada <- ifelse(telo_secuencia$extremo == "inicio", tel_r, tel_f)
telo_secuencia$orientacion_esperada_ok <- mejor$orientacion_esperada_ok[idx]
telo_secuencia$n_motivos_mejor <- mejor$n_motivos[idx]
telo_secuencia$dist_borde_pb <- mejor$dist_borde_pb[idx]
telo_secuencia$calidad <- mejor$calidad[idx]
telo_secuencia[[paste0("n_bloques_", tel_f)]] <- mapply(function(ch, ex) sum(telo_term$chrom == ch & telo_term$ubicacion == ex & telo_term$secuencia == tel_f),
                                                   telo_secuencia$chrom, telo_secuencia$extremo)
telo_secuencia[[paste0("n_bloques_", tel_r)]] <- mapply(function(ch, ex) sum(telo_term$chrom == ch & telo_term$ubicacion == ex & telo_term$secuencia == tel_r),
                                                   telo_secuencia$chrom, telo_secuencia$extremo)

rownames(telo_bloques) <- NULL

extremos <- carac_extremos[, c("chrom", "extremo", "tiene_telomero", "telo_start", "telo_end", "telo_long_pb",
                               "telo_copias", "telo_dist_borde_pb", "n_bloques_telo", "ST_detectada",
                               "ST_crit_longitud", "ST_crit_densidad", "ST_motivo_no", "ST_zona_densa_alcanzada",
                               "ST_start", "ST_end", "ST_long_pb", "n_genes_ST", "densidad_ST_genes_kb",
                               "GC_ST", "repeat_pct_ST")]
st_detalle <- carac_extremos[, c("chrom", "extremo", "dist_borde_primer_gen_pb", "gap_max_pb", "gap_medio_pb",
                                 "n_gaps_mayores_umbral", "dens_200kb_genes_kb", "dens_300kb_genes_kb",
                                 "dens_cromosoma_genes_kb", "genes_ST")]

wb <- createWorkbook()
indice <- data.frame(hoja = character(), descripcion = character(), stringsAsFactors = FALSE)
estilo_enc <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", halign = "center", border = "bottom")

agregar_hoja <- function(nombre, df, descripcion = "", ancho_max = 45) {
  if (is.null(df)) return(invisible(NULL))
  nombre <- substr(nombre, 1, 31)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df[] <- lapply(df, function(x) { if (is.factor(x)) x <- as.character(x); if (is.character(x)) substr(x, 1, 32000) else x })
  addWorksheet(wb, nombre)
  writeData(wb, nombre, df, headerStyle = estilo_enc, withFilter = TRUE)
  freezePane(wb, nombre, firstRow = TRUE)
  w <- vapply(seq_along(df), function(j) {
    x <- as.character(df[[j]]); x[is.na(x)] <- ""
    min(ancho_max, max(nchar(c(names(df)[j], head(x, 500)))) + 2)
  }, numeric(1))
  setColWidths(wb, nombre, cols = seq_along(df), widths = w)
  indice <<- rbind(indice, data.frame(hoja = nombre, descripcion = descripcion, stringsAsFactors = FALSE))
}
# -- Generales --
agregar_hoja("Parametros", parametros, "Parámetros con los que se generó el archivo")
agregar_hoja("Resumen_cromosomas", resumen_cromosomas, "Una fila por cromosoma: GC, genes, repetición, telómeros, ST, satélites")
agregar_hoja("Ventanas", tabla_out, "Todas las ventanas con todas las variables, filtros (valida/base/analisis) y banderas")
agregar_hoja("Umbrales_por_cromosoma", umbrales, "Umbrales (percentil) por cromosoma usados en las banderas bajo/alto")
agregar_hoja("Ventanas_excluidas", ventanas_revision, "Ventanas excluidas o marcadas como outlier, con motivo y genes")
# -- Telómeros, ST, satélites --

agregar_hoja("Extremos", extremos, "Telómero y subtélomero por extremo (ST exige longitud y, si mide >= dens_min_len, densidad génica < dens_st_max)")
agregar_hoja("Telomeros_por_cromosoma", telo_por_cromosoma,
             "Presencia de telómero en inicio/final de cada cromosoma, con coordenadas y longitud")
agregar_hoja("Telomeros_secuencia", telo_secuencia,
             "Secuencia (TTTAGGG o CCCTAAA), hebra y orientación del telómero en cada extremo")
agregar_hoja("ST_detalle", st_detalle, "Detalle de intervalos intergénicos y densidades en cada extremo")
agregar_hoja("ST_genes", ST_genes, "Genes dentro de cada tramo subtelomérico, con descripción")
agregar_hoja("Telomeros_bloques", telo_bloques, "Todos los bloques teloméricos encontrados")
agregar_hoja("Satelites_bloques", sat_bloques, "Bloques sat240/sat350 fusionados con su fuente")
agregar_hoja("TRF_consensos", trf_consensos, "Consensos de TRF para sat240/sat350")
agregar_hoja("GB_sat350_presencia", gb_presencia, "¿Están presentes las secuencias GenBank de sat350? pb cubiertos por cromosoma")
agregar_hoja("GB_sat350_hits", gb_hits[, c("query", "chrom", "acc", "start", "end", "strand", "longitud",
                                           "identidad_pct", "score", "cobertura_query_pct", "metodo")],
             "Cada hit de las secuencias GenBank en el genoma (score = bitscore en BLASTN, nº de semillas en k-mers)")
# -- Repeticiones --
agregar_hoja("Rep_top_familias", rep_top, "Top de repeticiones por bp con la secuencia de su copia más larga (genómica, orientada)")
agregar_hoja("Rep_top_ubicacion", rep_comp_top, "Dónde caen (exón/intrón/génico otro/intergénico) las repeticiones del top")
agregar_hoja("Rep_clase_familia", rep_familias, "Todas las clases/familias de RepeatMasker con bp y % del genoma")
agregar_hoja("Rep_clases_resumen", rep_clases_resumen, "Resumen por clase (LTR, LINE, SINE, DNA, Simple_repeat...)")
agregar_hoja("Rep_clases_pct_cromosoma", rep_clases_pct, "% de cada cromosoma ocupado por cada clase (las clases pueden solaparse levemente)")
agregar_hoja("Rep_clases_pb_cromosoma", rep_clases_pb, "pb por clase y cromosoma")
agregar_hoja("Rep_ubicacion_subconjuntos", rep_comp, "Repeticiones por compartimento en: genoma, ventanas alta/baja repetición, extremos, interior")
agregar_hoja("Rep_ubicacion_cromosoma", rep_comp_chr, "Repeticiones por compartimento en cada cromosoma")
agregar_hoja("Rep_ubicacion_clase", rep_comp_clase, "Repeticiones por compartimento para cada clase")
agregar_hoja("Rep_hits_resumen", rep_hits_res, "Nº de hits de RepeatMasker por clase y compartimento dominante")
if (exportar_hits_rep && !is.null(rep_hits) && nrow(rep_hits) < 1e6) agregar_hoja("Rep_hits", rep_hits, "Cada hit de RepeatMasker con pb en exón/intrón/génico otro/intergénico")
# -- Extremos --
agregar_hoja("Ext_vs_interior", ext_vs_int, "¿Más repetición y menos genes en extremos? comparación con varias definiciones de extremo")
agregar_hoja("Ext_vs_interior_chr", ext_chr, "Lo mismo por cromosoma, con FDR y veredicto del patrón")
agregar_hoja("Ext_correlacion_distancia", cor_dist, "Spearman de cada variable vs distancia al extremo")
agregar_hoja("Ext_perfil_distancia", perfil_dist, "Perfil por clases de distancia al extremo")
# -- Repetición alta vs baja --
agregar_hoja("Repeticion_alta_vs_baja", cmp_repeticion, "Ventanas de ALTA vs BAJA repetición (y vs resto): intrones, exones, genes, intergénico...")
agregar_hoja("Repeticion_perfil_grupos", perfil_rep, "Medianas por grupo de repetición")
# -- Intergénico --
agregar_hoja("Inter_orientacion", res_orient, "Regiones intergénicas por orientación de los genes flanqueantes")
agregar_hoja("Inter_distribucion", res_dist, "Distribución de longitudes intergénicas")
agregar_hoja("Inter_GC_por_elemento", resumen_elem, "GC de exones, intrones e intergénico")
agregar_hoja("Inter_concordancia", concordancia, "Baja densidad génica vs alto % intergénico (rho, Fisher, Jaccard)")
agregar_hoja("Inter_contingencia", as.data.frame(tab_ig), "Tabla de contingencia baja densidad vs alto intergénico")
if (!is.null(pw)) agregar_hoja("Inter_pairwise_p", data.frame(grupo = rownames(pw$p.value), round(pw$p.value, 4), check.names = FALSE),
                               "p ajustados (BH) de longitud intergénica entre orientaciones")
agregar_hoja("Inter_regiones", inter_all, "Cada región/par intergénico con orientación, GC, repetición")
# -- Estadística --
agregar_hoja("Est_cor_GC_genes_chr", cor_chr, "Spearman GC vs genes/Mb por cromosoma")
agregar_hoja("Est_mantel_haenszel", mh_res, "Mantel-Haenszel bajo GC vs baja densidad génica")
agregar_hoja("Est_bajoGC_vs_resto", cmp_bajoGC, "Bajo GC vs resto")
agregar_hoja("Est_bajaDens_vs_resto", cmp_bajoDen, "Baja densidad génica vs resto")
agregar_hoja("Est_altoIntron_vs_resto", cmp_intron, "Alto % intrones vs resto")
agregar_hoja("Est_perfil_bajoGC", perfil, "Medianas bajo GC vs resto")
agregar_hoja("Est_correlaciones", cor_tab, "Spearman de variables vs GC, genes/Mb y % repeticiones")
agregar_hoja("Est_heterogeneidad_GC", het_gc, "Heterogeneidad del GC por cromosoma")
agregar_hoja("Est_regresion_logistica", or_tab, "Regresión logística de bajo GC")
agregar_hoja("Est_composicion_bajoGC", res_perfil, "Composición de las ventanas de bajo GC")
agregar_hoja("Est_modelo_GC", coef_tab, "Modelo lineal del GC de la ventana")
# -- Zonas --
agregar_hoja("Zonas_criticas", zonas_df, "Zonas contiguas: bajo GC y/o baja densidad génica")
agregar_hoja("Zonas_genes", genes_zonas, "Genes dentro de cada zona crítica")
agregar_hoja("Zonas_baja_dens_alta_rep", zonas_rep_df, "Zonas con baja densidad génica Y alta repetición")
agregar_hoja("Zonas_baja_dens_alta_rep_genes", genes_zonas_rep, "Genes dentro de las zonas de baja densidad + alta repetición")

# -- SAT350 (GenBank/BLASTN/TRF) y microsatelites (antes: Sat350.R y microsatelites.R) --
agregar_hoja("Sat350_secuencias", sat350_refs_tab, "Secuencias de referencia SAT350: longitud en GenBank vs articulo")
agregar_hoja("Sat350_resumen_query", sat350_resumen_q, "Hits de BLASTN por secuencia de referencia")
agregar_hoja("Sat350_TRF_referencia", sat350_trf_ref, "TRF sobre las secuencias de referencia SAT350")
agregar_hoja("Sat350_regiones", sat350_regiones, "Regiones SAT350 fusionadas (gap_sat) con repeticion dominante de TRF y periodo del monomero")
agregar_hoja("Sat350_TRF_regiones", sat350_trf_regiones, "Todas las repeticiones TRF dentro de las regiones SAT350")
agregar_hoja("Micro_resumen_loci", micro_resumen, "Un locus por fila: amplicon, repeticion principal de TRF, motivo y comparacion con el articulo")
agregar_hoja("Micro_contexto", microsat_loci, "Microsatelites: posicion, distancia al borde, en subtelomero, dentro de SAT350")
agregar_hoja("Micro_TRF", micro_trf, "Todas las repeticiones TRF dentro de los amplicones")
agregar_hoja("Micro_amplicones", micro_amps, "Amplicones encontrados por seqkit")
agregar_hoja("Micro_primers", primers_micro, "Tabla de primers usada")

addWorksheet(wb, "Indice")
writeData(wb, "Indice", indice, headerStyle = estilo_enc)
setColWidths(wb, "Indice", cols = 1:2, widths = c(32, 110))
worksheetOrder(wb) <- c(length(names(wb)), seq_len(length(names(wb)) - 1))
archivo_xlsx <- paste0(carpeta_salida, "analisis_genomico_", ventana_label, ".xlsx")
saveWorkbook(wb, archivo_xlsx, overwrite = TRUE)
cat("\nExcel único guardado en:", archivo_xlsx, "(", nrow(indice), "hojas )\n")
