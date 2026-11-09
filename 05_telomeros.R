# ============================================================
# 05_telomeros.R  -  Detección de telómeros (llama a telomeros_v2.R)
# ------------------------------------------------------------
# REQUIERE   : 01 (unidad_telo, telo_*), 02 (fa, sl, armonizar, rep_red, restar),
#              04 (bins, tabla, cob_pct), main.R (dir_scripts)
#              + el archivo telomeros_v2.R en la misma carpeta de scripts
# CREA       : rotaciones(), rc()          utilidades de secuencia
#              rot, motivos                 todas las rotaciones + complemento inverso de unidad_telo
#              tel                          salida cruda de detectar_telomeros()
#              tel_df                       salida de a_tel_df(): seqnames, start, end, longitud, copias,
#                                           extremo (inicio/final/interno), dist_borde, es_telomero, hebra,
#                                           secuencia, n_motivos, orientacion_ok, calidad
#              tel_gr                       GRanges de tel_df
#              res_telo(acc, ext)           resumen del telómero de un extremo (1 fila)
#              tabla$telomere, tabla$telo_interno, tabla$repeat_sinTelo_pct
# FUNCIONES DE telomeros_v2.R (quedan disponibles): conteo_manual(), detectar_telomeros(),
#              a_tel_df(), resumen_extremos_telo()
# USADO POR  : 06_subtelomeros (res_telo, tel_df), 07_satelites (motivos, tel_gr),
#              10/11/12 (repeat_sinTelo_pct, tel_df), 14_excel (tel_df), 15_figuras (tel_df)
# ============================================================

# ============================================================
# 6. TELÓMEROS
# ============================================================

rotaciones <- function(s) {
  n <- nchar(s)
  vapply(0:(n - 1), function(i)
    paste0(substring(s, i + 1, n), substring(s, 1, i)),
    character(1))
}
rc <- function(s) as.character(reverseComplement(DNAString(s)))
rot <- unique(unlist(lapply(unidad_telo, rotaciones), use.names = FALSE))
motivos <- unique(c(rot, vapply(rot, rc, character(1))))

source(file.path(dir_scripts, "telomeros_v2.R"))

tel <- detectar_telomeros(fa, sl, unidad = unidad_telo[1], gap_max = telo_gap_max,
                          min_pureza = telo_min_pureza, niveles = telo_niveles,
                          min_interno = telo_min_interno, exigir_orientacion = telo_exigir_orientacion)

if (nrow(tel) == 0) {
  tel_df <- data.frame(seqnames = character(), start = integer(), end = integer(),
                       longitud = numeric(), copias = numeric(), extremo = character(),
                       dist_borde = numeric(), es_telomero = logical(),
                       hebra = character(), secuencia = character(), n_motivos = integer(),
                       orientacion_ok = logical(), calidad = character())
} else {
  tel_df <- a_tel_df(tel, unidad_telo[1])
}

tel_gr <- armonizar(GRanges(tel_df$seqnames, IRanges(tel_df$start, tel_df$end)))

tryCatch(print(resumen_extremos_telo(tel, fa, unidad = unidad_telo[1])),
         error = function(e) message("resumen_extremos_telo omitido: ", conditionMessage(e)))

tabla$telomere <- "No"
tabla$telo_interno <- "No"
if (nrow(tel_df) > 0) {
  t_term <- tel_gr[tel_df$es_telomero & tel_df$extremo != "interno"]
  t_int <- tel_gr[tel_df$es_telomero & tel_df$extremo == "interno"]
  tabla$telomere[unique(queryHits(findOverlaps(bins, t_term)))] <- "Si"
  tabla$telo_interno[unique(queryHits(findOverlaps(bins, t_int)))] <- "Si"
}
tabla$repeat_sinTelo_pct <- cob_pct(restar(rep_red, tel_gr))

res_telo <- function(acc, ext) {
  d <- tel_df[tel_df$seqnames == acc & tel_df$extremo == ext & tel_df$es_telomero, ]
  if (nrow(d) == 0)
    return(data.frame(tiene_telomero = "No", telo_start = NA_integer_, telo_end = NA_integer_,
                      telo_long_pb = NA_integer_, telo_copias = NA_real_,
                      telo_dist_borde_pb = NA_real_, n_bloques_telo = 0L, telo_total_pb = 0L,
                      telo_secuencia = NA_character_, telo_hebra = NA_character_,
                      telo_secuencias_todas = NA_character_, telo_orientacion_ok = NA_character_,
                      telo_calidad = NA_character_, stringsAsFactors = FALSE))
  d1 <- d[which.max(d$n_motivos), ]   # bloque con más copias exactas
  data.frame(tiene_telomero = "Si", telo_start = d1$start, telo_end = d1$end,
             telo_long_pb = d1$longitud, telo_copias = d1$copias,
             telo_dist_borde_pb = d1$dist_borde, n_bloques_telo = nrow(d),
             telo_total_pb = sum(d$longitud),
             telo_secuencia = d1$secuencia, telo_hebra = d1$hebra,
             telo_secuencias_todas = paste(sort(unique(d$secuencia)), collapse = " + "),
             telo_orientacion_ok = ifelse(d1$orientacion_ok, "Si", "No"),
             telo_calidad = d1$calidad, stringsAsFactors = FALSE)
}
###AZAR

# Esto te dirá cuántas veces aparece en cada cromosoma
table(tel_df$seqnames[tel_df$es_telomero])


