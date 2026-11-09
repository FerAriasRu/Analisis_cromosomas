# ============================================================
# telomeros_v2.R  -  FUNCIONES de detección de telómeros (no ejecuta nada al cargarse)
# ------------------------------------------------------------
# LLAMADO POR: 05_telomeros.R  (source)
# REQUIERE   : Biostrings, GenomicRanges; en uso: fa (DNAStringSet) y sl (longitudes por accession)
# FUNCIONES  :
#   conteo_manual(fa, w, unidad)          motivos sueltos en los primeros/últimos w pb
#   detectar_telomeros(fa, sl, unidad, gap_max, min_pureza, niveles, min_interno, exigir_orientacion)
#                                         agrupa motivos cercanos -> bloques con borde/orientación/calidad
#   a_tel_df(tel, unidad)                 formatea la salida para el resto del pipeline (tel_df);
#                                         `unidad` = motivo telomérico (el complemento inverso se deriva solo)
#   resumen_extremos_telo(tel, fa, w, unidad)   una fila por cromosoma y extremo (control de calidad)
# DIAGNÓSTICO: ver diagnostico_telomeros.R (se corre a mano, no forma parte del pipeline)
# ============================================================
# ============================================================
# telomeros_v2.R  -  Detección de telómeros (reemplaza la sección 6)
# Sin ejecutar todavía: probar y comparar contra la búsqueda manual.
#
# Ideas clave:
#  1. Un solo escaneo con vmatchPattern (C, vectorizado) en vez de 14 gregexpr.
#  2. Se agrupan motivos cercanos (gap_max) -> tolera repeticiones degeneradas
#     (1 mismatch entre motivos no rompe el bloque, a diferencia de la regex exacta).
#  3. Se exige: borde cercano + orientación biológica correcta + nº de motivos.
#     Orientación esperada (hebra + del FASTA):
#        inicio del cromosoma -> CCCTAAA (hebra "-" del motivo TTTAGGG)
#        final  del cromosoma -> TTTAGGG (hebra "+")
#  4. Se reporta cada candidato con su "calidad" para que veas los casos dudosos.
# ============================================================
suppressPackageStartupMessages({ library(Biostrings); library(GenomicRanges) })

# ---- 1. Replica programática de tu búsqueda manual (conteo de motivos sueltos) ----
conteo_manual <- function(fa, w = 8000, unidad = "TTTAGGG") {
  pf <- DNAString(unidad); pr <- reverseComplement(pf)
  L  <- width(fa)
  ini <- subseq(fa, 1, pmin(w, L))
  fin <- subseq(fa, pmax(1, L - w + 1), L)
  data.frame(
    accession = names(fa),
    TTTAGGG_inicio = vcountPattern(pf, ini), TTTAGGG_final = vcountPattern(pf, fin),
    CCCTAAA_inicio = vcountPattern(pr, ini), CCCTAAA_final = vcountPattern(pr, fin),
    row.names = NULL, stringsAsFactors = FALSE)
}

# ---- 2. Detección por agrupamiento de motivos ----

detectar_telomeros <- function(fa, sl,
                               unidad       = "TTTAGGG",
                               gap_max      = 30,
                               min_pureza   = 0.5,
                               niveles = data.frame(min_motivos = c(1, 1, 1),
                                                    dist_max    = c(10000, 1000, 6500)),
                               min_interno  = 5,      # bloques internos solo si tienen >= N motivos
                               exigir_orientacion = FALSE) {
  pf <- DNAString(unidad); pr <- reverseComplement(pf); k <- nchar(unidad)
  
  hits_de <- function(p, st) {
    mi <- vmatchPattern(p, fa)
    n  <- elementNROWS(mi)
    GRanges(rep(names(fa), n), unlist(mi, use.names = FALSE), strand = rep(st, sum(n)))
  }
  h <- c(hits_de(pf, "+"), hits_de(pr, "-"))
  seqlevels(h) <- names(sl); seqlengths(h) <- sl
  
  bl <- IRanges::reduce(h, min.gapwidth = gap_max + 1)
  if (length(bl) == 0) return(data.frame())
  n_mot  <- countOverlaps(bl, h)
  pureza <- n_mot * k / width(bl)
  
  acc   <- as.character(seqnames(bl)); L <- as.numeric(sl[acc])
  d_ini <- start(bl) - 1; d_fin <- L - end(bl)
  lado  <- ifelse(d_ini <= d_fin, "inicio", "final")
  dist  <- pmin(d_ini, d_fin)
  hebra <- as.character(strand(bl))
  orient_ok <- (lado == "inicio" & hebra == "-") | (lado == "final" & hebra == "+")
  
  # criterio escalonado: menos copias exigidas cuanto más cerca del borde
  pasa_term <- Reduce(`|`, lapply(seq_len(nrow(niveles)), function(i)
    n_mot >= niveles$min_motivos[i] & dist <= niveles$dist_max[i]))
  pasa_term <- pasa_term & pureza >= min_pureza
  
  keep <- pasa_term | (n_mot >= min_interno & pureza >= min_pureza)
  if (!any(keep)) return(data.frame())
  
  terminal <- pasa_term
  es_tel   <- terminal & (orient_ok | !exigir_orientacion)
  
  data.frame(
    seqnames = acc, start = start(bl), end = end(bl), longitud = width(bl),
    hebra = hebra, n_motivos = n_mot, pureza = round(pureza, 2),
    lado = lado, dist_borde = dist, terminal = terminal,
    orientacion_ok = orient_ok, es_telomero = es_tel,
    calidad = ifelse(terminal & orient_ok & n_mot >= 5, "Alta",
                     ifelse(terminal & orient_ok, "Media (pocas copias, cerca del borde)",
                            ifelse(terminal, "Dudosa (orientación invertida)", "Interna"))),
    stringsAsFactors = FALSE)[keep, ]
}

a_tel_df <- function(tel, unidad = "TTTAGGG") {
  k <- nchar(unidad)
  data.frame(seqnames = tel$seqnames, start = tel$start, end = tel$end,
             longitud = tel$longitud, copias = round(tel$longitud / k),
             extremo  = ifelse(tel$terminal, tel$lado, "interno"),
             dist_borde  = ifelse(tel$terminal, tel$dist_borde, NA_real_),
             es_telomero = TRUE,
             hebra       = tel$hebra,
             secuencia   = ifelse(tel$hebra == "+", unidad, as.character(reverseComplement(DNAString(unidad)))),
             n_motivos   = tel$n_motivos,
             orientacion_ok = tel$orientacion_ok,
             calidad     = tel$calidad,
             stringsAsFactors = FALSE)
}

# ---- 4. Resumen por extremo (una fila por cromosoma y extremo) ----
resumen_extremos_telo <- function(tel, fa, w = 8000, unidad = "TTTAGGG") {
  cm <- conteo_manual(fa, w, unidad)
  fila <- function(acc, lado) {
    t <- tel[tel$seqnames == acc & tel$lado == lado & tel$terminal, ]
    sueltos <- if (lado == "inicio") cm[cm$accession == acc, c("CCCTAAA_inicio", "TTTAGGG_inicio")]
               else                  cm[cm$accession == acc, c("TTTAGGG_final",  "CCCTAAA_final")]
    data.frame(accession = acc, extremo = lado,
               motivos_sueltos_esperada  = sueltos[[1]],   # orientación biológica correcta
               motivos_sueltos_invertida = sueltos[[2]],
               n_bloques = nrow(t),
               mejor_calidad = if (nrow(t)) t$calidad[which.max(t$n_motivos)] else "Sin bloque",
               n_motivos_mejor = if (nrow(t)) max(t$n_motivos) else 0L,
               stringsAsFactors = FALSE)
  }
  do.call(rbind, lapply(names(fa), function(a) rbind(fila(a, "inicio"), fila(a, "final"))))
}

# ---- Uso (en lugar de la sección 6 actual) ----
# tel      <- detectar_telomeros(fa, sl)
# tel_df   <- a_tel_df(tel)
# tel_gr   <- armonizar(GRanges(tel_df$seqnames, IRanges(tel_df$start, tel_df$end)))
# resumen  <- resumen_extremos_telo(tel, fa)
# print(resumen)
