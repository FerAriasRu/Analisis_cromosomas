# ============================================================
# 03_intergenico.R  -  Definición única de regiones intergénicas
# ------------------------------------------------------------
# REQUIERE   : 02_datos_base.R (sl, genes_gr, nombre_chr, armonizar)
# CREA       : construir_intergenico(g_gr)  función: tramos terminales y entre genes, con orientación
#              inter_all   data.frame: accession, tipo, start, end, dist_pb, orientacion, gen_izq, gen_der, chrom
#              ok_i        lógico: dist_pb > 0 (hay región intergénica real)
#              inter_gr    GRanges de las regiones intergénicas con dist_pb > 0
# MODIFICA   : inter_all se reordena/amplía en 08_intergenico_metricas.R
# USADO POR  : 04_ventanas, 06_subtelomeros, 08_intergenico_metricas, 10_repeticiones, 14_excel
# ============================================================

# ============================================================
# 4. INTERGÉNICO: DEFINICIÓN ÚNICA
#    dist_pb <= 0 -> genes solapados/contiguos: no hay región intergénica.
# ============================================================

construir_intergenico <- function(g_gr) {
  filas <- lapply(names(sl), function(acc) {
    g <- g_gr[as.character(seqnames(g_gr)) == acc]
    g <- g[order(start(g), end(g))]
    L <- as.integer(sl[acc]); n <- length(g)
    if (n == 0)
      return(data.frame(accession = acc, tipo = "sin_genes", start = 1L, end = L, dist_pb = L,
                        orientacion = "Terminal", gen_izq = NA_character_, gen_der = NA_character_,
                        stringsAsFactors = FALSE))
    s <- start(g); e <- end(g); st <- as.character(strand(g)); id <- g$gene_id
    emax <- cummax(e)
    imax <- integer(n); cur <- 1L
    for (i in seq_len(n)) { if (e[i] >= e[cur]) cur <- i; imax[i] <- cur }
    
    ini <- data.frame(accession = acc, tipo = "terminal_inicio", start = 1L, end = s[1] - 1L,
                      dist_pb = s[1] - 1L, orientacion = "Terminal",
                      gen_izq = NA_character_, gen_der = id[1], stringsAsFactors = FALSE)
    fin <- data.frame(accession = acc, tipo = "terminal_final", start = emax[n] + 1L, end = L,
                      dist_pb = L - emax[n], orientacion = "Terminal",
                      gen_izq = id[imax[n]], gen_der = NA_character_, stringsAsFactors = FALSE)
    if (n == 1) return(rbind(ini, fin))
    
    izq  <- imax[-n]; der <- 2:n
    dist <- s[-1] - emax[-n] - 1L
    sa <- st[izq]; sb <- st[der]
    orient <- ifelse(dist <= 0, "Solapados/contiguos",
                     ifelse(sa == "+" & sb == "-", "Convergentes (→←)",
                            ifelse(sa == "-" & sb == "+", "Divergentes (←→)",
                                   ifelse(sa == sb & sa %in% c("+", "-"), "Tándem (→→ / ←←)", "Indeterminada"))))
    pares <- data.frame(accession = acc, tipo = "entre_genes", start = emax[-n] + 1L, end = s[-1] - 1L,
                        dist_pb = dist, orientacion = orient,
                        gen_izq = id[izq], gen_der = id[der], stringsAsFactors = FALSE)
    rbind(ini, pares, fin)
  })
  do.call(rbind, filas)
}

inter_all <- construir_intergenico(genes_gr)
inter_all$chrom <- nombre_chr(inter_all$accession)
ok_i <- inter_all$dist_pb > 0
inter_gr <- armonizar(GRanges(inter_all$accession[ok_i], IRanges(inter_all$start[ok_i], inter_all$end[ok_i])))
stopifnot(length(inter_gr) == sum(ok_i))
cat("Regiones intergénicas:", length(inter_gr), "| pares de genes solapados/contiguos:",
    sum(inter_all$orientacion == "Solapados/contiguos"), "\n")

###############Informacion intergenico 

# ============================================================
# DISTANCIA DEL PRIMER GEN AL INICIO DE CADA CROMOSOMA
# ============================================================

primer_gen_chr <- do.call(rbind, lapply(names(sl), function(acc) {
  
  g <- genes_gr[as.character(seqnames(genes_gr)) == acc]
  
  if (length(g) == 0) {
    return(data.frame(
      accession = acc,
      primer_gen = NA_character_,
      inicio_gen = NA_integer_,
      distancia_inicio_pb = NA_integer_
    ))
  }
  
  # Primer gen según su posición genómica
  g <- g[order(start(g), end(g))]
  primero <- g[1]
  
  data.frame(
    accession = acc,
    primer_gen = as.character(primero$gene_id),
    inicio_gen = as.integer(start(primero)),
    distancia_inicio_pb = as.integer(start(primero) - 1L)
  )
}))

# Añadir nombres de cromosomas
primer_gen_chr$chrom <- nombre_chr(primer_gen_chr$accession)

# Mostrar resultados
primer_gen_chr <- primer_gen_chr[
  , c("accession", "chrom", "primer_gen",
      "inicio_gen", "distancia_inicio_pb")
]

print(primer_gen_chr)

