# ============================================================
# 01_parametros.R  -  TODOS los parametros editables del analisis (UNO o VARIOS genomas)
# ------------------------------------------------------------
# REQUIERE   : nada (es el primer modulo; main.R ya cargo las librerias y definio dir_scripts).
# CREA       : bin, dir_base, dir_salida, n_hilos, rehacer_externos, genomas, genomas_run, genomas_comparar,
#              herramientas externas, seleccion de cromosomas, lectura del GFF, filtro de ventanas, telomeros,
#              subtelomeros, satelites, SAT350 (sat350_refs), microsatelites (primers_micro), repeticiones, extremos
# EFECTO     : setwd(dir_base/...). Los datos del genoma activo los aplica 01b_genoma_activo.R.
# GENOMAS    : cada entrada de `genomas` necesita fna; con gff el genoma corre en modo "completo"; sin gff en
#              modo "ligero" (cromosomas + SAT350 + microsatelites, suficiente para compararlo).
#              fai/rep/trf son opcionales (NULL = se calculan / se usa softmask / se ejecuta TRF si esta instalado).
# ============================================================

# ---- 0. Corrida ----
bin        <- 50000                                    # tamano de ventana (pb): 20000, 50000, 100000...
dir_base   <- "/home/luisafernanda/Documentos/Tarea1.1"
dir_salida <- paste0(dir_base, "/Analisis_cromosomas6/")   # carpeta base de TODAS las salidas
setwd(paste0(dir_base, "/Analisis_cromosomas6"))
n_hilos          <- 4        # hilos para BLAST y seqkit
rehacer_externos <- FALSE    # TRUE: rehace BLAST, TRF y seqkit aunque exista cache en externos_<genoma>/

# ---- 1. Genomas ----
genomas <- list(
  RH = list(
    fna = file.path(dir_base, "RH/GCA_033216535.1_TgRH_pasteur.chromosomes.fna"),
    gff = file.path(dir_base, "RH/GCA_033216535.1_TgRH_pasteur.chromosomes.gff3"),
    fai = file.path(dir_base, "RH/GCA_033216535.1_TgRH_pasteur.chromosomes.fna.fai"),
    rep = file.path(dir_base, "RH/GCA_033216535.1_TgRH_pasteur.chromosomes.fna.out"),                      # RepeatMasker
    trf = file.path(dir_base, "RH/GCA_033216535.1_TgRH_pasteur.chromosomes.fna.2.7.7.80.10.50.500.dat"),   # TRF
    eq_manual = data.frame(                                 # nombres romanos y orden propio
      acc    = c( "CM065576.1", "CM065577.1", "CM065578.1", "CM065579.1",
                  "CM065580.1", "CM065581.1", "CM065582.1", "CM065583.1",
                  "CM065584.1", "CM065585.1", "CM065586.1", "CM065587.1",
                  "CM065588.1"),
      nombre = c(  "Ia", "Ib", "II", "III", "IV", "V", "VI", "VIIa",
                   "VIIb/VIII", "IX", "X", "XI", "XII"),
      stringsAsFactors = FALSE)),
  # Sin gff = modo ligero. Agrega gff = "..." (y eq_manual si los nombres no salen bien) para el analisis completo.
  RH88 = list(
    fna = file.path(dir_base, "RH88/ToxoDB-71_TgondiiRH88_GenomeRename.fasta"),
    gff = file.path(dir_base, "RH88/ToxoDB-71_TgondiiRH88_GenomeRename.gff3"),
    eq_manual = data.frame(
      acc = c(
        "RH88_chrIa", "RH88_chrIb", "RH88_chrII", "RH88_chrIII",
        "RH88_chrIV", "RH88_chrV", "RH88_chrVI", "RH88_chrVIIa",
        "RH88_chrVIIb_VIII", "RH88_chrXII", "RH88_chrXI",
        "RH88_chrIX", "RH88_chrX"
      ),
      nombre = c(
        "Ia", "Ib", "II", "III", "IV", "V", "VI", "VIIa",
        "VIIb/VIII", "XI", "X", "IX", "XII"
      ),
      stringsAsFactors = FALSE
    )
  ),
              
  `ME49-Architecture` = list(fna = file.path(dir_base, "ME49-Architecture-22/ME49_cromosomas_nucleares.fasta")),
  
  `ME49-Intrastrain` = list(
    fna = file.path(
      dir_base,
      "ME49-Intra-strain-25/GCA_049996585.1_ASM4999658v1_genomic.fna"
    ),
    eq_manual = data.frame(
      acc = c(
        "CM113603.1", "CM113604.1",
        "JBMUND010000005.1", "JBMUND010000006.1",
        "CM113605.1", "CM113606.1",
        "JBMUND010000013.1", "JBMUND010000014.1",
        "CM113607.1",
        "JBMUND010000017.1", "JBMUND010000018.1",
        "CM113608.1",
        "JBMUND010000009.1", "JBMUND010000010.1",
        "CM113609.1", "CM113610.1"
      ),
      nombre = c(
        "Ia", "Ib",
        "II", "II.1",
        "III", "IV",
        "V", "V.1",
        "VI",
        "VIIa", "VIIa.1",
        "VIIb/VIII",
        "IX", "IX.1",
        "X", "XI"
      ),
      stringsAsFactors = FALSE
    )
  ),
  
  `ME49-V24` = list(fna = file.path(dir_base, "ME49-V24/ToxoDB-24_TgondiiME49_Genome.fasta")),
  `ME49-V71` = list(fna = file.path(dir_base, "ME49-V71/ToxoDB-71_TgondiiME49_Genome.fasta"))
)
genomas_run      <- c("RH")    # genomas a CORRER ahora (vector de nombres de `genomas`; NULL = todos). Con uno solo no se compara
genomas_comparar <- c("RH", "RH88", "ME49-Intrastrain")    # genomas a COMPARAR: NULL = los de genomas_run; o c("RH","ME49-V71") si ya tienen objetos/*.rds

# ---- 2. Herramientas externas (si faltan archivos .out/.dat) ----
carpeta_externos    <- NULL   # la define 01b por genoma (externos_<genoma>/)
ejecutar_trf          <- TRUE
ejecutar_repeatmasker <- FALSE      # lento; necesita rm_species (nombre de especie de Dfam)
rm_species          <- NULL       # p. ej. "toxoplasma gondii"
rm_threads          <- 4
usar_softmask       <- TRUE       # sin RepeatMasker: tomar como repeticion las minusculas del FASTA

# ---- 3. Que secuencias se analizan como "cromosomas" (genomas sin eq_manual) ----
min_len_cromosoma <- 1e6          # se ignoran secuencias mas cortas (scaffolds, organelos). Genomas pequenos: bajalo
patron_incluir    <- NULL         # regex opcional sobre el nombre de secuencia (p. ej. "^chr")
patron_excluir    <- NULL         # regex opcional (p. ej. "mito|apico|Un|random")
orden_cromosomas  <- "archivo"    # "archivo" (orden del FASTA/.fai) o "longitud"
eq_manual         <- NULL         # lo fija 01b desde genomas[[...]]$eq_manual

# ---- 4. Cómo se lee el GFF ----
tipos_gen <- NULL     # NULL = autodetectar los tipos de gen (gene, protein_coding_gene, ncRNA_gene, pseudogene...)
                      # o fijarlos: c("protein_coding_gene", "ncRNA_gene")
attr_id     <- c("gene_id", "ID", "locus_tag", "Name")     # atributos candidatos para el identificador (en orden)
attr_nombre <- c("Name", "gene", "gene_name")
attr_desc   <- c("description", "product", "Note")
mapa_gff_a_fasta <- NULL   # si los nombres del GFF no coinciden con el FASTA: c("1" = "CM000001.1", ...)

# ---- 5. Estilo y parámetros del análisis ----
font_family <- "sans"
color_GC    <- "#2166AC"
color_genes <- "#B2182B"

# ---- Genes (conjunto único) ----
incluir_pseudogenes <- TRUE      # TRUE: los pseudogenes cuentan como genes
prom_up             <- 1000
incluir_primer_gen <- TRUE 

# ---- Filtro de ventanas ----
max_N_pct        <- 10             # ventanas con más de 10% de N se descartan
min_ancho_frac   <- 0.5            # ventanas más cortas que 50% de `bin` se descartan
k_mad            <- 4              # outlier = |x - mediana| / MAD > k_mad (por cromosoma)
mad_vars         <- c("GC_percent","genes_per_Mb")# variables donde se buscan outliers (puedes agregar "genes_per_Mb")
excluir_outliers_MAD      <- TRUE  # TRUE: los outliers se QUITAN del análisis
excluir_extremos_analisis <- TRUE  # TRUE: se quita la ventana terminal técnica de cada extremo
w_excluir_ext    <- 50000          # ventanas que tocan los primeros/últimos w_excluir_ext pb
umbral_pct       <- 0.15           # percentil (por cromosoma) de TODAS las banderas bajo/alto
min_bp_gen_ventana <- 5000
min_bp_gc_comp   <- 200

# ---- Telómeros ----
w_telo       <- 20000
gap_fusion   <- 100
min_copias   <- 2
unidad_telo <- "TTTAGGG"      # motivo telomérico (hebra G-rich). Otros: "TTAGGG" (vertebrados, tripanosomátidos)
                              # Se busca también su complemento inverso. Solo motivos sin ambigüedades IUPAC.
escanear_secuencia_telo <- TRUE
telo_gap_max   <- 30          # separación máxima entre motivos del mismo bloque (pb)
telo_min_pureza <- 0.5        # fracción mínima del bloque cubierta por motivos
telo_niveles   <- data.frame(min_motivos = c(1, 1, 1), dist_max = c(10000, 1000, 6500))  # nº motivos / distancia al borde
telo_min_interno <- 5         # bloques internos solo si tienen >= N motivos
telo_exigir_orientacion <- FALSE

# ---- Subtelómeros ----
w_st_max         <- 300000
gap_st           <- 4500
k_gaps_st        <- 5
min_st           <- 8000
dens_st_max      <- 0.06           # NUEVO: densidad génica máxima del ST (genes/kb). ST válido si densidad < 0.04
dens_min_len     <- 20000   # NUEVO: la densidad solo se exige si el ST mide >= dens_min_len pb
# (en ST cortos, 1 gen basta para superar el umbral y la densidad es inestable)


# ---- Satélites ----
patron_sat       <- c(sat240 = "sat240", sat350 = "sat350")
gap_sat          <- 1000
trf_min_copias   <- 3
buscar_denovo    <- TRUE
periodos_sat     <- list(sat240 = 230:250, sat350 = 335:365)
#dn_blk           <- 500
#dn_min_ident     <- 0.50
#dn_min_len       <- 1000

dn_win           <- 60      # ventana de resolución de los límites (pb)
dn_min_ident     <- 0.60    # identidad mínima a distancia p (el azar da ~0.27, GC ~52 %)
dn_min_len       <- 500     # baja de 1000 para no perder bloques cortos como el sat240 de Ia
usar_rm_trf_sat  <- FALSE   # TRUE: suma además los hits de RepeatMasker y TRF
# dn_blk ya no se usa

# ---- Secuencias de referencia de GenBank para sat350 ----
# Secuencias de referencia SAT350 (GenBank). `incluir = FALSE`: se descargan solo si pasan a TRUE.
sat350_refs <- data.frame(
  nombre = c("BGTg7","TGR1A","TGR1E","TGR2","TGR4","TG3RH","TgSat350-1","TgSat350-2","TgSat350-3",
             "TgSat350-4","TgSat350-5","TgSat350-6","sin_identificar"),
  acceso = c("X74557","M57916","M57917","M57918","M57919","AF534393","AY330208","AY330209","AY330210",
             "AY333784","AY333785","AY333786","AF022237"),
  longitud_paper = c(721,352,353,674,1032,353,689,736,700,1045,1043,650,NA),
  incluir = c(rep(TRUE, 12), FALSE),                       # AF022237: confirma que es SAT350 antes de incluirla
  stringsAsFactors = FALSE)
gb_accesiones  <- sat350_refs$acceso[sat350_refs$incluir]
gb_satelite    <- "sat350"                 # a qué satélite se suman estos hits
gb_fasta_local <- "sat350_genbank.fasta"   # si no hay internet: guarda aquí las secuencias (FASTA) en <dir_salida>/externos_referencias/
usar_blast     <- TRUE                     # TRUE: SAT350 por BLASTN (necesita blastn/makeblastdb; sin ellos se omite el modulo 06b)
gb_min_ident   <- 70                       # BLASTN: identidad mínima (%) del hit
gb_min_len     <- 100                      # BLASTN: longitud mínima del hit (pb)
gb_evalue      <- 1e-10                     # BLASTN: e-value máximo

# ---- Repeticiones ----
top_n_rep        <- 25             # nº de repeticiones (por nombre) en el top con secuencia
max_seq_chars    <- 2000           # longitud máxima de la secuencia mostrada en Excel (límite de celda: 32767)
exportar_hits_rep <- TRUE          # TRUE: hoja con cada hit de RepeatMasker clasificado (exón/intrón/intergénico)

# ---- Extremos ----
ext_umbrales_kb  <- c(100, 200, 300)  # definiciones de "extremo": distancia al borde <= X kb
ext_paso_kb      <- 100               # ancho de clase del perfil por distancia
ext_max_kb       <- 1000              # perfil hasta X kb (luego una clase "≥ X kb")


# ---- SAT350: regiones y TRF (BLASTN con gb_evalue / gb_min_ident / gb_min_len, arriba) ----
sat350_flanco  <- 0                           # pb extra a cada lado de la region antes de TRF
sat350_periodo <- c(300, 400)                 # periodo que se considera monomero SAT350
trf_params_sat <- c(2, 7, 7, 80, 10, 50, 500) # TRF: match mismatch indel PM PI minscore maxperiod

# ---- Microsatelites (primers -> seqkit amplicon -> TRF). primers_micro = NULL omite el modulo ----
micro_mismatches   <- 2                       # mismatches tolerados en primers
micro_max_amplicon <- 2000                    # descarta amplicones mas largos (pb)
micro_trf_params   <- c(2, 7, 7, 80, 10, 20, 6)
# motivo: esperados separados por ";" (las dos cadenas con "/"). cromosoma_paper y posiciones: ensamblaje ME49 del articulo
primers_micro <- data.frame(
  marcador = c("TUB2_multiplex","TUB2_intron","W35","TgM_A","B18","B17","M33","IV.1","XI.1",
               "M48","M102","N60","N82","AA","N61","N83"),
  cromosoma_paper = c("TGME49_chrIX","TGME49_chrIX","TGME49_chrII","TGME49_chrX","TGME49_chrVIIa","TGME49_chrXII",
                      "TGME49_chrIV","TGME49_chrIV","TGME49_chrXI","TGME49_chrIa","TGME49_chrVIIa","TGME49_chrIb",
                      "TGME49_chrXII","TGME49_chrVIII","TGME49_chrVIIb","TGME49_chrX"),
  fwd = c("GTCCGGGTGTTCCTACAAAA","CCAAGTTCTTCCGTCATTTC","GGTTCACTGGATCTTCTCCAA","GGCGTCGACATGAGTTTCTC",
          "TGGTCTTCACCCTTTCATCC","AACAGACACCCGATGCCTAC","TACGCTTCGCATTGTACCAG","GAAGTTCGGCCTGTTCCTC",
          "GCGTGTGACGAGTTCTGAAA","AACATGTCGCGTAAGATTCG","CAGTCCAGGCATACCTCACC","GAATCGTCGAGGTGCTATCC",
          "TGCGTGCTTGTCAGAGTTC","GATGTCCGGTCAATTTTGCT","ATCGGCGGTGGTTGTAGAT","ATGGGTGAACAGCGTAGACA"),
  rev = c("TTGGCCAAAGACGAAGTTGT","CCTCATTGTAGAACACATTGAT","AATGAACGTCGCTTGTTTCC","TGGGCATGTAAATGTAGAGATG",
          "AGGGATAAGTTTCTTCACAACGA","GGCAACAGGAGGTAGAGGAG","TCTTTTCTCCCCTTCGCTCT","TCTGCCTGGAAAAGGAAAGA",
          "AAGTCCCCTGAAAAGCCAAT","CTCTTCACTGAGCGCCTTTC","CAATCCCAAAATCCCAAACC","AACGGTTGACCTGTGGCGAGT",
          "GCGTCCTTGACATGCACAT","GACGGGAAGGACAGAAACAC","CCTGATGTTGATGTAAGGATGC","GCAGGACGAAGAGGATGAGA"),
  tipo = c(rep("tipificacion", 9), rep("fingerprinting", 7)),
  motivo = c("TG/AC","TG/AC","TC/AG;TG/AC","TG/AC","TG/AC","TC/AG","TC/AG","TG/AC","TG/AC",
             "TA/AT",NA,NA,"TA/AT","TA/AT","TA/AT","TA/AT"),
  pos_ini_paper = c(974608,NA,633241,4824879,2921536,6474746,672591,742419,189702,332951,3093491,1766079,1621472,5836880,4217145,1772898),
  pos_fin_paper = c(974896,NA,633482,4825083,2921693,6475079,672760,742693,190058,333166,3093664,1766221,1621585,5837144,4217238,1773209),
  stringsAsFactors = FALSE)

# ---- Validacion de la longitud del ST contra el articulo (solo T. gondii; ver 06_subtelomeros.R) ----
validar_st_articulo <- TRUE
