# ============================================================
# 17_preparar_proteinas.R
# Prepara las entradas de proteínas, calcula propiedades y PCA,
# y guarda datos_pca_proteinas.rds para el script de gráficas.
#
# Lee las rutas desde:
#   genomas[[genoma_activo]] en 01_parametros.R
#
# Campos por genoma:
#   fna              FASTA genómico
#   gff o gtf        anotación
#   protein_fasta    FASTA de proteínas existente o NULL
#   protein_products tabla de productos existente o NULL
# ============================================================

if (!exists("genomas", inherits = TRUE)) {
  source("01_parametros.R")
}

if (!exists("genoma_activo", inherits = TRUE)) {
  if (exists("genomas_run", inherits = TRUE) &&
      length(genomas_run) > 0 &&
      genomas_run[[1]] %in% names(genomas)) {
    genoma_activo <- genomas_run[[1]]
  } else {
    genoma_activo <- names(genomas)[[1]]
  }
}

if (!genoma_activo %in% names(genomas)) {
  stop("El genoma_activo '", genoma_activo,
       "' no está definido en genomas, en 01_parametros.R.")
}

cfg_genoma <- genomas[[genoma_activo]]

genome_file <- cfg_genoma[["fna"]]
annotation_file <- cfg_genoma[["gtf"]]

if (is.null(annotation_file) || !nzchar(annotation_file)) {
  annotation_file <- cfg_genoma[["gff"]]
}

protein_file <- cfg_genoma[["protein_fasta"]]
products_file <- cfg_genoma[["protein_products"]]

# En una corrida del main, carpeta_salida ya la define 01b_genoma_activo.R.
# Si se corre este script solo, crea una carpeta dentro de dir_salida.
if (exists("carpeta_salida", inherits = TRUE) &&
    length(carpeta_salida) == 1 &&
    nzchar(carpeta_salida)) {
  protein_outdir <- file.path(
    carpeta_salida,
    paste0("proteinas_", genoma_activo)
  )
} else {
  if (!exists("dir_salida", inherits = TRUE)) {
    stop("No existe dir_salida. Carga 01_parametros.R antes de este script.")
  }
  
  if (exists("bin", inherits = TRUE)) {
    ventana <- as.character(bin)
  } else {
    ventana <- "sin_ventana"
  }
  
  carpeta_genoma <- file.path(
    dir_salida,
    paste0("resultados_", genoma_activo, "_", ventana)
  )
  
  protein_outdir <- file.path(
    carpeta_genoma,
    paste0("proteinas_", genoma_activo)
  )
}

dir.create(protein_outdir, recursive = TRUE, showWarnings = FALSE)

if (is.null(protein_file) || !nzchar(protein_file)) {
  protein_file <- file.path(protein_outdir, "proteins.fasta")
}

if (is.null(products_file) || !nzchar(products_file)) {
  products_file <- file.path(protein_outdir, "products.tsv")
}

if (is.null(genome_file) || !file.exists(genome_file)) {
  stop("No existe el FASTA genómico configurado para '",
       genoma_activo, "': ", genome_file)
}

if (is.null(annotation_file) || !file.exists(annotation_file)) {
  stop("No existe una anotación GFF/GTF para '", genoma_activo,
       "'. Configura gff o gtf en 01_parametros.R.")
}

# ---- Dependencias ----
paquetes <- c(
  "Biostrings", "dplyr", "stringr", "ggplot2",
  "Peptides", "RColorBrewer"
)

faltan <- paquetes[
  !vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)
]

if (length(faltan)) {
  stop(
    "Faltan paquetes R requeridos: ",
    paste(faltan, collapse = ", "),
    "\nInstálalos antes de ejecutar este módulo."
  )
}

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(Peptides)
  library(RColorBrewer)
})

normalizar_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("\\s.*$", "", x)
  x <- sub(":pep$", "", x)
  x
}

# ---- Generar FASTA de proteínas si no existe ----
if (!file.exists(protein_file) || file.info(protein_file)$size == 0) {
  gffread_bin <- Sys.which("gffread")
  
  if (!nzchar(gffread_bin)) {
    stop(
      "No se encontró gffread en el PATH. Instálalo o configura ",
      "protein_fasta en 01_parametros.R con un FASTA existente."
    )
  }
  
  cat("Generando FASTA de proteínas para ", genoma_activo, "...\n", sep = "")
  
  mensajes <- system2(
    gffread_bin,
    args = c(
      shQuote(annotation_file),
      "-g", shQuote(genome_file),
      "-y", shQuote(protein_file)
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  
  estado <- attr(mensajes, "status")
  if (is.null(estado)) estado <- 0L
  
  if (estado != 0L ||
      !file.exists(protein_file) ||
      file.info(protein_file)$size == 0) {
    stop(
      "gffread no pudo generar el FASTA de proteínas.\n",
      paste(mensajes, collapse = "\n")
    )
  }
}

cat("Leyendo proteínas: ", protein_file, "\n", sep = "")
prot <- Biostrings::readAAStringSet(protein_file)

if (length(prot) == 0) {
  stop("El FASTA de proteínas está vacío: ", protein_file)
}

ids_prot <- names(prot)
ids_prot_norm <- normalizar_id(ids_prot)

# ---- Funciones para leer atributos GFF3 o GTF ----
extraer_atributo <- function(atributos, nombre) {
  # GFF3: nombre=valor
  patron_gff <- paste0("(?:^|;)", nombre, "=([^;]*)")
  hallazgo <- regexec(patron_gff, atributos, perl = TRUE)
  valor <- regmatches(atributos, hallazgo)[[1]]
  
  if (length(valor) >= 2) {
    return(utils::URLdecode(valor[[2]]))
  }
  
  # GTF: nombre "valor"
  patron_gtf <- paste0("(?:^|;\\s*)", nombre, "\\s+\"([^\"]*)\"")
  hallazgo <- regexec(patron_gtf, atributos, perl = TRUE)
  valor <- regmatches(atributos, hallazgo)[[1]]
  
  if (length(valor) >= 2) {
    return(valor[[2]])
  }
  
  NA_character_
}

leer_productos_existentes <- function(archivo) {
  x <- tryCatch(
    read.delim(
      archivo,
      header = FALSE,
      sep = "\t",
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      fill = TRUE
    ),
    error = function(e) NULL
  )
  
  if (is.null(x) || ncol(x) < 2) {
    return(NULL)
  }
  
  x <- x[, 1:2, drop = FALSE]
  names(x) <- c("protein_id", "product")
  x$protein_id <- normalizar_id(x$protein_id)
  x$product <- as.character(x$product)
  x
}

extraer_productos_anotacion <- function(archivo_anotacion) {
  anot <- tryCatch(
    read.delim(
      archivo_anotacion,
      header = FALSE,
      sep = "\t",
      comment.char = "#",
      quote = "",
      stringsAsFactors = FALSE,
      fill = TRUE
    ),
    error = function(e) NULL
  )
  
  if (is.null(anot) || ncol(anot) < 9) {
    return(data.frame(
      protein_id = character(),
      product = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  anot <- anot[, 1:9, drop = FALSE]
  names(anot) <- c(
    "seqid", "source", "type", "start", "end",
    "score", "strand", "phase", "attributes"
  )
  
  tipos_prioridad <- c(
    "polypeptide", "mRNA", "transcript",
    "CDS", "protein"
  )
  
  resultado <- list()
  
  for (tipo_actual in tipos_prioridad) {
    filas <- anot[anot$type == tipo_actual, , drop = FALSE]
    if (!nrow(filas)) next
    
    for (i in seq_len(nrow(filas))) {
      attrs <- filas$attributes[[i]]
      
      id <- extraer_atributo(attrs, "ID")
      
      if (is.na(id) || !nzchar(id)) {
        id <- extraer_atributo(attrs, "protein_id")
      }
      
      if (is.na(id) || !nzchar(id)) {
        id <- extraer_atributo(attrs, "transcript_id")
      }
      
      if (is.na(id) || !nzchar(id)) {
        id <- extraer_atributo(attrs, "Parent")
      }
      
      if (is.na(id) || !nzchar(id)) next
      
      producto <- extraer_atributo(attrs, "product")
      
      if (is.na(producto) || !nzchar(producto)) {
        producto <- extraer_atributo(attrs, "Name")
      }
      
      if (is.na(producto) || !nzchar(producto)) {
        producto <- extraer_atributo(attrs, "gene_name")
      }
      
      if (is.na(producto) || !nzchar(producto)) {
        producto <- extraer_atributo(attrs, "description")
      }
      
      if (is.na(producto) || !nzchar(producto)) next
      
      producto <- sub("^term=", "", producto)
      producto <- sub(";.*$", "", producto)
      id <- normalizar_id(id)
      
      resultado[[length(resultado) + 1L]] <- data.frame(
        protein_id = id,
        product = producto,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (!length(resultado)) {
    return(data.frame(
      protein_id = character(),
      product = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  bind_rows(resultado) %>%
    filter(
      !is.na(protein_id), protein_id != "",
      !is.na(product), product != ""
    ) %>%
    distinct(protein_id, .keep_all = TRUE)
}

# Usa products.tsv existente o lo extrae de la anotación.
if (file.exists(products_file) && file.info(products_file)$size > 0) {
  cat("Leyendo tabla de productos: ", products_file, "\n", sep = "")
  products <- leer_productos_existentes(products_file)
} else {
  cat("Extrayendo productos de la anotación...\n")
  products <- extraer_productos_anotacion(annotation_file)
  
  if (nrow(products)) {
    write.table(
      products,
      products_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    cat("Tabla de productos guardada en: ", products_file, "\n", sep = "")
  }
}

if (is.null(products)) {
  products <- data.frame(
    protein_id = character(),
    product = character(),
    stringsAsFactors = FALSE
  )
}

if (nrow(products)) {
  products$protein_id <- normalizar_id(products$protein_id)
  products$product <- as.character(products$product)
  
  products <- products %>%
    filter(!is.na(protein_id), protein_id != "") %>%
    distinct(protein_id, .keep_all = TRUE)
}

# Conserva proteínas sin descripción y las clasifica luego como "Otras".
ids_sin_producto <- setdiff(ids_prot_norm, products$protein_id)

if (length(ids_sin_producto)) {
  products <- bind_rows(
    products,
    data.frame(
      protein_id = ids_sin_producto,
      product = "Sin descripción funcional",
      stringsAsFactors = FALSE
    )
  )
}

# ---- Composición de aminoácidos ----
aa_freq <- Biostrings::letterFrequency(
  prot,
  letters = Biostrings::AA_STANDARD,
  as.prob = TRUE
)

df_aa <- as.data.frame(aa_freq)
df_aa$protein_id <- ids_prot_norm
df_aa <- df_aa %>%
  distinct(protein_id, .keep_all = TRUE)

# ---- Asignar familias funcionales ----
products <- products %>%
  mutate(
    product = utils::URLdecode(product),
    product = str_remove(product, "^term="),
    product = str_remove(product, ";.*$"),
    product = str_remove(product, ", putative$"),
    product = str_remove(product, " putative$"),
    product = str_trim(product),
    familia = case_when(
      str_detect(product, regex("SAG-related|\\bSRS[0-9]*\\b|surface antigen", ignore_case = TRUE)) ~ "SAG/SRS (superficie)",
      str_detect(product, regex("Toxoplasma gondii family", ignore_case = TRUE)) ~ "Familias Tg (A-E)",
      str_detect(product, regex("dense granule|\\bGRA[0-9]*\\b", ignore_case = TRUE)) ~ "GRA (gránulos densos)",
      str_detect(product, regex("rhoptry neck|\\bRON[0-9]*\\b", ignore_case = TRUE)) ~ "RON (cuello de roptria)",
      str_detect(product, regex("rhoptry|\\bROP[0-9]*\\b", ignore_case = TRUE)) ~ "ROP (roptrias)",
      str_detect(product, regex("microneme|\\bMIC[0-9]*\\b", ignore_case = TRUE)) ~ "MIC (micronemas)",
      str_detect(product, regex("inner membrane complex|\\bIMC[0-9]*\\b", ignore_case = TRUE)) ~ "IMC (complejo membrana)",
      str_detect(product, regex("apical cap|apical cone|apical neck", ignore_case = TRUE)) ~ "Apical (cap/cone)",
      str_detect(product, regex("alveolin", ignore_case = TRUE)) ~ "Alveolinas",
      str_detect(product, regex("\\bKRUF[0-9]*\\b", ignore_case = TRUE)) ~ "KRUF",
      str_detect(product, regex("SNARE|vesicle transport|clathrin|coatomer|dynamin|adaptin|\\bVPS[0-9]+\\b", ignore_case = TRUE)) ~ "Tráfico vesicular / SNARE",
      str_detect(product, regex("\\bRab\\b|Rho GTPase|Ras-related|ADP-ribosylation factor|\\bseptin", ignore_case = TRUE)) ~ "GTPasas pequeñas",
      str_detect(product, regex("apicoplast|\\bplastid\\b|ferredoxin|iron-sulfur|Fe-S cluster", ignore_case = TRUE)) ~ "Apicoplasto / Fe-S",
      str_detect(product, regex("AP2 domain transcription", ignore_case = TRUE)) ~ "AP2 factores transcripción",
      str_detect(product, regex("\\bcyclin\\b|cell cycle|\\bcentrin\\b|aurora kinase", ignore_case = TRUE)) ~ "Ciclo celular",
      str_detect(product, regex("\\bkinase\\b|kinase-like|\\bNEK[0-9]*\\b", ignore_case = TRUE)) ~ "Quinasas",
      str_detect(product, regex("phosphatase", ignore_case = TRUE)) ~ "Fosfatasas",
      str_detect(product, regex("subtilisin|protease|peptidase|proteasome", ignore_case = TRUE)) ~ "Proteasas",
      str_detect(product, regex("DNA polymerase|primase|\\bMCM[0-9]+\\b|mismatch repair|replication factor|replication protein", ignore_case = TRUE)) ~ "Replicación/reparación DNA",
      str_detect(product, regex("ribosomal|ribosome", ignore_case = TRUE)) ~ "Ribosomales",
      str_detect(product, regex("elongation factor|initiation factor|translation|tRNA", ignore_case = TRUE)) ~ "Traducción/tRNA",
      str_detect(product, regex("RNA polymerase|mediator complex|transcription", ignore_case = TRUE)) ~ "RNA polimerasa / transcripción",
      str_detect(product, regex("cytochrome|ATP synthase|oxidase|reductase|dehydrogenase", ignore_case = TRUE)) ~ "Metabolismo energético",
      str_detect(product, regex("transporter|permease|carrier|channel|ABC transporter|facilitator", ignore_case = TRUE)) ~ "Transportadores",
      str_detect(product, regex("glycosyltransferase|GPI anchor|glycosylphosphatidylinositol", ignore_case = TRUE)) ~ "Glicosilación / GPI",
      str_detect(product, regex("synthase|transferase|isomerase|lyase|ligase|mutase|hydrolase|phosphodiesterase|palmitoyltransferase", ignore_case = TRUE)) ~ "Metabolismo general",
      str_detect(product, regex("histone|methyltransferase.*histone|acetyltransferase|chromatin|SWI2/SNF2", ignore_case = TRUE)) ~ "Modificadores de cromatina",
      str_detect(product, regex("ubiquitin|HECT|OTU domain|carboxyl-terminal hydrolase", ignore_case = TRUE)) ~ "Ubiquitina / UPS",
      str_detect(product, regex("DnaJ|heat shock|Hsp|chaperone|thioredoxin", ignore_case = TRUE)) ~ "Chaperonas / Heat shock",
      str_detect(product, regex("tubulin|kinesin|dynein|myosin|actin|centrosomal|basal complex", ignore_case = TRUE)) ~ "Citoesqueleto / Motores",
      str_detect(product, regex("helicase|DEAD/DEAH", ignore_case = TRUE)) ~ "Helicasas",
      str_detect(product, regex("EF hand|calcium-dependent|calmodulin", ignore_case = TRUE)) ~ "Señalización calcio",
      str_detect(product, regex("MORN|kelch|RCC1|TBC domain|RAP domain|endonuclease|exonuclease", ignore_case = TRUE)) ~ "Otras enzimáticas",
      str_detect(product, regex("WD domain|RNA recognition motif|tetratricopeptide|leucine rich|HEAT repeat|ankyrin|zinc finger", ignore_case = TRUE)) ~ "Repeticiones (función desconocida)",
      str_detect(product, regex("hypothetical protein, conserved|conserved hypothetical", ignore_case = TRUE)) ~ "Hipotéticas conservadas",
      str_detect(product, regex("hypothetical", ignore_case = TRUE)) ~ "Hipotéticas",
      TRUE ~ "Otras"
    )
  )

# ---- Coordenadas cromosómicas de las proteínas ----
# Prioridad:
# 1. polypeptide: coordenada directa de la proteína en GFF3.
# 2. mRNA/transcript: intervalo del transcrito.
# 3. CDS: intervalo combinado de los CDS asociados al ID de proteína.

anot_loc <- read.delim(
  annotation_file,
  header = FALSE,
  sep = "\t",
  comment.char = "#",
  quote = "",
  stringsAsFactors = FALSE,
  fill = TRUE
)

if (ncol(anot_loc) < 9) {
  stop("La anotación no tiene las 9 columnas esperadas de GFF/GTF.")
}

anot_loc <- anot_loc[, 1:9, drop = FALSE]
names(anot_loc) <- c(
  "seqid", "source", "type", "start", "end",
  "score", "strand", "phase", "attributes"
)

coords_lista <- list()

for (i in seq_len(nrow(anot_loc))) {
  tipo <- anot_loc$type[[i]]
  attrs <- anot_loc$attributes[[i]]
  
  if (tipo == "polypeptide") {
    id <- extraer_atributo(attrs, "ID")
    prioridad <- 1L
    
  } else if (tipo %in% c("mRNA", "transcript")) {
    id <- extraer_atributo(attrs, "ID")
    if (is.na(id) || !nzchar(id)) {
      id <- extraer_atributo(attrs, "transcript_id")
    }
    prioridad <- 2L
    
  } else if (tipo == "CDS") {
    id <- extraer_atributo(attrs, "protein_id")
    if (is.na(id) || !nzchar(id)) {
      id <- extraer_atributo(attrs, "transcript_id")
    }
    if (is.na(id) || !nzchar(id)) {
      id <- extraer_atributo(attrs, "Parent")
    }
    prioridad <- 3L
    
  } else {
    next
  }
  
  if (is.na(id) || !nzchar(id)) next
  
  # Un CDS GFF3 puede tener más de un Parent.
  ids <- strsplit(id, ",", fixed = TRUE)[[1]]
  ids <- normalizar_id(ids)
  
  for (id_actual in ids) {
    coords_lista[[length(coords_lista) + 1L]] <- data.frame(
      protein_id = id_actual,
      cromosoma = anot_loc$seqid[[i]],
      start = as.numeric(anot_loc$start[[i]]),
      end = as.numeric(anot_loc$end[[i]]),
      strand = anot_loc$strand[[i]],
      prioridad = prioridad,
      stringsAsFactors = FALSE
    )
  }
}

if (length(coords_lista)) {
  coords_raw <- dplyr::bind_rows(coords_lista)
  
  coords_raw <- coords_raw %>%
    filter(
      !is.na(protein_id), protein_id != "",
      is.finite(start), is.finite(end)
    )
  
  # Por proteína conserva el tipo de coordenada de mayor prioridad.
  localizacion_proteinas <- coords_raw %>%
    group_by(protein_id) %>%
    filter(prioridad == min(prioridad)) %>%
    summarise(
      cromosoma = first(cromosoma),
      start = min(start),
      end = max(end),
      strand = first(strand),
      .groups = "drop"
    ) %>%
    inner_join(
      products %>% select(protein_id, familia),
      by = "protein_id"
    ) %>%
    distinct(protein_id, .keep_all = TRUE) %>%
    mutate(posicion = (start + end) / 2)
  
} else {
  localizacion_proteinas <- data.frame(
    protein_id = character(),
    cromosoma = character(),
    start = numeric(),
    end = numeric(),
    strand = character(),
    familia = character(),
    posicion = numeric(),
    stringsAsFactors = FALSE
  )
}

write.csv(
  localizacion_proteinas,
  file.path(protein_outdir, "localizacion_proteinas.csv"),
  row.names = FALSE
)

cat(
  "Proteínas con familia y coordenadas: ",
  nrow(localizacion_proteinas),
  "\n",
  sep = ""
)


# ---- Unir frecuencias, productos y secuencias ----
# Composición por proteína, expresada como porcentaje.
aa_cols <- setdiff(names(df_aa), "protein_id")

composicion_aa <- df_aa %>%
  inner_join(
    products %>% select(protein_id, familia),
    by = "protein_id"
  ) %>%
  mutate(
    across(all_of(aa_cols), ~ .x * 100)
  )

# ---- Calcular GRAVY por proteína ----
seqs_gravy <- as.character(prot)
names(seqs_gravy) <- ids_prot_norm
seqs_gravy <- seqs_gravy[!duplicated(names(seqs_gravy))]

composicion_aa$GRAVY <- vapply(
  composicion_aa$protein_id,
  function(id) {
    secuencia <- seqs_gravy[[id]]
    
    if (is.null(secuencia) || is.na(secuencia) || !nzchar(secuencia)) {
      return(NA_real_)
    }
    
    tryCatch(
      Peptides::hydrophobicity(
        secuencia,
        scale = "KyteDoolittle"
      ),
      error = function(e) NA_real_
    )
  },
  numeric(1)
)

df_full <- df_aa %>%
  inner_join(products, by = "protein_id")

if (nrow(df_full) == 0) {
  stop(
    "No coincidieron IDs entre proteínas y productos. ",
    "Revisa los encabezados del FASTA y los IDs de la anotación."
  )
}

cat("Proteínas leídas: ", length(prot), "\n", sep = "")
cat("IDs coincidentes: ", nrow(df_full), "\n", sep = "")

df_pca <- df_full %>%
  filter(!familia %in% c("Hipotéticas", "Hipotéticas conservadas"))

if (nrow(df_pca) < 3) {
  stop("Quedaron menos de tres proteínas para el PCA.")
}

# ---- Propiedades fisicoquímicas ----
seqs_chr <- as.character(prot)
names(seqs_chr) <- ids_prot_norm
seqs_chr <- seqs_chr[!duplicated(names(seqs_chr))]

seqs_pca <- seqs_chr[df_pca$protein_id]
seq_ok <- !is.na(seqs_pca) & nchar(seqs_pca) > 0

df_pca <- df_pca[seq_ok, , drop = FALSE]
seqs_pca <- seqs_pca[seq_ok]

if (length(seqs_pca) < 3) {
  stop("Quedaron menos de tres secuencias válidas para el PCA.")
}

physchem <- data.frame(
  protein_id = names(seqs_pca),
  mw = Peptides::mw(seqs_pca, monoisotopic = FALSE),
  pI = Peptides::pI(seqs_pca, pKscale = "EMBOSS"),
  charge_pH7 = Peptides::charge(seqs_pca, pH = 7, pKscale = "EMBOSS"),
  hydrophobicity = Peptides::hydrophobicity(
    seqs_pca, scale = "KyteDoolittle"
  ),
  instability = Peptides::instaIndex(seqs_pca),
  aliphatic_index = Peptides::aIndex(seqs_pca),
  boman_index = Peptides::boman(seqs_pca),
  stringsAsFactors = FALSE
)

df_arom <- df_aa %>%
  transmute(protein_id, aromaticity = F + Y + W)

physchem <- physchem %>%
  left_join(df_arom, by = "protein_id")

df_pca_full <- df_pca %>%
  select(protein_id, familia) %>%
  inner_join(physchem, by = "protein_id")

vars_pca <- c(
  "mw", "pI", "charge_pH7", "hydrophobicity",
  "instability", "aliphatic_index", "boman_index", "aromaticity"
)

mat_pca <- as.matrix(df_pca_full[, vars_pca, drop = FALSE])
filas_ok <- complete.cases(mat_pca)
mat_pca <- mat_pca[filas_ok, , drop = FALSE]
df_pca_full <- df_pca_full[filas_ok, , drop = FALSE]

sd_vars <- apply(mat_pca, 2, sd, na.rm = TRUE)
vars_validas <- names(sd_vars)[is.finite(sd_vars) & sd_vars > 0]

if (length(vars_validas) < 2) {
  stop("Hay menos de dos propiedades variables para calcular el PCA.")
}

mat_pca <- mat_pca[, vars_validas, drop = FALSE]

if (nrow(mat_pca) < 3) {
  stop("Hay menos de tres proteínas completas para el PCA.")
}

pca_res <- prcomp(mat_pca, center = TRUE, scale. = TRUE)
var_exp <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100

df_scores <- data.frame(
  protein_id = df_pca_full$protein_id,
  familia = df_pca_full$familia,
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  stringsAsFactors = FALSE
)

# Estadísticas descriptivas por familia
df_stats <- df_pca_full %>%
  inner_join(df_scores %>% select(protein_id, PC1, PC2), by = "protein_id")

tabla_stats <- df_stats %>%
  group_by(familia) %>%
  summarise(
    n = n(),
    across(
      c(PC1, PC2, all_of(vars_validas)),
      list(
        media = ~mean(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(n)) %>%
  mutate(across(-c(familia, n), ~round(.x, 2)))

write.csv(
  tabla_stats,
  file.path(protein_outdir, "tabla_estadisticas_por_familia.csv"),
  row.names = FALSE
)

loadings <- as.data.frame(pca_res$rotation[, 1:2, drop = FALSE])
loadings$variable <- rownames(loadings)
loadings <- loadings[
  order(-abs(loadings$PC1)),
  c("variable", "PC1", "PC2")
]

write.csv(
  loadings,
  file.path(protein_outdir, "PCA_loadings.csv"),
  row.names = FALSE
)

# El módulo 18 lee este archivo para dibujar las gráficas.
datos_pca_proteinas <- list(
  genoma = genoma_activo,
  scores = df_scores,
  pca = pca_res,
  var_exp = var_exp,
  vars_pca = vars_validas,
  localizacion = localizacion_proteinas,
  composicion_aa = composicion_aa
)

saveRDS(
  datos_pca_proteinas,
  file.path(protein_outdir, "datos_pca_proteinas.rds")
)

cat("\nPreparación de proteínas terminada para ", genoma_activo, ".\n", sep = "")
cat("Datos y tablas guardados en: ", protein_outdir, "\n", sep = "")