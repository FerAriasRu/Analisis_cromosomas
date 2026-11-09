# ============================================================
# 17_pca_familias.R - PCA de propiedades fisicoquimicas por familia
#
# Ejecucion (desde la carpeta del proyecto):
# Rscript 17_pca_familias.R --fasta /ruta/proteins_CM.fasta \
#   --products /ruta/products.tsv [--outdir /ruta/salida]
#
# El analisis solo se ejecuta si se proporcionan ambos archivos requeridos.
# No forma parte de 00_main.R: se corre bajo demanda.
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(args) {
  res <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!key %in% c("--fasta", "--products", "--outdir")) {
      stop("Argumento desconocido: ", key,
           "\nUso: Rscript 17_pca_familias.R --fasta archivo.fasta --products products.tsv [--outdir carpeta]")
    }
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      stop("Falta el valor para ", key)
    }
    res[[substring(key, 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  res
}

cfg <- parse_args(args)
if (is.null(cfg$fasta) || is.null(cfg$products)) {
  stop("Este analisis requiere que indiques --fasta y --products.\n",
       "Ejemplo: Rscript 17_pca_familias.R --fasta proteins_CM.fasta --products products.tsv")
}
if (!file.exists(cfg$fasta)) stop("No existe el FASTA: ", cfg$fasta)
if (!file.exists(cfg$products)) stop("No existe el archivo products: ", cfg$products)

paquetes <- c("Biostrings", "dplyr", "tidyr", "stringr", "ggplot2", "ggrepel",
              "Peptides", "RColorBrewer")
faltan <- paquetes[!vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltan)) stop("Faltan paquetes: ", paste(faltan, collapse = ", "),
                         "\nInstala los paquetes requeridos antes de correr este analisis.")
suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(Peptides)
  library(RColorBrewer)
})

fasta <- normalizePath(cfg$fasta)
products_file <- normalizePath(cfg$products)
outdir <- if (is.null(cfg$outdir)) dirname(fasta) else cfg$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outdir <- normalizePath(outdir)
setwd(outdir)

# ============================================================
# FUNCIÓN: Normalizar IDs de proteína
# ============================================================
normalizar_id <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- sub("\\..*", "", x)   # quita ".1:pep" o ".1"
  x <- sub(":.*", "", x)     # quita ":pep" si quedó
  x <- str_trim(x)
  return(x)
}

# ============================================================
# 1. LEER PROTEÍNAS Y CALCULAR COMPOSICIÓN
# ============================================================
cat("=== 1. LEYENDO PROTEÍNAS ===\n")
prot <- readAAStringSet(fasta)
cat("Secuencias leídas:", length(prot), "\n")

ids_prot <- names(prot)
cat("Ejemplo names(prot):\n")
print(head(ids_prot, 5))

# Frecuencias de aminoácidos (se reutilizan más abajo para aromaticidad)
aa_freq <- letterFrequency(prot, letters = AA_STANDARD, as.prob = TRUE)

df_aa <- as.data.frame(aa_freq)
df_aa$protein_id <- normalizar_id(ids_prot)

df_aa <- df_aa %>%
  select(protein_id, everything())

cat("\nProteínas en df_aa:", nrow(df_aa), "\n")
cat("Ejemplo IDs df_aa:\n")
print(head(df_aa$protein_id, 5))
cat("\n")

# ============================================================
# 2. LEER PRODUCTS DEL GFF3
# ============================================================
cat("=== 2. LEYENDO PRODUCTS ===\n")
products <- read.delim(products_file, header = FALSE,
                       col.names = c("protein_id", "product"),
                       stringsAsFactors = FALSE) %>%
  mutate(
    protein_id = normalizar_id(protein_id),
    product = URLdecode(product),
    product = str_remove(product, "^term="),
    product = str_split_i(product, ";", 1),
    product = str_remove(product, ", putative$"),
    product = str_remove(product, " putative$"),
    product = str_trim(product)
  )

cat("Products en df:", nrow(products), "\n")
cat("Ejemplo IDs df:\n")
print(head(products$protein_id, 5))
cat("\n")

# ============================================================
# 3. DIAGNÓSTICO DE COINCIDENCIAS
# ============================================================
cat("=== 3. DIAGNÓSTICO DE IDs ===\n")
comunes <- intersect(df_aa$protein_id, products$protein_id)
cat("IDs únicos en df_aa:", length(unique(df_aa$protein_id)), "\n")
cat("IDs únicos en df:", length(unique(products$protein_id)), "\n")
cat("IDs en común:", length(comunes), "\n\n")

if (length(comunes) == 0) {
  cat("⚠️  NO HAY COINCIDENCIAS. Ejemplos de cada lado:\n")
  cat("df_aa:", head(df_aa$protein_id, 3), "\n")
  cat("df:   ", head(products$protein_id, 3), "\n")
  stop("Revisa el formato de IDs. Los regex no están normalizando bien.")
} else {
  cat("✓ IDs normalizados correctamente\n\n")
}

# ============================================================
# 4. ASIGNAR FAMILIAS (v2: límites de palabra + categorías nuevas)
# ============================================================
cat("=== 4. ASIGNANDO FAMILIAS ===\n")
df <- products %>%
  mutate(familia = case_when(
    # ===== APICOMPLEXA-ESPECÍFICAS =====
    # OJO: \\b evita falsos positivos tipo "genomic"->MIC, "iron"/"environment"->RON,
    # "tropism"->ROP. El [0-9]* cubre nomenclatura tipo ROP18, GRA1, MIC2.
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
    
    # ===== NUEVO: tráfico / señalización de membrana =====
    str_detect(product, regex("SNARE|vesicle transport|clathrin|coatomer|dynamin|adaptin|\\bVPS[0-9]+\\b", ignore_case = TRUE)) ~ "Tráfico vesicular / SNARE",
    str_detect(product, regex("\\bRab\\b|Rho GTPase|Ras-related|ADP-ribosylation factor|\\bseptin", ignore_case = TRUE)) ~ "GTPasas pequeñas",
    str_detect(product, regex("apicoplast|\\bplastid\\b|ferredoxin|iron-sulfur|Fe-S cluster", ignore_case = TRUE)) ~ "Apicoplasto / Fe-S",
    
    # ===== REGULADORES =====
    str_detect(product, regex("AP2 domain transcription", ignore_case = TRUE)) ~ "AP2 factores transcripción",
    # Ciclo celular antes que quinasas genéricas, para no perder "cyclin-dependent kinase"
    str_detect(product, regex("\\bcyclin\\b|cell cycle|\\bcentrin\\b|aurora kinase", ignore_case = TRUE)) ~ "Ciclo celular",
    str_detect(product, regex("\\bkinase\\b|kinase-like|\\bNEK[0-9]*\\b", ignore_case = TRUE)) ~ "Quinasas",
    str_detect(product, regex("phosphatase", ignore_case = TRUE)) ~ "Fosfatasas",
    str_detect(product, regex("subtilisin|protease|peptidase|proteasome", ignore_case = TRUE)) ~ "Proteasas",
    
    # ===== NUEVO: replicación / reparación de DNA =====
    str_detect(product, regex("DNA polymerase|primase|\\bMCM[0-9]+\\b|mismatch repair|replication factor|replication protein", ignore_case = TRUE)) ~ "Replicación/reparación DNA",
    
    # ===== METABOLISMO =====
    str_detect(product, regex("ribosomal|ribosome", ignore_case = TRUE)) ~ "Ribosomales",
    str_detect(product, regex("elongation factor|initiation factor|translation|tRNA", ignore_case = TRUE)) ~ "Traducción/tRNA",
    str_detect(product, regex("RNA polymerase|mediator complex|transcription", ignore_case = TRUE)) ~ "RNA polimerasa / transcripción",
    str_detect(product, regex("cytochrome|ATP synthase|oxidase|reductase|dehydrogenase", ignore_case = TRUE)) ~ "Metabolismo energético",
    str_detect(product, regex("transporter|permease|carrier|channel|ABC transporter|facilitator", ignore_case = TRUE)) ~ "Transportadores",
    # ===== NUEVO: glicosilación / anclaje GPI =====
    str_detect(product, regex("glycosyltransferase|GPI anchor|glycosylphosphatidylinositol", ignore_case = TRUE)) ~ "Glicosilación / GPI",
    str_detect(product, regex("synthase|transferase|isomerase|lyase|ligase|mutase|hydrolase|phosphodiesterase|palmitoyltransferase", ignore_case = TRUE)) ~ "Metabolismo general",
    
    # ===== OTROS =====
    str_detect(product, regex("histone|methyltransferase.*histone|acetyltransferase|chromatin|SWI2/SNF2", ignore_case = TRUE)) ~ "Modificadores de cromatina",
    str_detect(product, regex("ubiquitin|HECT|OTU domain|carboxyl-terminal hydrolase", ignore_case = TRUE)) ~ "Ubiquitina / UPS",
    str_detect(product, regex("DnaJ|heat shock|Hsp|chaperone|thioredoxin", ignore_case = TRUE)) ~ "Chaperonas / Heat shock",
    str_detect(product, regex("tubulin|kinesin|dynein|myosin|actin|centrosomal|basal complex", ignore_case = TRUE)) ~ "Citoesqueleto / Motores",
    str_detect(product, regex("helicase|DEAD/DEAH", ignore_case = TRUE)) ~ "Helicasas",
    str_detect(product, regex("EF hand|calcium-dependent|calmodulin", ignore_case = TRUE)) ~ "Señalización calcio",
    str_detect(product, regex("MORN|kelch|RCC1|TBC domain|RAP domain|endonuclease|exonuclease", ignore_case = TRUE)) ~ "Otras enzimáticas",
    
    # ===== REPETICIONES (función desconocida) =====
    str_detect(product, regex("WD domain|RNA recognition motif|tetratricopeptide|leucine rich|HEAT repeat|ankyrin|zinc finger", ignore_case = TRUE)) ~ "Repeticiones (función desconocida)",
    
    # ===== HIPOTÉTICAS =====
    str_detect(product, regex("hypothetical protein, conserved|conserved hypothetical", ignore_case = TRUE)) ~ "Hipotéticas conservadas",
    str_detect(product, regex("hypothetical", ignore_case = TRUE)) ~ "Hipotéticas",
    
    # ===== RESTO =====
    TRUE ~ "Otras"
  ))

cat("\n=== DISTRIBUCIÓN DE FAMILIAS ===\n")
resumen_familias <- df %>%
  count(familia, sort = TRUE) %>%
  as.data.frame()
print(resumen_familias)
cat("\n")

# ============================================================
# 5. UNIR CON COMPOSICIÓN AMINOACÍDICA
# ============================================================
cat("=== 5. UNIENDO df_aa CON df ===\n")
df_full <- df_aa %>%
  inner_join(df, by = "protein_id")

cat("Filas en df_full:", nrow(df_full), "\n\n")

if (nrow(df_full) == 0) {
  stop("El join sigue dando 0. Revisa el diagnóstico de IDs del paso 3.")
}

# ============================================================
# 6. FILTRAR PARA PCA (excluir hipotéticas)
# ============================================================
df_pca <- df_full %>%
  filter(!familia %in% c("Hipotéticas", "Hipotéticas conservadas"))

cat("=== 6. PROTEÍNAS PARA PCA ===\n")
cat("Total (sin hipotéticas):", nrow(df_pca), "\n\n")

# ============================================================
# CONTINUACIÓN DEL SCRIPT: PCA de propiedades fisicoquímicas
# Requiere haber corrido antes: pasos 1-6 (lectura de fasta,
# products, asignación de familias, df_pca)
# ============================================================

# --- LIBRERÍA ADICIONAL (agregar junto a las demás al inicio) ---
library(RColorBrewer)

# ============================================================
# 7. CALCULAR PROPIEDADES FISICOQUÍMICAS (paquete Peptides)
# ============================================================
cat("=== 7. CALCULANDO PROPIEDADES FISICOQUÍMICAS ===\n")

# Vector de secuencias como character, nombrado por protein_id normalizado
seqs_chr <- as.character(prot)
names(seqs_chr) <- normalizar_id(names(prot))

# Nos quedamos solo con las proteínas que van al PCA (sin hipotéticas)
seqs_pca <- seqs_chr[df_pca$protein_id]

# Chequeo de secuencias válidas
seq_ok <- !is.na(seqs_pca) & nchar(seqs_pca) > 0
if (any(!seq_ok)) {
  cat("⚠️  Proteínas sin secuencia válida, se excluyen:", sum(!seq_ok), "\n")
}
df_pca   <- df_pca[seq_ok, ]
seqs_pca <- seqs_pca[seq_ok]

# Peptides::* está vectorizado: acepta el vector completo de secuencias
physchem <- data.frame(
  protein_id       = names(seqs_pca),
  mw               = mw(seqs_pca, monoisotopic = FALSE),
  pI               = pI(seqs_pca, pKscale = "EMBOSS"),
  charge_pH7       = charge(seqs_pca, pH = 7, pKscale = "EMBOSS"),
  hydrophobicity   = hydrophobicity(seqs_pca, scale = "KyteDoolittle"),
  instability      = instaIndex(seqs_pca),
  aliphatic_index  = aIndex(seqs_pca),
  boman_index      = boman(seqs_pca),
  stringsAsFactors = FALSE
)
rownames(physchem) <- NULL

# Aromaticidad (Phe + Tyr + Trp) a partir de las frecuencias ya calculadas
df_arom <- df_aa %>%
  transmute(protein_id, aromaticity = F + Y + W)

physchem <- physchem %>%
  left_join(df_arom, by = "protein_id")

cat("Propiedades calculadas para", nrow(physchem), "proteínas\n\n")

# ============================================================
# 8. UNIR PROPIEDADES FISICOQUÍMICAS CON FAMILIAS
# ============================================================
df_pca_full <- df_pca %>%
  select(protein_id, familia) %>%
  inner_join(physchem, by = "protein_id")

cat("=== 8. DATOS FINALES PARA PCA ===\n")
cat("Filas:", nrow(df_pca_full), "\n\n")

# ============================================================
# 9. PCA SOBRE PROPIEDADES FISICOQUÍMICAS
# ============================================================
cat("=== 9. CALCULANDO PCA ===\n")

vars_pca <- c("mw", "pI", "charge_pH7", "hydrophobicity",
              "instability", "aliphatic_index", "boman_index", "aromaticity")

mat_pca <- df_pca_full %>%
  select(all_of(vars_pca)) %>%
  as.matrix()

# Quitar filas con NA (por si algún cálculo dio NA)
filas_ok <- complete.cases(mat_pca)
if (any(!filas_ok)) {
  cat("⚠️  Se excluyen", sum(!filas_ok), "proteínas con NA en propiedades\n")
}
df_pca_full <- df_pca_full[filas_ok, ]
mat_pca     <- mat_pca[filas_ok, ]

pca_res <- prcomp(mat_pca, center = TRUE, scale. = TRUE)

# Varianza explicada
var_exp <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
pc1_lab <- paste0("PC1 (", round(var_exp[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(var_exp[2], 1), "%)")

cat("Varianza explicada PC1:", round(var_exp[1], 1), "%\n")
cat("Varianza explicada PC2:", round(var_exp[2], 1), "%\n\n")

# Data frame con los scores del PCA + metadatos
df_scores <- df_pca_full %>%
  select(protein_id, familia) %>%
  bind_cols(as.data.frame(pca_res$x[, 1:2]))

# Rangos comunes de ejes (para que las 3 gráficas sean comparables)
rango_x <- range(df_scores$PC1) * 1.05
rango_y <- range(df_scores$PC2) * 1.05

# Paleta con un color por familia (se reutiliza en gráficas 2 y 3)
familias_orden <- sort(unique(df_scores$familia))
n_fam <- length(familias_orden)
paleta_fam <- colorRampPalette(brewer.pal(12, "Paired"))(n_fam)
names(paleta_fam) <- familias_orden

# ============================================================
# 10. GRÁFICA 1: PCA sin colorear (nube gris claro)
# ============================================================
cat("=== 10. GRÁFICA 1: nube gris ===\n")

p1 <- ggplot(df_scores, aes(x = PC1, y = PC2)) +
  geom_point(color = "grey65", size = 1.4, alpha = 0.75, stroke = 0) +
  coord_equal(xlim = rango_x, ylim = rango_y) +
  labs(
    title = "PCA de propiedades fisicoquímicas - T. gondii TgRH_pasteur",
    x = pc1_lab,
    y = pc2_lab
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title   = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.title   = element_text(face = "bold")
  )

ggsave("PCA_01_nube_gris.png", p1, width = 6, height = 5.5, dpi = 600)
print(p1)

# ============================================================
# 11. GRÁFICA 2: FACET por familia (fondo gris + familia coloreada)
# ============================================================
cat("=== 11. GRÁFICA 2: facet por familia ===\n")

# df sin la columna 'familia' -> esta capa se repite en TODOS los paneles
df_bg <- df_scores %>% select(-familia)

p2 <- ggplot(df_scores, aes(x = PC1, y = PC2)) +
  geom_point(data = df_bg, color = "grey85", size = 0.8, alpha = 0.5, stroke = 0) +
  geom_point(aes(color = familia), size = 1.1, alpha = 0.9, stroke = 0) +
  facet_wrap(~ familia, ncol = 5) +           # scales = "fixed" por defecto -> mismos rangos
  scale_color_manual(values = paleta_fam, guide = "none") +
  coord_equal(xlim = rango_x, ylim = rango_y) +
  labs(
    title = "PCA de propiedades fisicoquímicas por familia funcional",
    x = pc1_lab,
    y = pc2_lab
  ) +
  theme_bw(base_size = 9) +
  theme(
    strip.text        = element_text(face = "bold", size = 7),
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.text         = element_text(size = 6)
  )

n_filas_facet <- ceiling(n_fam / 5)
ggsave("PCA_02_facet_familias.pdf", p2,
       width = 14, height = n_filas_facet * 2.6, dpi = 600, limitsize = FALSE)
ggsave("PCA_02_facet_familias.png", p2,
       width = 14, height = n_filas_facet * 2.6, dpi = 600, limitsize = FALSE)
print(p2)

# ============================================================
# 12. GRÁFICA 3: PCA completo coloreado por familia
# ============================================================
cat("=== 12. GRÁFICA 3: PCA completo coloreado ===\n")

p3 <- ggplot(df_scores, aes(x = PC1, y = PC2, color = familia)) +
  geom_point(size = 1.6, alpha = 0.85, stroke = 0) +
  scale_color_manual(values = paleta_fam, name = "Familia funcional") +
  coord_equal(xlim = rango_x, ylim = rango_y) +
  labs(
    title = "PCA de propiedades fisicoquímicas coloreado por familia",
    x = pc1_lab,
    y = pc2_lab
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.title    = element_text(face = "bold"),
    legend.text   = element_text(size = 7),
    legend.title  = element_text(face = "bold", size = 9),
    legend.key.size = unit(0.35, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))

ggsave("PCA_03_completo_por_familia.png", p3, width = 9, height = 6.5, dpi = 600)
ggsave("PCA_03_completo_por_familia.pdf", p3, width = 9, height = 6.5, dpi = 600)
print(p3)

cat("\n✓ Listo. Archivos generados:\n",
    " - PCA_01_nube_gris.png\n",
    " - PCA_02_facet_familias.pdf / .png\n",
    " - PCA_03_completo_por_familia.png / .pdf\n",
    " - PCA_04_elipses_confianza.png / .pdf\n",
    " - tabla_estadisticas_por_familia.csv\n")

# ============================================================
# CONTINUACIÓN DEL SCRIPT: PCA de propiedades fisicoquímicas
# Requiere haber corrido antes: pasos 1-6 (lectura de fasta,
# products, asignación de familias, df_pca)
# ============================================================

# --- LIBRERÍA ADICIONAL (agregar junto a las demás al inicio) ---
library(RColorBrewer)


# ============================================================
# 13. GRÁFICA 4: PCA completo con elipses de confianza (95%) por familia
# ============================================================
cat("=== 13. GRÁFICA 4: PCA con elipses de confianza ===\n")

# stat_ellipse (tipo "t") necesita al menos 4 puntos por grupo para
# poder estimarse; identificamos qué familias cumplen esa condición
n_por_familia <- df_scores %>% count(familia, name = "n")
familias_con_elipse <- n_por_familia %>% filter(n >= 4) %>% pull(familia)
familias_sin_elipse <- n_por_familia %>% filter(n < 4) %>% pull(familia)

if (length(familias_sin_elipse) > 0) {
  cat("⚠️  Familias con n < 4, sin elipse (solo se grafican sus puntos):\n")
  print(n_por_familia %>% filter(n < 4))
}

df_elipse <- df_scores %>% filter(familia %in% familias_con_elipse)

p4 <- ggplot(df_scores, aes(x = PC1, y = PC2, color = familia)) +
  geom_point(size = 1.4, alpha = 0.55, stroke = 0) +
  stat_ellipse(
    data = df_elipse,
    aes(group = familia),
    type = "t", level = 0.95,
    linewidth = 0.5, alpha = 0.9,
    show.legend = FALSE
  ) +
  scale_color_manual(values = paleta_fam, name = "Familia funcional") +
  coord_equal(xlim = rango_x, ylim = rango_y) +
  labs(
    title    = "PCA de propiedades fisicoquímicas con elipses de confianza (95%) por familia",
    subtitle = paste0("Elipses calculadas solo para familias con n \u2265 4 (", 
                      length(familias_con_elipse), " de ", nrow(n_por_familia), " familias)"),
    x = pc1_lab,
    y = pc2_lab
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey30"),
    axis.title    = element_text(face = "bold"),
    legend.text   = element_text(size = 7),
    legend.title  = element_text(face = "bold", size = 9),
    legend.key.size = unit(0.35, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))

ggsave("PCA_04_elipses_confianza.png", p4, width = 9.5, height = 6.5, dpi = 600)
ggsave("PCA_04_elipses_confianza.pdf", p4, width = 9.5, height = 6.5, dpi = 600)
print(p4)

# ============================================================
# 14. TABLA FINAL: ESTADÍSTICOS DESCRIPTIVOS POR FAMILIA
# ============================================================
cat("\n=== 14. TABLA DE ESTADÍSTICOS POR FAMILIA ===\n")

# Unimos scores de PCA (PC1, PC2) con las propiedades fisicoquímicas originales
df_stats_full <- df_pca_full %>%
  inner_join(df_scores %>% select(protein_id, PC1, PC2), by = "protein_id")

tabla_stats <- df_stats_full %>%
  group_by(familia) %>%
  summarise(
    n = n(),
    across(
      c(PC1, PC2, all_of(vars_pca)),
      list(media = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(n)) %>%
  mutate(across(-c(familia, n), ~round(.x, 2)))

# Nota: familias con n = 1 tendrán sd = NA (no calculable con un solo dato)

cat("Tabla generada con", nrow(tabla_stats), "familias y", ncol(tabla_stats), "columnas\n")
print(tabla_stats)

write.csv(tabla_stats, "tabla_estadisticas_por_familia.csv", row.names = FALSE)

cat("\n✓ Listo. Archivos generados:\n",
    " - PCA_01_nube_gris.png\n",
    " - PCA_02_facet_familias.pdf / .png\n",
    " - PCA_03_completo_por_familia.png / .pdf\n",
    " - PCA_04_elipses_confianza.png / .pdf\n",
    " - tabla_estadisticas_por_familia.csv\n")


# PERMANOVA y betadisper fueron retirados: este script genera solo PCA y descriptivos.
loadings <- as.data.frame(pca_res$rotation[, 1:2])
loadings$variable <- rownames(loadings)
loadings <- loadings[order(-abs(loadings$PC1)), c("variable", "PC1", "PC2")]
print(loadings)

