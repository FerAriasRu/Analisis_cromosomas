# Análisis de cromosomas y proteínas

Pipeline en R para analizar la composición, anotación y organización de uno o varios genomas en ventanas cromosómicas. También incluye módulos complementarios para preparar y analizar secuencias de proteínas y familias proteicas.

El análisis genómico se coordina desde `00_main.R`. El análisis de proteínas se ejecuta mediante los scripts `17_preparar_proteinas.R` y `18_graficas_pca_proteinas.R`, según estén configurados en el proyecto.

## Qué analiza

### Análisis genómico

- Composición de bases por ventana: %GC, %AT y contenido de N.
- Densidad génica (genes/Mb), cobertura de genes, exones, intrones, promotores e intergénico.
- Repeticiones, satélites y microsatélites cuando hay datos o herramientas disponibles.
- Telómeros y regiones subteloméricas.
- Estadísticas por cromosoma, extremos, zonas con valores extremos y genes que caen en esas zonas.
- Comparaciones entre genomas configurados.

El análisis genómico completo requiere un FASTA y su anotación GFF3. El modo ligero admite FASTA sin GFF3 y sirve para análisis comparativos que no requieren anotación génica.

### Análisis de proteínas

Según los datos y análisis habilitados en los scripts, puede incluir:

- Análisis PCA de las proteínas y gráficos por familia.
- Localización de proteínas seleccionadas a lo largo de los cromosomas.
- Composición de aminoácidos por proteína y familia.
- Boxplots de composición de aminoácidos por familia y por genoma.
- Comparación estadística de aminoácidos entre familias mediante Kruskal-Wallis y pruebas Dunn, con corrección por comparaciones múltiples.
- Resúmenes del tamaño del efecto.
- Heatmaps de la mediana de composición de aminoácidos por familia.

Las familias y los genomas incluidos se seleccionan en los parámetros o en los scripts correspondientes.

## Cómo ejecutar

1. Instala R y RStudio.
2. Instala los paquetes R requeridos (ver [Dependencias](#dependencias)).
3. Edita `01_parametros.R`: define las rutas de entrada en `genomas`, selecciona los genomas en `genomas_run` y ajusta los parámetros necesarios.
4. Abre el proyecto `Cromosomas4.Rproj` en RStudio y ejecuta `00_main.R` con **Source** para correr el pipeline genómico.

Para correr solo hasta un módulo, asigna su nombre a `hasta` en `00_main.R`, por ejemplo `"13_zonas_criticas.R"`. El valor `NULL` corre el pipeline completo. `saltar` permite omitir módulos, pero solo debe usarse si los módulos siguientes no dependen de sus resultados.

### Ejecutar el análisis de proteínas por separado

El análisis de proteínas se puede ejecutar sin correr todo `00_main.R`. En una sesión de R, configura los parámetros y el genoma activo y ejecuta los módulos de proteínas en este orden:

```r
source("01_parametros.R")

# Selecciona los genomas que se analizarán.
genomas_run_main <- c("RH", "RH88")

# Selecciona el genoma activo para la preparación.
genoma_activo <- "RH88"
source("01b_genoma_activo.R")

# Prepara los datos de proteínas del genoma activo.
source("17_preparar_proteinas.R")

# Genera las figuras y tablas a partir de los datos preparados.
source("18_graficas_pca_proteinas.R")
