# Análisis de cromosomas y proteínas

Pipeline en R para analizar la composición, anotación y organización de uno o varios genomas en ventanas cromosómicas. Incluye módulos complementarios para preparar y analizar secuencias de proteínas y familias proteicas.

El análisis genómico se coordina desde `00_main.R`. El análisis de proteínas se ejecuta mediante `17_preparar_proteinas.R` y `18_graficas_pca_proteinas.R`.

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

- Análisis PCA de proteínas y gráficos por familia.
- Localización de proteínas seleccionadas a lo largo de los cromosomas.
- Composición de aminoácidos por proteína y familia.
- Boxplots de composición de aminoácidos por familia y por genoma.
- Comparación estadística entre familias mediante Kruskal-Wallis y pruebas Dunn, con corrección por comparaciones múltiples.
- Resúmenes del tamaño del efecto.
- Heatmaps de la mediana de composición de aminoácidos por familia.

Las familias y los genomas incluidos se seleccionan en los parámetros o en los scripts correspondientes.

## Cómo ejecutar

1. Instala R y RStudio.
2. Instala los paquetes R requeridos (ver [Dependencias](#dependencias)).
3. Edita `01_parametros.R`: define las rutas de entrada en `genomas`, selecciona los genomas en `genomas_run` y ajusta los parámetros necesarios.
4. Abre `Cromosomas4.Rproj` en RStudio y ejecuta `00_main.R` con **Source** para correr el pipeline genómico.

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

# Prepara los datos del genoma activo.
source("17_preparar_proteinas.R")

# Genera las figuras y tablas a partir de los datos preparados.
source("18_graficas_pca_proteinas.R")
```

Si cambias los archivos de proteínas, las anotaciones o el cálculo de métricas en el módulo 17, vuelve a ejecutar `17_preparar_proteinas.R` antes de `18_graficas_pca_proteinas.R`.

Si el módulo 17 procesa solamente el genoma activo, repite la preparación para cada genoma que quieras analizar. La lista de genomas que se grafica debe coincidir con los genomas cuyos datos fueron preparados.

## Archivos de entrada

Las rutas genómicas se especifican por genoma en la lista `genomas` de `01_parametros.R`:

- `fna`: FASTA del genoma (requerido para el análisis genómico).
- `gff`: anotación GFF3; al incluirla se usa el modo completo.
- `fai`: índice FASTA opcional; si no está disponible, se calculan las longitudes desde el FASTA.
- `rep`: salida de RepeatMasker opcional.
- `trf`: salida de Tandem Repeats Finder opcional.
- `eq_manual`: tabla opcional con las columnas `acc` y `nombre` para asignar nombres y orden cromosómico explícitos.

Sin `eq_manual`, las secuencias se seleccionan mediante `min_len_cromosoma`, `patron_incluir` y `patron_excluir`. Revisa que los nombres de secuencia del GFF3 coincidan con los del FASTA; `mapa_gff_a_fasta` permite traducirlos cuando sea necesario.

Los datos genómicos no forman parte de este repositorio: configura rutas a archivos disponibles localmente.

### Archivos de proteínas

Los archivos de proteínas se configuran por genoma en `01_parametros.R`. Según la configuración de los scripts, pueden incluir:

- `protein_fasta`: FASTA de secuencias de proteínas. Si se deja como `NULL`, el módulo de preparación puede generar el FASTA a partir del genoma y la anotación, siempre que estén disponibles los archivos y herramientas necesarios.
- `product_tsv`: tabla opcional de productos o descripciones proteicas, si está configurada para complementar la asignación de familias.

Cuando las proteínas se generan desde la anotación, el identificador del FASTA debe poder relacionarse con el identificador correspondiente del GFF3. Estos identificadores pueden variar entre genomas. Por ejemplo, un GFF3 puede utilizar el identificador del transcrito en `mRNA`, mientras que otro puede presentar atributos distintos. El módulo de preparación debe extraer el identificador apropiado para cada formato.

Si el archivo de proteínas ya existe, se puede configurar directamente como `protein_fasta`. Si no existe, se deben proporcionar los archivos de entrada necesarios para que el script 17 lo genere.

## Ventanas y zonas críticas

El tamaño predeterminado de ventana es 50.000 pb (`bin <- 50000`). Cada ventana se evalúa con filtros de validez, contenido de N y, opcionalmente, exclusión de outliers MAD y ventanas terminales técnicas.

El parámetro `umbral_pct <- 0.15` de `01_parametros.R` es una fracción para calcular percentiles **por cromosoma** sobre las ventanas que quedan en el conjunto de análisis:

- **P15 de GC:** una ventana se marca como `Bajo GC` si su GC es menor o igual al valor GC del percentil 15 de ese cromosoma.
- **P15 de genes/Mb:** se marca como `Baja densidad génica` si su densidad es menor o igual al valor correspondiente al percentil 15.
- Los criterios se calculan por separado. Una ventana puede cumplir uno o ambos.
- El módulo `09_filtro_banderas.R` usa también P85 (el percentil `1 - umbral_pct`) para banderas altas de repeticiones, intrones e intergénico. No crea una bandera de GC alto.

P15 **no significa 15% de GC**: significa el valor de corte que deja aproximadamente el 15% de las ventanas del cromosoma en ese extremo bajo. El corte numérico de GC puede ser, por ejemplo, 52% en un cromosoma y otro valor en otro cromosoma.

El módulo `13_zonas_criticas.R` combina ventanas contiguas con las mismas banderas y crea las categorías `Bajo GC`, `Baja densidad génica` y `Bajo GC + baja densidad génica`. Una zona etiquetada `Baja densidad génica` puede tener GC alto: su clasificación se debe a genes/Mb. Las franjas de zonas críticas en las figuras proceden de esas categorías; la línea superpuesta representa la variable indicada en el eje y y puede mostrar valores altos dentro de una franja clasificada por otra variable. Las tablas `Umbrales_por_cromosoma` y `Zonas_*` del Excel permiten consultar los cortes y resúmenes numéricos.

## Composición de aminoácidos y comparaciones entre familias

La composición de aminoácidos se calcula a partir de las secuencias proteicas asignadas a cada familia. Los resultados se pueden resumir por genoma y familia mediante boxplots y tablas.

Las comparaciones entre familias pueden usar Kruskal-Wallis para cada aminoácido y pruebas Dunn para comparaciones por pares. La corrección por comparaciones múltiples —por ejemplo, Bonferroni o FDR— debe indicarse en los resultados. Los valores de significación se interpretan junto con una medida de tamaño del efecto, porque muestras grandes pueden producir valores p pequeños incluso ante diferencias reducidas.

El heatmap de medianas resume la composición de aminoácidos de las familias. Si se aplica un z-score por aminoácido, los colores muestran diferencias relativas entre familias para cada aminoácido, no porcentajes absolutos.

## Entropía y regiones de baja complejidad

La entropía de Shannon puede utilizarse para resumir la diversidad de aminoácidos de una secuencia o de ventanas dentro de una proteína:

\[
H = -\sum_i p_i \log_2(p_i)
\]

donde \(p_i\) es la frecuencia del aminoácido \(i\). Una entropía baja indica que la región está dominada por pocos aminoácidos; por ejemplo, puede reflejar repeticiones o enriquecimiento de N, Q, S o E. Para comparar secuencias, la entropía puede normalizarse con respecto a la entropía máxima de los aminoácidos considerados.

SEG y fLPS2 son herramientas para detectar regiones de baja complejidad. Sus resultados dependen del método y de los parámetros usados, por lo que esos parámetros deben registrarse junto con las tablas y figuras.

La entropía de Shannon y la detección mediante SEG o fLPS2 deben describirse como análisis disponibles solo si están incorporados y activados en los scripts. Son métodos relacionados, pero distintos; si se usan ambos, conserva y reporta sus resultados por separado.

## Estructura de los scripts

`00_main.R` coordina los módulos genómicos. En el modo completo, el flujo principal es:

1. `01_parametros.R` y `01b_genoma_activo.R`: parámetros y selección del genoma.
2. `02_datos_base.R` y `03_intergenico.R`: lectura de FASTA/GFF3 y objetos genómicos base.
3. `04_ventanas.R` a `08_intergenico_metricas.R`: ventanas, extremos, telómeros, subtelómeros, elementos repetitivos y métricas.
4. `09_filtro_banderas.R` a `13_zonas_criticas.R`: filtros, análisis, estadística y zonas.
5. `14_excel.R`, `90_guardar_genoma.R` y `15_figuras.R`: exportación, guardado de objetos y figuras.
6. `16_comparacion.R`: comparación entre genomas.

El análisis de proteínas se organiza en módulos separados:

- `17_preparar_proteinas.R`: carga o generación de datos proteicos y preparación de tablas de análisis.
- `18_graficas_pca_proteinas.R`: generación de las figuras y resúmenes a partir de los datos preparados por el módulo 17.

El modo ligero usa `02l_base_ligera.R`, `06b_sat350.R`, `07b_microsatelites.R` y `15F08_idiograma_ligero.R`, sin requerir anotación GFF3.

## Resultados

Los resultados genómicos se guardan en carpetas `resultados_<genoma>_<ventana>/`; el nombre y la ubicación dependen de la configuración. En el modo completo pueden incluir:

- `analisis_genomico_<ventana>.xlsx`: tablas de ventanas, umbrales, zonas, genes asociados y resultados de los módulos.
- Figuras en PNG y PDF, además de `leyendas_figuras.md`.
- Objetos R serializados para reutilizar resultados en comparaciones posteriores.
- Archivos auxiliares de herramientas externas en carpetas `externos_<genoma>/`.

Las comparaciones genómicas se guardan en `comparacion_genomas/` cuando hay al menos dos genomas disponibles según la configuración.

Los resultados de proteínas se guardan dentro de la carpeta de resultados del genoma, en una subcarpeta como `proteinas_RH` o `proteinas_RH88`, según la configuración. Pueden incluir:

- Objetos R con datos preparados para las figuras.
- Figuras PCA y figuras por familia.
- Gráficos de localización cromosómica.
- Boxplots y heatmaps de composición de aminoácidos.
- Tablas de resultados estadísticos y tamaños del efecto.
- Tablas de entropía o regiones de baja complejidad, si esos análisis están implementados y activados.

## Dependencias

### Paquetes de R para el análisis genómico

Los paquetes requeridos por `00_main.R` incluyen:

- `GenomicFeatures`
- `GenomicRanges`
- `GenomeInfoDb`
- `Biostrings`
- `rtracklayer`
- `txdbmaker`
- `ggplot2`
- `patchwork`
- `dplyr`
- `openxlsx`
- `BSgenome`
- `scales`

### Paquetes de R para proteínas y gráficos

Los módulos de proteínas pueden requerir:

- `Biostrings`, para leer secuencias de proteínas.
- `ggplot2`, para generar figuras.
- `dplyr` y `tidyr`, para manipular y reorganizar tablas.
- `RColorBrewer`, si se usa para las paletas de colores.
- `stats`, que forma parte de R y proporciona funciones estadísticas como `kruskal.test`, `pairwise.wilcox.test`, `dist` y `hclust`.

Instala los paquetes que aparecen en las llamadas `library()` o `requireNamespace()` de los módulos que ejecutes. Los paquetes de CRAN se pueden instalar con `install.packages()`. Para paquetes de Bioconductor, usa `BiocManager::install()`.

Ejemplo:

```r
install.packages(c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "RColorBrewer",
  "patchwork",
  "openxlsx",
  "scales"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "GenomicFeatures",
  "GenomicRanges",
  "GenomeInfoDb",
  "Biostrings",
  "rtracklayer",
  "txdbmaker",
  "BSgenome"
))
```

### Herramientas externas

Herramientas externas como TRF, BLASTN, seqkit y RepeatMasker se usan en módulos concretos y dependen de los parámetros y archivos configurados. Las repeticiones también pueden obtenerse de regiones en minúscula (softmask) del FASTA si `usar_softmask <- TRUE`.

`gffread` puede ser necesario para generar secuencias de proteínas a partir del FASTA genómico y la anotación. SEG y fLPS2 son opcionales y se requieren únicamente si los scripts invocan esas herramientas para detectar regiones de baja complejidad.

Las herramientas externas deben estar instaladas y disponibles desde el sistema si el módulo correspondiente las ejecuta.

## Reproducibilidad

Guarda junto a los resultados la configuración utilizada, en especial:

- Genomas y familias seleccionados.
- Rutas y versiones de los FASTA, GFF3 y archivos de proteínas.
- Tamaño de ventana y selección de cromosomas.
- Filtros y valor de `umbral_pct`.
- Parámetros de composición y análisis estadísticos.
- Tamaño y paso de ventana y umbral de entropía, si se calculan regiones de baja entropía.
- Herramienta y parámetros de detección, si se utiliza SEG o fLPS2.

El Excel generado registra los parámetros principales y las fuentes de datos detectadas durante una corrida genómica.
