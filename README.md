# Análisis de cromosomas

Pipeline en R para analizar composición, anotación y organización de uno o varios genomas en ventanas cromosómicas. Genera tablas, un libro Excel y figuras por genoma; si se seleccionan varios genomas, también calcula comparaciones.

## Qué analiza

- Composición de bases por ventana: %GC, %AT y contenido de N.
- Densidad génica (genes/Mb), cobertura de genes, exones, intrones, promotores e intergénico.
- Repeticiones, satélites y microsatélites cuando hay datos o herramientas disponibles.
- Telómeros y regiones subteloméricas.
- Estadísticas por cromosoma, extremos, zonas con valores extremos y genes que caen en esas zonas.
- Comparaciones entre genomas configurados.

El análisis completo requiere un FASTA y su anotación GFF3. El modo ligero admite FASTA sin GFF3 y sirve para los análisis comparativos que no requieren anotación génica.

## Cómo ejecutar

1. Instala R y RStudio.
2. Instala los paquetes R requeridos (ver [Dependencias](#dependencias)).
3. Edita `01_parametros.R`: define las rutas de entrada en `genomas`, selecciona los genomas en `genomas_run` y ajusta los parámetros necesarios.
4. Abre el proyecto `Cromosomas4.Rproj` en RStudio y ejecuta `00_main.R` con **Source**.

Para correr solo hasta un módulo, asigna su nombre a `hasta` en `00_main.R`, por ejemplo `"13_zonas_criticas.R"`. El valor `NULL` corre el pipeline completo. `saltar` permite omitir módulos, pero solo debe usarse si los módulos siguientes no dependen de sus resultados.

## Archivos de entrada

Las rutas se especifican por genoma en la lista `genomas` de `01_parametros.R`:

- `fna`: FASTA del genoma (requerido).
- `gff`: anotación GFF3; al incluirla se usa el modo completo.
- `fai`: índice FASTA opcional; si no está disponible se calculan las longitudes desde el FASTA.
- `rep`: salida de RepeatMasker opcional.
- `trf`: salida de Tandem Repeats Finder opcional.
- `eq_manual`: tabla opcional con las columnas `acc` y `nombre` para asignar nombres y orden cromosómico explícitos.

Sin `eq_manual`, las secuencias se seleccionan mediante `min_len_cromosoma`, `patron_incluir` y `patron_excluir`. Revisa que los nombres de secuencia del GFF3 coincidan con los del FASTA; `mapa_gff_a_fasta` permite traducirlos cuando sea necesario.

Los datos genómicos no forman parte de este repositorio: configura rutas a archivos que tengas disponibles localmente.

## Ventanas y zonas críticas

El tamaño predeterminado de ventana es 50.000 pb (`bin <- 50000`). Cada ventana se evalúa con filtros de validez, contenido de N y, opcionalmente, exclusión de outliers MAD y ventanas terminales técnicas.

El parámetro `umbral_pct <- 0.15` de `01_parametros.R` es una fracción para calcular percentiles **por cromosoma** sobre las ventanas que quedan en el conjunto de análisis:

- **P15 de GC:** una ventana se marca como `Bajo GC` si su GC es menor o igual al valor GC del percentil 15 de ese cromosoma.
- **P15 de genes/Mb:** se marca como `Baja densidad génica` si su densidad es menor o igual al valor correspondiente al percentil 15.
- Los criterios se calculan por separado. Una ventana puede cumplir uno o ambos.
- El módulo `09_filtro_banderas.R` usa también P85 (el percentil `1 - umbral_pct`) para banderas altas de repeticiones, intrones e intergénico. No crea una bandera de GC alto.

P15 **no significa 15% de GC**: significa el valor de corte que deja aproximadamente el 15% de las ventanas del cromosoma en ese extremo bajo. El corte numérico de GC puede ser, por ejemplo, 52% en un cromosoma y otro valor en otro cromosoma.

El módulo `13_zonas_criticas.R` combina ventanas contiguas con las mismas banderas y crea las categorías `Bajo GC`, `Baja densidad génica` y `Bajo GC + baja densidad génica`. Por tanto, una zona etiquetada `Baja densidad génica` puede tener GC alto: su clasificación se debe a genes/Mb. Las franjas de zonas críticas que aparecen en las figuras proceden de esas categorías; la línea superpuesta representa la variable indicada en el eje y puede mostrar valores altos dentro de una franja clasificada por otra variable. Las tablas `Umbrales_por_cromosoma` y `Zonas_*` del Excel permiten consultar los cortes y resúmenes numéricos.

## Estructura de los scripts

`00_main.R` coordina los módulos. En el modo completo, el flujo principal es:

1. `01_parametros.R` y `01b_genoma_activo.R`: parámetros y selección del genoma.
2. `02_datos_base.R` y `03_intergenico.R`: lectura de FASTA/GFF3 y objetos genómicos base.
3. `04_ventanas.R` a `08_intergenico_metricas.R`: ventanas, extremos, telómeros, subtelómeros, elementos repetitivos y métricas.
4. `09_filtro_banderas.R` a `13_zonas_criticas.R`: filtros, análisis, estadística y zonas.
5. `14_excel.R`, `90_guardar_genoma.R` y `15_figuras.R`: exportación, guardado de objetos y figuras.
6. `16_comparacion.R` y `17_pca_familias.R`: módulos de comparación y análisis complementarios.

El modo ligero usa `02l_base_ligera.R`, `06b_sat350.R`, `07b_microsatelites.R` y `15F08_idiograma_ligero.R`, sin requerir anotación GFF3.

## Resultados

Los resultados se guardan en carpetas `resultados_<genoma>_<ventana>/` (el nombre y la ubicación dependen de la configuración). En el modo completo pueden incluir:

- `analisis_genomico_<ventana>.xlsx`: tablas de ventanas, umbrales, zonas, genes asociados y resultados de los módulos.
- Figuras en PNG y PDF, más `leyendas_figuras.md`.
- Objetos R serializados para reutilizar resultados en comparaciones posteriores.
- Archivos auxiliares de herramientas externas en carpetas `externos_<genoma>/`.

Las comparaciones se guardan en `comparacion_genomas/` cuando hay al menos dos genomas disponibles según la configuración.

## Dependencias

Paquetes R requeridos por `00_main.R`: `GenomicFeatures`, `GenomicRanges`, `GenomeInfoDb`, `Biostrings`, `rtracklayer`, `txdbmaker`, `ggplot2`, `patchwork`, `dplyr`, `openxlsx`, `BSgenome` y `scales`. Se pueden instalar desde CRAN o Bioconductor según corresponda; para paquetes Bioconductor se puede usar `BiocManager::install()`.

Herramientas externas como TRF, BLASTN, seqkit y RepeatMasker se usan en módulos concretos y dependen de los parámetros y archivos configurados. Las repeticiones también pueden obtenerse de regiones en minúscula (softmask) del FASTA si `usar_softmask <- TRUE`.

## Reproducibilidad

Guarda junto a los resultados la configuración utilizada, en especial tamaño de ventana, selección de cromosomas, filtros y `umbral_pct`. El Excel generado registra los parámetros principales y las fuentes de datos detectadas durante la corrida.
