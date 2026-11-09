# ============================================================
# 18_graficas_pca_proteinas.R
# Lee los RDS preparados por 17_preparar_proteinas.R y genera:
#   - PCA por genoma
#   - Boxplots de composición de aminoácidos por genoma/familia
#   - Mapas de ubicación cromosómica de las familias seleccionadas
#
# Rutas configuradas en 01_parametros.R.
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

if (!exists("familias_localizacion", inherits = TRUE) ||
    length(familias_localizacion) == 0) {
  stop(
    "Define familias_localizacion en 01_parametros.R. ",
    "Ejemplo: c(\"SAG/SRS (superficie)\", \"ROP (roptrias)\")."
  )
}

paquetes <- c("ggplot2", "dplyr", "tidyr", "RColorBrewer", "Biostrings")
faltan <- paquetes[
  !vapply(paquetes, requireNamespace, logical(1), quietly = TRUE)
]

if (length(faltan)) {
  stop(
    "Faltan paquetes R requeridos: ",
    paste(faltan, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
})

# Genomas que se procesarán. Para comparar RH y RH88,
# genomas_run_main o genomas_run debe contener ambos.
if (exists("genomas_run_main", inherits = TRUE)) {
  genomas_grafica <- genomas_run_main
} else if (exists("genomas_run", inherits = TRUE) &&
           length(genomas_run) > 0) {
  genomas_grafica <- genomas_run
} else {
  genomas_grafica <- genoma_activo
}

# Devuelve la carpeta de resultados de cada genoma.
obtener_carpeta_genoma <- function(gen) {
  if (identical(gen, genoma_activo) &&
      exists("carpeta_salida", inherits = TRUE) &&
      length(carpeta_salida) == 1 &&
      !is.na(carpeta_salida)) {
    return(as.character(carpeta_salida))
  }
  
  if (exists("resumen_corridas", inherits = TRUE) &&
      !is.null(resumen_corridas[[gen]]) &&
      "salida" %in% names(resumen_corridas[[gen]]) &&
      !is.na(resumen_corridas[[gen]][["salida"]])) {
    return(unname(as.character(resumen_corridas[[gen]][["salida"]])))
  }
  
  if (exists("salidas_proteinas", inherits = TRUE) &&
      !is.null(salidas_proteinas[[gen]]) &&
      !is.na(salidas_proteinas[[gen]])) {
    return(as.character(salidas_proteinas[[gen]]))
  }
  
  if (exists("dir_salida", inherits = TRUE) &&
      exists("bin", inherits = TRUE)) {
    return(file.path(
      dir_salida,
      paste0("resultados_", gen, "_", bin)
    ))
  }
  
  NA_character_
}

# Convierte accesiones del FASTA a etiquetas cromosómicas de eq_manual.
normalizar_cromosomas <- function(x, eq) {
  x <- as.character(x)
  
  if (is.data.frame(eq) &&
      all(c("acc", "nombre") %in% names(eq))) {
    nombre <- eq$nombre[match(x, eq$acc)]
    x[!is.na(nombre)] <- nombre[!is.na(nombre)]
  }
  
  x
}

# Extiende un rango para evitar recortar puntos en los PCA.
ampliar_rango <- function(x, fraccion = 0.05) {
  r <- range(x, na.rm = TRUE)
  
  if (any(!is.finite(r))) {
    return(c(-1, 1))
  }
  
  if (diff(r) == 0) {
    return(r + c(-1, 1))
  }
  
  r + c(-1, 1) * diff(r) * fraccion
}

datos_por_genoma <- list()
composiciones_aa_lista <- list()
carpetas_proteinas_lista <- list()
datos_localizacion <- list()
datos_longitudes <- list()

for (gen in genomas_grafica) {
  carpeta_gen <- obtener_carpeta_genoma(gen)
  
  if (is.na(carpeta_gen) || !nzchar(carpeta_gen)) {
    message("No se pudo determinar la carpeta de resultados de ", gen, ".")
    next
  }
  
  protein_outdir <- file.path(
    carpeta_gen,
    paste0("proteinas_", gen)
  )
  
  archivo_rds <- file.path(
    protein_outdir,
    "datos_pca_proteinas.rds"
  )
  
  if (!file.exists(archivo_rds)) {
    message("No existe el RDS de ", gen, ": ", archivo_rds)
    next
  }
  
  datos <- readRDS(archivo_rds)
  datos_por_genoma[[gen]] <- datos
  
  # ----------------------------------------------------------
  # PCA: genera una serie de figuras para cada genoma
  # ----------------------------------------------------------
  df_scores <- datos$scores
  var_exp <- datos$var_exp
  
  if (is.null(df_scores) ||
      !all(c("PC1", "PC2", "familia") %in% names(df_scores)) ||
      nrow(df_scores) < 2) {
    message("El RDS de ", gen, " no tiene scores suficientes para el PCA.")
  } else {
    familias_orden <- sort(unique(as.character(df_scores$familia)))
    n_fam <- length(familias_orden)
    
    paleta_fam <- grDevices::hcl.colors(
      max(1, n_fam),
      palette = "Dark 3"
    )
    names(paleta_fam) <- familias_orden
    
    rango_x <- ampliar_rango(df_scores$PC1)
    rango_y <- ampliar_rango(df_scores$PC2)
    
    pc1_lab <- paste0("PC1 (", round(var_exp[1], 1), "%)")
    pc2_lab <- paste0("PC2 (", round(var_exp[2], 1), "%)")
    
    tema_pca <- ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        panel.border = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
          face = "bold",
          hjust = 0.5
        ),
        axis.title = ggplot2::element_text(face = "bold")
      )
    
    # PCA 1: puntos grises
    p1 <- ggplot2::ggplot(df_scores, ggplot2::aes(PC1, PC2)) +
      ggplot2::geom_point(
        colour = "grey65",
        size = 1.4,
        alpha = 0.75
      ) +
      ggplot2::coord_equal(
        xlim = rango_x,
        ylim = rango_y
      ) +
      ggplot2::labs(
        title = paste("PCA de proteínas —", gen),
        x = pc1_lab,
        y = pc2_lab
      ) +
      tema_pca
    
    ggplot2::ggsave(
      file.path(protein_outdir, "PCA_01_nube_gris.png"),
      p1,
      width = 6,
      height = 5.5,
      dpi = 300,
      bg = "white"
    )
    
    # PCA 2: un facet por familia
    df_bg <- df_scores %>% dplyr::select(-familia)
    n_filas_facet <- ceiling(n_fam / 5)
    
    p2 <- ggplot2::ggplot(df_scores, ggplot2::aes(PC1, PC2)) +
      ggplot2::geom_point(
        data = df_bg,
        colour = "grey85",
        size = 0.8,
        alpha = 0.5
      ) +
      ggplot2::geom_point(
        ggplot2::aes(colour = familia),
        size = 1.1,
        alpha = 0.9
      ) +
      ggplot2::facet_wrap(~familia, ncol = 5) +
      ggplot2::scale_colour_manual(
        values = paleta_fam,
        guide = "none"
      ) +
      ggplot2::coord_equal(
        xlim = rango_x,
        ylim = rango_y
      ) +
      ggplot2::labs(
        title = paste("PCA por familia —", gen),
        x = pc1_lab,
        y = pc2_lab
      ) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(
        panel.border = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(
          face = "bold",
          size = 7
        ),
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
          face = "bold",
          hjust = 0.5
        ),
        axis.text = ggplot2::element_text(size = 6)
      )
    
    ggplot2::ggsave(
      file.path(protein_outdir, "PCA_02_facet_familias.pdf"),
      p2,
      width = 14,
      height = max(2.6, n_filas_facet * 2.6),
      limitsize = FALSE
    )
    
    ggplot2::ggsave(
      file.path(protein_outdir, "PCA_02_facet_familias.png"),
      p2,
      width = 14,
      height = max(2.6, n_filas_facet * 2.6),
      dpi = 300,
      bg = "white",
      limitsize = FALSE
    )
    
    # PCA 3: coloreado por familia
    p3 <- ggplot2::ggplot(
      df_scores,
      ggplot2::aes(PC1, PC2, colour = familia)
    ) +
      ggplot2::geom_point(size = 1.6, alpha = 0.85) +
      ggplot2::scale_colour_manual(
        values = paleta_fam,
        name = "Familia"
      ) +
      ggplot2::coord_equal(
        xlim = rango_x,
        ylim = rango_y
      ) +
      ggplot2::labs(
        title = paste("PCA coloreado por familia —", gen),
        x = pc1_lab,
        y = pc2_lab
      ) +
      tema_pca
    
    ggplot2::ggsave(
      file.path(protein_outdir, "PCA_03_completo_por_familia.png"),
      p3,
      width = 9,
      height = 6.5,
      dpi = 300,
      bg = "white"
    )
    
    ggplot2::ggsave(
      file.path(protein_outdir, "PCA_03_completo_por_familia.pdf"),
      p3,
      width = 9,
      height = 6.5
    )
    
    # PCA 4: elipses para familias con al menos cuatro proteínas
    n_por_familia <- df_scores %>%
      dplyr::count(familia, name = "n")
    
    familias_con_elipse <- n_por_familia %>%
      dplyr::filter(n >= 4) %>%
      dplyr::pull(familia)
    
    df_elipse <- df_scores %>%
      dplyr::filter(familia %in% familias_con_elipse)
    
    p4 <- ggplot2::ggplot(
      df_scores,
      ggplot2::aes(PC1, PC2, colour = familia)
    ) +
      ggplot2::geom_point(size = 1.4, alpha = 0.55)
    
    if (length(familias_con_elipse) > 0) {
      p4 <- p4 +
        ggplot2::stat_ellipse(
          data = df_elipse,
          ggplot2::aes(group = familia),
          type = "t",
          level = 0.95,
          linewidth = 0.5,
          show.legend = FALSE
        )
    }
    
    p4 <- p4 +
      ggplot2::scale_colour_manual(
        values = paleta_fam,
        name = "Familia"
      ) +
      ggplot2::coord_equal(
        xlim = rango_x,
        ylim = rango_y
      ) +
      ggplot2::labs(
        title = paste("PCA con elipses —", gen),
        subtitle = "Elipses solo para familias con al menos cuatro proteínas",
        x = pc1_lab,
        y = pc2_lab
      ) +
      tema_pca
    
    ggplot2::ggsave(
      file.path(protein_outdir, "PCA_04_elipses_confianza.png"),
      p4,
      width = 9.5,
      height = 6.5,
      dpi = 300,
      bg = "white"
    )
    
    ggplot2::ggsave(
      file.path(protein_outdir, "PCA_04_elipses_confianza.pdf"),
      p4,
      width = 9.5,
      height = 6.5
    )
    
    message("PCA guardados para ", gen, " en ", protein_outdir)
  }
  
  # ----------------------------------------------------------
  # Composición de aminoácidos para boxplots por genoma
  # Requiere composicion_aa guardada dentro del RDS del módulo 17.
  # ----------------------------------------------------------
  if (!is.null(datos$composicion_aa) &&
      nrow(datos$composicion_aa) > 0) {
    composiciones_aa_lista[[gen]] <- datos$composicion_aa
    carpetas_proteinas_lista[[gen]] <- protein_outdir
  } else {
    message(
      "El RDS de ", gen,
      " no contiene composicion_aa. Vuelve a correr el módulo 17 para ese genoma."
    )
  }
  
  # ----------------------------------------------------------
  # Coordenadas cromosómicas de proteínas
  # ----------------------------------------------------------
  if (is.null(genomas[[gen]][["fna"]]) ||
      !file.exists(genomas[[gen]][["fna"]])) {
    message("No existe el FASTA genómico de ", gen, "; se omite su mapa.")
    next
  }
  
  if (is.null(datos$localizacion) ||
      nrow(datos$localizacion) == 0) {
    message("El RDS de ", gen, " no contiene localizaciones; se omite su mapa.")
    next
  }
  
  eq_gen <- genomas[[gen]][["eq_manual"]]
  
  loc <- datos$localizacion
  loc$genoma <- gen
  loc$cromosoma <- normalizar_cromosomas(loc$cromosoma, eq_gen)
  loc$posicion_mb <- ((loc$start + loc$end) / 2) / 1e6
  
  datos_localizacion[[gen]] <- loc
  
  # Longitudes de las secuencias del FASTA genómico.
  secuencias_genoma <- Biostrings::readDNAStringSet(
    genomas[[gen]][["fna"]]
  )
  
  nombres_todos <- sub("\\s.*$", "", names(secuencias_genoma))
  longitudes_todas <- as.numeric(
    Biostrings::width(secuencias_genoma)
  )
  
  # Si eq_manual está definido, incluye solo sus accesiones.
  # En RH88 esto excluye apicoplasto y scaffolds.
  if (is.data.frame(eq_gen) && "acc" %in% names(eq_gen)) {
    conservar <- nombres_todos %in% eq_gen$acc
  } else {
    conservar <- rep(TRUE, length(nombres_todos))
    
    if (exists("min_len_cromosoma", inherits = TRUE) &&
        !is.null(min_len_cromosoma)) {
      conservar <- conservar &
        longitudes_todas >= min_len_cromosoma
    }
    
    if (exists("patron_incluir", inherits = TRUE) &&
        !is.null(patron_incluir)) {
      conservar <- conservar &
        grepl(patron_incluir, nombres_todos)
    }
    
    if (exists("patron_excluir", inherits = TRUE) &&
        !is.null(patron_excluir)) {
      conservar <- conservar &
        !grepl(patron_excluir, nombres_todos)
    }
  }
  
  nombres_fasta <- nombres_todos[conservar]
  longitud_pb <- longitudes_todas[conservar]
  
  long <- data.frame(
    genoma = rep(gen, length(nombres_fasta)),
    cromosoma = normalizar_cromosomas(nombres_fasta, eq_gen),
    longitud_pb = longitud_pb,
    stringsAsFactors = FALSE
  )
  
  datos_longitudes[[gen]] <- long
}

# ============================================================
# Boxplots de composición: un archivo por genoma, un panel por familia
# ============================================================
# Clasificación fisicoquímica corregida y estandarizada
clase_aa <- c(
  # Alifáticos / Hidrofóbicos pequeños
  A = "Alifáticos", G = "Alifáticos", 
  
  # Alifáticos / Hidrofóbicos de cadena lineal o ramificada
  V = "Hidrofóbicos", I = "Hidrofóbicos", L = "Hidrofóbicos", 
  M = "Hidrofóbicos", P = "Hidrofóbicos",
  
  # Aromáticos (también son hidrofóbicos, pero prima su estructura)
  F = "Aromáticos", W = "Aromáticos", Y = "Aromáticos",
  
  # Polares sin carga (se incluye Serina aquí por reactividad)
  S = "Polares", N = "Polares", Q = "Polares", T = "Polares", C = "Polares",
  
  # Cargados Positivos / Básicos
  K = "Básicos (+)", R = "Básicos (+)", H = "Básicos (+)",
  
  # Cargados Negativos / Ácidos
  D = "Ácidos (-)", E = "Ácidos (-)"
)

# Paleta de colores consistente con la nueva nomenclatura
colores_clase_aa <- c(
  "Alifáticos"   = "#56B4E9", # Tu color para "Pequeños"
  "Hidrofóbicos" = "#0072B2",
  "Aromáticos"   = "#E69F00",
  "Polares"      = "#009E73",
  "Básicos (+)"  = "#D55E00", # Tu color para "Cargados (+)"
  "Ácidos (-)"   = "#CC79A7"  # Tu color para "Cargados (-)"
)

ajuste_pvalores <- "BH"  # Cambia a "bonferroni" si prefieres Bonferroni.
for (gen in names(composiciones_aa_lista)) {
  comp_gen <- composiciones_aa_lista[[gen]] %>%
    dplyr::filter(familia %in% familias_localizacion)
  
  if (nrow(comp_gen) == 0) {
    message("Sin composición para las familias seleccionadas en ", gen)
    next
  }
  
  aa_presentes <- intersect(
    names(clase_aa),
    names(comp_gen)
  )
  
  if (length(aa_presentes) == 0) {
    message("No hay columnas de aminoácidos en el RDS de ", gen)
    next
  }
  
  comp_larga <- comp_gen %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(aa_presentes),
      names_to = "aminoacido",
      values_to = "porcentaje"
    ) %>%
    dplyr::mutate(
      clase = unname(clase_aa[aminoacido]),
      aminoacido = factor(
        aminoacido,
        levels = names(clase_aa)
      )
    ) %>%
    dplyr::filter(
      !is.na(porcentaje),
      !is.na(clase)
    )
  
  # ---- Kruskal-Wallis por aminoácido y tamaño del efecto epsilon² ----
  
  kw_lista <- lapply(
    split(comp_larga, comp_larga$aminoacido),
    function(d) {
      d <- d[is.finite(d$porcentaje) & !is.na(d$familia), , drop = FALSE]
      d$familia <- droplevels(factor(d$familia))
      
      n <- nrow(d)
      k <- nlevels(d$familia)
      aminoacido_actual <- as.character(d$aminoacido[1])
      
      if (k < 2 || n <= k) {
        return(data.frame(
          aminoacido = aminoacido_actual,
          n = n,
          n_familias = k,
          estadistico_H = NA_real_,
          gl = NA_real_,
          p = NA_real_,
          epsilon2 = NA_real_
        ))
      }
      
      prueba <- kruskal.test(porcentaje ~ familia, data = d)
      H <- unname(prueba$statistic)
      
      # Epsilon cuadrado para Kruskal-Wallis.
      epsilon2 <- max(0, (H - k + 1) / (n - k))
      
      data.frame(
        aminoacido = aminoacido_actual,
        n = n,
        n_familias = k,
        estadistico_H = H,
        gl = unname(prueba$parameter),
        p = prueba$p.value,
        epsilon2 = epsilon2
      )
    }
  )
  
  tabla_kw <- dplyr::bind_rows(kw_lista) %>%
    dplyr::mutate(
      p_ajustado = p.adjust(p, method = ajuste_pvalores),
      significativo = !is.na(p_ajustado) & p_ajustado < 0.05
    ) %>%
    dplyr::arrange(p_ajustado)
  
  write.csv(
    tabla_kw,
    file.path(
      carpetas_proteinas_lista[[gen]],
      paste0("Kruskal_Wallis_aminoacidos_", gen, ".csv")
    ),
    row.names = FALSE
  )
  
  # ---- Dunn por aminoácido: comparaciones entre pares de familias ----
  # Se ajustan conjuntamente todas las comparaciones de todos los aminoácidos.
  
  dunn_un_aminoacido <- function(d) {
    d <- d[is.finite(d$porcentaje) & !is.na(d$familia), , drop = FALSE]
    d$familia <- droplevels(factor(d$familia))
    
    k <- nlevels(d$familia)
    n <- nrow(d)
    aminoacido_actual <- as.character(d$aminoacido[1])
    
    if (k < 2 || n <= k || n < 2) {
      return(data.frame())
    }
    
    rangos <- rank(d$porcentaje, ties.method = "average")
    medias_rango <- tapply(rangos, d$familia, mean)
    n_por_familia <- table(d$familia)
    
    tamanos_empate <- as.numeric(table(d$porcentaje))
    correccion_empates <- sum(tamanos_empate^3 - tamanos_empate)
    
    varianza_rangos <- n * (n + 1) / 12 -
      correccion_empates / (12 * (n - 1))
    
    if (!is.finite(varianza_rangos) || varianza_rangos <= 0) {
      return(data.frame())
    }
    
    pares <- utils::combn(names(medias_rango), 2, simplify = FALSE)
    
    dunn_lista <- lapply(pares, function(par) {
      f1 <- par[1]
      f2 <- par[2]
      
      error_estandar <- sqrt(
        varianza_rangos *
          (1 / n_por_familia[[f1]] + 1 / n_por_familia[[f2]])
      )
      
      z <- (medias_rango[[f1]] - medias_rango[[f2]]) / error_estandar
      p <- 2 * stats::pnorm(-abs(z))
      
      mediana1 <- median(d$porcentaje[d$familia == f1], na.rm = TRUE)
      mediana2 <- median(d$porcentaje[d$familia == f2], na.rm = TRUE)
      
      data.frame(
        aminoacido = aminoacido_actual,
        familia1 = f1,
        familia2 = f2,
        n1 = as.integer(n_por_familia[[f1]]),
        n2 = as.integer(n_por_familia[[f2]]),
        mediana1 = mediana1,
        mediana2 = mediana2,
        diferencia_medianas = mediana1 - mediana2,
        z = z,
        p = p,
        stringsAsFactors = FALSE
      )
    })
    
    dplyr::bind_rows(dunn_lista)
  }
  
  dunn_lista <- lapply(
    split(comp_larga, comp_larga$aminoacido),
    dunn_un_aminoacido
  )
  
  tabla_dunn <- dplyr::bind_rows(dunn_lista)
  
  if (nrow(tabla_dunn) > 0) {
    tabla_dunn <- tabla_dunn %>%
      dplyr::mutate(
        p_ajustado = p.adjust(p, method = ajuste_pvalores),
        significativo = p_ajustado < 0.05
      ) %>%
      dplyr::arrange(p_ajustado)
  }
  
  write.csv(
    tabla_dunn,
    file.path(
      carpetas_proteinas_lista[[gen]],
      paste0("Dunn_aminoacidos_", gen, ".csv")
    ),
    row.names = FALSE
  )
  
  p_aa_genoma <- ggplot2::ggplot(
    comp_larga,
    ggplot2::aes(
      x = aminoacido,
      y = porcentaje,
      fill = clase
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.75,
      outlier.size = 0.35,
      outlier.alpha = 0.35
    ) +
    ggplot2::facet_wrap(
      ~familia,
      ncol = 2
    ) +
    ggplot2::scale_fill_manual(
      values = colores_clase_aa,
      name = "Clase"
    ) +
    ggplot2::labs(
      title = paste("Composición de aminoácidos —", gen),
      subtitle = paste(
        "Familias:",
        paste(familias_localizacion, collapse = ", ")
      ),
      x = "Aminoácido",
      y = "Composición por proteína (%)"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "bottom"
    )
  
  salida_gen <- carpetas_proteinas_lista[[gen]]
  dir.create(salida_gen, recursive = TRUE, showWarnings = FALSE)
  
  ggplot2::ggsave(
    file.path(
      salida_gen,
      paste0("boxplots_composicion_aminoacidos_", gen, ".png")
    ),
    p_aa_genoma,
    width = 12,
    height = max(6, 3 * length(unique(comp_larga$familia))),
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
  
  ggplot2::ggsave(
    file.path(
      salida_gen,
      paste0("boxplots_composicion_aminoacidos_", gen, ".pdf")
    ),
    p_aa_genoma,
    width = 12,
    height = max(6, 3 * length(unique(comp_larga$familia))),
    bg = "white",
    limitsize = FALSE
  )
  
  # ---- Heatmap de medianas: familias × aminoácidos ----
  
  medianas_aa <- comp_larga %>%
    dplyr::group_by(familia, aminoacido) %>%
    dplyr::summarise(
      mediana = median(porcentaje, na.rm = TRUE),
      .groups = "drop"
    )
  
  medianas_wide <- tidyr::pivot_wider(
    medianas_aa,
    names_from = aminoacido,
    values_from = mediana
  )
  
  mat_medianas <- as.matrix(
    medianas_wide[, setdiff(names(medianas_wide), "familia"), drop = FALSE]
  )
  rownames(mat_medianas) <- medianas_wide$familia
  
  # Z-score por aminoácido (columna).
  mat_z <- mat_medianas
  
  for (j in seq_len(ncol(mat_medianas))) {
    valores <- mat_medianas[, j]
    sd_col <- stats::sd(valores, na.rm = TRUE)
    
    if (length(valores) < 2 || !is.finite(sd_col) || sd_col == 0) {
      mat_z[, j] <- 0
    } else {
      mat_z[, j] <- as.numeric(scale(valores))
    }
  }
  
  # Ordena las familias por similitud de perfil cuando hay al menos dos.
  if (nrow(mat_z) >= 2) {
    orden_familias <- rownames(mat_z)[
      stats::hclust(stats::dist(mat_z))$order
    ]
  } else {
    orden_familias <- rownames(mat_z)
  }
  
  heatmap_df <- as.data.frame(as.table(mat_z))
  names(heatmap_df) <- c("familia", "aminoacido", "z_score")
  
  heatmap_df$familia <- factor(
    heatmap_df$familia,
    levels = orden_familias
  )
  
  heatmap_df$aminoacido <- factor(
    heatmap_df$aminoacido,
    levels = names(clase_aa)
  )
  
  p_heatmap_aa <- ggplot2::ggplot(
    heatmap_df,
    ggplot2::aes(
      x = aminoacido,
      y = familia,
      fill = z_score
    )
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "Z-score"
    ) +
    ggplot2::labs(
      title = paste("Mediana de composición de aminoácidos —", gen),
      subtitle = "Z-score por aminoácido; familias agrupadas por similitud",
      x = "Aminoácido",
      y = "Familia"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )
  
  ggplot2::ggsave(
    file.path(
      carpetas_proteinas_lista[[gen]],
      paste0("heatmap_mediana_aminoacidos_", gen, ".png")
    ),
    p_heatmap_aa,
    width = 10,
    height = max(5, 0.45 * nrow(mat_z) + 2),
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
  
  ggplot2::ggsave(
    file.path(
      carpetas_proteinas_lista[[gen]],
      paste0("heatmap_mediana_aminoacidos_", gen, ".pdf")
    ),
    p_heatmap_aa,
    width = 10,
    height = max(5, 0.45 * nrow(mat_z) + 2),
    bg = "white",
    limitsize = FALSE
  )
  
  print(p_aa_genoma)
  message("Boxplots de composición guardados para ", gen, " en ", salida_gen)
}

# ============================================================
# Mapas de ubicación cromosómica
# Una figura por cromosoma y otra con todos los cromosomas.
# ============================================================

if (!length(datos_longitudes)) {
  message("No hay longitudes cromosómicas para dibujar mapas.")
} else {
  longitudes_mapa <- dplyr::bind_rows(datos_longitudes)
  genomas_mapa_validos <- unique(longitudes_mapa$genoma)
  cromosomas_mapa <- unique(longitudes_mapa$cromosoma)
  
  if (length(datos_localizacion)) {
    loc_proteinas <- dplyr::bind_rows(datos_localizacion) %>%
      dplyr::filter(familia %in% familias_localizacion)
  } else {
    loc_proteinas <- data.frame()
  }
  
  colores_familias <- grDevices::hcl.colors(
    max(1, length(familias_localizacion)),
    palette = "Dark 3"
  )
  names(colores_familias) <- familias_localizacion
  
  dir_mapa_proteinas <- file.path(
    dir_salida,
    "comparacion_proteinas"
  )
  dir.create(
    dir_mapa_proteinas,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dibujar_cromosoma_proteinas <- function(chr) {
    pistas <- data.frame(
      genoma = genomas_mapa_validos,
      y = rev(seq_along(genomas_mapa_validos)) * 2,
      stringsAsFactors = FALSE
    )
    
    long_chr <- longitudes_mapa[
      longitudes_mapa$cromosoma == chr,
      c("genoma", "longitud_pb"),
      drop = FALSE
    ]
    
    pistas <- merge(
      pistas,
      long_chr,
      by = "genoma",
      all.x = TRUE
    )
    
    validas <- pistas[
      !is.na(pistas$longitud_pb),
      ,
      drop = FALSE
    ]
    
    if (nrow(validas) == 0) return(NULL)
    
    puntos <- loc_proteinas[
      loc_proteinas$cromosoma == chr,
      ,
      drop = FALSE
    ]
    
    if (nrow(puntos) > 0) {
      puntos <- merge(
        puntos,
        pistas[, c("genoma", "y")],
        by = "genoma",
        all = FALSE
      )
    }
    
    xmax <- max(validas$longitud_pb, na.rm = TRUE) / 1e6
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_segment(
        data = validas,
        ggplot2::aes(
          x = 0,
          xend = longitud_pb / 1e6,
          y = y,
          yend = y
        ),
        linewidth = 7,
        colour = "#A7B0BC",
        lineend = "round"
      ) +
      ggplot2::geom_segment(
        data = validas,
        ggplot2::aes(
          x = 0,
          xend = longitud_pb / 1e6,
          y = y,
          yend = y
        ),
        linewidth = 6,
        colour = "#E9EDF2",
        lineend = "round"
      )
    
    if (nrow(puntos) > 0) {
      p <- p +
        ggplot2::geom_point(
          data = puntos,
          ggplot2::aes(
            x = posicion_mb,
            y = y,
            fill = familia
          ),
          shape = 21,
          size = 2.8,
          colour = "white",
          stroke = 0.45
        )
    }
    
    p +
      ggplot2::scale_fill_manual(
        values = colores_familias,
        name = "Familia"
      ) +
      ggplot2::scale_x_continuous(
        name = "Posición cromosómica (Mb)",
        expand = ggplot2::expansion(mult = c(0, 0.03))
      ) +
      ggplot2::scale_y_continuous(
        breaks = pistas$y,
        labels = pistas$genoma,
        limits = c(
          min(pistas$y) - 0.8,
          max(pistas$y) + 0.8
        ),
        expand = c(0, 0)
      ) +
      ggplot2::labs(
        title = paste("Ubicación de proteínas — cromosoma", chr),
        subtitle = paste(
          "Familias:",
          paste(familias_localizacion, collapse = ", ")
        ),
        y = NULL
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        panel.border = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(
          fill = "white",
          colour = NA
        ),
        plot.background = ggplot2::element_rect(
          fill = "white",
          colour = NA
        ),
        legend.position = "top",
        axis.line.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(face = "bold"),
        panel.grid.major.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
          face = "bold",
          hjust = 0.5
        )
      ) +
      ggplot2::coord_cartesian(
        xlim = c(0, xmax * 1.03),
        clip = "off"
      )
  }
  
  figuras_proteinas <- list()
  
  for (chr in cromosomas_mapa) {
    p <- dibujar_cromosoma_proteinas(chr)
    if (!is.null(p)) {
      figuras_proteinas[[chr]] <- p
    }
  }
  
  if (length(figuras_proteinas) > 0) {
    alto_pdf <- 2.5 + 0.9 * length(genomas_mapa_validos)
    
    grDevices::pdf(
      file.path(
        dir_mapa_proteinas,
        "mapa_familias_proteinas_por_cromosoma.pdf"
      ),
      width = 12,
      height = alto_pdf,
      onefile = TRUE
    )
    
    for (p in figuras_proteinas) {
      print(p)
    }
    
    grDevices::dev.off()
    
    for (chr in names(figuras_proteinas)) {
      nombre_seguro <- gsub("[^A-Za-z0-9_.-]", "_", chr)
      
      ggplot2::ggsave(
        file.path(
          dir_mapa_proteinas,
          paste0("cromosoma_", nombre_seguro, ".png")
        ),
        figuras_proteinas[[chr]],
        width = 12,
        height = alto_pdf,
        dpi = 300,
        bg = "white"
      )
    }
  }
  
  # ---- Figura única con todos los cromosomas ----
  pistas_todas <- longitudes_mapa %>%
    dplyr::mutate(
      etiqueta = paste(cromosoma, genoma, sep = " — ")
    )
  
  pistas_todas <- pistas_todas[
    order(
      match(pistas_todas$cromosoma, cromosomas_mapa),
      match(pistas_todas$genoma, genomas_mapa_validos)
    ),
    ,
    drop = FALSE
  ]
  
  pistas_todas$y <- rev(seq_len(nrow(pistas_todas)))
  
  if (nrow(loc_proteinas) > 0) {
    puntos_todos <- merge(
      loc_proteinas,
      pistas_todas[
        ,
        c("genoma", "cromosoma", "longitud_pb", "y"),
        drop = FALSE
      ],
      by = c("genoma", "cromosoma"),
      all = FALSE
    )
    
    puntos_todos$posicion_mb <-
      (puntos_todos$start + puntos_todos$end) / 2e6
    
    puntos_todos <- puntos_todos[
      is.finite(puntos_todos$posicion_mb) &
        puntos_todos$posicion_mb >= 0 &
        puntos_todos$posicion_mb <= puntos_todos$longitud_pb / 1e6,
      ,
      drop = FALSE
    ]
  } else {
    puntos_todos <- data.frame()
  }
  
  xmax_todos <- max(pistas_todas$longitud_pb, na.rm = TRUE) / 1e6
  
  grafica_todos_cromosomas <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = pistas_todas,
      ggplot2::aes(
        x = 0,
        xend = longitud_pb / 1e6,
        y = y,
        yend = y
      ),
      linewidth = 3.5,
      colour = "#C9CED4",
      lineend = "round"
    )
  
  if (nrow(puntos_todos) > 0) {
    grafica_todos_cromosomas <- grafica_todos_cromosomas +
      ggplot2::geom_point(
        data = puntos_todos,
        ggplot2::aes(
          x = posicion_mb,
          y = y,
          fill = familia
        ),
        shape = 21,
        size = 2.6,
        colour = "white",
        stroke = 0.35,
        alpha = 0.95
      )
  }
  
  grafica_todos_cromosomas <- grafica_todos_cromosomas +
    ggplot2::scale_fill_manual(
      values = colores_familias,
      name = "Familia"
    ) +
    ggplot2::scale_x_continuous(
      name = "Posición cromosómica (Mb)",
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::scale_y_continuous(
      breaks = pistas_todas$y,
      labels = pistas_todas$etiqueta,
      limits = range(pistas_todas$y) + c(-0.7, 0.7),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      title = "Familias de proteínas en todos los cromosomas",
      subtitle = paste(
        "Familias seleccionadas:",
        paste(familias_localizacion, collapse = ", ")
      ),
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 7, face = "bold"),
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      legend.position = "top"
    ) +
    ggplot2::coord_cartesian(clip = "off")
  
  ggplot2::ggsave(
    file.path(
      dir_mapa_proteinas,
      "mapa_todos_los_cromosomas.png"
    ),
    grafica_todos_cromosomas,
    width = 14,
    height = max(8, 0.28 * nrow(pistas_todas) + 2),
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
  
  ggplot2::ggsave(
    file.path(
      dir_mapa_proteinas,
      "mapa_todos_los_cromosomas.pdf"
    ),
    grafica_todos_cromosomas,
    width = 14,
    height = max(8, 0.28 * nrow(pistas_todas) + 2),
    bg = "white",
    limitsize = FALSE
  )
  
  write.csv(
    loc_proteinas,
    file.path(
      dir_mapa_proteinas,
      "proteinas_familias_seleccionadas.csv"
    ),
    row.names = FALSE
  )
  
  print(grafica_todos_cromosomas)
  
  cat(
    "Mapas cromosómicos guardados en: ",
    dir_mapa_proteinas,
    "\n",
    sep = ""
  )
}