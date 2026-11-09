# ============================================================
# 15_figuras.R  -  UNICO script de graficas (genoma actual y comparacion entre genomas)
# ------------------------------------------------------------
# MODOS      : fig_modo = "genoma"      (lo fija el main tras cada genoma completo): figuras F01-F20
#              fig_modo = "comparacion" (lo fija el main tras 16_comparacion.R)    : figuras G1-G4
# REGLA CLAVE: las "zonas criticas" salen SOLO de zonas_df (script 13). Cada figura se dibuja UNA vez.
# REQUIERE   : modo genoma: casi todos los modulos 01-13 (ver cabeceras), 06b (gb_mat, sat350_regiones),
#              07b (microsat_loci). Modo comparacion: 16 (cmp_*, niveles_chr, col_gen, carpeta_cmp)
# FIGURAS    : F01 panorama | F02 perfiles con zonas | F03 zoom a zonas | F04 zonas vs resto | F05 GC vs genes/Mb
#              F06 extremos vs interior | F07 clases de repeticion | F08 repeticiones ST vs interior
#              F09 alta vs baja repeticion | F10 intergenico por orientacion | F11 satelites | F12 GenBank sat350
#              F13 idiograma | F14 tabla-mapa de extremos | F15 validacion ST vs articulo | F16 genes en extremos
#              F17 mapa de subtelomeros | F18 satelites + microsatelites | F19 microsatelites y contexto
#              F20 periodo TRF en regiones SAT350 | G1 mapas por cromosoma y genoma | G2 mapa de calor SAT350
#              G3 resumen por genoma | G4 subtelomeros y telomeros por extremo
# SALIDA     : PDF + PNG (300 dpi) y leyendas_figuras.md (carpeta_salida) / comparacion_genomas/ (comparacion)
# ============================================================
# Modo: "comparacion" solo si hay datos de comparacion (16_comparacion.R los crea); si no, se dibuja el genoma actual
if (!exists("fig_modo") || (fig_modo == "comparacion" && !isTRUE(get0("cmp_ok", ifnotfound = FALSE)))) fig_modo <- "genoma"

# ---------- 0. Utilidades ----------
# Cada figura corre en su propio entorno: las variables temporales no se pisan entre figuras.
intentar <- function(nombre, expr, cond = TRUE) {
  ok <- tryCatch(isTRUE(cond), error = function(e) FALSE)
  if (!ok) { message("Figura omitida (sin datos): ", nombre); return(invisible(NULL)) }
  entorno <- new.env(parent = parent.frame())
  tryCatch(eval(substitute(expr), entorno), error = function(e) {
    message("Figura '", nombre, "' falló: ", conditionMessage(e)); invisible(NULL)
  })
}

tema_publicacion <- theme_classic(base_size = 10, base_family = font_family) +
  theme(axis.line = element_line(color = "grey20", linewidth = 0.4),
        axis.ticks = element_line(color = "grey40", linewidth = 0.3),
        axis.text = element_text(color = "black", size = 8),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 11, margin = margin(b = 3)),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(), legend.position = "none")

tema_art <- tema_publicacion +
  theme(plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey25"),
        plot.caption  = element_text(size = 7, hjust = 0, color = "grey30"))

# PDF + PNG (+ TIFF opcional)
guardar <- function(p, nombre, ancho, alto, tiff = FALSE) {
  base <- paste0(carpeta_salida, nombre, "_", ventana_label)
  ggsave(paste0(base, ".pdf"), p, width = ancho, height = alto, device = cairo_pdf)
  ggsave(paste0(base, ".png"), p, width = ancho, height = alto, dpi = 300, bg = "white")
  if (tiff) ggsave(paste0(base, ".tiff"), p, width = ancho, height = alto, dpi = 600, compression = "lzw")
}

pal <- list(telomero = "#1D2D44", st_valido = "#3A86FF", st_no = "#8D99AE",
            gen = "#D90429", otras = "grey65", zona = "#D55E00")
col_sat <- c(sat240 = "#8338EC", sat350 = "#00F5D4")

# Leyendas de figura
leyendas_file <- if (fig_modo == "genoma") paste0(carpeta_salida, "leyendas_figuras.md") else paste0(carpeta_cmp, "leyendas_comparacion.md")
cat("# Leyendas de figuras (", if (fig_modo == "genoma") paste0(genoma_id, ", ventanas de ", ventana_label) else paste(cmp_nombres, collapse = ", "),
    ")\n\n", sep = "", file = leyendas_file)
leyenda <- function(titulo, texto)
  cat("## ", titulo, "\n\n", texto, "\n\n", sep = "", file = leyendas_file, append = TRUE)

# Flechas de gen (polígonos con punta según la hebra)
poli_gen <- function(x0, x1, str, y = 0.5, h = 0.5, hl_max = 1) {
  if (length(x0) == 0) return(NULL)
  do.call(rbind, lapply(seq_along(x0), function(k) {
    hl <- min(hl_max, (x1[k] - x0[k]) * 0.6)
    xs <- c(x0[k], x1[k] - hl, x1[k] - hl, x1[k], x1[k] - hl, x1[k] - hl, x0[k])
    ys <- y + h * c(-0.25, -0.25, -0.5, 0, 0.5, 0.25, 0.25)
    if (str[k] == "-") xs <- x0[k] + x1[k] - xs
    data.frame(id = k, x = xs, y = ys)
  }))
}

# Estadística por faceta (2 grupos: Wilcoxon + delta de Cliff; >2: Kruskal-Wallis)
stats_grupos <- function(s) {
  s <- s[!is.na(s$valor), ]; s$grupo <- droplevels(s$grupo)
  n <- table(s$grupo); k <- length(n)
  if (k < 2 || any(n < 2)) return(NA_character_)
  if (k == 2) {
    x1 <- s$valor[s$grupo == levels(s$grupo)[1]]; x2 <- s$valor[s$grupo == levels(s$grupo)[2]]
    w <- suppressWarnings(wilcox.test(x1, x2))
    delta <- 2 * unname(w$statistic) / (length(x1) * length(x2)) - 1
    sprintf("Wilcoxon p = %s | delta de Cliff = %.2f\nn = %s", format.pval(w$p.value, digits = 2, eps = 1e-4),
            delta, paste(n, collapse = " / "))
  } else {
    kw <- kruskal.test(valor ~ grupo, data = s)
    sprintf("Kruskal-Wallis p = %s\nn = %s", format.pval(kw$p.value, digits = 2, eps = 1e-4), paste(n, collapse = " / "))
  }
}

# Boxplots por faceta con p, tamaño de efecto y n
box_anot <- function(dl, titulo, subtitulo, fills, ncol = 4) {
  dl$grupo <- factor(dl$grupo); dl$variable <- factor(dl$variable, levels = unique(dl$variable))
  st <- do.call(rbind, lapply(split(dl, as.character(dl$variable)), function(s)
    data.frame(variable = as.character(s$variable[1]), label = stats_grupos(s), stringsAsFactors = FALSE)))
  st$variable <- factor(st$variable, levels = levels(dl$variable)); st$x <- (nlevels(dl$grupo) + 1) / 2
  ggplot(dl, aes(grupo, valor, fill = grupo)) +
    geom_boxplot(outlier.size = 0.4, linewidth = 0.4, width = 0.6) +
    geom_text(data = st, aes(x = x, y = Inf, label = label), inherit.aes = FALSE,
              vjust = 1.15, size = 2.2, lineheight = 0.9) +
    facet_wrap(~ variable, scales = "free_y", ncol = ncol) +
    scale_fill_manual(values = fills) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.28))) +
    labs(x = NULL, y = NULL, title = titulo, subtitle = subtitulo) +
    tema_art + theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 7))
}

# Agrupa microsatelites situados a < gap_mb en una sola etiqueta (columnas: genoma [opcional], chrom, pos_mb, nombre, tipo)
agrupar_microsat <- function(m, gap_mb = 0.001) {
  if (is.null(m) || nrow(m) == 0) return(NULL)
  if (!"genoma" %in% names(m)) m$genoma <- "g"
  m <- m[order(m$genoma, m$chrom, m$pos_mb), ]
  m$grupo <- ave(seq_len(nrow(m)), m$genoma, m$chrom, FUN = function(i) cumsum(c(TRUE, diff(m$pos_mb[i]) > gap_mb)))
  key <- paste(m$genoma, m$chrom, m$grupo, sep = "|")
  out <- do.call(rbind, lapply(split(m, key), function(x)
    data.frame(genoma = x$genoma[1], chrom = x$chrom[1], pos_mb = mean(x$pos_mb),
               etiqueta = paste(unique(x$nombre), collapse = " / "), tipo = x$tipo[1], stringsAsFactors = FALSE)))
  rownames(out) <- NULL
  out
}

if (fig_modo == "genoma") {
  # ---------- 1. Datos de dibujo compartidos ----------
  # Ventanas: NA fuera del conjunto de análisis (la línea se corta, no cruza ventanas excluidas)
  cols_p <- c("GC_percent", "genes_per_Mb", "repeat_pct", "intron_pct", "intergenic_pct")
  tp <- tabla; tp[!tp$analisis, cols_p] <- NA
  #tp <- tabla
  # Zonas críticas: UNA sola definición (zonas_df) para todas las figuras
  hay_zonas <- exists("zonas_df") && is.data.frame(zonas_df) && nrow(zonas_df) > 0
  if (hay_zonas) {
    stopifnot(all(c("chrom", "start_Mb", "end_Mb", "tipo", "zona_id") %in% names(zonas_df)))
    zonas_df$chrom <- as.character(zonas_df$chrom)
    zonas_df$tipo  <- as.character(zonas_df$tipo)
    zona_de <- function(df) {
      out <- rep("Fuera de zonas", nrow(df))
      for (i in seq_len(nrow(zonas_df))) {
        z <- zonas_df[i, ]
        h <- as.character(df$chrom) == z$chrom & df$end_Mb > z$start_Mb & df$start_Mb < z$end_Mb
        out[h] <- z$tipo
      }
      out
    }
    tp$zona_tipo <- factor(zona_de(tp), levels = c(unique(c(niveles_zona, zonas_df$tipo)), "Fuera de zonas"))
  }
  
  # Genes (coordenadas, hebra y etiqueta corta; "hp" si es hipotético)
  mc_g <- mcols(genes_gr)
  col_o <- function(n) if (n %in% names(mc_g)) as.character(mc_g[[n]]) else rep(NA_character_, length(genes_gr))
  g_lab <- col_o("gene_nombre"); g_des <- col_o("descripcion"); g_id <- col_o("gene_id")
  sin_nom <- is.na(g_lab) | !nzchar(g_lab) | g_lab == g_id
  g_lab[sin_nom] <- ifelse(grepl("hypothetical|putative|unspecified|uncharacterized", g_des[sin_nom], ignore.case = TRUE),
                           "hp", g_des[sin_nom])
  g_lab <- ifelse(!is.na(g_lab) & nchar(g_lab) > 14, paste0(substr(g_lab, 1, 12), "..."), g_lab)
  g_str <- as.character(strand(genes_gr)); g_str[g_str == "*"] <- "+"
  g_all <- data.frame(acc = as.character(seqnames(genes_gr)), s = start(genes_gr), e = end(genes_gr),
                      str = g_str, lab = g_lab, stringsAsFactors = FALSE)
  genes_rug <- data.frame(chrom = factor(nombre_chr(as.character(seqnames(genes_mid))), levels = eq$nombre),
                          pos = start(genes_mid) / 1e6)
  st_plot <- if (length(st_gr) > 0) data.frame(chrom = factor(nombre_chr(as.character(seqnames(st_gr))), levels = eq$nombre),
                                               xmin = start(st_gr) / 1e6, xmax = end(st_gr) / 1e6) else NULL
  
  # RepeatMasker con clase simplificada (satélites con su nombre)
  tiene_rm <- exists("rm_df") && is.data.frame(rm_df) && nrow(rm_df) > 0 && ncol(rm_df) >= 11
  if (tiene_rm) {
    rm_g <- rm_df[as.character(rm_df$V5) %in% names(sl), ]
    rm_g$V5 <- as.character(rm_g$V5)
    rm_g$clase <- sub("/.*", "", as.character(rm_g$V11))
    rm_g$clase[is.na(rm_g$clase) | !nzchar(rm_g$clase)] <- "Otras"
    for (sn in names(patron_sat)) rm_g$clase[grepl(patron_sat[[sn]], rm_g$V10, ignore.case = TRUE)] <- sn
    bp_cl <- sort(tapply(rm_g$V7 - rm_g$V6 + 1, rm_g$clase, sum), decreasing = TRUE)
    top_cl_art <- head(setdiff(names(bp_cl), c(names(patron_sat), "Otras")), 8)
    rm_g$clase_f <- ifelse(rm_g$clase %in% c(top_cl_art, names(patron_sat)), rm_g$clase, "Otras")
    niv_cl <- c(top_cl_art, names(patron_sat), "Otras")
    col_clase <- setNames(grDevices::hcl.colors(max(3, length(niv_cl)), "Set 2")[seq_along(niv_cl)], niv_cl)
    usar_s <- intersect(names(patron_sat), names(col_sat)); col_clase[usar_s] <- col_sat[usar_s]
    col_clase["Otras"] <- "grey65"
  }
  
  # ============================================================
  # F01. Panorama genómico (repetición por ventana + intrones + promotores + telómeros + zonas críticas)
  # ============================================================
  intentar("F01. Panorama genómico", {
    intrones_df <- data.frame(seqnames = factor(as.character(seqnames(intrones_red)), levels = names(sl)),
                              start = start(intrones_red) / 1e6, end = end(intrones_red) / 1e6)
    prom_df     <- data.frame(seqnames = factor(as.character(seqnames(prom_red)), levels = names(sl)),
                              start = start(prom_red) / 1e6, end = end(prom_red) / 1e6)
    df_bins     <- data.frame(seqnames = factor(tabla$chrom_acc, levels = names(sl)),
                              start = tabla$start, end = tabla$end, pct = tabla$repeat_pct)
    tel_plot <- tel_df[tel_df$es_telomero, ]
    tel_plot <- if (nrow(tel_plot) > 0) data.frame(seqnames = factor(tel_plot$seqnames, levels = names(sl)),
                                                   pos = tel_plot$start / 1e6) else NULL
    zon_pan <- NULL
    if (hay_zonas) {
      acc_de <- setNames(eq$acc, eq$nombre)
      zon_pan <- data.frame(seqnames = factor(acc_de[zonas_df$chrom], levels = names(sl)),
                            start = zonas_df$start_Mb, end = zonas_df$end_Mb)
    }
    p <- ggplot(df_bins, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0, ymax = 1, fill = pct)) +
      geom_rect() +
      geom_segment(data = intrones_df, aes(x = start, xend = end, y = 0.75, yend = 0.75),
                   inherit.aes = FALSE, linewidth = 0.5, colour = "#6a040f", alpha = 0.7) +
      geom_segment(data = prom_df, aes(x = start, xend = end, y = 0.25, yend = 0.25),
                   inherit.aes = FALSE, linewidth = 0.8, colour = "#457b9d", alpha = 0.8) +
      
      ########ZONAS CRITICAS######################################################################################
    
      #{ if (!is.null(zon_pan)) geom_segment(data = zon_pan, aes(x = start, xend = end, y = 0.04, yend = 0.04),
      #                                      inherit.aes = FALSE, linewidth = 1.6, colour = pal$zona) } +
      #{ if (!is.null(tel_plot)) geom_point(data = tel_plot, aes(x = pos, y = 0.5), inherit.aes = FALSE,
      #                                     shape = "|", size = 3, colour = "#1d3557") } +
    
      #############################################################################################################  
      facet_grid(seqnames ~ ., switch = "y", labeller = as_labeller(lab)) +
      scale_fill_gradientn(colours = c("white", "#ffd166", "#e63946", "#6a040f"), name = "% repeticiones", limits = c(0, 100)) +
      scale_x_continuous(expand = c(0, 0), breaks = seq(0, ceiling(max(sl) / 1e6), 1)) +
      labs(x = "Posición en el cromosoma (Mb)", y = NULL,
           title = paste0("Panorama genómico por cromosoma (ventanas de ", ventana_label, ")"),
           subtitle = "Rojo oscuro: intrones | azul: promotores | | : telómero | naranja inferior: zonas críticas") +
      theme_minimal(base_family = font_family) +
      theme(axis.text.y = element_blank(), panel.grid = element_blank(),
            strip.text.y.left = element_text(angle = 0, hjust = 1), panel.spacing = unit(2, "pt"))
    print(p); guardar(p, "panorama_genomico", 10, 1.5 + 0.55 * n_filas)
    leyenda("F01. Panorama genómico", "Color = % de repeticiones por ventana. Segmentos: intrones (arriba), promotores (abajo), telómeros (|) y zonas críticas (naranja, de zonas_df).")
  })
  
  # ============================================================
  # F02. Perfiles por cromosoma con las zonas críticas sombreadas
  #      Línea de color = variable (eje izq.), línea negra = genes/Mb (eje der.)
  #   
  # ============================================================
  intentar("F02. Perfiles por cromosoma con zonas", {
    tp <- tabla
    cap <- quantile(d$genes_per_Mb, 0.99, na.rm = TRUE)
    lab_excl <- "Ventana excluida"; lab_st <- "Subtélomero (ST)"
    excl_d <- tabla[!tabla$analisis, ]; excl_d$chrom <- factor(excl_d$chrom, levels = eq$nombre)
    #####################################################
    ###Dibujar todas las zonas, se debe comentar las proximas 2 lineas
    #zon_p <- zonas_df; zon_p$chrom <- factor(zon_p$chrom, levels = eq$nombre)
    #####################################################
    zon_p <- zonas_df[as.character(zonas_df$tipo) == "Bajo GC + baja densidad génica", ]
    zon_p$chrom <- factor(zon_p$chrom, levels = eq$nombre)
    fills <- c(colores_zona, setNames(c("grey60", "black"), c(lab_excl, lab_st)))
    brk_fill <- names(fills)[names(fills) %in% c(unique(zon_p$tipo), lab_st, if (nrow(excl_d) > 0) lab_excl)]
    
    #################################
    ##########Por cromosomas#########
    #################################
    tp <- tabla
    
    variables <- c(
      "GC_percent", "repeat_pct", "intergenic_pct", "intron_pct"
    )
    
    etiquetas <- c(
      GC_percent = "GC (%)",
      repeat_pct = "Repeticiones (%)",
      intergenic_pct = "Intergénico (%)",
      intron_pct = "Intrones (%)",
      genes_per_Mb = "Genes/Mb"
    )
    
    for (chr in eq$nombre) {
      datos_chr <- tp[as.character(tp$chrom) == chr, ]
      
      datos_largos <- tidyr::pivot_longer(
        datos_chr,
        cols = dplyr::all_of(variables),
        names_to = "variable",
        values_to = "valor"
      )
      
      p <- ggplot2::ggplot(
        datos_largos,
        ggplot2::aes(x = position_Mb, y = valor, color = variable)
      ) +
        ggplot2::geom_line(linewidth = 0.6, na.rm = TRUE) +
        ggplot2::scale_color_discrete(labels = etiquetas) +
        ggplot2::labs(
          title = paste("Perfiles del cromosoma", chr),
          x = "Posición genómica (Mb)",
          y = "Valor",
          color = NULL
        ) +
        tema_art +
        ggplot2::theme(legend.position = "bottom")
      
      print(p)
      guardar(p, paste0("perfil_todas_lineas_", gsub("[^A-Za-z0-9]+", "_", chr)),
              9, 5)
    }
    
    # Figura general: todos los cromosomas en facetas
    datos_largos <- tidyr::pivot_longer(
      tp,
      cols = dplyr::all_of(variables),
      names_to = "variable",
      values_to = "valor"
    )
    
    datos_largos$chrom <- factor(
      as.character(datos_largos$chrom),
      levels = eq$nombre
    )
    
    p_facet <- ggplot2::ggplot(
      datos_largos,
      ggplot2::aes(
        x = position_Mb,
        y = valor,
        color = variable,
        group = variable
      )
    ) +
      ggplot2::geom_line(linewidth = 0.5, na.rm = TRUE) +
      ggplot2::facet_wrap(~ chrom, scales = "free_x", ncol = 3) +
      ggplot2::scale_color_discrete(labels = etiquetas) +
      ggplot2::labs(
        title = "Perfiles por cromosoma",
        x = "Posición genómica (Mb)",
        y = "Valor",
        color = NULL
      ) +
      tema_art +
      ggplot2::theme(legend.position = "bottom")
    
    print(p_facet)
    guardar(p_facet, "perfil_todos_cromosomas_facet", 12, 9)

    #################################################################################################
    
    fig_perfil <- function(ycol, serie, color, titulo, lim = c(0, 100)) {
      lo <- lim[1]; hi <- lim[2]
      esc <- function(g) lo + pmin(g, cap) / cap * (hi - lo)
      d_ln <- rbind(data.frame(chrom = factor(tp$chrom, levels = eq$nombre), pos = tp$position_Mb, y = tp[[ycol]], serie = serie),
                    data.frame(chrom = factor(tp$chrom, levels = eq$nombre), pos = tp$position_Mb, y = esc(tp$genes_per_Mb),
                               serie = "Genes/Mb (eje der.)"))
      ggplot() +
        ##############################
        #ventanas excluidas, se eliminan de la leyenda tambien 
        #geom_rect(data = excl_d, aes(xmin = start_Mb, xmax = end_Mb, ymin = -Inf, ymax = Inf, fill = lab_excl), alpha = 0.35) +
        ###############################
        geom_rect(data = zon_p, aes(xmin = start_Mb, xmax = end_Mb, ymin = -Inf, ymax = Inf, fill = tipo), alpha = 0.40) +
        geom_line(data = d_ln, aes(pos, y, color = serie), linewidth = 0.5, na.rm = TRUE) +
        geom_rug(data = genes_rug, aes(x = pos), sides = "b", inherit.aes = FALSE, alpha = 0.35, linewidth = 0.2, color = "grey20") +
        ########################
        #Subtelomero
        #{ if (!is.null(st_plot)) geom_rect(data = st_plot, aes(xmin = xmin, xmax = xmax, ymin = lo, ymax = lo + 0.06 * (hi - lo), fill = lab_st),
        #                                   alpha = 0.9) } +
        ################################
        #Nombre de las zonas
        #geom_text(data = zon_p, aes(x = (start_Mb + end_Mb) / 2, y = Inf, label = sub(".*_", "", zona_id)),
        #          inherit.aes = FALSE, vjust = 1.4, size = 2, fontface = "bold") +
        ##############################################
        facet_grid(chrom ~ ., switch = "y") +
        scale_color_manual(values = setNames(c(color, "black"), c(serie, "Genes/Mb (eje der.)")), name = NULL) +
        scale_fill_manual(values = fills, breaks = brk_fill, name = NULL) +
        scale_x_continuous(limits = c(0, max_length_Mb), expand = c(0, 0)) +
        scale_y_continuous(limits = c(lo, hi), breaks = pretty(c(lo, hi), n = 3),
                           sec.axis = sec_axis(~ (. - lo) / (hi - lo) * cap, name = NULL, breaks = pretty(c(0, cap), n = 3))) +
        guides(color = guide_legend(order = 1, override.aes = list(linewidth = 1)),
               fill = guide_legend(order = 2, override.aes = list(alpha = 0.6))) +
        labs(x = "Posición genómica (Mb)", y = NULL, title = titulo,
             subtitle = paste0("Zonas críticas = zonas_df (percentil ", 100 * umbral_pct, " por cromosoma) ")) +
        tema_art + theme(legend.position = "bottom", legend.box = "vertical", strip.background = element_blank(),
                         strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold", size = 8),
                         panel.spacing = unit(3, "pt"))
    }
    
    lim_gc <- c(floor(min(d$GC_percent, na.rm = TRUE)), ceiling(max(d$GC_percent, na.rm = TRUE)))
    defs <- list(
      list("GC_percent",     "GC (%)",                 color_GC,  "GC y densidad génica por cromosoma",            lim_gc, "perfil_GC"),
      list("repeat_pct",     "% repeticiones",         "#6a040f", "Repeticiones y densidad génica por cromosoma",   c(0, 100), "perfil_repeticiones"),
      list("intergenic_pct", "% intergénico",          "#6a040f", "Contenido intergénico y densidad génica",        c(0, 100), "perfil_intergenico"),
      list("intron_pct",     "% de intrones",          "#6a040f", "Contenido intrónico y densidad génica",          c(0, 100), "perfil_intrones"))
    for (df_ in defs) {
      p <- fig_perfil(df_[[1]], df_[[2]], df_[[3]], df_[[4]], df_[[5]])
      print(p); guardar(p, df_[[6]], 9, 1.4 * n_filas + 2, tiff = (df_[[6]] == "perfil_GC"))
    }
    leyenda("F02. Perfiles por cromosoma",
            paste0("Línea de color: variable (eje izquierdo); línea negra: genes/Mb (eje derecho, tope en el percentil 99). Franjas: zonas críticas de zonas_df ",
                   "(percentil ", 100 * umbral_pct, " por cromosoma), coloreadas por tipo. Una zona de 'baja densidad génica' NO tiene por qué mostrar GC bajo."))
  }, cond = hay_zonas)
  
  # ============================================================
  # F04. Zonas críticas vs resto + perfiles genómicos de las zonas
  # ============================================================
  intentar("F04. Zonas críticas vs resto y perfiles", {
    flanco_kb <- 100
    
    # Mantener la tabla resumen original: cada zona frente a su cromosoma
    resumen <- do.call(rbind, lapply(seq_len(nrow(zonas_df)), function(i) {
      z <- zonas_df[i, ]
      chr <- as.character(z$chrom)
      
      w <- tp[
        as.character(tp$chrom) == chr &
          tp$analisis &
          tp$end_Mb > z$start_Mb &
          tp$start_Mb < z$end_Mb,
      ]
      
      chr_d <- tp[
        as.character(tp$chrom) == chr & tp$analisis,
      ]
      
      data.frame(
        zona_id = z$zona_id,
        chrom = chr,
        tipo = as.character(z$tipo),
        inicio_Mb = z$start_Mb,
        fin_Mb = z$end_Mb,
        n_ventanas = nrow(w),
        GC_zona = mean(w$GC_percent, na.rm = TRUE),
        GC_mediana_chr = median(chr_d$GC_percent, na.rm = TRUE),
        genes_Mb_zona = mean(w$genes_per_Mb, na.rm = TRUE),
        genes_Mb_mediana_chr = median(chr_d$genes_per_Mb, na.rm = TRUE),
        repeat_zona = mean(w$repeat_pct, na.rm = TRUE),
        repeat_mediana_chr = median(chr_d$repeat_pct, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
    
    write.csv(
      resumen,
      paste0(carpeta_salida, "zonas_criticas_resumen_", ventana_label, ".csv"),
      row.names = FALSE
    )
    
    # Mantener el gráfico estadístico original de F04
    vv <- vars[names(vars) %in% names(tp)]
    dz <- tp[tp$analisis & !is.na(tp$zona_tipo), ]
    
    if (nrow(dz) > 0 && length(vv) > 0) {
      dl_lista <- lapply(names(vv), function(v) {
        keep <- !is.na(dz$zona_tipo) & !is.na(dz[[v]])
        if (!any(keep)) return(NULL)
        
        data.frame(
          grupo = dz$zona_tipo[keep],
          valor = dz[[v]][keep],
          variable = rep(unname(vv[[v]]), sum(keep))
        )
      })
      
      dl_lista <- Filter(Negate(is.null), dl_lista)
      
      if (length(dl_lista) > 0) {
        dl <- do.call(rbind, dl_lista)
        fills <- c(colores_zona, "Fuera de zonas" = "grey80")
        
        p_stat <- box_anot(
          dl,
          paste0("Ventanas dentro de zonas críticas vs resto (", ventana_label, ")"),
          "Zonas definidas en zonas_df; se muestran las ventanas incluidas en el análisis",
          fills
        )
        
        print(p_stat)
        guardar(p_stat, "zonas_criticas_vs_resto", 11, 7)
      } else {
        message("F04: no hay valores disponibles para el gráfico estadístico.")
      }
    } else {
      message("F04: no hay ventanas incluidas con valores para comparar.")
    }
    # Para los perfiles se usa tabla_zonas si existe;
    # si no, se muestran las zonas combinadas.
    if (exists("tabla_zonas") &&
        is.data.frame(tabla_zonas) &&
        nrow(tabla_zonas) > 0) {
      zz <- tabla_zonas
    } else {
      zz <- zonas_df[
        as.character(zonas_df$tipo) == "Bajo GC + baja densidad génica",
      ]
    }
    
    if (nrow(zz) == 0) {
      message("F04: no hay zonas para dibujar perfiles.")
    } else {
      perfiles <- vector("list", nrow(zz))
      names(perfiles) <- as.character(zz$zona_id)
      
      for (i in seq_len(nrow(zz))) {
        z <- zz[i, ]
        chr <- as.character(z$chrom)
        acc <- eq$acc[match(chr, eq$nombre)]
        
        if (length(acc) == 0 || is.na(acc)) next
        
        largo_chr <- as.numeric(sl[acc])
        zona_ini <- as.numeric(z$start_pb)
        zona_fin <- as.numeric(z$end_pb)
        
        ini <- max(1, zona_ini - flanco_kb * 1000)
        fin <- min(largo_chr, zona_fin + flanco_kb * 1000)
        xlim <- c(ini, fin) / 1e6
        
        banda_zona <- data.frame(
          xmin = zona_ini / 1e6,
          xmax = zona_fin / 1e6
        )
        
        # GC por ventana, incluyendo todos los valores disponibles

        # GC de alta resolución: ventana deslizante de 2 kb, paso de 500 pb
        win_gc <- 2000
        paso_gc <- 500
        
        secuencia_contexto <- Biostrings::subseq(
          fa[[acc]],
          start = ini,
          end = fin
        )
        
        gc_fino <- as.numeric(
          Biostrings::letterFrequencyInSlidingView(
            secuencia_contexto,
            win_gc,
            "GC",
            as.prob = TRUE
          )
        ) * 100
        idx <- seq(1, length(gc_fino), by = paso_gc)
        
        dgc <- data.frame(
          pos = (ini - 1 + idx + win_gc / 2) / 1e6,
          gc = gc_fino[idx]
        )
        
        # GC de cada ventana original: segmentos negros
        gc_win <- tabla[
          as.character(tabla$chrom) == chr &
            tabla$end >= ini &
            tabla$start <= fin &
            !is.na(tabla$GC_percent),
        ]
        
        # Referencia del cromosoma y umbral de GC
        gc_chr <- tp$GC_percent[
          as.character(tp$chrom) == chr & !is.na(tp$GC_percent)
        ]
        
        gc_mediana_chr <- if (length(gc_chr) > 0) {
          median(gc_chr, na.rm = TRUE)
        } else {
          NA_real_
        }
        
        gc_umbral_chr <- if (length(gc_chr) > 0) {
          as.numeric(quantile(gc_chr, umbral_pct, na.rm = TRUE))
        } else {
          NA_real_
        }
        
        # Punto destacado: mayor desviación de la mediana dentro de la zona
        gc_en_zona <- dgc[
          dgc$pos >= zona_ini / 1e6 &
            dgc$pos <= zona_fin / 1e6,
        ]
        
        if (nrow(gc_en_zona) > 0 && is.finite(gc_mediana_chr)) {
          pk <- gc_en_zona[
            which.max(abs(gc_en_zona$gc - gc_mediana_chr)),
          ]
        } else {
          pk <- data.frame(pos = numeric(), gc = numeric())
        }
        
        # ---- Paleta (apta para daltonismo) ----
        col_zona   <- "#B8C4D6"   # zona (gris azulado)
        col_zona_b <- "#7F93AD"   # borde de la zona
        col_curva  <- "#1B3A57"
        col_bajo   <- "#E07A5F"   # por debajo de la mediana
        col_alto   <- "#3FA796"   # por encima de la mediana
        col_med    <- "#E9A23B"   # mediana
        col_umb    <- "#6A4C93"   # umbral
        col_vent   <- "#4A5568"   # ventanas originales
        col_rep    <- "#6C757D"
        col_plus   <- "#176B87"
        col_minus  <- "#B23A48"
        
        x_zi <- zona_ini / 1e6
        x_zf <- zona_fin / 1e6
        rango_x <- xlim[2] - xlim[1]
        
        # ---- Repeticiones (estaba entre p_gc y p_tracks en el original) ----
        rep_plot <- data.frame(xmin = numeric(), xmax = numeric())
        
        if (exists("rep_red") && length(rep_red) > 0) {
          rr <- rep_red[
            as.character(seqnames(rep_red)) == acc &
              end(rep_red) >= ini &
              start(rep_red) <= fin
          ]
          
          if (length(rr) > 0) {
            rep_plot <- data.frame(
              xmin = pmax(start(rr), ini) / 1e6,
              xmax = pmin(end(rr), fin) / 1e6
            )
            # Ancho mínimo para que se vean como bloques
            rep_plot$xmax <- pmax(rep_plot$xmax, rep_plot$xmin + rango_x * 0.003)
          }
        }
        
        # ---- Genes ----
        genes_plot <- g_all[
          g_all$acc == acc &
            g_all$e >= ini &
            g_all$s <= fin,
        ]
        
        if (nrow(genes_plot) > 0) {
          genes_plot$x0 <- pmax(genes_plot$s, ini) / 1e6
          genes_plot$x1 <- pmin(genes_plot$e, fin) / 1e6
          genes_plot$x    <- ifelse(genes_plot$str == "-", genes_plot$x1, genes_plot$x0)
          genes_plot$xend <- ifelse(genes_plot$str == "-", genes_plot$x0, genes_plot$x1)
        }
        # Etiquetar los genes dentro de la zona y el más próximo a cada lado
        genes_etiqueta <- genes_plot
        
        if (nrow(genes_plot) > 0) {
          dentro <- which(
            genes_plot$s <= zona_fin &
              genes_plot$e >= zona_ini
          )
          
          antes <- which(genes_plot$e < zona_ini)
          despues <- which(genes_plot$s > zona_fin)
          
          vecino_antes <- if (length(antes) > 0) {
            antes[which.max(genes_plot$e[antes])]
          } else {
            integer(0)
          }
          
          vecino_despues <- if (length(despues) > 0) {
            despues[which.min(genes_plot$s[despues])]
          } else {
            integer(0)
          }
          
          idx_etiqueta <- unique(c(dentro, vecino_antes, vecino_despues))
          
          genes_etiqueta <- genes_plot[idx_etiqueta, ]
          genes_etiqueta <- genes_etiqueta[
            !is.na(genes_etiqueta$lab) & nzchar(genes_etiqueta$lab),
          ]
        }
        
        
        # ---- Insertar los cruces exactos con la mediana (evita franjas falsas) ----
        cruce_mediana <- function(d, m) {
          out <- vector("list", nrow(d) * 2)
          n <- 0
          for (k in seq_len(nrow(d) - 1)) {
            n <- n + 1; out[[n]] <- d[k, c("pos", "gc")]
            if ((d$gc[k] - m) * (d$gc[k + 1] - m) < 0) {
              t <- (m - d$gc[k]) / (d$gc[k + 1] - d$gc[k])
              n <- n + 1
              out[[n]] <- data.frame(pos = d$pos[k] + t * (d$pos[k + 1] - d$pos[k]),
                                     gc  = m)
            }
          }
          n <- n + 1; out[[n]] <- d[nrow(d), c("pos", "gc")]
          do.call(rbind, out[seq_len(n)])
        }
        
        dgc2 <- cruce_mediana(dgc, gc_mediana_chr)
        dgc2$y_bajo <- pmin(dgc2$gc, gc_mediana_chr)
        dgc2$y_alto <- pmax(dgc2$gc, gc_mediana_chr)
        
        # Límites de las ventanas originales recortados al contexto
        lim_vent <- sort(unique(c(gc_win$start_Mb, gc_win$end_Mb)))
        lim_vent <- lim_vent[lim_vent > xlim[1] & lim_vent < xlim[2]]
        
        ylim_gc <- range(c(dgc$gc, gc_win$GC_percent), na.rm = TRUE) + c(-1.2, 0.8)
        
        # ---- Panel A: GC ----
        # Límites de las ventanas originales dentro del contexto
        lim_vent <- sort(unique(c(gc_win$start_Mb, gc_win$end_Mb)))
        lim_vent <- lim_vent[lim_vent > xlim[1] & lim_vent < xlim[2]]
        
        ylim_gc <- range(dgc$gc, na.rm = TRUE) + c(-1.2, 1.0)
        
        p_gc <- ggplot() +
          # Zona crítica
          annotate("rect", xmin = x_zi, xmax = x_zf, ymin = -Inf, ymax = Inf,
                   fill = col_zona, alpha = 0.35) +
          annotate("segment", x = c(x_zi, x_zf), xend = c(x_zi, x_zf),
                   y = -Inf, yend = Inf, colour = col_zona_b, linewidth = 0.4) +
          # Límites de ventana: verticales y tenues
          geom_vline(xintercept = lim_vent, colour = "grey70",
                     linewidth = 0.3, linetype = "dashed", alpha = 0.8) +
          # Referencias
          geom_hline(yintercept = gc_mediana_chr, colour = col_med,
                     linetype = "dashed", linewidth = 0.7) +
          geom_hline(yintercept = gc_umbral_chr, colour = col_umb,
                     linetype = "dotted", linewidth = 0.8) +
          # Curva
          geom_line(data = dgc, aes(x = pos, y = gc),
                    colour = col_curva, linewidth = 0.6) +
          # Rótulos a la derecha, separados para que no se pisen
          annotate("text", x = xlim[2], y = gc_mediana_chr,
                   label = paste0("Mediana ", sprintf("%.1f%%", gc_mediana_chr)),
                   hjust = -0.05, vjust = -0.4, size = 2.9,
                   colour = col_med, fontface = "bold") +
          annotate("text", x = xlim[2], y = gc_umbral_chr,
                   label = paste0("Umbral P", umbral_pct * 100, " ",
                                  sprintf("%.1f%%", gc_umbral_chr)),
                   hjust = -0.05, vjust = 1.3, size = 2.9,
                   colour = col_umb, fontface = "bold") +
          # Etiqueta de la zona
          annotate("label", x = (x_zi + x_zf) / 2, y = Inf,
                   label = as.character(z$zona_id),
                   vjust = 1.3, size = 3.2, fontface = "bold",
                   fill = col_zona_b, colour = "white", linewidth = 0) +
          # Mínimo de GC en la zona
          { if (nrow(pk) > 0) list(
            geom_point(data = pk, aes(x = pos, y = gc),
                       shape = 21, fill = "white", colour = col_bajo,
                       size = 2.8, stroke = 1.2),
            geom_label(data = pk,
                       aes(x = pos, y = gc, label = sprintf("%.1f%% GC", gc)),
                       hjust = -0.2, vjust = 0.5, size = 3, fontface = "bold",
                       fill = "white", colour = col_bajo, linewidth = 0.25)
          ) } +
          labs(y = "Contenido de GC (%)", x = NULL) +
          coord_cartesian(xlim = xlim, ylim = ylim_gc, expand = FALSE, clip = "off") +
          tema_art +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                panel.grid = element_blank(),
                plot.margin = margin(5, 90, 5, 5))
        
        
        # ---- Panel B: genes + repeticiones ----
        p_tracks <- ggplot() +
          annotate("rect", xmin = x_zi, xmax = x_zf, ymin = 0, ymax = 2,
                   fill = col_zona, alpha = 0.35) +
          annotate("segment", x = c(x_zi, x_zf), xend = c(x_zi, x_zf),
                   y = 0, yend = 2, colour = col_zona_b, linewidth = 0.4) +
          geom_vline(xintercept = lim_vent, colour = "grey70",
                     linewidth = 0.3, linetype = "dashed", alpha = 0.8) +
          { if (nrow(rep_plot) > 0)
            geom_rect(data = rep_plot,
                      aes(xmin = xmin, xmax = xmax, ymin = 0.2, ymax = 0.8),
                      fill = col_rep, colour = NA) } +
          { if (nrow(genes_plot) > 0)
            geom_segment(data = genes_plot,
                         aes(x = x, xend = xend, y = 1.45, yend = 1.45, colour = str),
                         linewidth = 1.6,
                         arrow = grid::arrow(length = grid::unit(2.2, "mm"),
                                             type = "closed"),
                         lineend = "butt", linejoin = "mitre") } +
          {
            if (nrow(genes_etiqueta) > 0) {
              if (requireNamespace("ggrepel", quietly = TRUE)) {
                ggrepel::geom_text_repel(
                  data = genes_etiqueta,
                  aes(
                    x = (s + e) / 2e6,
                    y = 1.78,
                    label = lab
                  ),
                  inherit.aes = FALSE,
                  direction = "x",
                  size = 2.6,
                  colour = "#263746",
                  box.padding = 0.22,
                  point.padding = 0.12,
                  min.segment.length = 0,
                  segment.colour = "grey50",
                  segment.size = 0.25,
                  max.overlaps = Inf,
                  seed = 42
                )
              } else {
                geom_text(
                  data = genes_etiqueta,
                  aes(
                    x = (s + e) / 2e6,
                    y = 1.78,
                    label = lab
                  ),
                  inherit.aes = FALSE,
                  size = 2.4,
                  colour = "#263746",
                  check_overlap = TRUE,
                  vjust = 0
                )
              }
            }
          } +
          scale_colour_manual(values = c("+" = col_plus, "-" = col_minus),
                              breaks = c("+", "-"),
                              labels = c("Gen en hebra +", "Gen en hebra −"),
                              name = NULL) +
          scale_y_continuous(breaks = c(0.5, 1.45),
                             labels = c("Repeticiones", "Genes"),
                             limits = c(0, 2.1), expand = c(0, 0)) +
          scale_x_continuous(labels = function(x) sprintf("%.2f", x)) +
          labs(x = paste0("Posición en el cromosoma ", chr, " (Mb)"), y = NULL) +
          coord_cartesian(xlim = xlim, expand = FALSE) +
          tema_art +
          theme(axis.text.y = element_text(face = "bold", size = 9),
                panel.grid = element_blank(),
                plot.margin = margin(5, 70, 5, 5))
        
        perfiles[[i]] <- (p_gc / p_tracks) +
          patchwork::plot_layout(heights = c(2.4, 0.8), guides = "collect") +
          patchwork::plot_annotation(
            title = paste0("Zona crítica ", z$zona_id, " (cromosoma ", chr, ")"),
            subtitle = paste0("Curva: GC con ventana de ", win_gc / 1000, " kb | ",
                              "líneas verticales discontinuas: límites de las ventanas originales | ",
                              "banda gris: zona ± ", flanco_kb, " kb de contexto"),
            tag_levels = "A"
          ) &
          theme(legend.position = "bottom",
                plot.title = element_text(face = "bold", size = 13, hjust = 0),
                plot.subtitle = element_text(size = 8.5, colour = "grey30", hjust = 0),
                plot.tag = element_text(face = "bold", size = 14))
      }
      perfiles <- Filter(Negate(is.null), perfiles)
      
      if (length(perfiles) > 0) {
        # PDF multipágina: una página por zona
        archivo_pdf <- paste0(
          carpeta_salida, "F04_perfiles_zonas_", ventana_label, ".pdf"
        )
        
        grDevices::cairo_pdf(
          archivo_pdf,
          width = 10,
          height = 6,
          onefile = TRUE
        )
        for (p in perfiles) print(p)
        grDevices::dev.off()
        
        # PNG individual por zona
        for (nm in names(perfiles)) {
          nombre_seguro <- gsub("[^A-Za-z0-9_-]", "-", nm)
          
          ggsave(
            filename = paste0(
              carpeta_salida, "F04_perfil_", nombre_seguro,
              "_", ventana_label, ".png"
            ),
            plot = perfiles[[nm]],
            width = 10,
            height = 6,
            dpi = 300,
            bg = "white"
          )
        }
        
        # Mostrar el primer perfil en RStudio
        print(perfiles[[1]])
      }
    }
    
    leyenda(
      "F04. Zonas críticas vs resto y perfiles",
      paste0(
        "Se conserva el resumen estadístico zonas vs resto y el CSV por zona. ",
        "Perfiles: GC por ventana; genes con flecha según hebra; repeticiones como bloques ",
        "en una pista única. La franja naranja marca la zona y se incluyen ",
        flanco_kb, " kb de contexto."
      )
    )
  }, cond = hay_zonas && exists("vars"))
  
  # ============================================================
  # F05. GC vs genes/Mb por cromosoma (con rho de Spearman)
  # ============================================================
  intentar("F05. GC vs genes/Mb", {
    rho_c <- sapply(split(d, d$chrom), function(s) suppressWarnings(cor(s$GC_percent, s$genes_per_Mb, method = "spearman", use = "complete.obs")))
    n_c <- table(d$chrom)
    dd <- d; ch <- as.character(dd$chrom)
    etq <- sprintf("%s (rho = %.2f, n = %d)", eq$nombre, rho_c[eq$nombre], as.integer(n_c[eq$nombre]))
    dd$chrom_lab <- factor(sprintf("%s (rho = %.2f, n = %d)", ch, rho_c[ch], as.integer(n_c[ch])), levels = etq)
    p <- ggplot(dd, aes(GC_percent, genes_per_Mb)) + geom_point(alpha = 0.5, size = 0.8, color = "grey30") +
      geom_smooth(method = "lm", se = FALSE, color = color_genes, linewidth = 0.6, formula = y ~ x) +
      facet_wrap(~ chrom_lab, scales = "free", ncol = 4) +
      labs(x = "GC (%)", y = "Genes/Mb", title = "Relación GC vs densidad génica por cromosoma",
           subtitle = sprintf("Mediana de rho (Spearman) entre cromosomas = %.2f", median(rho_c, na.rm = TRUE))) +
      tema_art + theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 7))
    print(p); guardar(p, "scatter_GC_vs_genes_rho", 10, 8)
    leyenda("F05. GC vs genes/Mb", "Cada punto es una ventana del conjunto de análisis; recta = ajuste lineal; rho = Spearman por cromosoma.")
  })
  
  # ============================================================
  # F06. Extremos: perfil por distancia al borde + extremos vs interior (con estadística)
  # ============================================================
  intentar("F06. Extremos vs interior", {
    etq_v <- c(
      repeat_pct = "% repeticiones",
      genes_per_Mb = "Genes/Mb",
      GC_percent = "GC (%)"
    )
    etq_v <- etq_v[names(etq_v) %in% names(db)]
    
    gr_e <- factor(
      ifelse(
        db$dist_extremo_pb <= ext_ref_kb * 1000,
        "Extremo",
        "Interior"
      ),
      levels = c("Extremo", "Interior")
    )
    
    dl <- do.call(rbind, lapply(names(etq_v), function(v) {
      data.frame(
        grupo = gr_e,
        valor = db[[v]],
        variable = etq_v[[v]]
      )
    }))
    
    p2 <- box_anot(
      dl,
      paste0("Extremos (≤ ", ext_ref_kb, " kb del borde) vs interior"),
      "Delta de Cliff > 0: valores mayores en el extremo",
      c(Extremo = "#B2182B", Interior = "grey80"),
      ncol = 3
    )
    
    print(p2)
    guardar(p2, "extremos_vs_interior_stats", 8, 4)
    
    leyenda(
      "F06. Extremos",
      paste0(
        "Boxplots de GC, repeticiones y Genes/Mb en extremos e interior. ",
        "Wilcoxon de dos colas y delta de Cliff (-1 a 1); extremo = distancia ≤ ",
        ext_ref_kb, " kb."
      )
    )
  }, cond = exists("db") && "dist_extremo_pb" %in% names(db))

  # ============================================================
  # F07. Alta vs baja repetición (con estadística)
  # ============================================================
  intentar("F7. Alta vs baja repetición", {
    vb <- c(intron_pct = "% intrones", exon_pct = "% exones", gene_pct = "% cobertura génica", intergenic_pct = "% intergénico",
            genes_per_Mb = "Genes/Mb", GC_percent = "GC (%)", introns_por_gen = "Intrones por gen")
    vb <- vb[names(vb) %in% names(d)]
    dg <- d[d$alto_repeat | d$bajo_repeat, ]
    gr <- factor(ifelse(dg$alto_repeat, "Alta repetición", "Baja repetición"), levels = c("Alta repetición", "Baja repetición"))
    dl <- do.call(rbind, lapply(names(vb), function(v) data.frame(grupo = gr, valor = dg[[v]], variable = vb[[v]])))
    p <- box_anot(dl, "Ventanas de alta vs baja repetición",
                  paste0("Percentiles ", 100 * umbral_pct, "/", 100 * (1 - umbral_pct), " por cromosoma"),
                  c("Alta repetición" = "#B2182B", "Baja repetición" = "grey80"))
    print(p); guardar(p, "alta_vs_baja_repeticion_stats", 11, 6)
    leyenda("F09. Alta vs baja repetición", "Wilcoxon de dos colas y delta de Cliff por variable; n = ventanas por grupo.")
  }, cond = tiene_rep)
  
  # ============================================================
  # F08 Ideograma de satelites y microsatelites. 
  # ============================================================
  intentar("F08. Idiograma", {
    ce <- carac_extremos
    
    # ---- Opciones ----
    mostrar_nombres_micro <- TRUE   # TRUE: añade el nombre de cada microsatélite
    
    # ---- Cromosomas ----
    chr_df <- data.frame(
      chrom = eq$nombre,
      L = as.numeric(sl[eq$acc]) / 1e6,
      stringsAsFactors = FALSE
    )
    chr_df$y <- rev(seq_len(nrow(chr_df))) * 2.2
    yi    <- setNames(chr_df$y, chr_df$chrom)
    Li    <- setNames(chr_df$L, chr_df$chrom)
    xmax  <- max(chr_df$L, na.rm = TRUE)
    alto_chr <- 0.30
    rx       <- xmax * 0.014        # radio horizontal de las puntas redondeadas
    
    # ---- Helper: polígono en forma de cápsula (recortado al cromosoma) ----
    capsula <- function(x0, x1, L, y, h, rx, id, n = 60) {
      r  <- min(rx, L / 2)
      x0 <- max(0, x0); x1 <- min(L, x1)
      x  <- seq(x0, x1, length.out = n)
      d  <- pmax(r - x, x - (L - r), 0)
      off <- h * sqrt(pmax(0, 1 - (d / r)^2))
      data.frame(x = c(x, rev(x)), y = c(y + off, rev(y - off)), id = id)
    }
    armar <- function(d) {
      if (is.null(d) || nrow(d) == 0) return(NULL)
      do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
        p <- capsula(d$x0[i], d$x1[i], d$L[i], d$y[i], alto_chr, rx, i)
        p$tipo <- d$tipo[i]
        p
      }))
    }
    
    # Contorno/relleno de cada cromosoma completo
    chr_poly <- do.call(rbind, lapply(seq_len(nrow(chr_df)), function(i) {
      p <- capsula(0, chr_df$L[i], chr_df$L[i], chr_df$y[i], alto_chr, rx, i)
      p
    }))
    
    # ---- Subtelómeros ----
    #ok_st <- !is.na(ce$ST_start) &
    #  ce$ST_long_pb > 0 &
    #  ce$ST_detectada == "Si"
    ok_st <- !is.na(ce$ST_start) &
      !is.na(ce$ST_end) &
      ce$ST_long_pb > 0
    st_d <- data.frame(
      chrom = as.character(ce$chrom[ok_st]),
      x0 = ce$ST_start[ok_st] / 1e6,
      x1 = ce$ST_end[ok_st] / 1e6,
      tipo = "Subtelómero",
      stringsAsFactors = FALSE
    )
    st_d$y <- unname(yi[st_d$chrom])
    st_d$L <- unname(Li[st_d$chrom])
    st_d <- st_d[!is.na(st_d$y), ]
    st_poly <- armar(st_d)
    
    # ---- Familias de satélites (dentro del cromosoma) ----
    sat_d <- data.frame()
    if (exists("sat_bloques") && nrow(sat_bloques) > 0) {
      sat_d <- data.frame(
        chrom = as.character(sat_bloques$chrom),
        x0 = sat_bloques$start / 1e6,
        x1 = sat_bloques$end / 1e6,
        tipo = as.character(sat_bloques$satelite),
        stringsAsFactors = FALSE
      )
      sat_d$x1 <- pmax(sat_d$x1, sat_d$x0 + xmax * 0.003)  # ancho mínimo visible
      sat_d$y <- unname(yi[sat_d$chrom])
      sat_d$L <- unname(Li[sat_d$chrom])
      sat_d <- sat_d[!is.na(sat_d$y), ]
    }
    sat_poly <- armar(sat_d)
    
    # ---- Microsatélites ----
    micro_d <- data.frame()
    if (exists("microsat_loci") && !is.null(microsat_loci) &&
        nrow(microsat_loci) > 0) {
      micro_d <- data.frame(
        chrom = as.character(microsat_loci$chrom),
        x = microsat_loci$pos_mb,
        marcador = ifelse(
          !is.na(microsat_loci$nombre) & nzchar(microsat_loci$nombre),
          microsat_loci$nombre,
          microsat_loci$marcador
        ),
        tipo = "Microsatélite",
        stringsAsFactors = FALSE
      )
      micro_d$y_chr <- unname(yi[micro_d$chrom])
      micro_d <- micro_d[!is.na(micro_d$y_chr), ]
      micro_d$y_ini <- micro_d$y_chr - alto_chr
      micro_d$y_pt  <- micro_d$y_chr - alto_chr - 0.32
    }
    
    # ---- Paleta (sobria, tipo revista) ----
    col_sat <- c(sat240 = "#7B5EA7", sat350 = "#E08E45")
    if (nrow(sat_d) > 0) {
      faltan <- setdiff(unique(sat_d$tipo), names(col_sat))
      if (length(faltan) > 0) {
        extra <- grDevices::hcl.colors(length(faltan), "Dark 3")
        names(extra) <- faltan
        col_sat <- c(col_sat, extra)
      }
    }
    cf <- c(
      "Subtelómero"   = "#FF6B6B",
      col_sat,
      "Microsatélite" = "#0F7B8A"
    )
    
    tipos_presentes <- unique(c(
      st_d$tipo,
      if (nrow(sat_d) > 0) sat_d$tipo,
      if (nrow(micro_d) > 0) micro_d$tipo
    ))
    breaks_fill <- names(cf)[names(cf) %in% tipos_presentes]
    
    hay_repel <- requireNamespace("ggrepel", quietly = TRUE)
    
    # Etiquetas de longitud (Mb) al final de cada cromosoma
    chr_df$lab_L <- sprintf("%.2f Mb", chr_df$L)
    
    # ---- Gráfico ----
    p <- ggplot() +
      # Cuerpo del cromosoma (relleno)
      geom_polygon(
        data = chr_poly,
        aes(x = x, y = y, group = id),
        fill = "#EEF1F4", colour = NA
      ) +
      
      # Subtelómeros
      {
        if (!is.null(st_poly))
          geom_polygon(
            data = st_poly,
            aes(x = x, y = y, group = id, fill = tipo),
            colour = NA
          )
      } +
      
      # Satélites
      {
        if (!is.null(sat_poly))
          geom_polygon(
            data = sat_poly,
            aes(x = x, y = y, group = id, fill = tipo),
            colour = NA
          )
      } +
      
      # Contorno fino del cromosoma
      geom_polygon(
        data = chr_poly,
        aes(x = x, y = y, group = id),
        fill = NA, colour = "#4A5568", linewidth = 0.3
      ) +
      
      # Longitud de cada cromosoma
      geom_text(
        data = chr_df,
        aes(x = L + xmax * 0.012, y = y, label = lab_L),
        hjust = 0, size = 2.5, colour = "#6B7280"
      ) +
      
      # Microsatélites: línea fina + triángulo bajo el cromosoma
      {
        if (nrow(micro_d) > 0)
          geom_segment(
            data = micro_d,
            aes(x = x, xend = x, y = y_ini, yend = y_pt),
            colour = "#0F7B8A", linewidth = 0.25
          )
      } +
      {
        if (nrow(micro_d) > 0)
          geom_point(
            data = micro_d,
            aes(x = x, y = y_pt, fill = tipo),
            shape = 24, size = 1.9, stroke = 0.25, colour = "white"
          )
      } +
      
      # Nombres de microsatélites
      {
        if (mostrar_nombres_micro && nrow(micro_d) > 0) {
          if (hay_repel) {
            ggrepel::geom_text_repel(
              data = micro_d,
              aes(x = x, y = y_pt, label = marcador),
              inherit.aes = FALSE,
              nudge_y = -0.38, direction = "x",
              size = 2.3, colour = "#14535D", fontface = "italic",
              min.segment.length = 0, segment.size = 0.15,
              segment.colour = "grey70", max.overlaps = Inf, seed = 42
            )
          } else {
            geom_text(
              data = micro_d,
              aes(x = x, y = y_pt - 0.38, label = marcador),
              inherit.aes = FALSE,
              size = 2.3, colour = "#14535D", fontface = "italic",
              check_overlap = TRUE
            )
          }
        }
      } +
      
      scale_fill_manual(
        values = cf, breaks = breaks_fill, name = NULL, drop = FALSE,
        guide = guide_legend(
          nrow = 1,
          override.aes = list(shape = 22, size = 3.5, colour = "#4A5568",
                              stroke = 0.2, alpha = 1)
        )
      ) +
      scale_x_continuous(
        name = "Posición cromosómica (Mb)",
        limits = c(0, xmax * 1.12),
        breaks = scales::pretty_breaks(n = 6),
        expand = c(0, 0)
      ) +
      scale_y_continuous(
        breaks = chr_df$y,
        labels = chr_df$chrom,
        limits = c(min(chr_df$y) - 1.1, max(chr_df$y) + 0.7),
        expand = c(0, 0)
      ) +
      labs(title = NULL, y = NULL) +
      tema_art +
      theme(
        panel.background   = element_rect(fill = "white", colour = NA),
        plot.background    = element_rect(fill = "white", colour = NA),
        legend.position    = "top",
        legend.justification = "left",
        legend.key.size    = grid::unit(3, "mm"),
        legend.key         = element_blank(),
        legend.text        = element_text(size = 8, colour = "#374151"),
        legend.margin      = margin(0, 0, 2, 0),
        axis.title.x       = element_text(size = 8.5, colour = "#374151",
                                          margin = margin(t = 5)),
        axis.text.x        = element_text(size = 8, colour = "#4B5563"),
        axis.line.x        = element_line(colour = "#4A5568", linewidth = 0.3),
        axis.ticks.x       = element_line(colour = "#4A5568", linewidth = 0.3),
        axis.ticks.length.x = grid::unit(1.5, "mm"),
        axis.title.y       = element_blank(),
        axis.line.y        = element_blank(),
        axis.ticks.y       = element_blank(),
        axis.text.y        = element_text(face = "bold", size = 8.5,
                                          colour = "#1F2937"),
        panel.grid.major   = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.border       = element_blank(),
        plot.margin        = margin(6, 10, 6, 6)
      )
    
    print(p)
    # 7.1 in ≈ 180 mm (doble columna Nature)
    guardar(p, "idiograma_resumen", 7.1, 1.4 + 0.55 * nrow(eq))
    
    leyenda(
      "F13. Idiograma",
      paste0(
        "Cada barra representa un cromosoma, con su longitud (Mb) indicada a la derecha. ",
        "Los subtelómeros y las familias de satélites se muestran como bloques de color ",
        "dentro del cromosoma. Los triángulos bajo cada barra indican la posición de los ",
        "microsatélites detectados con primers."
      )
    )
  }, cond = exists("carac_extremos"))
  
  # ============================================================
  # F09. Presencia/ausencia de genes en los extremos (un panel por cromosoma y extremo)
  # ============================================================
  intentar("F09. Presencia/ausencia de genes en los extremos", {
    ventana_kb <- 200; min_gen_kb <- 0.4; paso_kb <- 20
    col_ST <- "#7ac943"; col_telo <- "#00a0e3"
    ce <- carac_extremos
    ce$acc <- eq$acc[match(ce$chrom, eq$nombre)]
    ce <- ce[!is.na(ce$acc), ]
    genes_df <- g_all[, c("acc", "s", "e")]
    
    paneles <- vector("list", nrow(ce)); gen_l <- paneles; st_l <- paneles; tel_l <- paneles
    for (i in seq_len(nrow(ce))) {
      r <- ce[i, ]; L <- as.numeric(r$long_cromosoma_pb); lado <- as.character(r$extremo)
      id <- paste0("ST_", r$chrom, "_", lado)
      st_len <- ifelse(is.na(r$ST_long_pb), 0, r$ST_long_pb)
      #W <- min(max(ventana_kb * 1000, ceiling(st_len / 1e5) * 1e5), floor(L / 2))
      flanco_st_kb <- 100  # contexto adicional hacia el interior del cromosoma
      
      W_base <- max(
        ventana_kb * 1000,
        ceiling(st_len / 1e5) * 1e5
      )
      
      W <- min(
        W_base + flanco_st_kb * 1000,
        floor(L / 2)
      )
      a <- if (lado == "inicio") 0 else L - W
      b <- if (lado == "inicio") W else L
      paneles[[i]] <- data.frame(panel = id, chrom = r$chrom, lado = lado, a = a / 1e3, b = b / 1e3, L = L / 1e3,
                                 tiene_telo = r$tiene_telomero == "Si", stringsAsFactors = FALSE)
      g <- genes_df[genes_df$acc == r$acc & genes_df$e >= a & genes_df$s <= b, ]
      if (nrow(g) > 0)
        gen_l[[i]] <- data.frame(panel = id, xmin = pmax(g$s, a) / 1e3, xmax = pmin(g$e, b) / 1e3, stringsAsFactors = FALSE)
      if (!is.na(r$ST_start) && st_len > 0) {
        borde_int <- if (lado == "inicio") r$ST_end else r$ST_start
        st_l[[i]] <- data.frame(panel = id, xmin = max(r$ST_start, a) / 1e3, xmax = min(r$ST_end, b) / 1e3,
                                x_int = borde_int / 1e3, valido = r$ST_detectada == "Si", stringsAsFactors = FALSE)
      }
      if (r$tiene_telomero == "Si")
        tel_l[[i]] <- data.frame(panel = id, x = if (lado == "inicio") 0 else L / 1e3, stringsAsFactors = FALSE)
    }
    pan <- do.call(rbind, paneles); gens <- do.call(rbind, gen_l)
    sts <- do.call(rbind, st_l);    tels <- do.call(rbind, tel_l)
    
    orden_chr <- eq$nombre[eq$nombre %in% pan$chrom]
    niv_panel <- as.vector(rbind(paste0("ST_", orden_chr, "_inicio"), paste0("ST_", orden_chr, "_final")))
    niv_panel <- niv_panel[niv_panel %in% pan$panel]
    fac <- function(x) { if (!is.null(x)) x$panel <- factor(x$panel, levels = niv_panel); x }
    pan <- fac(pan); gens <- fac(gens); sts <- fac(sts); tels <- fac(tels)
    
    regla <- pan; regla$lab <- paste0(round(regla$b - regla$a), " kb"); regla$xm <- (regla$a + regla$b) / 2
    sin_telo <- pan[!pan$tiene_telo, ]
    sin_telo$x  <- ifelse(sin_telo$lado == "inicio", sin_telo$a, sin_telo$b)
    sin_telo$hj <- ifelse(sin_telo$lado == "inicio", 1.3, -0.3)
    
    p <- ggplot() +
      geom_segment(data = regla, aes(x = a, xend = b, y = 1.0, yend = 1.0),
                   arrow = arrow(ends = "both", length = unit(0.06, "in"), type = "open"), linewidth = 0.3) +
      geom_text(data = regla, aes(x = xm, y = 1.0, label = lab), vjust = -0.5, size = 2.6) +
      geom_segment(data = pan, aes(x = a, xend = b, y = 0.55, yend = 0.55), color = "grey70", linewidth = 0.2) +
      { if (!is.null(gens)) geom_rect(data = gens, aes(xmin = xmin, xmax = pmax(xmax, xmin + min_gen_kb), ymin = 0.42, ymax = 0.68),
                                      fill = "black", color = NA) } +
      { if (!is.null(sts)) geom_segment(data = sts, aes(x = xmin, xend = xmax, y = 0.2, yend = 0.2, linetype = valido),
                                        color = col_ST, linewidth = 1.1) } +
      { if (!is.null(sts)) geom_segment(data = sts, aes(x = x_int, xend = x_int, y = 0.2, yend = 0.38), color = col_ST, linewidth = 1.1) } +
      { if (!is.null(tels)) geom_segment(data = tels, aes(x = x, xend = x, y = 0.05, yend = 0.8), color = col_telo, linewidth = 1.6) } +
      geom_text(data = sin_telo, aes(x = x, y = 0.2, label = "?", hjust = hj), color = "red", size = 3, fontface = "bold") +
      facet_wrap(~ panel, ncol = 2, scales = "free_x", strip.position = "left") +
      scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22"), guide = "none") +
      scale_y_continuous(limits = c(0, 1.25), breaks = NULL, expand = c(0, 0)) +
      scale_x_continuous(breaks = scales::breaks_width(paso_kb), expand = expansion(mult = 0.01)) +
      labs(x = "Posición en el cromosoma (kb)", y = NULL,
           title = "Presencia y ausencia de genes en los extremos cromosómicos",
           subtitle = "Negro: genes | Verde: subtélomero (discontinuo = tramo no válido) | Azul: telómero | ? rojo: extremo sin telómero") +
      tema_art + theme(strip.background = element_blank(), strip.placement = "outside",
                       strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7),
                       axis.line.y = element_blank(), panel.grid = element_blank(),
                       panel.spacing.y = unit(6, "pt"), panel.spacing.x = unit(14, "pt"), axis.text.x = element_text(size = 6))
    print(p); guardar(p, "presencia_ausencia_genes_ST", 11, 0.9 * length(orden_chr) + 1.5, tiff = TRUE)
    leyenda("F16. Genes en los extremos", "Cada panel es un extremo (ventana de 200 kb, o más si el ST es mayor). Negro: genes; verde: ST; azul: telómero; '?': extremo sin telómero.")
  }, cond = length(genes_gr) > 0 && exists("carac_extremos") && nrow(carac_extremos) > 0)
  
} # cierra el bloque if (fig_modo == "genoma")

if (fig_modo == "comparacion") {
   # ============================================================
  # MODO COMPARACION: figuras entre genomas (datos de 16_comparacion.R)
  # ============================================================
  guardar_cmp <- function(p, nombre, ancho, alto) {
    ggsave(paste0(carpeta_cmp, nombre, ".pdf"), p, width = ancho, height = alto, device = cairo_pdf, limitsize = FALSE)
    ggsave(paste0(carpeta_cmp, nombre, ".png"), p, width = ancho, height = alto, dpi = 300, bg = "white", limitsize = FALSE)
  }
  col_sat350 <- "#B13C67"
  
  # G1. Un panel por cromosoma, una pista por genoma: SAT350 + microsatélites (reemplaza microysat.R)
  intentar("G1. Mapas por cromosoma", {
    mic_g <- agrupar_microsat(cmp_micro)
    dibujar_chr <- function(chr) {
      pistas <- data.frame(genoma = cmp_nombres, y = rev(seq_along(cmp_nombres)) * 2, stringsAsFactors = FALSE)
      pistas <- merge(pistas, cmp_long[cmp_long$chrom == chr, c("genoma", "long_pb")], by = "genoma", all.x = TRUE)
      validas <- pistas[!is.na(pistas$long_pb), ]
      if (nrow(validas) == 0) return(NULL)
      sat <- if (!is.null(cmp_sat350)) merge(cmp_sat350[cmp_sat350$chrom == chr, ], pistas[, c("genoma", "y")], by = "genoma") else NULL
      mic <- if (!is.null(mic_g)) merge(mic_g[mic_g$chrom == chr, ], pistas[, c("genoma", "y")], by = "genoma") else NULL
      xmax <- max(pistas$long_pb, na.rm = TRUE) / 1e6
      p <- ggplot() +
        geom_segment(data = validas, aes(x = 0, xend = long_pb / 1e6, y = y, yend = y), linewidth = 7, colour = "#A7B0BC", lineend = "round") +
        geom_segment(data = validas, aes(x = 0, xend = long_pb / 1e6, y = y, yend = y), linewidth = 6, colour = "#E9EDF2", lineend = "round")
      if (!is.null(sat) && nrow(sat) > 0)
        p <- p + geom_segment(data = sat, aes(x = ini / 1e6, xend = fin / 1e6, y = y, yend = y, colour = "SAT350"), linewidth = 4.8)
      if (!is.null(mic) && nrow(mic) > 0) {
        p <- p + geom_point(data = mic, aes(pos_mb, y, fill = "Microsatélites"), shape = 21, size = 2.5, colour = "white", stroke = 0.45)
        p <- p + if (requireNamespace("ggrepel", quietly = TRUE))
          ggrepel::geom_text_repel(data = mic, aes(x = pos_mb, y = y + 0.25, label = etiqueta), direction = "x", nudge_y = 0.12, seed = 42,
                                   size = 2.6, box.padding = 0.28, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.25, show.legend = FALSE)
        else geom_text(data = mic, aes(x = pos_mb, y = y + 0.4, label = etiqueta), size = 2.6, check_overlap = TRUE)
      }
      p + scale_colour_manual(values = c(SAT350 = col_sat350), name = NULL) +
        scale_fill_manual(values = c("Microsatélites" = "#187C8E"), name = NULL) +
        scale_x_continuous(name = "Posición cromosómica (Mb)", limits = c(0, xmax * 1.03), expand = c(0, 0)) +
        scale_y_continuous(breaks = pistas$y, labels = pistas$genoma, limits = c(min(pistas$y) - 0.75, max(pistas$y) + 1.8), expand = c(0, 0)) +
        labs(title = paste("Cromosoma", chr), subtitle = paste0(length(cmp_nombres), " genomas  |  longitud máxima: ", round(xmax, 2), " Mb")) +
        tema_art + theme(legend.position = "top", axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                         axis.text.y = element_text(face = "bold", size = 9), panel.grid.major.y = element_blank()) +
        coord_cartesian(clip = "off")
    }
    figs <- list()
    for (chr in niveles_chr) { p <- dibujar_chr(chr); if (!is.null(p)) figs[[chr]] <- p }
    if (length(figs) == 0) stop("Sin cromosomas que dibujar")
    alto <- 2.5 + 0.9 * length(cmp_nombres)
    dir.create(paste0(carpeta_cmp, "mapas_cromosoma"), showWarnings = FALSE)
    cairo_pdf(paste0(carpeta_cmp, "mapas_cromosoma.pdf"), width = 12, height = alto, onefile = TRUE)
    for (p in figs) print(p)
    dev.off()
    for (nm in names(figs)) ggsave(paste0(carpeta_cmp, "mapas_cromosoma/cromosoma_", nm, ".png"), figs[[nm]], width = 12, height = alto, dpi = 300, bg = "white")
    print(figs[[1]])
    leyenda("G1. Mapas por cromosoma", "Un panel por cromosoma (PDF multipágina + un PNG por cromosoma); una pista por genoma a escala física. Rojo: regiones SAT350 (BLASTN); puntos: microsatélites.")
  }, cond = !is.null(cmp_sat350) || !is.null(cmp_micro))
  
  # G2. Mapa de calor: kb de SAT350 por cromosoma y genoma
  intentar("G2. Mapa de calor SAT350", {
    s <- cmp_sat350; s$kb <- (s$fin - s$ini + 1) / 1e3
    a <- aggregate(kb ~ genoma + chrom, data = s, FUN = sum)
    hm <- merge(cmp_long[, c("genoma", "chrom")], a, by = c("genoma", "chrom"), all.x = TRUE); hm$kb[is.na(hm$kb)] <- 0
    hm$chrom <- factor(hm$chrom, levels = niveles_chr); hm$genoma <- factor(hm$genoma, levels = rev(cmp_nombres))
    hm$etq <- ifelse(hm$kb > 0, formatC(hm$kb, format = "f", digits = 1), "–")
    p <- ggplot(hm, aes(chrom, genoma, fill = ifelse(kb > 0, log10(kb), NA_real_))) +
      geom_tile(colour = "white", linewidth = 0.6) + geom_text(aes(label = etq), size = 2.4) +
      scale_fill_gradient(low = "#fddbc7", high = "#b2182b", na.value = "grey92", name = "log10(kb)") +
      labs(x = NULL, y = NULL, title = "SAT350 (kb) por cromosoma y genoma", subtitle = "Regiones BLASTN fusionadas. '–' = ausente") +
      tema_art + theme(legend.position = "right", axis.line = element_blank())
    print(p); guardar_cmp(p, "sat350_kb_cromosoma_genoma", 3 + 0.6 * length(niveles_chr), 1.8 + 0.5 * length(cmp_nombres))
    leyenda("G2. SAT350 por cromosoma", "kb de SAT350 (regiones BLASTN fusionadas) por cromosoma y genoma.")
  }, cond = !is.null(cmp_sat350))
  
  # G3. Resumen por genoma
  intentar("G3. Resumen por genoma", {
    r <- cmp_resumen
    dl <- rbind(
      data.frame(genoma = r$genoma, valor = r$n_regiones_SAT350, variable = "Regiones SAT350 (n)"),
      data.frame(genoma = r$genoma, valor = r$kb_SAT350, variable = "SAT350 (kb)"),
      data.frame(genoma = r$genoma, valor = r$pct_genoma_SAT350, variable = "SAT350 (% del genoma)"),
      data.frame(genoma = r$genoma, valor = r$SAT350_en_extremos, variable = paste0("SAT350 a <= ", w_st_max / 1e3, " kb del borde (n)")),
      data.frame(genoma = r$genoma, valor = r$n_microsat, variable = "Microsatélites ubicados (n)"),
      data.frame(genoma = r$genoma, valor = r$extremos_con_ST_valido, variable = "Extremos con ST válido (n)"),
      data.frame(genoma = r$genoma, valor = r$extremos_con_telomero, variable = "Extremos con telómero (n)"))
    dl$genoma <- factor(dl$genoma, levels = cmp_nombres); dl$variable <- factor(dl$variable, levels = unique(dl$variable))
    p <- ggplot(dl, aes(genoma, valor, fill = genoma)) + geom_col(width = 0.7, na.rm = TRUE) +
      facet_wrap(~ variable, scales = "free_y", nrow = 2) + scale_fill_manual(values = col_gen) +
      labs(x = NULL, y = NULL, title = "Comparación entre genomas") +
      tema_art + theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 7),
                       axis.text.x = element_text(angle = 40, hjust = 1))
    print(p); guardar_cmp(p, "resumen_genomas", 11, 6)
    leyenda("G3. Resumen por genoma", "SAT350, microsatélites, subtélomeros y telómeros por genoma (los dos últimos solo para genomas en modo completo).")
  })
  
  # G4. Subtélomeros y telómeros por extremo
  intentar("G4. Subtélomeros y telómeros por extremo", {
    ce <- cmp_extremos
    ce$chrom <- factor(ce$chrom, levels = niveles_chr); ce$genoma <- factor(ce$genoma, levels = cmp_nombres)
    ce$extremo <- factor(ce$extremo, levels = c("inicio", "final")); ce$st_kb <- ce$ST_long_pb / 1e3
    p1 <- ggplot(ce, aes(chrom, st_kb, colour = genoma, shape = ST_detectada)) +
      geom_point(position = position_dodge(width = 0.6), size = 2, na.rm = TRUE) + facet_grid(extremo ~ .) +
      scale_shape_manual(values = c(Si = 16, No = 1), name = "ST válido") + scale_colour_manual(values = col_gen, name = "Genoma") +
      labs(x = NULL, y = "Longitud del subtélomero (kb)", title = "Subtélomeros y telómeros por extremo y genoma") +
      tema_art + theme(legend.position = "right", strip.background = element_blank(), strip.text = element_text(face = "bold"),
                       axis.text.x = element_blank())
    p2 <- ggplot(ce, aes(chrom, genoma, fill = tiene_telomero)) + geom_tile(colour = "white") + facet_grid(extremo ~ .) +
      scale_fill_manual(values = c(Si = "#1D2D44", No = "#d73027"), name = "Telómero") + labs(x = NULL, y = NULL) +
      tema_art + theme(legend.position = "right", axis.line = element_blank(), strip.background = element_blank(),
                       strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 30, hjust = 1))
    p <- p1 / p2 + plot_layout(heights = c(1.2, 1))
    print(p); guardar_cmp(p, "subtelomeros_telomeros_genomas", 2 + 0.55 * length(niveles_chr), 4 + 0.5 * length(cmp_nombres))
    leyenda("G4. Subtélomeros y telómeros", "Arriba: longitud del subtélomero por extremo (círculo lleno = válido). Abajo: presencia de telómero (azul oscuro = sí, rojo = no).")
  }, cond = !is.null(cmp_extremos))
  
  
  # ============================================================
  # G5. Idiogramas F08 comparativos: una gráfica por cromosoma
  #     (una cápsula por genoma; estilo F08, formato tipo Nature)
  # REQUIERE (de 16_comparacion.R): cmp_long, cmp_extremos, cmp_sat_bloques,
  #   cmp_micro, cmp_nombres, niveles_chr, carpeta_cmp, agrupar_microsat()
  #   y de 15_figuras: intentar(), leyenda()
  # PAQUETES: ggplot2, ragg, ggrepel (opcional), scales
  # SALIDA : idiogramas_F08_comparativos.pdf (multipágina)
  #          idiogramas_F08_por_cromosoma/idiograma_<chr>.png / .tiff / .pdf
  # ============================================================
  intentar("G5. Idiogramas F08 comparativos por cromosoma", {
    
    # ---- Opciones ----
    solo_st_validos        <- FALSE   # TRUE: solo subtélomeros con ST_detectada == "Si"
    guardar_pdf_individual <- TRUE    # TRUE: además de PNG/TIFF, un PDF por cromosoma
    guardar_tiff           <- TRUE    # TRUE: TIFF con compresión LZW
    dpi_png                <- 600     # Nature: 600 dpi para line art
    ancho_fig              <- 7.2     # 183 mm (doble columna). Usa 3.5 para 89 mm
    alto_fig               <- 1.0 + 0.9 * length(cmp_nombres)
    alto_chr               <- 0.30    # semialtura de cada cápsula
    separacion             <- 2.7     # distancia vertical entre genomas
    
    # ---- Tipografía ----
    # Si Arial no existe: "Liberation Sans" o "Helvetica".
    # Ver disponibles: systemfonts::system_fonts()$family
    fuente    <- "Arial"
    pt_texto  <- 7                    # texto general (pt)
    pt_titulo <- 8                    # título (pt)
    pt_micro  <- 8.5                    # etiquetas de microsatélites (pt)
    pt2mm     <- function(pt) pt / ggplot2::.pt   # geoms usan mm
    
    # ---- Comprobaciones ----
    if (is.null(cmp_long) || nrow(cmp_long) == 0)
      stop("No hay longitudes cromosómicas para dibujar.")
    if (is.null(cmp_extremos) || !all(c("ST_start", "ST_end") %in% names(cmp_extremos)))
      stop("Faltan ST_start/ST_end en cmp_extremos. Actualiza 16_comparacion.R.")
    if (!requireNamespace("ragg", quietly = TRUE))
      stop("Instala ragg: install.packages('ragg')")
    
    micro_all <- if (exists("cmp_micro") && !is.null(cmp_micro) && nrow(cmp_micro) > 0)
      agrupar_microsat(cmp_micro) else NULL
    
    y_genoma <- setNames(rev(seq_along(cmp_nombres)) * separacion, cmp_nombres)
    
    # ---- Helpers ----
    # Polígono en forma de cápsula, recortado al cromosoma [0, L]
    capsula <- function(x0, x1, L, y, h, rx, id, n = 60) {
      r  <- min(rx, L / 2)
      x0 <- max(0, min(L, x0)); x1 <- max(x0, min(L, x1))
      xx <- seq(x0, x1, length.out = n)
      d  <- pmax(r - xx, xx - (L - r), 0)
      off <- h * sqrt(pmax(0, 1 - (d / r)^2))
      data.frame(x = c(xx, rev(xx)), y = c(y + off, rev(y - off)), id = id)
    }
    
    poligonos_bloque <- function(d, rx) {
      if (is.null(d) || nrow(d) == 0) return(NULL)
      out <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
        z <- capsula(d$x0[i], d$x1[i], d$L[i], d$y[i], alto_chr, rx, i)
        z$tipo <- as.character(d$tipo[i])
        z
      }))
      rownames(out) <- NULL
      out
    }
    
    # Añade L (Mb) e y (altura de la pista) según genoma y cromosoma
    anadir_pista <- function(d, base) {
      if (is.null(d) || nrow(d) == 0) return(NULL)
      k <- match(d$genoma, base$genoma)
      d$L <- base$L[k]
      d$y <- base$y[k]
      d[!is.na(d$y) & !is.na(d$L), , drop = FALSE]
    }
    
    # ---- Paleta (Okabe-Ito / Wong, apta para daltonismo) ----
    pal_sat_base <- c(sat240 = "#0072B2", sat350 = "#E69F00")   # azul / naranja
    extra_wong   <- c("#56B4E9", "#D55E00", "#F0E442", "#999999")
    col_st       <- "#CC79A7"   # subtelómero
    col_micro    <- "#009E73"   # microsatélite
    col_micro_tx <- "#006B4F"   # texto de microsatélite (más oscuro)
    
    dibujar_f08_chr <- function(chr) {
      base <- cmp_long[cmp_long$chrom == chr, , drop = FALSE]
      if (nrow(base) == 0) return(NULL)
      base$L <- base$long_pb / 1e6
      base$y <- unname(y_genoma[as.character(base$genoma)])
      base <- base[!is.na(base$y) & is.finite(base$L) & base$L > 0, , drop = FALSE]
      if (nrow(base) == 0) return(NULL)
      
      xmax <- max(base$L)
      rx   <- xmax * 0.014
      
      # Cuerpo de cada cromosoma
      cuerpo <- do.call(rbind, lapply(seq_len(nrow(base)), function(i) {
        z <- capsula(0, base$L[i], base$L[i], base$y[i], alto_chr, rx, i)
        z$genoma <- as.character(base$genoma[i])
        z
      }))
      
      # Subtelómeros
      st <- cmp_extremos[cmp_extremos$chrom == chr &
                           !is.na(cmp_extremos$ST_start) & !is.na(cmp_extremos$ST_end) &
                           cmp_extremos$ST_end > cmp_extremos$ST_start, , drop = FALSE]
      if (solo_st_validos && "ST_detectada" %in% names(st))
        st <- st[st$ST_detectada == "Si", , drop = FALSE]
      st <- anadir_pista(st, base)
      st_poly <- NULL
      if (!is.null(st) && nrow(st) > 0) {
        st$x0 <- st$ST_start / 1e6; st$x1 <- st$ST_end / 1e6; st$tipo <- "Subtelómero"
        st_poly <- poligonos_bloque(st, rx)
      }
      
      # Satélites
      sat <- NULL; sat_poly <- NULL
      if (exists("cmp_sat_bloques") && !is.null(cmp_sat_bloques) && nrow(cmp_sat_bloques) > 0) {
        sat <- anadir_pista(cmp_sat_bloques[cmp_sat_bloques$chrom == chr, , drop = FALSE], base)
        if (!is.null(sat) && nrow(sat) > 0) {
          sat$x0 <- sat$start / 1e6
          sat$x1 <- pmax(sat$end / 1e6, sat$x0 + xmax * 0.003)   # ancho mínimo visible
          sat$tipo <- as.character(sat$satelite)
          sat_poly <- poligonos_bloque(sat, rx)
        }
      }
      
      # Microsatélites
      micro <- NULL
      if (!is.null(micro_all) && nrow(micro_all) > 0) {
        micro <- anadir_pista(micro_all[micro_all$chrom == chr, , drop = FALSE], base)
        if (!is.null(micro) && nrow(micro) > 0) {
          micro$tipo  <- "Microsatélite"
          micro$y_ini <- micro$y - alto_chr
          micro$y_pt  <- micro$y - alto_chr - 0.32
          micro$y_lab <- micro$y_pt - 0.38
        } else micro <- NULL
      }
      
      # Paleta final
      tipos_sat <- if (!is.null(sat) && nrow(sat) > 0) unique(as.character(sat$tipo)) else character(0)
      pal_sat <- pal_sat_base
      extra <- setdiff(tipos_sat, names(pal_sat))
      if (length(extra) > 0)
        pal_sat <- c(pal_sat, setNames(extra_wong[(seq_along(extra) - 1) %% length(extra_wong) + 1], extra))
      cf <- c("Subtelómero" = col_st, pal_sat, "Microsatélite" = col_micro)
      presentes <- unique(c(if (!is.null(st_poly)) "Subtelómero", tipos_sat,
                            if (!is.null(micro)) "Microsatélite"))
      
      # ---- Gráfico ----
      p <- ggplot() +
        geom_polygon(data = cuerpo, aes(x, y, group = interaction(id, genoma)),
                     fill = "#F3F4F6", colour = NA)
      if (!is.null(st_poly))
        p <- p + geom_polygon(data = st_poly, aes(x, y, group = id, fill = tipo), colour = NA)
      if (!is.null(sat_poly))
        p <- p + geom_polygon(data = sat_poly, aes(x, y, group = id, fill = tipo), colour = NA)
      
      p <- p +
        geom_polygon(data = cuerpo, aes(x, y, group = interaction(id, genoma)),
                     fill = NA, colour = "#222222", linewidth = 0.25) +
        geom_text(data = base, aes(x = L + xmax * 0.012, y = y, label = sprintf("%.2f Mb", L)),
                  hjust = 0, size = pt2mm(pt_texto), family = fuente, colour = "#333333")
      
      if (!is.null(micro)) {
        p <- p +
          geom_segment(data = micro, aes(x = pos_mb, xend = pos_mb, y = y_ini, yend = y_pt),
                       colour = col_micro, linewidth = 0.2) +
          geom_point(data = micro, aes(x = pos_mb, y = y_pt, fill = tipo),
                     shape = 24, size = 1.6, stroke = 0.3, colour = "white")
        p <- p + if (requireNamespace("ggrepel", quietly = TRUE)) {
          ggrepel::geom_text_repel(data = micro, aes(x = pos_mb, y = y_lab, label = etiqueta),
                                   direction = "x", nudge_y = -0.10, size = pt2mm(pt_micro),
                                   family = fuente, colour = col_micro_tx, fontface = "bold.italic",
                                   min.segment.length = 0, segment.size = 0.2,
                                   segment.colour = "grey60", max.overlaps = Inf, seed = 42)
        } else {
          geom_text(data = micro, aes(x = pos_mb, y = y_lab, label = etiqueta),
                    size = pt2mm(pt_micro), family = fuente, colour = col_micro_tx,
                    fontface = "bold.italic", check_overlap = TRUE)
        }
      }
      
      p +
        scale_fill_manual(values = cf, breaks = names(cf)[names(cf) %in% presentes], name = NULL,
                          guide = guide_legend(nrow = 1,
                                               override.aes = list(shape = 22, size = 3.2,
                                                                   colour = "#222222", stroke = 0.2))) +
        scale_x_continuous(name = "Posición cromosómica (Mb)", limits = c(0, xmax * 1.14),
                           breaks = scales::pretty_breaks(n = 6), expand = c(0, 0)) +
        scale_y_continuous(breaks = base$y, labels = as.character(base$genoma),
                           limits = c(min(base$y) - 1.25, max(base$y) + 0.7), expand = c(0, 0)) +
        labs(title = paste("Cromosoma", chr), y = NULL) +
        theme_classic(base_size = pt_texto, base_family = fuente) +
        theme(
          plot.title           = element_text(face = "bold", size = pt_titulo, hjust = 0,
                                              margin = margin(b = 4)),
          panel.background     = element_rect(fill = "white", colour = NA),
          plot.background      = element_rect(fill = "white", colour = NA),
          legend.position      = "bottom",
          legend.justification = "left",
          legend.key.size      = grid::unit(3, "mm"),
          legend.key           = element_blank(),
          legend.text          = element_text(size = pt_texto, colour = "black"),
          legend.margin        = margin(0, 0, 2, 0),
          axis.title.x         = element_text(size = pt_texto, colour = "black", margin = margin(t = 4)),
          axis.text.x          = element_text(size = pt_texto, colour = "black"),
          axis.text.y          = element_text(face = "bold", size = pt_texto, colour = "black"),
          axis.line.x          = element_line(colour = "black", linewidth = 0.25),
          axis.ticks.x         = element_line(colour = "black", linewidth = 0.25),
          axis.ticks.length.x  = grid::unit(1.2, "mm"),
          axis.line.y          = element_blank(),
          axis.ticks.y         = element_blank(),
          panel.grid           = element_blank(),
          panel.border         = element_blank(),
          plot.margin          = margin(4, 10, 4, 4)
        ) +
        coord_cartesian(clip = "off")
    }
    
    # ---- Construir todas las figuras ----
    figs_f08 <- list()
    for (chr in niveles_chr) {
      p_chr <- tryCatch(dibujar_f08_chr(chr),
                        error = function(e) { message("  ", chr, ": ", conditionMessage(e)); NULL })
      if (!is.null(p_chr)) figs_f08[[chr]] <- p_chr
    }
    if (length(figs_f08) == 0) stop("No hubo cromosomas para graficar.")
    
    # ---- PDF multipágina vectorial (una página por cromosoma, fuente embebida) ----
    dir.create(carpeta_cmp, showWarnings = FALSE, recursive = TRUE)
    pdf_multi <- file.path(carpeta_cmp, "idiogramas_F08_comparativos.pdf")
    grDevices::cairo_pdf(pdf_multi, width = ancho_fig, height = alto_fig,
                         onefile = TRUE, family = fuente)
    for (p_chr in figs_f08) print(p_chr)
    grDevices::dev.off()
    
    # ---- Imágenes por separado ----
    dir_f08 <- file.path(carpeta_cmp, "idiogramas_F08_por_cromosoma")
    dir.create(dir_f08, showWarnings = FALSE, recursive = TRUE)
    for (chr in names(figs_f08)) {
      nom       <- gsub("[^A-Za-z0-9]+", "_", chr)
      base_ruta <- file.path(dir_f08, paste0("idiograma_", nom))
      
      # PNG (ragg: mejor antialiasing)
      ggsave(paste0(base_ruta, ".png"), figs_f08[[chr]], device = ragg::agg_png,
             width = ancho_fig, height = alto_fig, units = "in",
             dpi = dpi_png, bg = "white", limitsize = FALSE)
      
      # TIFF con compresión LZW (formato pedido por muchas revistas)
      if (guardar_tiff)
        ggsave(paste0(base_ruta, ".tiff"), figs_f08[[chr]], device = ragg::agg_tiff,
               width = ancho_fig, height = alto_fig, units = "in",
               dpi = dpi_png, compression = "lzw", bg = "white", limitsize = FALSE)
      
      # PDF vectorial
      if (guardar_pdf_individual)
        ggsave(paste0(base_ruta, ".pdf"), figs_f08[[chr]], device = grDevices::cairo_pdf,
               width = ancho_fig, height = alto_fig, bg = "white", limitsize = FALSE)
    }
    
    message("G5: ", length(figs_f08), " cromosomas | PDF: ", pdf_multi, " | imágenes: ", dir_f08)
    print(figs_f08[[1]])
    
    leyenda("G5. Idiogramas F08 comparativos",
            paste0("Un gráfico por cromosoma; una cápsula por genoma, con su longitud (Mb) a la derecha. ",
                   "Subtelómeros y familias de satélites se muestran como segmentos de color dentro del cromosoma; ",
                   "los triángulos y etiquetas bajo cada cápsula indican los microsatélites."))
  }, cond = exists("cmp_long") && !is.null(cmp_long) && nrow(cmp_long) > 0)
 
  
#
  
  cat("Leyendas de la comparación en:", leyendas_file, "\n")
}

if (fig_modo == "genoma") {
  cat("Leyendas de figuras escritas en:", leyendas_file, "\n")
  cat("\nFiguras guardadas en:", carpeta_salida, "\n")
} else cat("\nFiguras de comparacion guardadas en:", carpeta_cmp, "\n")
