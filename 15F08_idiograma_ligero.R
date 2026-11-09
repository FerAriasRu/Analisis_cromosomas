# F08 individual para genomas en modo ligero.
# Muestra cromosomas, regiones SAT350 y microsatélites con primers.
# En modo ligero no se calculan subtelómeros.

chr_df <- data.frame(
  chrom = as.character(eq$nombre),
  L = as.numeric(sl[eq$acc]) / 1e6,
  stringsAsFactors = FALSE
)
chr_df$y <- rev(seq_len(nrow(chr_df))) * 2.2
chr_df$lab_L <- sprintf("%.2f Mb", chr_df$L)

xmax <- max(chr_df$L, na.rm = TRUE)
alto_chr <- 0.30
rx <- xmax * 0.014

capsula <- function(x0, x1, L, y, h, rx, id, n = 60) {
  r <- min(rx, L / 2)
  x0 <- max(0, min(L, x0))
  x1 <- max(x0, min(L, x1))
  xx <- seq(x0, x1, length.out = n)
  d <- pmax(r - xx, xx - (L - r), 0)
  off <- h * sqrt(pmax(0, 1 - (d / r)^2))
  data.frame(
    x = c(xx, rev(xx)),
    y = c(y + off, rev(y - off)),
    id = id
  )
}

cuerpo <- do.call(rbind, lapply(seq_len(nrow(chr_df)), function(i) {
  p <- capsula(0, chr_df$L[i], chr_df$L[i],
               chr_df$y[i], alto_chr, rx, i)
  p$chrom <- chr_df$chrom[i]
  p
}))

# Regiones de satélites detectadas, si existen
sat_d <- NULL

if (exists("sat_bloques") && !is.null(sat_bloques) &&
    nrow(sat_bloques) > 0) {
  sat_d <- data.frame(
    chrom = as.character(sat_bloques$chrom),
    x0 = sat_bloques$start / 1e6,
    x1 = sat_bloques$end / 1e6,
    tipo = as.character(sat_bloques$satelite),
    stringsAsFactors = FALSE
  )
}

# Añade las regiones SAT350 de GenBank/BLAST del módulo 06b
if (exists("sat350_regiones") && !is.null(sat350_regiones) &&
    nrow(sat350_regiones) > 0) {
  sat350_d <- data.frame(
    chrom = as.character(sat350_regiones$chrom),
    x0 = sat350_regiones$start / 1e6,
    x1 = sat350_regiones$end / 1e6,
    tipo = "sat350",
    stringsAsFactors = FALSE
  )
  sat_d <- if (is.null(sat_d)) sat350_d else rbind(sat_d, sat350_d)
}

sat_poly <- NULL
if (!is.null(sat_d) && nrow(sat_d) > 0) {
  sat_d$y <- chr_df$y[match(sat_d$chrom, chr_df$chrom)]
  sat_d$L <- chr_df$L[match(sat_d$chrom, chr_df$chrom)]
  sat_d <- sat_d[!is.na(sat_d$y) & !is.na(sat_d$L), ]
  
  if (nrow(sat_d) > 0) {
    sat_d$x1 <- pmax(sat_d$x1, sat_d$x0 + xmax * 0.003)
    sat_poly <- do.call(rbind, lapply(seq_len(nrow(sat_d)), function(i) {
      p <- capsula(sat_d$x0[i], sat_d$x1[i], sat_d$L[i],
                   sat_d$y[i], alto_chr, rx, i)
      p$tipo <- sat_d$tipo[i]
      p
    }))
  }
}

# Microsatélites encontrados con primers
micro_d <- NULL
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
    stringsAsFactors = FALSE
  )
  micro_d$y_chr <- chr_df$y[match(micro_d$chrom, chr_df$chrom)]
  micro_d <- micro_d[!is.na(micro_d$y_chr), ]
  micro_d$y_ini <- micro_d$y_chr - alto_chr
  micro_d$y_pt <- micro_d$y_chr - alto_chr - 0.32
  micro_d$y_lab <- micro_d$y_pt - 0.38
}

pal_sat <- c(sat240 = "#7B5EA7", sat350 = "#E08E45")
tipos_sat <- if (!is.null(sat_d) && nrow(sat_d) > 0) {
  unique(sat_d$tipo)
} else character(0)

extras <- setdiff(tipos_sat, names(pal_sat))
if (length(extras) > 0) {
  pal_sat <- c(
    pal_sat,
    setNames(grDevices::hcl.colors(length(extras), "Dark 3"), extras)
  )
}

colores <- c(
  pal_sat,
  "Microsatélite" = "#0F7B8A"
)

p_f08 <- ggplot2::ggplot() +
  ggplot2::geom_polygon(
    data = cuerpo,
    ggplot2::aes(x = x, y = y, group = interaction(id, chrom)),
    fill = "#EEF1F4", colour = NA
  )

if (!is.null(sat_poly) && nrow(sat_poly) > 0) {
  p_f08 <- p_f08 +
    ggplot2::geom_polygon(
      data = sat_poly,
      ggplot2::aes(x = x, y = y, group = id, fill = tipo),
      colour = NA
    )
}

p_f08 <- p_f08 +
  ggplot2::geom_polygon(
    data = cuerpo,
    ggplot2::aes(x = x, y = y, group = interaction(id, chrom)),
    fill = NA, colour = "#4A5568", linewidth = 0.3
  ) +
  ggplot2::geom_text(
    data = chr_df,
    ggplot2::aes(x = L + xmax * 0.012, y = y, label = lab_L),
    hjust = 0, size = 2.5, colour = "#6B7280"
  )

if (!is.null(micro_d) && nrow(micro_d) > 0) {
  p_f08 <- p_f08 +
    ggplot2::geom_segment(
      data = micro_d,
      ggplot2::aes(x = x, xend = x, y = y_ini, yend = y_pt),
      colour = "#0F7B8A", linewidth = 0.25
    ) +
    ggplot2::geom_point(
      data = micro_d,
      ggplot2::aes(x = x, y = y_pt),
      shape = 24, size = 1.9, fill = "#0F7B8A",
      stroke = 0.25, colour = "white"
    )
  
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_f08 <- p_f08 +
      ggrepel::geom_text_repel(
        data = micro_d,
        ggplot2::aes(x = x, y = y_lab, label = marcador),
        direction = "x", size = 2.3,
        colour = "#14535D", fontface = "italic",
        max.overlaps = Inf, seed = 42
      )
  } else {
    p_f08 <- p_f08 +
      ggplot2::geom_text(
        data = micro_d,
        ggplot2::aes(x = x, y = y_lab, label = marcador),
        size = 2.3, colour = "#14535D",
        fontface = "italic", check_overlap = TRUE
      )
  }
}

presentes <- unique(c(
  tipos_sat,
  if (!is.null(micro_d) && nrow(micro_d) > 0) "Microsatélite"
))

p_f08 <- p_f08 +
  ggplot2::scale_fill_manual(
    values = colores,
    breaks = names(colores)[names(colores) %in% presentes],
    name = NULL
  ) +
  ggplot2::scale_x_continuous(
    name = "Posición cromosómica (Mb)",
    limits = c(0, xmax * 1.14),
    breaks = scales::pretty_breaks(n = 6),
    expand = c(0, 0)
  ) +
  ggplot2::scale_y_continuous(
    breaks = chr_df$y,
    labels = chr_df$chrom,
    limits = c(min(chr_df$y) - 1.25, max(chr_df$y) + 0.7),
    expand = c(0, 0)
  ) +
  ggplot2::labs(
    title = paste("Idiograma de", genoma_id),
    subtitle = "Cápsulas: cromosomas | bloques: satélites | triángulos: microsatélites",
    y = NULL
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.line.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
  ) +
  ggplot2::coord_cartesian(clip = "off")

print(p_f08)

ggplot2::ggsave(
  file.path(carpeta_salida, paste0("F08_idiograma_", ventana_label, ".png")),
  p_f08, width = 7.2, height = 1.4 + 0.55 * nrow(chr_df),
  dpi = 300, bg = "white", limitsize = FALSE
)
ggplot2::ggsave(
  file.path(carpeta_salida, paste0("F08_idiograma_", ventana_label, ".pdf")),
  p_f08, width = 7.2, height = 1.4 + 0.55 * nrow(chr_df),
  device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE
)