# ============================================================
# FIGURAS A1-A4 DEL ANEXO - VERSION SIN PATCHWORK
# Hallmark ssGSEA: M. bovis vs Lineage 5
#
# Esta versión evita por completo patchwork y wrap_dims().
#
# Archivos de entrada:
#   All_limma_results.tsv
#   All_Hallmark_ssGSEA_scores_long.tsv
#
# Paquetes:
#   tidyverse, grid, gridExtra
# ============================================================

library(tidyverse)
library(grid)
library(gridExtra)

# -----------------------------
# 1. ARCHIVOS
# -----------------------------
limma_file  <- "All_limma_results.tsv"
scores_file <- "All_Hallmark_ssGSEA_scores_long.tsv"

outdir <- "Anexo_A1_A4"
dir.create(outdir, showWarnings = FALSE)

limma_res <- read.delim(limma_file, check.names = FALSE)
scores    <- read.delim(scores_file, check.names = FALSE)

# -----------------------------
# 2. COLORES
# -----------------------------
col_lineage5 <- "#4C78A8"
col_mbovis   <- "#E45756"

condition_cols <- c(
  "Lineage_5" = col_lineage5,
  "M_bovis"   = col_mbovis
)

heat_cols <- c(
  "#4575B4",
  "#ABD9E9",
  "#FFFFBF",
  "#FDAE61",
  "#D73027"
)

# -----------------------------
# 3. FUNCIONES AUXILIARES
# -----------------------------
pretty_pathway <- function(x) {
  x %>%
    str_remove("^HALLMARK_") %>%
    str_replace_all("_", " ") %>%
    str_to_title()
}

priority <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_HYPOXIA",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_APOPTOSIS"
)

# -----------------------------
# 4. PANEL A: HEATMAP + TIRA DE CONDICION
# -----------------------------
make_heatmap_panel <- function(sc, heat_paths, model, day) {

  sample_info <- sc %>%
    distinct(Sample, Condition) %>%
    mutate(
      Condition = factor(
        Condition,
        levels = c("Lineage_5", "M_bovis")
      )
    ) %>%
    arrange(Condition, Sample)

  sample_order <- sample_info$Sample

  hm_wide <- sc %>%
    filter(Pathway %in% heat_paths) %>%
    select(Pathway, Sample, Pathway_score) %>%
    pivot_wider(
      names_from = Sample,
      values_from = Pathway_score
    )

  hm_mat <- hm_wide %>%
    column_to_rownames("Pathway") %>%
    as.matrix()

  hm_mat <- hm_mat[heat_paths, sample_order, drop = FALSE]

  # Z-score por vía
  hm_z <- t(scale(t(hm_mat)))
  hm_z[is.na(hm_z)] <- 0

  pathway_levels <- rev(heat_paths)
  pathway_labels <- pretty_pathway(pathway_levels)

  heat_df <- as.data.frame(hm_z) %>%
    rownames_to_column("Pathway") %>%
    pivot_longer(
      -Pathway,
      names_to = "Sample",
      values_to = "Z"
    ) %>%
    mutate(
      Sample = factor(Sample, levels = sample_order),
      Pathway = factor(Pathway, levels = pathway_levels)
    )

  # Heatmap principal
  p_heat <- ggplot(
    heat_df,
    aes(x = Sample, y = Pathway, fill = Z)
  ) +
    geom_tile(
      colour = "grey55",
      linewidth = 0.25
    ) +
    scale_fill_gradientn(
      colours = heat_cols,
      limits = c(-1.5, 1.5),
      oob = scales::squish,
      name = "Z-score",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        direction = "horizontal",
        barwidth = unit(5.0, "cm"),
        barheight = unit(0.35, "cm")
      )
    ) +
    scale_y_discrete(
      position = "right",
      breaks = pathway_levels,
      labels = pathway_labels
    ) +
    labs(
      title = paste0(
        model, " ",
        str_replace(day, "Day", "Day "),
        " : M_bovis vs Lineage_5"
      ),
      x = NULL,
      y = NULL,
      tag = "A"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 8
      ),
      axis.text.y = element_text(
        size = 9,
        margin = margin(l = 8)
      ),
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 13,
        margin = margin(b = 8)
      ),
      plot.tag = element_text(
        size = 16,
        face = "plain"
      ),
      legend.position = "bottom",
      plot.margin = margin(
        t = 5, r = 25, b = 5, l = 5
      )
    )

  # Tira superior de condición
  annotation_df <- sample_info %>%
    mutate(
      Sample = factor(Sample, levels = sample_order),
      y = 1
    )

  p_ann <- ggplot(
    annotation_df,
    aes(x = Sample, y = y, fill = Condition)
  ) +
    geom_tile() +
    scale_fill_manual(
      values = condition_cols,
      name = "Condition"
    ) +
    scale_x_discrete(drop = FALSE) +
    guides(
      fill = guide_legend(
        title.position = "top"
      )
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      plot.margin = margin(
        t = 0, r = 25, b = 0, l = 5
      )
    )

  # Convertir a grobs y apilar SIN patchwork
  ann_grob  <- ggplotGrob(p_ann)
  heat_grob <- ggplotGrob(p_heat)

  arrangeGrob(
    ann_grob,
    heat_grob,
    ncol = 1,
    heights = c(0.16, 1)
  )
}

# -----------------------------
# 5. BOXPLOTS B-F
# -----------------------------
make_boxplot <- function(sc, stat_row, pathway, letter) {

  plot_df <- sc %>%
    filter(Pathway == pathway) %>%
    mutate(
      Condition = factor(
        Condition,
        levels = c("Lineage_5", "M_bovis")
      )
    )

  ggplot(
    plot_df,
    aes(
      x = Condition,
      y = Pathway_score,
      fill = Condition
    )
  ) +
    geom_boxplot(
      width = 0.55,
      outlier.shape = NA,
      alpha = 0.70,
      colour = "grey30"
    ) +
    geom_jitter(
      width = 0.05,
      height = 0,
      size = 2.3,
      colour = "black"
    ) +
    scale_fill_manual(
      values = condition_cols
    ) +
    labs(
      title = pretty_pathway(pathway),
      subtitle = paste0(
        "logFC = ",
        sprintf("%.3f", stat_row$logFC),
        "; FDR = ",
        format(
          stat_row$adj.P.Val,
          scientific = TRUE,
          digits = 3
        )
      ),
      x = NULL,
      y = "ssGSEA score",
      tag = letter
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 11
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 9,
        colour = "#607080"
      ),
      plot.tag = element_text(
        size = 16,
        face = "plain"
      ),
      axis.text.x = element_text(
        angle = 20,
        hjust = 1
      ),
      panel.grid.minor = element_blank()
    )
}

# -----------------------------
# 6. GENERAR UNA FIGURA A1-A4
# -----------------------------
make_annex_figure <- function(model, day, code) {

  lr <- limma_res %>%
    filter(
      Model == model,
      Day == day,
      Contrast == "M_bovis_vs_Lineage_5"
    )

  sig <- lr %>%
    filter(adj.P.Val < 0.05) %>%
    mutate(abs_logFC = abs(logFC)) %>%
    arrange(adj.P.Val, desc(abs_logFC))

  if (nrow(sig) == 0) {
    stop(
      paste(
        "No hay vías significativas para",
        model, day
      )
    )
  }

  # 12 vías para el heatmap
  heat_paths <- head(sig$Pathway, 12)

  # 5 vías para B-F
  selected <- c(
    priority[priority %in% sig$Pathway],
    sig$Pathway[!sig$Pathway %in% priority]
  )

  selected_box <- head(unique(selected), 5)

  if (length(selected_box) < 5) {
    stop(
      paste(
        "Hay menos de 5 vías significativas para",
        model, day
      )
    )
  }

  sc <- scores %>%
    filter(
      Model == model,
      Day == day,
      Condition %in% c(
        "Lineage_5",
        "M_bovis"
      )
    )

  # Panel A
  panel_A <- make_heatmap_panel(
    sc = sc,
    heat_paths = heat_paths,
    model = model,
    day = day
  )

  # Paneles B-F
  letters <- c("B", "C", "D", "E", "F")
  box_grobs <- vector("list", 5)

  for (i in seq_along(selected_box)) {

    pathway <- selected_box[i]

    stat_row <- lr %>%
      filter(Pathway == pathway) %>%
      slice(1)

    p <- make_boxplot(
      sc = sc,
      stat_row = stat_row,
      pathway = pathway,
      letter = letters[i]
    )

    box_grobs[[i]] <- ggplotGrob(p)
  }

  # Filas B-C, D-E, F-vacío
  row_BC <- arrangeGrob(
    box_grobs[[1]],
    box_grobs[[2]],
    ncol = 2
  )

  row_DE <- arrangeGrob(
    box_grobs[[3]],
    box_grobs[[4]],
    ncol = 2
  )

  row_F <- arrangeGrob(
    box_grobs[[5]],
    nullGrob(),
    ncol = 2
  )

  # Cuerpo completo
  body <- arrangeGrob(
    panel_A,
    row_BC,
    row_DE,
    row_F,
    ncol = 1,
    heights = c(1.45, 1, 1, 1)
  )

  # Título general
  title_grob <- textGrob(
    paste0(
      "Figure ", code,
      ". Hallmark pathway activity in ",
      model,
      " at ",
      day,
      " (M_bovis vs Lineage_5)"
    ),
    gp = gpar(
      fontsize = 16,
      fontface = "bold"
    )
  )

  final_grob <- arrangeGrob(
    title_grob,
    body,
    ncol = 1,
    heights = c(0.07, 1)
  )

  # -------------------------
  # GUARDAR
  # -------------------------
  output_file <- file.path(
    outdir,
    paste0(
      "Figura_",
      code, "_",
      str_replace_all(model, "-", ""),
      "_",
      day,
      ".png"
    )
  )

  png(
    filename = output_file,
    width = 12,
    height = 16,
    units = "in",
    res = 300,
    bg = "white"
  )

  grid.newpage()
  grid.draw(final_grob)

  dev.off()

  message("Generada: ", output_file)
}

# -----------------------------
# 7. GENERAR A1-A4
# -----------------------------
make_annex_figure("THP-1", "Day1", "A1")
make_annex_figure("THP-1", "Day3", "A2")
make_annex_figure("BoMac", "Day1", "A3")
make_annex_figure("BoMac", "Day3", "A4")

message("Figuras A1-A4 completadas.")
