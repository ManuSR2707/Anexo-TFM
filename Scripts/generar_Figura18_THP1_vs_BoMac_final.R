###############################################################################
# FIGURA 18 FINAL:
# COMPARACIÓN FUNCIONAL ENTRE LOS MODELOS THP-1 Y BoMac
###############################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

required_packages <- c("dplyr", "tidyr", "ggplot2", "readr", "tibble")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

results_directory <- "Pathway_analysis_results"

limma_file <- file.path(
  results_directory,
  "04_Statistical_tables",
  "All_limma_results.tsv"
)

output_directory <- file.path(
  results_directory,
  "09_Model_comparison"
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(limma_file)) {
  stop(
    paste0(
      "\nNo se encuentra el archivo:\n",
      limma_file,
      "\n\nCarpeta de trabajo actual:\n",
      getwd()
    )
  )
}

limma_results <- read_tsv(
  limma_file,
  show_col_types = FALSE,
  progress = FALSE
)

selected_pathways <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_HYPOXIA",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "HALLMARK_ANGIOGENESIS"
)

pathway_labels <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE" = "Interferón α",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE" = "Interferón γ",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB" = "TNFα/NFκB",
  "HALLMARK_INFLAMMATORY_RESPONSE" = "Respuesta inflamatoria",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING" = "IL6/JAK/STAT3",
  "HALLMARK_HYPOXIA" = "Hipoxia",
  "HALLMARK_GLYCOLYSIS" = "Glucólisis",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS" = "Homeostasis del colesterol",
  "HALLMARK_ANGIOGENESIS" = "Angiogénesis"
)

pathway_order <- unname(pathway_labels[selected_pathways])

plot_data <- limma_results |>
  filter(
    Contrast == "M_bovis_vs_Control",
    Model %in% c("THP-1", "BoMac"),
    Day %in% c("Day1", "Day3"),
    Pathway %in% selected_pathways
  ) |>
  mutate(
    Day_label = recode(
      Day,
      Day1 = "Día 1",
      Day3 = "Día 3"
    ),
    Pathway_label = unname(pathway_labels[Pathway]),
    Significance_symbol = case_when(
      adj.P.Val < 0.001 ~ "***",
      adj.P.Val < 0.01  ~ "**",
      adj.P.Val < 0.05  ~ "*",
      TRUE ~ ""
    ),
    Cell_label = paste0(
      sprintf("%.2f", logFC),
      ifelse(
        Significance_symbol == "",
        "",
        paste0("\n", Significance_symbol)
      )
    ),
    Model = factor(
      Model,
      levels = c("THP-1", "BoMac")
    ),
    Day_label = factor(
      Day_label,
      levels = c("Día 1", "Día 3")
    ),
    Pathway_label = factor(
      Pathway_label,
      levels = rev(pathway_order)
    )
  )

write_tsv(
  plot_data |>
    select(
      Model,
      Day,
      Contrast,
      Pathway,
      Pathway_label,
      Day_label,
      logFC,
      adj.P.Val,
      Significance_symbol
    ),
  file.path(
    output_directory,
    "Figure18_data_final.tsv"
  )
)

maximum_absolute_logFC <- max(
  abs(plot_data$logFC),
  na.rm = TRUE
)

figure_18 <- ggplot(
  plot_data,
  aes(
    x = Day_label,
    y = Pathway_label,
    fill = logFC
  )
) +
  geom_tile(
    color = "white",
    linewidth = 1.1
  ) +
  geom_text(
    aes(label = Cell_label),
    size = 3.1,
    lineheight = 0.82,
    color = "black"
  ) +
  facet_grid(
    cols = vars(Model),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_gradient2(
    low = "#4C72B0",
    mid = "white",
    high = "#C44E52",
    midpoint = 0,
    limits = c(
      -maximum_absolute_logFC,
      maximum_absolute_logFC
    ),
    name = NULL
  ) +
  labs(
    title = "Comparación funcional entre los modelos THP-1 y BoMac",
    subtitle = "* FDR < 0,05; ** FDR < 0,01; *** FDR < 0,001",
    x = NULL,
    y = NULL,
    caption = paste0(
      "Valores positivos indican mayor actividad en M. bovis; ",
      "valores negativos indican mayor actividad en el control."
    )
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 9.5,
      color = "#4B5563",
      margin = margin(b = 10)
    ),
    strip.text = element_text(
      face = "bold",
      size = 11,
      color = "black"
    ),
    strip.background = element_blank(),
    axis.text.x = element_text(
      face = "bold",
      size = 10,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 10,
      color = "black"
    ),
    panel.spacing.x = unit(0.35, "cm"),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.caption = element_text(
      hjust = 0,
      size = 8.8,
      color = "#5B6670",
      margin = margin(t = 10)
    ),
    plot.margin = margin(
      15,
      20,
      15,
      15
    )
  )

ggsave(
  filename = file.path(
    output_directory,
    "Figure18_THP1_vs_BoMac_Mbovis_vs_Control_final.png"
  ),
  plot = figure_18,
  width = 9.5,
  height = 6.8,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_directory,
    "Figure18_THP1_vs_BoMac_Mbovis_vs_Control_final.pdf"
  ),
  plot = figure_18,
  width = 9.5,
  height = 6.8,
  units = "in",
  bg = "white"
)

cat(
  "\n============================================================\n",
  "FIGURA 18 FINAL GENERADA CORRECTAMENTE\n",
  "Resultados guardados en:\n",
  normalizePath(
    output_directory,
    mustWork = FALSE
  ),
  "\n============================================================\n",
  sep = ""
)
