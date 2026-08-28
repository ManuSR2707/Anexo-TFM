###############################################################################
# BUBBLE PLOT DE INTERACCIONES LIGANDO-RECEPTOR RECURRENTES
#
# Entrada:
#   Archivos BulkSignalR_*_significant.tsv en la carpeta de trabajo.
#
# Salida:
#   1. Bubble_plot_BulkSignalR_top_interactions.png
#   2. Bubble_plot_BulkSignalR_top_interactions.pdf
#   3. BulkSignalR_interaction_frequencies.tsv
#
# Idea:
#   Cada interacción L-R se cuenta una sola vez por comparación, aunque aparezca
#   repetida en varias filas por estar anotada en más de una ruta de Reactome.
###############################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

# =============================================================================
# 1. PAQUETES
# =============================================================================

required_packages <- c("dplyr", "ggplot2", "stringr", "readr", "forcats")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(readr)
  library(forcats)
})

# =============================================================================
# 2. CONFIGURACIÓN
# =============================================================================

# Número de interacciones que se mostrarán en la figura.
top_n_interactions <- 30

# Número mínimo de comparaciones en las que debe aparecer una interacción.
# Usa 1 para no excluir ninguna antes de seleccionar el Top N.
minimum_comparisons <- 1

# Carpeta que contiene los archivos.
input_directory <- "."

# Carpeta de salida.
output_directory <- "Bubble_plot_results"

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

# =============================================================================
# 3. LOCALIZAR ARCHIVOS DE BULKSIGNALR
# =============================================================================

bulk_files <- list.files(
  path = input_directory,
  pattern = "^BulkSignalR_.*_significant(?:\\([0-9]+\\))?\\.tsv$",
  full.names = TRUE
)

if (length(bulk_files) == 0) {
  stop(
    paste0(
      "No se encontraron archivos BulkSignalR significativos en:\n",
      normalizePath(input_directory),
      "\n\nPatrón buscado:\n",
      "^BulkSignalR_.*_significant(?:\\([0-9]+\\))?\\.tsv$"
    )
  )
}

message("Archivos detectados: ", length(bulk_files))

# =============================================================================
# 4. LEER Y RESUMIR CADA COMPARACIÓN
# =============================================================================

read_one_comparison <- function(file_name) {

  data <- read_tsv(
    file_name,
    show_col_types = FALSE,
    progress = FALSE
  )

  required_columns <- c("L", "R")
  missing_columns <- setdiff(required_columns, colnames(data))

  if (length(missing_columns) > 0) {
    stop(
      "Faltan columnas en ",
      basename(file_name),
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  comparison_name <- basename(file_name) |>
    str_remove("^BulkSignalR_") |>
    str_remove("_significant(?:\\([0-9]+\\))?\\.tsv$")

  # Una interacción se cuenta solo una vez por comparación.
  data |>
    select(L, R) |>
    filter(
      !is.na(L),
      !is.na(R),
      L != "",
      R != ""
    ) |>
    distinct(L, R) |>
    mutate(
      Comparison = comparison_name,
      Interaction = paste(L, R, sep = " \u2192 ")
    )
}

all_unique_interactions <- lapply(
  bulk_files,
  read_one_comparison
) |>
  bind_rows()

# =============================================================================
# 5. CALCULAR FRECUENCIA ENTRE COMPARACIONES
# =============================================================================

number_of_comparisons <- n_distinct(
  all_unique_interactions$Comparison
)

interaction_frequency <- all_unique_interactions |>
  group_by(L, R, Interaction) |>
  summarise(
    Number_comparisons = n_distinct(Comparison),
    Percentage_comparisons =
      100 * Number_comparisons / number_of_comparisons,
    Comparisons = paste(
      sort(unique(Comparison)),
      collapse = "; "
    ),
    .groups = "drop"
  ) |>
  arrange(
    desc(Number_comparisons),
    Interaction
  )

write_tsv(
  interaction_frequency,
  file.path(
    output_directory,
    "BulkSignalR_interaction_frequencies.tsv"
  )
)

# =============================================================================
# 6. SELECCIONAR INTERACCIONES PARA LA FIGURA
# =============================================================================

plot_data <- interaction_frequency |>
  filter(
    Number_comparisons >= minimum_comparisons
  ) |>
  slice_head(
    n = top_n_interactions
  ) |>
  mutate(
    Interaction = fct_reorder(
      Interaction,
      Number_comparisons
    )
  )

if (nrow(plot_data) == 0) {
  stop(
    "Ninguna interacción cumple minimum_comparisons = ",
    minimum_comparisons
  )
}

# =============================================================================
# 7. GENERAR BUBBLE PLOT
# =============================================================================

bubble_plot <- ggplot(
  plot_data,
  aes(
    x = Number_comparisons,
    y = Interaction,
    size = Number_comparisons,
    color = Percentage_comparisons
  )
) +
  geom_point(
    alpha = 0.88
  ) +
  geom_text(
    aes(
      label = Number_comparisons
    ),
    color = "black",
    size = 3.3,
    fontface = "bold"
  ) +
  scale_size_continuous(
    range = c(5, 14),
    name = "N.º de comparaciones"
  ) +
  scale_color_gradient(
    low = "#90CAF9",
    high = "#C62828",
    name = "% de comparaciones"
  ) +
  scale_x_continuous(
    breaks = seq(
      0,
      number_of_comparisons,
      by = 2
    ),
    limits = c(
      0,
      max(plot_data$Number_comparisons) + 1
    ),
    expand = expansion(
      mult = c(0.01, 0.05)
    )
  ) +
  labs(
    title = "Interacciones ligando-receptor más recurrentes",
    subtitle = paste0(
      "Frecuencia calculada en ",
      number_of_comparisons,
      " comparaciones de BulkSignalR"
    ),
    x = "Número de comparaciones en las que aparece la interacción",
    y = NULL,
    caption = paste0(
      "Cada interacción se contabilizó una sola vez por comparación. ",
      "Se muestran las ",
      nrow(plot_data),
      " interacciones más frecuentes."
    )
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 17
    ),
    plot.subtitle = element_text(
      color = "#4B5563",
      margin = margin(
        b = 14
      )
    ),
    axis.text.y = element_text(
      size = 9
    ),
    panel.grid.major.y = element_line(
      color = "#E5E7EB"
    ),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.caption = element_text(
      hjust = 0,
      color = "#6B7280",
      size = 9
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
    "Bubble_plot_BulkSignalR_top_interactions.png"
  ),
  plot = bubble_plot,
  width = 12,
  height = 11,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_directory,
    "Bubble_plot_BulkSignalR_top_interactions.pdf"
  ),
  plot = bubble_plot,
  width = 12,
  height = 11,
  units = "in",
  bg = "white"
)

cat(
  "\n============================================================\n",
  "BUBBLE PLOT GENERADO CORRECTAMENTE\n",
  "Comparaciones analizadas: ",
  number_of_comparisons,
  "\nResultados guardados en:\n",
  normalizePath(
    output_directory,
    mustWork = FALSE
  ),
  "\n============================================================\n",
  sep = ""
)
