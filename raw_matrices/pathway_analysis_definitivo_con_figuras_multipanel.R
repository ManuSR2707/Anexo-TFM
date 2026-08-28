###############################################################################
# PATHWAY ANALYSIS DEFINITIVO
# RNA-seq THP-1 (humano) y BoMac (bovino)
#
# Archivos esperados en la carpeta de trabajo:
#   humanRawCounts.txt
#   bosRawCounts.txt
#   metadata_rnaseq.tsv
#
# El script:
#   1) comprueba y carga los archivos;
#   2) separa THP-1 y BoMac;
#   3) filtra y normaliza cada modelo por separado con edgeR;
#   4) obtiene los 50 pathways Hallmark mediante msigdbr;
#   5) calcula scores ssGSEA por muestra;
#   6) compara las condiciones dentro de cada día mediante Wilcoxon y limma;
#   7) genera tablas, boxplots, heatmaps, figuras multipanel A-F y un Excel resumen.
###############################################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(1234)

# =============================================================================
# 1. PAQUETES
# =============================================================================

cran_packages <- c(
  "dplyr", "tidyr", "tibble", "ggplot2",
  "msigdbr", "pheatmap", "writexl", "patchwork"
)

bioconductor_packages <- c(
  "edgeR", "limma", "GSVA"
)

install_packages_if_needed <- function(packages, bioconductor = FALSE) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) == 0) {
    return(invisible(NULL))
  }

  if (bioconductor) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }

    BiocManager::install(
      missing_packages,
      ask = FALSE,
      update = FALSE
    )
  } else {
    install.packages(
      missing_packages,
      dependencies = TRUE
    )
  }
}

install_packages_if_needed(cran_packages)
install_packages_if_needed(
  bioconductor_packages,
  bioconductor = TRUE
)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(msigdbr)
  library(pheatmap)
  library(writexl)
  library(patchwork)
  library(edgeR)
  library(limma)
  library(GSVA)
})

# =============================================================================
# 2. CONFIGURACIÓN
# =============================================================================

human_counts_file <- "humanRawCounts.txt"
bovine_counts_file <- "bosRawCounts.txt"
metadata_file <- "metadata_rnaseq.tsv"

results_directory <- "Pathway_analysis_results"

minimum_genes_per_pathway <- 10
maximum_genes_per_pathway <- 500
number_pathways_heatmap <- 20
number_pathways_boxplot <- 6
fdr_threshold <- 0.05

result_subdirectories <- c(
  "01_QC",
  "02_Normalized_expression",
  "03_Pathway_scores",
  "04_Statistical_tables",
  "05_Boxplots",
  "06_Heatmaps",
  "07_Individual_boxplots",
  "08_Publication_figures"
)

dir.create(
  results_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

for (directory_name in result_subdirectories) {
  dir.create(
    file.path(results_directory, directory_name),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

required_files <- c(
  human_counts_file,
  bovine_counts_file,
  metadata_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "\nNo se encuentran los siguientes archivos:\n",
      paste(missing_files, collapse = "\n"),
      "\n\nCarpeta de trabajo actual:\n",
      getwd(),
      "\n\nArchivos disponibles:\n",
      paste(list.files(), collapse = "\n")
    )
  )
}

# =============================================================================
# 3. FUNCIONES AUXILIARES
# =============================================================================

safe_filename <- function(text) {
  gsub(
    pattern = "[^A-Za-z0-9_-]+",
    replacement = "_",
    x = text
  )
}

write_tsv <- function(data, filename, row_names = FALSE) {
  write.table(
    data,
    file = filename,
    sep = "\t",
    quote = FALSE,
    row.names = row_names,
    col.names = if (row_names) NA else TRUE
  )
}

read_raw_counts <- function(filename) {
  message("Leyendo raw counts: ", filename)

  # Los raw counts proporcionados están separados por uno o más espacios.
  # La cabecera no tiene nombre para la primera columna de identificadores,
  # por lo que row.names = 1 es la forma correcta de leerlos.
  count_table <- read.table(
    file = filename,
    header = TRUE,
    sep = "",
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    fill = FALSE
  )

  if (ncol(count_table) < 1) {
    stop("No se detectaron columnas de muestras en ", filename)
  }

  count_matrix <- as.matrix(count_table)

  suppressWarnings(
    storage.mode(count_matrix) <- "numeric"
  )

  if (anyNA(count_matrix)) {
    stop(
      "Se detectaron valores no numéricos o NA en ",
      filename,
      "."
    )
  }

  if (any(count_matrix < 0)) {
    stop(
      "Se detectaron valores negativos en ",
      filename,
      ". El archivo debe contener raw counts."
    )
  }

  if (any(abs(count_matrix - round(count_matrix)) > 1e-8)) {
    warning(
      filename,
      " contiene valores no enteros. Se continuará, pero conviene confirmar ",
      "que realmente son raw counts."
    )
  }

  gene_ids <- sub(
    pattern = "\\.[0-9]+$",
    replacement = "",
    x = rownames(count_matrix)
  )

  rownames(count_matrix) <- gene_ids

  if (anyDuplicated(rownames(count_matrix))) {
    count_matrix <- rowsum(
      count_matrix,
      group = rownames(count_matrix),
      reorder = FALSE
    )
  }

  if (anyDuplicated(colnames(count_matrix))) {
    stop("Hay nombres de muestras duplicados en ", filename)
  }

  message(
    "  Cargados ",
    nrow(count_matrix),
    " genes y ",
    ncol(count_matrix),
    " muestras."
  )

  count_matrix
}

prepare_metadata <- function(filename) {
  message("Leyendo metadata: ", filename)

  metadata_original <- read.delim(
    file = filename,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )

  required_columns <- c(
    "RNA_ID",
    "LINEAGE/ECOTYPE",
    "CELL",
    "DAYS_POST_INFECTION"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(metadata_original)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Faltan estas columnas en ",
      filename,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  metadata <- metadata_original |>
    transmute(
      Sample = trimws(as.character(RNA_ID)),
      Model = trimws(as.character(CELL)),
      Day_number = as.integer(DAYS_POST_INFECTION),
      Day = paste0("Day", Day_number),
      Original_condition = trimws(
        as.character(`LINEAGE/ECOTYPE`)
      ),
      Condition = case_when(
        Original_condition == "Control" ~ "Control",
        Original_condition == "Lineage 5" ~ "Lineage_5",
        Original_condition == "M. bovis" ~ "M_bovis",
        TRUE ~ NA_character_
      )
    )

  if (anyNA(metadata$Sample) || any(metadata$Sample == "")) {
    stop("Hay identificadores de muestra vacíos en la metadata.")
  }

  if (anyDuplicated(metadata$Sample)) {
    duplicated_samples <- unique(
      metadata$Sample[duplicated(metadata$Sample)]
    )

    stop(
      "Hay muestras duplicadas en la metadata: ",
      paste(duplicated_samples, collapse = ", ")
    )
  }

  if (anyNA(metadata$Model)) {
    stop("Hay valores vacíos en la columna CELL.")
  }

  if (anyNA(metadata$Day_number)) {
    stop(
      "DAYS_POST_INFECTION contiene valores no numéricos o vacíos."
    )
  }

  if (anyNA(metadata$Condition)) {
    unexpected_conditions <- unique(
      metadata$Original_condition[is.na(metadata$Condition)]
    )

    stop(
      "Se encontraron condiciones no reconocidas: ",
      paste(unexpected_conditions, collapse = ", ")
    )
  }

  expected_models <- c("THP-1", "BoMac")
  unexpected_models <- setdiff(
    unique(metadata$Model),
    expected_models
  )

  if (length(unexpected_models) > 0) {
    stop(
      "Se encontraron modelos celulares no reconocidos: ",
      paste(unexpected_models, collapse = ", ")
    )
  }

  expected_days <- c(1L, 3L)
  unexpected_days <- setdiff(
    unique(metadata$Day_number),
    expected_days
  )

  if (length(unexpected_days) > 0) {
    stop(
      "Se encontraron días no esperados: ",
      paste(unexpected_days, collapse = ", ")
    )
  }

  metadata |>
    select(
      Sample,
      Model,
      Day_number,
      Day,
      Condition
    )
}

match_counts_and_metadata <- function(
    counts,
    metadata,
    model_name
) {
  model_metadata <- metadata |>
    filter(Model == model_name)

  missing_in_counts <- setdiff(
    model_metadata$Sample,
    colnames(counts)
  )

  missing_in_metadata <- setdiff(
    colnames(counts),
    model_metadata$Sample
  )

  if (length(missing_in_counts) > 0) {
    stop(
      "\nMuestras de ",
      model_name,
      " presentes en la metadata pero ausentes en counts:\n",
      paste(missing_in_counts, collapse = "\n")
    )
  }

  if (length(missing_in_metadata) > 0) {
    stop(
      "\nMuestras de ",
      model_name,
      " presentes en counts pero ausentes en la metadata:\n",
      paste(missing_in_metadata, collapse = "\n")
    )
  }

  ordered_counts <- counts[
    ,
    model_metadata$Sample,
    drop = FALSE
  ]

  if (!identical(
    colnames(ordered_counts),
    model_metadata$Sample
  )) {
    stop(
      "No se pudo ordenar correctamente las muestras de ",
      model_name,
      "."
    )
  }

  list(
    counts = ordered_counts,
    metadata = model_metadata
  )
}

normalize_expression <- function(
    counts,
    metadata,
    model_name
) {
  group <- interaction(
    metadata$Day,
    metadata$Condition,
    drop = TRUE
  )

  dge <- DGEList(
    counts = counts,
    group = group
  )

  keep_genes <- filterByExpr(
    dge,
    group = group
  )

  if (sum(keep_genes) == 0) {
    stop(
      "Ningún gen superó filterByExpr en ",
      model_name,
      "."
    )
  }

  dge <- dge[
    keep_genes,
    ,
    keep.lib.sizes = FALSE
  ]

  dge <- calcNormFactors(
    dge,
    method = "TMM"
  )

  log_cpm <- cpm(
    dge,
    log = TRUE,
    prior.count = 1,
    normalized.lib.sizes = TRUE
  )

  normalization_summary <- tibble(
    Model = model_name,
    Initial_genes = nrow(counts),
    Retained_genes = nrow(log_cpm),
    Removed_genes = nrow(counts) - nrow(log_cpm)
  )

  library_summary <- tibble(
    Sample = colnames(counts),
    Model = model_name,
    Raw_library_size = colSums(counts),
    Normalization_factor = dge$samples$norm.factors,
    Effective_library_size =
      dge$samples$lib.size *
      dge$samples$norm.factors
  )

  list(
    dge = dge,
    log_cpm = log_cpm,
    normalization_summary = normalization_summary,
    library_summary = library_summary
  )
}

get_hallmark_gene_sets <- function(
    species_name,
    expression_gene_ids,
    model_name
) {
  message(
    "Obteniendo Hallmark para ",
    species_name,
    "..."
  )

  hallmark_data <- tryCatch(
    {
      msigdbr(
        db_species = "HS",
        species = species_name,
        collection = "H"
      )
    },
    error = function(error_message) {
      # Compatibilidad con versiones antiguas de msigdbr.
      msigdbr(
        species = species_name,
        category = "H"
      )
    }
  )

  if (nrow(hallmark_data) == 0) {
    stop(
      "msigdbr no devolvió pathways Hallmark para ",
      species_name,
      "."
    )
  }

  possible_ensembl_columns <- c(
    "ensembl_gene",
    "ortholog_ensembl_gene",
    "db_ensembl_gene"
  )

  ensembl_column <- possible_ensembl_columns[
    possible_ensembl_columns %in% colnames(hallmark_data)
  ][1]

  if (is.na(ensembl_column)) {
    stop(
      "No se encontró una columna de identificadores Ensembl ",
      "en la salida de msigdbr para ",
      species_name,
      ". Columnas disponibles: ",
      paste(colnames(hallmark_data), collapse = ", ")
    )
  }

  pathway_column <- if (
    "gs_name" %in% colnames(hallmark_data)
  ) {
    "gs_name"
  } else {
    stop("msigdbr no devolvió la columna gs_name.")
  }

  hallmark_data$Selected_Ensembl_ID <- sub(
    "\\.[0-9]+$",
    "",
    as.character(hallmark_data[[ensembl_column]])
  )

  hallmark_data <- hallmark_data |>
    filter(
      !is.na(Selected_Ensembl_ID),
      Selected_Ensembl_ID != ""
    )

  gene_sets <- split(
    hallmark_data$Selected_Ensembl_ID,
    hallmark_data[[pathway_column]]
  )

  gene_sets <- lapply(
    gene_sets,
    unique
  )

  genes_present <- vapply(
    gene_sets,
    function(pathway_genes) {
      sum(pathway_genes %in% expression_gene_ids)
    },
    integer(1)
  )

  coverage_table <- tibble(
    Model = model_name,
    Species = species_name,
    Pathway = names(gene_sets),
    Pathway_genes_in_database = vapply(
      gene_sets,
      length,
      integer(1)
    ),
    Genes_present_in_expression = as.integer(
      genes_present
    ),
    Included = (
      genes_present >= minimum_genes_per_pathway &
      genes_present <= maximum_genes_per_pathway
    )
  ) |>
    arrange(
      desc(Included),
      desc(Genes_present_in_expression)
    )

  filtered_gene_sets <- gene_sets[
    coverage_table$Included[
      match(
        names(gene_sets),
        coverage_table$Pathway
      )
    ]
  ]

  if (length(filtered_gene_sets) == 0) {
    stop(
      "Ningún pathway Hallmark alcanzó al menos ",
      minimum_genes_per_pathway,
      " genes presentes en ",
      model_name,
      "."
    )
  }

  message(
    "  ",
    length(filtered_gene_sets),
    " pathways Hallmark conservados para ",
    model_name,
    "."
  )

  list(
    gene_sets = filtered_gene_sets,
    coverage = coverage_table,
    msigdbr_table = hallmark_data
  )
}

calculate_ssgsea <- function(
    expression_matrix,
    gene_sets,
    model_name
) {
  message("Calculando ssGSEA para ", model_name, "...")

  # API actual de GSVA.
  if (
    exists(
      "ssgseaParam",
      where = asNamespace("GSVA"),
      mode = "function"
    )
  ) {
    ssgsea_parameters <- GSVA::ssgseaParam(
      exprData = expression_matrix,
      geneSets = gene_sets,
      minSize = minimum_genes_per_pathway,
      maxSize = maximum_genes_per_pathway,
      normalize = TRUE
    )

    score_matrix <- GSVA::gsva(
      ssgsea_parameters,
      verbose = TRUE
    )
  } else {
    # Compatibilidad con versiones antiguas de GSVA.
    score_matrix <- GSVA::gsva(
      expr = expression_matrix,
      gset.idx.list = gene_sets,
      method = "ssgsea",
      kcdf = "Gaussian",
      min.sz = minimum_genes_per_pathway,
      max.sz = maximum_genes_per_pathway,
      ssgsea.norm = TRUE,
      verbose = TRUE
    )
  }

  score_matrix <- as.matrix(score_matrix)

  if (
    nrow(score_matrix) == 0 ||
    ncol(score_matrix) == 0
  ) {
    stop(
      "GSVA no produjo scores para ",
      model_name,
      "."
    )
  }

  if (anyNA(score_matrix)) {
    stop(
      "Los scores ssGSEA contienen NA para ",
      model_name,
      "."
    )
  }

  score_matrix
}

scores_to_long <- function(
    score_matrix,
    metadata
) {
  score_table <- as.data.frame(
    score_matrix,
    check.names = FALSE
  )

  score_table$Pathway <- rownames(score_table)

  score_long <- score_table |>
    pivot_longer(
      cols = -Pathway,
      names_to = "Sample",
      values_to = "Pathway_score"
    ) |>
    left_join(
      metadata,
      by = "Sample"
    )

  if (anyNA(score_long$Model)) {
    samples_without_metadata <- unique(
      score_long$Sample[is.na(score_long$Model)]
    )

    stop(
      "No se encontró metadata para estas muestras: ",
      paste(samples_without_metadata, collapse = ", ")
    )
  }

  score_long
}

rank_biserial_effect <- function(
    group_2_values,
    group_1_values
) {
  pairwise_differences <- outer(
    group_2_values,
    group_1_values,
    FUN = "-"
  )

  (
    sum(pairwise_differences > 0) -
    sum(pairwise_differences < 0)
  ) / length(pairwise_differences)
}

run_wilcoxon_comparisons <- function(
    score_matrix,
    score_long,
    metadata,
    model_name
) {
  comparisons <- tribble(
    ~Group1,     ~Group2,
    "Control",   "Lineage_5",
    "Control",   "M_bovis",
    "Lineage_5", "M_bovis"
  )

  all_results <- list()
  result_index <- 1L

  for (current_day in sort(unique(metadata$Day))) {
    for (
      comparison_index in seq_len(nrow(comparisons))
    ) {
      group_1 <- comparisons$Group1[comparison_index]
      group_2 <- comparisons$Group2[comparison_index]

      comparison_data <- score_long |>
        filter(
          Day == current_day,
          Condition %in% c(group_1, group_2)
        )

      number_group_1 <- sum(
        metadata$Day == current_day &
        metadata$Condition == group_1
      )

      number_group_2 <- sum(
        metadata$Day == current_day &
        metadata$Condition == group_2
      )

      if (
        number_group_1 < 2 ||
        number_group_2 < 2
      ) {
        warning(
          "Comparación omitida por falta de muestras: ",
          model_name,
          ", ",
          current_day,
          ", ",
          group_2,
          " vs ",
          group_1
        )

        next
      }

      pathway_results <- lapply(
        unique(comparison_data$Pathway),
        function(current_pathway) {
          current_data <- comparison_data |>
            filter(Pathway == current_pathway)

          group_1_values <- current_data |>
            filter(Condition == group_1) |>
            pull(Pathway_score)

          group_2_values <- current_data |>
            filter(Condition == group_2) |>
            pull(Pathway_score)

          wilcoxon_test <- suppressWarnings(
            wilcox.test(
              x = group_2_values,
              y = group_1_values,
              alternative = "two.sided",
              exact = FALSE,
              conf.int = FALSE
            )
          )

          tibble(
            Model = model_name,
            Day = current_day,
            Contrast = paste0(
              group_2,
              "_vs_",
              group_1
            ),
            Group1 = group_1,
            Group2 = group_2,
            Pathway = current_pathway,
            N_Group1 = length(group_1_values),
            N_Group2 = length(group_2_values),
            Mean_Group1 = mean(group_1_values),
            Mean_Group2 = mean(group_2_values),
            Median_Group1 = median(group_1_values),
            Median_Group2 = median(group_2_values),
            Delta_mean = (
              mean(group_2_values) -
              mean(group_1_values)
            ),
            Delta_median = (
              median(group_2_values) -
              median(group_1_values)
            ),
            Rank_biserial = rank_biserial_effect(
              group_2_values,
              group_1_values
            ),
            P_value = wilcoxon_test$p.value
          )
        }
      ) |>
        bind_rows() |>
        mutate(
          FDR = p.adjust(
            P_value,
            method = "BH"
          ),
          Absolute_effect_size = abs(Rank_biserial),
          Absolute_delta_median = abs(Delta_median),
          Direction = case_when(
            Delta_median > 0 ~ paste0(
              "Higher_in_",
              group_2
            ),
            Delta_median < 0 ~ paste0(
              "Higher_in_",
              group_1
            ),
            TRUE ~ "No_median_difference"
          ),
          Significant_FDR_0_05 = FDR < fdr_threshold
        ) |>
        arrange(
          FDR,
          desc(Absolute_effect_size),
          desc(Absolute_delta_median)
        )

      all_results[[result_index]] <- pathway_results
      result_index <- result_index + 1L

      comparison_name <- safe_filename(
        paste(
          model_name,
          current_day,
          group_2,
          "vs",
          group_1,
          sep = "_"
        )
      )

      write_tsv(
        pathway_results,
        file.path(
          results_directory,
          "04_Statistical_tables",
          paste0(
            "Wilcoxon_",
            comparison_name,
            ".tsv"
          )
        )
      )

      # -----------------------------------------------------------------------
      # Boxplots
      # -----------------------------------------------------------------------

      pathways_for_boxplot <- pathway_results |>
        arrange(
          FDR,
          desc(Absolute_effect_size)
        ) |>
        slice_head(
          n = number_pathways_boxplot
        ) |>
        pull(Pathway)

      boxplot_data <- comparison_data |>
        filter(
          Pathway %in% pathways_for_boxplot
        ) |>
        mutate(
          Condition = factor(
            Condition,
            levels = c(group_1, group_2)
          ),
          Pathway = factor(
            Pathway,
            levels = pathways_for_boxplot
          )
        )

      boxplot_figure <- ggplot(
        boxplot_data,
        aes(
          x = Condition,
          y = Pathway_score
        )
      ) +
        geom_boxplot(
          outlier.shape = NA,
          width = 0.60
        ) +
        geom_jitter(
          width = 0.08,
          height = 0,
          size = 2.5
        ) +
        facet_wrap(
          ~Pathway,
          scales = "free_y",
          ncol = 2
        ) +
        theme_bw(
          base_size = 12
        ) +
        labs(
          title = paste(
            model_name,
            current_day,
            ":",
            group_2,
            "vs",
            group_1
          ),
          x = NULL,
          y = "ssGSEA score"
        ) +
        theme(
          axis.text.x = element_text(
            angle = 30,
            hjust = 1
          ),
          strip.text = element_text(
            size = 8
          )
        )

      ggsave(
        filename = file.path(
          results_directory,
          "05_Boxplots",
          paste0(
            "Boxplots_",
            comparison_name,
            ".pdf"
          )
        ),
        plot = boxplot_figure,
        width = 9,
        height = 9,
        units = "in"
      )

      # -----------------------------------------------------------------------
      # Heatmap
      # -----------------------------------------------------------------------

      pathways_for_heatmap <- pathway_results |>
        arrange(
          FDR,
          desc(Absolute_effect_size)
        ) |>
        slice_head(
          n = number_pathways_heatmap
        ) |>
        pull(Pathway)

      selected_samples <- metadata |>
        filter(
          Day == current_day,
          Condition %in% c(
            group_1,
            group_2
          )
        ) |>
        pull(Sample)

      heatmap_matrix <- score_matrix[
        pathways_for_heatmap,
        selected_samples,
        drop = FALSE
      ]

      heatmap_annotation <- metadata |>
        filter(
          Sample %in% selected_samples
        ) |>
        select(
          Sample,
          Condition
        ) |>
        column_to_rownames("Sample")

      heatmap_annotation <- heatmap_annotation[
        colnames(heatmap_matrix),
        ,
        drop = FALSE
      ]

      pdf(
        file = file.path(
          results_directory,
          "06_Heatmaps",
          paste0(
            "Heatmap_",
            comparison_name,
            ".pdf"
          )
        ),
        width = 9,
        height = 10
      )

      pheatmap(
        mat = heatmap_matrix,
        scale = "row",
        annotation_col = heatmap_annotation,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_colnames = TRUE,
        show_rownames = TRUE,
        fontsize_row = 7,
        main = paste(
          model_name,
          current_day,
          ":",
          group_2,
          "vs",
          group_1
        )
      )

      dev.off()
    }
  }

  bind_rows(all_results)
}

run_limma_comparisons <- function(
    score_matrix,
    metadata,
    model_name
) {
  all_results <- list()
  result_index <- 1L

  for (current_day in sort(unique(metadata$Day))) {
    day_metadata <- metadata |>
      filter(Day == current_day)

    condition_counts <- table(
      factor(
        day_metadata$Condition,
        levels = c(
          "Control",
          "Lineage_5",
          "M_bovis"
        )
      )
    )

    if (any(condition_counts < 2)) {
      warning(
        "limma omitido para ",
        model_name,
        " ",
        current_day,
        " porque falta alguna condición."
      )

      next
    }

    day_scores <- score_matrix[
      ,
      day_metadata$Sample,
      drop = FALSE
    ]

    condition <- factor(
      day_metadata$Condition,
      levels = c(
        "Control",
        "Lineage_5",
        "M_bovis"
      )
    )

    design <- model.matrix(
      ~0 + condition
    )

    colnames(design) <- levels(condition)

    contrast_matrix <- makeContrasts(
      Lineage_5_vs_Control =
        Lineage_5 - Control,
      M_bovis_vs_Control =
        M_bovis - Control,
      M_bovis_vs_Lineage_5 =
        M_bovis - Lineage_5,
      levels = design
    )

    fit <- lmFit(
      day_scores,
      design
    )

    fit <- contrasts.fit(
      fit,
      contrast_matrix
    )

    fit <- eBayes(
      fit,
      trend = TRUE,
      robust = TRUE
    )

    for (
      contrast_name in colnames(
        contrast_matrix
      )
    ) {
      contrast_results <- topTable(
        fit,
        coef = contrast_name,
        number = Inf,
        adjust.method = "BH",
        sort.by = "P"
      ) |>
        rownames_to_column(
          "Pathway"
        ) |>
        mutate(
          Model = model_name,
          Day = current_day,
          Contrast = contrast_name,
          Significant_FDR_0_05 =
            adj.P.Val < fdr_threshold,
          Direction = case_when(
            logFC > 0 ~ "Higher_in_numerator",
            logFC < 0 ~ "Higher_in_denominator",
            TRUE ~ "No_change"
          )
        ) |>
        select(
          Model,
          Day,
          Contrast,
          Pathway,
          everything()
        )

      all_results[[result_index]] <- contrast_results
      result_index <- result_index + 1L

      output_name <- safe_filename(
        paste(
          model_name,
          current_day,
          contrast_name,
          sep = "_"
        )
      )

      write_tsv(
        contrast_results,
        file.path(
          results_directory,
          "04_Statistical_tables",
          paste0(
            "limma_",
            output_name,
            ".tsv"
          )
        )
      )
    }
  }

  bind_rows(all_results)
}


# =============================================================================
# FUNCIONES PARA BOXPLOTS INDIVIDUALES Y FIGURAS MULTIPANEL A-F
# =============================================================================

pretty_pathway_label <- function(pathway_name) {
  cleaned_name <- gsub(
    pattern = "^HALLMARK_",
    replacement = "",
    x = pathway_name
  )

  cleaned_name <- gsub(
    pattern = "_",
    replacement = " ",
    x = cleaned_name
  )

  tools::toTitleCase(tolower(cleaned_name))
}

select_publication_pathways <- function(
    contrast_results,
    number_pathways = 5
) {
  priority_pathways <- c(
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

  significant_results <- contrast_results |>
    filter(adj.P.Val < fdr_threshold) |>
    arrange(adj.P.Val, desc(abs(logFC)))

  candidate_results <- if (nrow(significant_results) >= number_pathways) {
    significant_results
  } else {
    contrast_results |>
      arrange(adj.P.Val, desc(abs(logFC)))
  }

  priority_selected <- priority_pathways[
    priority_pathways %in% candidate_results$Pathway
  ]

  remaining_pathways <- candidate_results$Pathway[
    !candidate_results$Pathway %in% priority_selected
  ]

  selected_pathways <- unique(
    c(priority_selected, remaining_pathways)
  )

  head(selected_pathways, number_pathways)
}

create_one_limma_boxplot <- function(
    score_long,
    contrast_results,
    pathway_name,
    current_day,
    group_1,
    group_2,
    model_name,
    comparison_name
) {
  pathway_statistics <- contrast_results |>
    filter(Pathway == pathway_name) |>
    slice_head(n = 1)

  plot_data <- score_long |>
    filter(
      Day == current_day,
      Condition %in% c(group_1, group_2),
      Pathway == pathway_name
    ) |>
    mutate(
      Condition = factor(
        Condition,
        levels = c(group_1, group_2)
      )
    )

  pathway_title <- pretty_pathway_label(pathway_name)

  subtitle_text <- if (nrow(pathway_statistics) == 1) {
    paste0(
      "logFC = ",
      formatC(pathway_statistics$logFC, format = "f", digits = 3),
      "; FDR = ",
      formatC(pathway_statistics$adj.P.Val, format = "g", digits = 3)
    )
  } else {
    NULL
  }

  current_plot <- ggplot(
    plot_data,
    aes(
      x = Condition,
      y = Pathway_score,
      fill = Condition
    )
  ) +
    geom_boxplot(
      width = 0.60,
      outlier.shape = NA,
      alpha = 0.70
    ) +
    geom_jitter(
      width = 0.08,
      height = 0,
      size = 2.7,
      color = "black"
    ) +
    scale_fill_manual(
      values = c(
        "Control" = "#BDBDBD",
        "Lineage_5" = "#4C78A8",
        "M_bovis" = "#E45756"
      ),
      drop = FALSE
    ) +
    labs(
      title = pathway_title,
      subtitle = subtitle_text,
      x = NULL,
      y = "ssGSEA score"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        face = "bold",
        size = 10,
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 8,
        hjust = 0.5,
        color = "#4B5563"
      ),
      axis.text.x = element_text(
        angle = 20,
        hjust = 1
      ),
      panel.grid.minor = element_blank()
    )

  individual_base_name <- safe_filename(
    paste(
      comparison_name,
      pathway_name,
      sep = "_"
    )
  )

  ggsave(
    filename = file.path(
      results_directory,
      "07_Individual_boxplots",
      paste0(individual_base_name, ".pdf")
    ),
    plot = current_plot,
    width = 4.5,
    height = 4.0,
    units = "in",
    bg = "white"
  )

  ggsave(
    filename = file.path(
      results_directory,
      "07_Individual_boxplots",
      paste0(individual_base_name, ".png")
    ),
    plot = current_plot,
    width = 4.5,
    height = 4.0,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  current_plot
}

create_publication_multipanel <- function(
    score_matrix,
    score_long,
    metadata,
    limma_results,
    model_name,
    current_day,
    group_1 = "Control",
    group_2 = "M_bovis",
    figure_number
) {
  contrast_name <- paste0(group_2, "_vs_", group_1)

  contrast_results <- limma_results |>
    filter(
      Model == model_name,
      Day == current_day,
      Contrast == contrast_name
    ) |>
    arrange(adj.P.Val, desc(abs(logFC)))

  if (nrow(contrast_results) == 0) {
    warning(
      "No hay resultados limma para ",
      model_name, " ", current_day, " ", contrast_name,
      ". No se generará la figura multipanel."
    )
    return(NULL)
  }

  selected_boxplot_pathways <- select_publication_pathways(
    contrast_results,
    number_pathways = 5
  )

  heatmap_pathways <- contrast_results |>
    filter(adj.P.Val < fdr_threshold) |>
    arrange(adj.P.Val, desc(abs(logFC))) |>
    slice_head(n = 12) |>
    pull(Pathway)

  if (length(heatmap_pathways) < 6) {
    heatmap_pathways <- contrast_results |>
      arrange(adj.P.Val, desc(abs(logFC))) |>
      slice_head(n = 12) |>
      pull(Pathway)
  }

  selected_samples <- metadata |>
    filter(
      Day == current_day,
      Condition %in% c(group_1, group_2)
    ) |>
    pull(Sample)

  heatmap_matrix <- score_matrix[
    heatmap_pathways,
    selected_samples,
    drop = FALSE
  ]

  heatmap_annotation <- metadata |>
    filter(Sample %in% selected_samples) |>
    select(Sample, Condition) |>
    column_to_rownames("Sample")

  heatmap_annotation <- heatmap_annotation[
    colnames(heatmap_matrix),
    ,
    drop = FALSE
  ]

  row_labels <- vapply(
    rownames(heatmap_matrix),
    pretty_pathway_label,
    character(1)
  )

  heatmap_object <- pheatmap(
    mat = heatmap_matrix,
    scale = "row",
    annotation_col = heatmap_annotation,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_colnames = TRUE,
    show_rownames = TRUE,
    labels_row = row_labels,
    fontsize_row = 8,
    fontsize_col = 8,
    main = paste(
      model_name,
      current_day,
      ":",
      group_2,
      "vs",
      group_1
    ),
    silent = TRUE
  )

  heatmap_panel <- patchwork::wrap_elements(
    full = heatmap_object$gtable
  )

  comparison_name <- safe_filename(
    paste(
      model_name,
      current_day,
      group_2,
      "vs",
      group_1,
      sep = "_"
    )
  )

  boxplot_list <- lapply(
    selected_boxplot_pathways,
    function(pathway_name) {
      create_one_limma_boxplot(
        score_long = score_long,
        contrast_results = contrast_results,
        pathway_name = pathway_name,
        current_day = current_day,
        group_1 = group_1,
        group_2 = group_2,
        model_name = model_name,
        comparison_name = comparison_name
      )
    }
  )

  multipanel_design <- "
AA
BC
DE
F#
"

  combined_figure <- heatmap_panel +
    boxplot_list[[1]] +
    boxplot_list[[2]] +
    boxplot_list[[3]] +
    boxplot_list[[4]] +
    boxplot_list[[5]] +
    plot_layout(
      design = multipanel_design,
      heights = c(1.35, 1, 1, 1)
    ) +
    plot_annotation(
      title = paste0(
        "Figure ",
        figure_number,
        ". Hallmark pathway activity in ",
        model_name,
        " at ",
        current_day,
        " (",
        group_2,
        " vs ",
        group_1,
        ")"
      ),
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(
          face = "bold",
          size = 15,
          hjust = 0.5
        ),
        plot.tag = element_text(
          face = "bold",
          size = 15
        )
      )
    )

  output_base_name <- paste0(
    "Figure",
    figure_number,
    "_",
    comparison_name,
    "_multipanel"
  )

  ggsave(
    filename = file.path(
      results_directory,
      "08_Publication_figures",
      paste0(output_base_name, ".pdf")
    ),
    plot = combined_figure,
    width = 11,
    height = 15,
    units = "in",
    bg = "white"
  )

  ggsave(
    filename = file.path(
      results_directory,
      "08_Publication_figures",
      paste0(output_base_name, ".png")
    ),
    plot = combined_figure,
    width = 11,
    height = 15,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  selection_table <- contrast_results |>
    filter(Pathway %in% selected_boxplot_pathways) |>
    mutate(
      Figure = paste0("Figure ", figure_number),
      Panel = LETTERS[match(Pathway, selected_boxplot_pathways) + 1]
    ) |>
    arrange(match(Pathway, selected_boxplot_pathways)) |>
    select(
      Figure,
      Panel,
      Model,
      Day,
      Contrast,
      Pathway,
      logFC,
      P.Value,
      adj.P.Val
    )

  write_tsv(
    selection_table,
    file.path(
      results_directory,
      "08_Publication_figures",
      paste0(output_base_name, "_selected_pathways.tsv")
    )
  )

  selection_table
}

# =============================================================================
# 4. CARGA Y COMPROBACIÓN DE DATOS
# =============================================================================

human_raw_counts <- read_raw_counts(
  human_counts_file
)

bovine_raw_counts <- read_raw_counts(
  bovine_counts_file
)

metadata <- prepare_metadata(
  metadata_file
)

human_data <- match_counts_and_metadata(
  counts = human_raw_counts,
  metadata = metadata,
  model_name = "THP-1"
)

bovine_data <- match_counts_and_metadata(
  counts = bovine_raw_counts,
  metadata = metadata,
  model_name = "BoMac"
)

human_raw_counts <- human_data$counts
human_metadata <- human_data$metadata

bovine_raw_counts <- bovine_data$counts
bovine_metadata <- bovine_data$metadata

sample_summary <- metadata |>
  count(
    Model,
    Day,
    Condition,
    name = "Number_of_samples"
  ) |>
  arrange(
    Model,
    Day,
    Condition
  )

print(sample_summary)

write_tsv(
  sample_summary,
  file.path(
    results_directory,
    "01_QC",
    "Sample_summary.tsv"
  )
)

# =============================================================================
# 5. NORMALIZACIÓN
# =============================================================================

human_normalization <- normalize_expression(
  counts = human_raw_counts,
  metadata = human_metadata,
  model_name = "THP-1"
)

bovine_normalization <- normalize_expression(
  counts = bovine_raw_counts,
  metadata = bovine_metadata,
  model_name = "BoMac"
)

human_log_cpm <- human_normalization$log_cpm
bovine_log_cpm <- bovine_normalization$log_cpm

normalization_summary <- bind_rows(
  human_normalization$normalization_summary,
  bovine_normalization$normalization_summary
)

library_summary <- bind_rows(
  human_normalization$library_summary,
  bovine_normalization$library_summary
)

write_tsv(
  normalization_summary,
  file.path(
    results_directory,
    "01_QC",
    "Gene_filtering_summary.tsv"
  )
)

write_tsv(
  library_summary,
  file.path(
    results_directory,
    "01_QC",
    "Library_size_and_TMM_factors.tsv"
  )
)

write_tsv(
  human_log_cpm,
  file.path(
    results_directory,
    "02_Normalized_expression",
    "THP1_TMM_logCPM.tsv"
  ),
  row_names = TRUE
)

write_tsv(
  bovine_log_cpm,
  file.path(
    results_directory,
    "02_Normalized_expression",
    "BoMac_TMM_logCPM.tsv"
  ),
  row_names = TRUE
)

# =============================================================================
# 6. PATHWAYS HALLMARK
# =============================================================================

human_hallmark <- get_hallmark_gene_sets(
  species_name = "Homo sapiens",
  expression_gene_ids = rownames(
    human_log_cpm
  ),
  model_name = "THP-1"
)

bovine_hallmark <- get_hallmark_gene_sets(
  species_name = "Bos taurus",
  expression_gene_ids = rownames(
    bovine_log_cpm
  ),
  model_name = "BoMac"
)

pathway_coverage <- bind_rows(
  human_hallmark$coverage,
  bovine_hallmark$coverage
)

write_tsv(
  pathway_coverage,
  file.path(
    results_directory,
    "01_QC",
    "Hallmark_pathway_gene_coverage.tsv"
  )
)

write_tsv(
  human_hallmark$msigdbr_table,
  file.path(
    results_directory,
    "01_QC",
    "msigdbr_Hallmark_Homo_sapiens.tsv"
  )
)

write_tsv(
  bovine_hallmark$msigdbr_table,
  file.path(
    results_directory,
    "01_QC",
    "msigdbr_Hallmark_Bos_taurus.tsv"
  )
)

# =============================================================================
# 7. ssGSEA
# =============================================================================

human_scores <- calculate_ssgsea(
  expression_matrix = human_log_cpm,
  gene_sets = human_hallmark$gene_sets,
  model_name = "THP-1"
)

bovine_scores <- calculate_ssgsea(
  expression_matrix = bovine_log_cpm,
  gene_sets = bovine_hallmark$gene_sets,
  model_name = "BoMac"
)

write_tsv(
  human_scores,
  file.path(
    results_directory,
    "03_Pathway_scores",
    "THP1_Hallmark_ssGSEA_scores.tsv"
  ),
  row_names = TRUE
)

write_tsv(
  bovine_scores,
  file.path(
    results_directory,
    "03_Pathway_scores",
    "BoMac_Hallmark_ssGSEA_scores.tsv"
  ),
  row_names = TRUE
)

human_scores_long <- scores_to_long(
  score_matrix = human_scores,
  metadata = human_metadata
)

bovine_scores_long <- scores_to_long(
  score_matrix = bovine_scores,
  metadata = bovine_metadata
)

all_scores_long <- bind_rows(
  human_scores_long,
  bovine_scores_long
)

write_tsv(
  all_scores_long,
  file.path(
    results_directory,
    "03_Pathway_scores",
    "All_Hallmark_ssGSEA_scores_long.tsv"
  )
)

# =============================================================================
# 8. WILCOXON
# =============================================================================

human_wilcoxon_results <- run_wilcoxon_comparisons(
  score_matrix = human_scores,
  score_long = human_scores_long,
  metadata = human_metadata,
  model_name = "THP-1"
)

bovine_wilcoxon_results <- run_wilcoxon_comparisons(
  score_matrix = bovine_scores,
  score_long = bovine_scores_long,
  metadata = bovine_metadata,
  model_name = "BoMac"
)

all_wilcoxon_results <- bind_rows(
  human_wilcoxon_results,
  bovine_wilcoxon_results
)

write_tsv(
  all_wilcoxon_results,
  file.path(
    results_directory,
    "04_Statistical_tables",
    "All_Wilcoxon_results.tsv"
  )
)

top_wilcoxon_results <- all_wilcoxon_results |>
  group_by(
    Model,
    Day,
    Contrast
  ) |>
  arrange(
    FDR,
    desc(Absolute_effect_size),
    desc(Absolute_delta_median),
    .by_group = TRUE
  ) |>
  slice_head(
    n = 20
  ) |>
  ungroup()

write_tsv(
  top_wilcoxon_results,
  file.path(
    results_directory,
    "04_Statistical_tables",
    "Top20_Wilcoxon_results_per_comparison.tsv"
  )
)

# =============================================================================
# 9. LIMMA SOBRE LOS SCORES ssGSEA
# =============================================================================

human_limma_results <- run_limma_comparisons(
  score_matrix = human_scores,
  metadata = human_metadata,
  model_name = "THP-1"
)

bovine_limma_results <- run_limma_comparisons(
  score_matrix = bovine_scores,
  metadata = bovine_metadata,
  model_name = "BoMac"
)

all_limma_results <- bind_rows(
  human_limma_results,
  bovine_limma_results
)

write_tsv(
  all_limma_results,
  file.path(
    results_directory,
    "04_Statistical_tables",
    "All_limma_results.tsv"
  )
)

top_limma_results <- all_limma_results |>
  group_by(
    Model,
    Day,
    Contrast
  ) |>
  arrange(
    adj.P.Val,
    desc(abs(logFC)),
    .by_group = TRUE
  ) |>
  slice_head(
    n = 20
  ) |>
  ungroup()

write_tsv(
  top_limma_results,
  file.path(
    results_directory,
    "04_Statistical_tables",
    "Top20_limma_results_per_comparison.tsv"
  )
)


# =============================================================================
# 10. BOXPLOTS INDIVIDUALES Y FIGURAS MULTIPANEL A-F
# =============================================================================

publication_figure_selections <- bind_rows(
  create_publication_multipanel(
    score_matrix = human_scores,
    score_long = human_scores_long,
    metadata = human_metadata,
    limma_results = all_limma_results,
    model_name = "THP-1",
    current_day = "Day1",
    group_1 = "Control",
    group_2 = "M_bovis",
    figure_number = 14
  ),
  create_publication_multipanel(
    score_matrix = human_scores,
    score_long = human_scores_long,
    metadata = human_metadata,
    limma_results = all_limma_results,
    model_name = "THP-1",
    current_day = "Day3",
    group_1 = "Control",
    group_2 = "M_bovis",
    figure_number = 15
  ),
  create_publication_multipanel(
    score_matrix = bovine_scores,
    score_long = bovine_scores_long,
    metadata = bovine_metadata,
    limma_results = all_limma_results,
    model_name = "BoMac",
    current_day = "Day1",
    group_1 = "Control",
    group_2 = "M_bovis",
    figure_number = 16
  ),
  create_publication_multipanel(
    score_matrix = bovine_scores,
    score_long = bovine_scores_long,
    metadata = bovine_metadata,
    limma_results = all_limma_results,
    model_name = "BoMac",
    current_day = "Day3",
    group_1 = "Control",
    group_2 = "M_bovis",
    figure_number = 17
  )
)

write_tsv(
  publication_figure_selections,
  file.path(
    results_directory,
    "08_Publication_figures",
    "All_selected_pathways_for_multipanel_figures.tsv"
  )
)

# =============================================================================
# 11. HEATMAP GLOBAL DE CAMBIOS DE PATHWAYS
# =============================================================================

if (nrow(all_limma_results) > 0) {
  global_pathways <- all_limma_results |>
    group_by(Pathway) |>
    summarise(
      Maximum_absolute_logFC = max(
        abs(logFC),
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    arrange(
      desc(Maximum_absolute_logFC)
    ) |>
    slice_head(
      n = 30
    ) |>
    pull(Pathway)

  global_heatmap_matrix <- all_limma_results |>
    filter(
      Pathway %in% global_pathways
    ) |>
    mutate(
      Full_comparison = paste(
        Model,
        Day,
        Contrast,
        sep = " | "
      )
    ) |>
    select(
      Pathway,
      Full_comparison,
      logFC
    ) |>
    distinct() |>
    pivot_wider(
      names_from = Full_comparison,
      values_from = logFC
    ) |>
    column_to_rownames(
      "Pathway"
    ) |>
    as.matrix()

  pdf(
    file = file.path(
      results_directory,
      "06_Heatmaps",
      "Global_Hallmark_logFC_heatmap.pdf"
    ),
    width = 14,
    height = 11
  )

  pheatmap(
    mat = global_heatmap_matrix,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize_row = 7,
    fontsize_col = 7,
    main = "Hallmark pathway changes across all comparisons"
  )

  dev.off()
}

# =============================================================================
# 12. EXCEL RESUMEN Y REPRODUCIBILIDAD
# =============================================================================

excel_sheets <- list(
  Sample_summary = sample_summary,
  Filtering_summary = normalization_summary,
  Library_summary = library_summary,
  Pathway_coverage = pathway_coverage,
  Wilcoxon_all = all_wilcoxon_results,
  Wilcoxon_top20 = top_wilcoxon_results,
  Limma_all = all_limma_results,
  Limma_top20 = top_limma_results,
  Figure_pathways = publication_figure_selections
)

write_xlsx(
  excel_sheets,
  path = file.path(
    results_directory,
    "Pathway_analysis_summary.xlsx"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    results_directory,
    "sessionInfo.txt"
  )
)

writeLines(
  c(
    "Pathway analysis completed successfully.",
    paste0(
      "Working directory: ",
      getwd()
    ),
    paste0(
      "Results directory: ",
      normalizePath(
        results_directory,
        mustWork = FALSE
      )
    ),
    paste0(
      "Date: ",
      Sys.time()
    )
  ),
  con = file.path(
    results_directory,
    "ANALYSIS_COMPLETED.txt"
  )
)

cat(
  "\n============================================================\n",
  "ANÁLISIS COMPLETADO CORRECTAMENTE\n",
  "Resultados guardados en:\n",
  normalizePath(
    results_directory,
    mustWork = FALSE
  ),
  "\n============================================================\n",
  sep = ""
)
