library(dplyr)

crear_red_global <- function(top_n_por_comparacion = 40) {
  
  # =========================
  # 1. Cargar resultados
  # =========================
  
  infected <- read.delim("BulkSignalR_results_immune.tsv",
                         header = TRUE,
                         sep = "\t",
                         check.names = FALSE)
  
  day <- read.delim("BulkSignalR_Day3_vs_Day1_results.tsv",
                    header = TRUE,
                    sep = "\t",
                    check.names = FALSE)
  
  lineage <- read.delim("BulkSignalR_M_bovis_vs_Lineage_5_results.tsv",
                        header = TRUE,
                        sep = "\t",
                        check.names = FALSE)
  
  
  # =========================
  # 2. Añadir nombre de comparación
  # =========================
  
  infected$comparison <- "Infected_vs_Control"
  day$comparison <- "Day3_vs_Day1"
  lineage$comparison <- "M_bovis_vs_Lineage_5"
  
  
  # =========================
  # 3. Ordenar y seleccionar top N de cada comparación
  # =========================
  
  infected <- infected[order(infected$qval), ]
  day <- day[order(day$qval), ]
  lineage <- lineage[order(lineage$qval), ]
  
  infected_top <- head(infected, top_n_por_comparacion)
  day_top <- head(day, top_n_por_comparacion)
  lineage_top <- head(lineage, top_n_por_comparacion)
  
  
  # =========================
  # 4. Unir resultados
  # =========================
  
  res_all <- rbind(infected_top, day_top, lineage_top)
  
  
  # =========================
  # 5. Eliminar duplicados L-R
  # =========================
  
  res_clean <- res_all %>%
    group_by(L, R) %>%
    summarise(
      qval = min(qval),
      pval = min(pval),
      L.logFC = mean(L.logFC, na.rm = TRUE),
      R.logFC = mean(R.logFC, na.rm = TRUE),
      L.expr = mean(L.expr, na.rm = TRUE),
      R.expr = mean(R.expr, na.rm = TRUE),
      comparisons = paste(unique(comparison), collapse = "; "),
      n_comparisons = n_distinct(comparison),
      pathways = paste(unique(pw.name), collapse = "; "),
      n_pathways = n_distinct(pw.name),
      .groups = "drop"
    ) %>%
    arrange(qval)
  
  
  # =========================
  # 6. Crear tabla de aristas
  # =========================
  
  edges <- data.frame(
    Source = res_clean$L,
    Target = res_clean$R,
    interaction = paste(res_clean$L, res_clean$R, sep = "_"),
    qval = res_clean$qval,
    pval = res_clean$pval,
    minus_log10_qval = -log10(res_clean$qval),
    comparisons = res_clean$comparisons,
    n_comparisons = res_clean$n_comparisons,
    n_pathways = res_clean$n_pathways,
    pathways = res_clean$pathways
  )
  
  
  # =========================
  # 7. Crear tabla de nodos
  # =========================
  
  ligands <- data.frame(
    id = res_clean$L,
    type = "Ligand",
    logFC = res_clean$L.logFC,
    expr = res_clean$L.expr
  )
  
  receptors <- data.frame(
    id = res_clean$R,
    type = "Receptor",
    logFC = res_clean$R.logFC,
    expr = res_clean$R.expr
  )
  
  nodes <- rbind(ligands, receptors)
  
  nodes <- nodes %>%
    group_by(id) %>%
    summarise(
      type = paste(unique(type), collapse = "/"),
      logFC = mean(logFC, na.rm = TRUE),
      expr = mean(expr, na.rm = TRUE),
      .groups = "drop"
    )
  
  degree_table <- table(c(edges$Source, edges$Target))
  nodes$degree <- as.numeric(degree_table[nodes$id])
  
  
  # =========================
  # 8. Exportar
  # =========================
  
  write.table(edges,
              "Network_Global_edges_clean.tsv",
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  write.table(nodes,
              "Network_Global_nodes_clean.tsv",
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
}

crear_red_global(top_n_por_comparacion = 40)

