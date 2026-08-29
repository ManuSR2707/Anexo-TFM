# =========================
# REDES LIMPIAS PARA CYTOSCAPE
# =========================

library(dplyr)

crear_red_limpia <- function(input_file, output_prefix, top_n = 40) {
  
  res <- read.delim(input_file,
                    header = TRUE,
                    sep = "\t",
                    check.names = FALSE)
  
  res <- res[order(res$qval), ]
  
  # Eliminar duplicados ligando-receptor
  res_clean <- res %>%
    group_by(L, R) %>%
    summarise(
      qval = min(qval),
      pval = min(pval),
      L.logFC = first(L.logFC),
      R.logFC = first(R.logFC),
      L.expr = first(L.expr),
      R.expr = first(R.expr),
      pathways = paste(unique(pw.name), collapse = "; "),
      n_pathways = n_distinct(pw.name),
      .groups = "drop"
    ) %>%
    arrange(qval)
  
  # Quedarse con top N
  res_top <- head(res_clean, top_n)
  
  # Tabla de aristas
  edges <- data.frame(
    Source = res_top$L,
    Target = res_top$R,
    interaction = paste(res_top$L, res_top$R, sep = "_"),
    qval = res_top$qval,
    pval = res_top$pval,
    minus_log10_qval = -log10(res_top$qval),
    n_pathways = res_top$n_pathways,
    pathways = res_top$pathways
  )
  
  # Tabla de nodos
  ligands <- data.frame(
    id = res_top$L,
    type = "Ligand",
    logFC = res_top$L.logFC,
    expr = res_top$L.expr
  )
  
  receptors <- data.frame(
    id = res_top$R,
    type = "Receptor",
    logFC = res_top$R.logFC,
    expr = res_top$R.expr
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
  
  # Calcular degree
  degree_table <- table(c(edges$Source, edges$Target))
  nodes$degree <- as.numeric(degree_table[nodes$id])
  
  # Exportar
  write.table(edges,
              paste0(output_prefix, "_edges_clean.tsv"),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
  
  write.table(nodes,
              paste0(output_prefix, "_nodes_clean.tsv"),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
}


# =========================
# 1. INFECTED VS CONTROL
# =========================

crear_red_limpia(
  input_file = "BulkSignalR_results_immune.tsv",
  output_prefix = "Network_Infected_vs_Control",
  top_n = 40
)


# =========================
# 2. DAY 3 VS DAY 1
# =========================

crear_red_limpia(
  input_file = "BulkSignalR_Day3_vs_Day1_results.tsv",
  output_prefix = "Network_Day3_vs_Day1",
  top_n = 40
)


# =========================
# 3. M. BOVIS VS LINEAGE 5
# =========================

crear_red_limpia(
  input_file = "BulkSignalR_M_bovis_vs_Lineage_5_results.tsv",
  output_prefix = "Network_M_bovis_vs_Lineage_5",
  top_n = 40
)
