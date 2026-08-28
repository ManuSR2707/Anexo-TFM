# =========================
# 1. LIBRERÍAS
# =========================

library(BulkSignalR)
library(edgeR)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(orthogene)


# =========================
# 2. CARGA DE DATOS
# =========================

human_raw <- read.table("humanRawCounts.txt",
                        header = TRUE,
                        row.names = 1,
                        check.names = FALSE)

bos_raw <- read.table("bosRawCounts.txt",
                      header = TRUE,
                      row.names = 1,
                      check.names = FALSE)

meta <- read.delim("metadata_rnaseq.tsv",
                   header = TRUE,
                   check.names = FALSE)


# =========================
# 3. HUMANO: ENSEMBL → HUGO
# =========================

rownames(human_raw) <- sub("\\..*", "", rownames(human_raw))

symbols_human <- mapIds(org.Hs.eg.db,
                        keys = rownames(human_raw),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")

keep_human <- !is.na(symbols_human)

human_hugo <- rowsum(as.matrix(human_raw[keep_human, ]),
                     group = symbols_human[keep_human])


# =========================
# 4. BOVINO: ORTÓLOGOS → HUMANO/HUGO
# =========================

rownames(bos_raw) <- sub("\\..*", "", rownames(bos_raw))

ortholog_dict <- findOrthoGenes(
  from_organism = "btaurus",
  from_values = rownames(bos_raw)
)

bos_hugo <- convertToHuman(
  counts = bos_raw,
  dictionary = ortholog_dict
)

bos_hugo <- as.matrix(bos_hugo)


# =========================
# 5. UNIR MATRICES POR GENES COMUNES
# =========================

common_genes <- intersect(rownames(human_hugo),
                          rownames(bos_hugo))

human_hugo <- human_hugo[common_genes, ]
bos_hugo   <- bos_hugo[common_genes, ]

matriz_combined <- cbind(human_hugo, bos_hugo)


# =========================
# 6. ALINEAR METADATA
# =========================

meta <- meta[match(colnames(matriz_combined), meta$RNA_ID), ]

stopifnot(all(colnames(matriz_combined) == meta$RNA_ID))


# =========================
# 7. CREAR GRUPO COMBINADO
# =========================

meta$LINEAGE_CELL <- paste(meta$`LINEAGE/ECOTYPE`,
                           meta$CELL,
                           sep = "_")

meta$LINEAGE_CELL <- gsub(" ", "_", meta$LINEAGE_CELL)
meta$LINEAGE_CELL <- gsub("\\.", "", meta$LINEAGE_CELL)

group_all <- factor(meta$LINEAGE_CELL)

table(group_all)


# =========================
# 8. CREAR COMPARACIONES 1 A 1
# =========================

comparisons <- combn(levels(group_all), 2, simplify = FALSE)


# =========================
# 9. FUNCIÓN PRINCIPAL
# =========================

run_bulksignalr_comparison <- function(group_A, group_B) {
  
  message("====================================")
  message("Analizando: ", group_A, " vs ", group_B)
  message("====================================")
  
  idx <- group_all %in% c(group_A, group_B)
  
  matriz_sub <- matriz_combined[, idx]
  group_sub <- droplevels(group_all[idx])
  
  print(table(group_sub))
  
  if (any(table(group_sub) < 2)) {
    message("Saltando comparación: algún grupo tiene menos de 2 muestras")
    return(NULL)
  }
  
  dge <- DGEList(counts = matriz_sub)
  dge <- calcNormFactors(dge)
  
  design <- model.matrix(~ group_sub)
  
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  res <- glmQLFTest(fit, coef = 2)
  
  dge_table <- topTags(res, n = Inf)$table
  
  stats <- data.frame(
    pval = dge_table$PValue,
    logFC = dge_table$logFC,
    expr = rowMeans(cpm(dge, log = TRUE))
  )
  
  rownames(stats) <- rownames(dge_table)
  
  expr_norm <- cpm(dge, log = FALSE)
  
  bsrdm <- BSRDataModel(expr_norm)
  bsrdm.comp <- as(bsrdm, "BSRDataModelComp")
  
  common <- intersect(rownames(ncounts(bsrdm.comp)),
                      rownames(stats))
  
  stats2 <- stats[common, ]
  expr2  <- expr_norm[common, ]
  
  bsrdm <- BSRDataModel(expr2)
  bsrdm.comp <- as(bsrdm, "BSRDataModelComp")
  
  stats2 <- stats2[rownames(ncounts(bsrdm.comp)), ]
  
  stopifnot(nrow(stats2) == nrow(ncounts(bsrdm.comp)))
  stopifnot(all(rownames(stats2) == rownames(ncounts(bsrdm.comp))))
  
  colA <- which(group_sub == group_A)
  colB <- which(group_sub == group_B)
  
  comparison_name <- paste(group_A, "vs", group_B, sep = "_")
  
  bsrcc <- BSRClusterComp(bsrdm.comp,
                          colA,
                          colB,
                          stats2)
  
  bsrdm.comp <- addClusterComp(bsrdm.comp,
                               bsrcc,
                               comparison_name)
  
  bsrinf.comp <- tryCatch(
    BSRInferenceComp(
      bsrdm.comp,
      reference = "REACTOME",
      max.pval = 0.05,
      comparison_name
    ),
    error = function(e) {
      message("No se encontraron interacciones para: ", comparison_name)
      message("Error: ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(bsrinf.comp)) {
    return(NULL)
  }
  
  resLR <- LRinter(bsrinf.comp)
  resLR <- resLR[order(resLR$qval), ]
  
  resLR_sig <- resLR[resLR$qval < 0.05, ]
  
  immune_keywords <- "IL|TNF|CXCL|CCL|IFN|TLR|CD40|CCR|CXCR"
  
  resLR_immune <- resLR_sig[
    grepl(immune_keywords, resLR_sig$L) |
      grepl(immune_keywords, resLR_sig$R) |
      grepl(immune_keywords, resLR_sig$pw.name),
  ]
  
  resLR_top20 <- head(resLR_sig, 20)
  

  write.table(resLR_sig,
              paste0("BulkSignalR_", comparison_name, "_significant.tsv"),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)

  message("Finalizado: ", comparison_name)
}


# =========================
# 10. EJECUTAR TODAS LAS COMPARACIONES
# =========================

for (comp in comparisons) {
  run_bulksignalr_comparison(comp[1], comp[2])
}
