# =========================
# 1. LIBRERÍAS
# =========================

library(BulkSignalR)
library(edgeR)
library(org.Hs.eg.db)
library(AnnotationDbi)


# =========================
# 2. CARGA DE DATOS
# =========================

matriz_humano <- read.table("humanRawCounts.txt",
                            header = TRUE,
                            row.names = 1)

meta <- read.table("metadata_rnaseq.tsv",
                   header = TRUE,
                   sep = "\t")


# =========================
# 3. CONVERSIÓN ENSEMBL → HUGO
# =========================

rownames(matriz_humano) <- sub("\\..*", "", rownames(matriz_humano))

symbols <- mapIds(org.Hs.eg.db,
                  keys = rownames(matriz_humano),
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first")

keep <- !is.na(symbols)

matriz_humano <- rowsum(as.matrix(matriz_humano[keep, ]),
                        group = symbols[keep])

dim(matriz_humano)
head(rownames(matriz_humano))
head(matriz_humano[, 1:3])


# =========================
# 4. ALINEAR METADATA
# =========================

meta <- read.delim("metadata_rnaseq.tsv",
                   header = TRUE,
                   check.names = FALSE)

meta <- meta[match(colnames(matriz_humano), meta$RNA_ID), ]

all(colnames(matriz_humano) == meta$RNA_ID)


# =========================
# 5. FILTRAR SOLO LINAJES DE INTERÉS
# =========================

linaje_A <- "Lineage 5"
linaje_B <- "M. bovis"

idx_lineage <- meta$`LINEAGE/ECOTYPE` %in% c(linaje_A, linaje_B)

matriz_humano <- matriz_humano[, idx_lineage]
meta <- meta[idx_lineage, ]

group <- factor(meta$`LINEAGE/ECOTYPE`,
                levels = c(linaje_A, linaje_B))

table(group)


# =========================
# 6. NORMALIZACIÓN Y ANÁLISIS DIFERENCIAL
# =========================

dge <- DGEList(counts = matriz_humano, group = group)

keep_genes <- filterByExpr(dge, group = group)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)

expr_norm <- cpm(dge, log = FALSE)

design <- model.matrix(~ group)

dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

res <- glmQLFTest(fit, coef = 2)

dge_table <- topTags(res, n = Inf)$table

head(dge_table)
summary(dge_table$PValue)


# =========================
# 7. PREPARAR STATS
# =========================

stats <- data.frame(
  pval = dge_table$PValue,
  logFC = dge_table$logFC,
  expr = rowMeans(cpm(dge, log = TRUE))
)

rownames(stats) <- rownames(dge_table)


# =========================
# 8. PREPARACIÓN BULKSIGNALR
# =========================

bsrdm <- BSRDataModel(expr_norm)
bsrdm.comp <- as(bsrdm, "BSRDataModelComp")

nrow(stats)
nrow(ncounts(bsrdm.comp))

common_genes <- intersect(rownames(ncounts(bsrdm.comp)), rownames(stats))

stats2 <- stats[common_genes, ]
expr2  <- expr_norm[common_genes, ]

bsrdm <- BSRDataModel(expr2)
bsrdm.comp <- as(bsrdm, "BSRDataModelComp")

stats2 <- stats2[rownames(ncounts(bsrdm.comp)), ]

nrow(stats2) == nrow(ncounts(bsrdm.comp))
all(rownames(stats2) == rownames(ncounts(bsrdm.comp)))


# =========================
# 9. DEFINIR COMPARACIÓN M. BOVIS VS LINEAGE 5
# =========================

colLineage5 <- which(group == linaje_A)
colBovis    <- which(group == linaje_B)

bsrcc <- BSRClusterComp(bsrdm.comp,
                        colBovis,
                        colLineage5,
                        stats2)

comparison_name <- "M_bovis_vs_Lineage_5"

bsrdm.comp <- addClusterComp(bsrdm.comp,
                             bsrcc,
                             comparison_name)


# =========================
# 10. INFERENCIA L-R M. BOVIS VS LINEAGE 5
# =========================

bsrinf.comp <- BSRInferenceComp(
  bsrdm.comp,
  reference = "REACTOME",
  max.pval = 0.05,
  comparison_name
)

resLR <- LRinter(bsrinf.comp)

head(resLR)


# =========================
# 11. FILTRADO DE RESULTADOS
# =========================

resLR <- resLR[order(resLR$qval), ]

resLR_sig <- resLR[resLR$qval < 0.05, ]

immune_keywords <- "IL|TNF|CXCL|CCL|IFN|TLR|CD40|CCR|CXCR"

resLR_immune <- resLR_sig[
  grepl(immune_keywords, resLR_sig$L) |
    grepl(immune_keywords, resLR_sig$R) |
    grepl(immune_keywords, resLR_sig$pw.name),
]

resLR_top20 <- head(resLR_sig, 20)


# =========================
# 12. EXPORTAR RESULTADOS
# =========================

write.table(resLR,
            "BulkSignalR_M_bovis_vs_Lineage_5_results.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

write.table(resLR_sig,
            "BulkSignalR_M_bovis_vs_Lineage_5_significant.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

write.table(resLR_immune,
            "BulkSignalR_M_bovis_vs_Lineage_5_immune.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

write.table(resLR_top20,
            "BulkSignalR_M_bovis_vs_Lineage_5_top20.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
