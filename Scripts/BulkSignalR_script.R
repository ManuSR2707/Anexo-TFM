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
# 3. CONVERSIÓN ENSMBL → HUGO
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

# Comprobación
dim(matriz_humano)
head(rownames(matriz_humano))
head(matriz_humano[,1:3])


# =========================
# 4. NORMALIZACIÓN (edgeR)
# =========================

dge <- DGEList(counts = matriz_humano)
dge <- calcNormFactors(dge)

expr_norm <- cpm(dge, log = FALSE)


# =========================
# 5. ALINEAR METADATA
# =========================

meta <- read.delim("metadata_rnaseq.tsv",
                   header = TRUE,
                   check.names = FALSE)

meta <- meta[match(colnames(matriz_humano), meta$RNA_ID), ]

all(colnames(matriz_humano) == meta$RNA_ID)


# =========================
# 6. DEFINIR GRUPOS
# =========================

group <- ifelse(meta$STRAIN == "Uninfected", "Control", "Infected")
group <- factor(group)

table(group)


# =========================
# 7. ANALISIS DIFERENCIAL
# =========================

design <- model.matrix(~ group)

dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

res <- glmQLFTest(fit, coef = 2)

dge_table <- topTags(res, n = Inf)$table

# Comprobación
head(dge_table)
summary(dge_table$PValue)


# =========================
# 8. PREPARAR STATS
# =========================

stats <- data.frame(
  pval = dge_table$PValue,
  logFC = dge_table$logFC,
  expr = rowMeans(cpm(dge, log = TRUE))
)

rownames(stats) <- rownames(dge_table)


# =========================
# 9. PREPARACIÓN BULKSIGNALR
# =========================

# Crear objeto inicial
bsrdm <- BSRDataModel(expr_norm)
bsrdm.comp <- as(bsrdm, "BSRDataModelComp")

# Comprobar dimensiones
nrow(stats)
nrow(ncounts(bsrdm.comp))

# Intersección de genes
common_genes <- intersect(rownames(ncounts(bsrdm.comp)), rownames(stats))

stats2 <- stats[common_genes, ]
expr2  <- expr_norm[common_genes, ]

# Recrear modelo con genes filtrados
bsrdm <- BSRDataModel(expr2)
bsrdm.comp <- as(bsrdm, "BSRDataModelComp")

# Reordenar stats
stats2 <- stats2[rownames(ncounts(bsrdm.comp)), ]

# Comprobaciones
nrow(stats2) == nrow(ncounts(bsrdm.comp))
all(rownames(stats2) == rownames(ncounts(bsrdm.comp)))


# =========================
# 10. DEFINIR COMPARACIÓN CONTROL VS INFECTED
# =========================

colControl  <- which(group == "Control")
colInfected <- which(group == "Infected")

bsrcc <- BSRClusterComp(bsrdm.comp,
                        colInfected,
                        colControl,
                        stats2)

bsrdm.comp <- addClusterComp(bsrdm.comp,
                             bsrcc,
                             "Infected_vs_Control")


# ===================
# 11. INFERENCIA L-R CONTROL VS INFECTED
# ===================

bsrinf.comp <- BSRInferenceComp(
  bsrdm.comp,
  reference = "REACTOME",
  max.pval = 0.05,
  "Infected_vs_Control"
)

resLR <- LRinter(bsrinf.comp)

# Comprobación
head(resLR)


# =========================
# 12. FILTRADO DE RESULTADOS
# =========================

# Ordenar por significancia
resLR <- resLR[order(resLR$qval), ]

# Significativos
resLR_sig <- resLR[resLR$qval < 0.05, ]

# Interacciones inmunes
immune_keywords <- "IL|TNF|CXCL|CCL|IFN|TLR|CD40|CCR|CXCR"

resLR_immune <- resLR_sig[
  grepl(immune_keywords, resLR_sig$L) |
    grepl(immune_keywords, resLR_sig$R) |
    grepl(immune_keywords, resLR_sig$pw.name),
]

# Top 20
resLR_top20 <- head(resLR_sig, 20)


# =========================
# 13. EXPORTAR RESULTADOS
# =========================

write.table(resLR,
            "BulkSignalR_results.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

write.table(resLR_sig,
            "BulkSignalR_results_significant.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

write.table(resLR_immune,
            "BulkSignalR_results_immune.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

write.table(resLR_top20,
            "BulkSignalR_results_top20.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
