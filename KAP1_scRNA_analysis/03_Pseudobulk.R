## 03_Pseudobulk.R
##
## Pseudobulk (cell-type level) limma-voom differential expression between
## KO and WT for each of the 5 annotated lineages, plus GO/KEGG
## over-representation for each comparison.
##
## Input : RDS/SeuratObj.rds (from 02_Integration.R)
##         Annotation/Mus_musculus.gene_info.gz (NCBI gene info table)
## Output: CPM_afterFiltering_KAP1_scRNA.csv
##         GeneList/DE_<comparison>.csv, GeneList/KEGG-<comparison>.csv,
##         GeneList/GO-<comparison>.csv, GeneList/glimmaMAplot-<comparison>.html

source("setup.R")
source("config.R")
library(limma)
library(Glimma)

data_seurat <- readRDS(file.path(ROBJECT_DIR, "SeuratObj.rds"))
DefaultAssay(data_seurat) <- "integrated"
data_seurat <- subset(data_seurat, subset = cell_type != "Cluster6")

## ---- Build pseudobulk samples ---------------------------------------------
y <- Seurat2PB(data_seurat, sample = "sample_name", cluster = "cell_type")
condition_map <- setNames(samples$condition, samples$sample_id)
y$samples$condition <- unname(condition_map[sub("_.*", "", rownames(y$samples))])
if (anyNA(y$samples$condition)) {
  stop("Condition is missing for one or more pseudobulk samples. Check the sample sheet in config.R.")
}
y$samples$group <- paste(y$samples$condition, y$samples$cluster, sep = "_")

## Shorten cluster labels for sample names (e.g. "SK-K1_clusterBasal" -> "SK-K1_Ba").
for (long in names(pb_cluster_labels)) {
  short <- pb_cluster_labels[long]
  rownames(y$samples) <- gsub(paste0("cluster", long), short, rownames(y$samples), fixed = TRUE)
  colnames(y$counts)  <- gsub(paste0("cluster", long), short, colnames(y$counts), fixed = TRUE)
}

lib_sizes <- y$samples$lib.size
bar_cols <- bar_colors_for(rownames(y$samples))

par(mar = c(12, 5, 4, 2) + 0.1)
bp <- barplot(lib_sizes / 1e6, ylim = c(0, 70), col = bar_cols,
              ylab = "Library Size (in millions)", main = "Pseudobulk Library Sizes",
              xaxt = "n", yaxt = "n")
axis(1, at = bp, labels = rownames(y$samples), las = 2)
axis(2, at = seq(0, 70, by = 5), las = 2)

print(table(data_seurat$sample_name, data_seurat$cell_type))

## ---- Gene annotation + filtering ------------------------------------------
gene_info <- read.delim(file.path(annotation_dir, "Mus_musculus.gene_info.gz"), header = TRUE)[, -1]
m <- match(rownames(y), gene_info$Symbol)
genes <- cbind(y$genes$gene, gene_info[m, c(2, 6, 9)])
colnames(genes) <- c("gene", "Symbol", "Chr", "Type")

keep_annotated <- !is.na(genes$Symbol) & !genes$Chr %in% c("MT", "Un") &
  !genes$Type %in% c("pseudo", "unknown", "other")
y$genes <- genes
y <- y[keep_annotated, , keep.lib.size = FALSE]

y$genes$Entrez <- mapIds(org.Mm.eg.db, keys = y$genes$Symbol, column = "ENTREZID",
                          keytype = "SYMBOL", multiVals = "first")

y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
y <- normLibSizes(y)

cpm_df <- data.frame(GeneID = rownames(y), cpm(y, log = FALSE), check.names = FALSE)
cpm_df <- merge(y$genes, cpm_df, by.x = "Symbol", by.y = "GeneID")
write.csv(cpm_df, file = "CPM_afterFiltering_KAP1_scRNA.csv", row.names = FALSE)

## ---- MDS plots ---------------------------------------------------------
group <- factor(y$samples$group)
condition <- factor(y$samples$condition)
cluster <- factor(y$samples$cluster, levels = c("Basal", "Cycling.Basal", "LP", "Cycling.LP", "ML"))
donor <- factor(y$samples$sample)

mds <- plotMDS(y, plot = FALSE)
par(mfrow = c(2, 2))
plotMDS(mds, col = col.p2[group], main = "MDS by group")
legend("topleft", legend = levels(group), pch = 16, col = col.p2, ncol = 2)
plotMDS(mds, col = col.p2[21:length(col.p2)][condition], main = "MDS by condition")
legend("topleft", legend = levels(condition), pch = 16, col = col.p2[21:length(col.p2)])
plotMDS(mds, col = col.p2[23:length(col.p2)][cluster], main = "MDS by cell type")
legend("topleft", legend = levels(cluster), pch = 16, col = col.p2[23:length(col.p2)])
plotMDS(mds, col = col.pLight[donor], main = "MDS by donor")
legend("topleft", legend = levels(donor), pch = 16, col = col.pLight)

## ---- Design, contrasts, voom fit -------------------------------------------
design <- model.matrix(~0 + group)
colnames(design) <- gsub("group", "", colnames(design))

contrast.matrix <- makeContrasts(
  Basal_KOvWT    = KO_Basal - WT_Basal,
  CycBasal_KOvWT = KO_Cycling.Basal - WT_Cycling.Basal,
  LP_KOvWT       = KO_LP - WT_LP,
  CycLP_KOvWT    = KO_Cycling.LP - WT_Cycling.LP,
  ML_KOvWT       = KO_ML - WT_ML,
  levels = design
)

donor <- factor(y$samples$sample)
v <- voomLmFit(y, design, block = donor, sample.weights = TRUE, plot = TRUE, save.plot = TRUE)

par(mar = c(10, 4, 4, 2))
barplot(v$targets$sample.weight, names = rownames(v$targets), main = "Sample-specific weights",
        ylab = "Weight", xlab = "", las = 2, ylim = c(0, 3), yaxt = "n",
        col = bar_colors_for(rownames(y$samples)))
axis(side = 2, at = seq(0, 3, by = 0.2))
abline(h = 1, col = 2, lty = 2)
mtext("Sample", side = 1, line = 8, adj = 0.5)

cfit <- contrasts.fit(v, contrast.matrix)
efit <- eBayes(cfit)

print(plotSA(efit, main = "Final model: Mean-variance trend"))
dt <- decideTests(efit)
print(summary(dt))

## ---- Per-comparison DE + GO/KEGG -------------------------------------------
comparison_titles <- c(
  Basal_KOvWT    = "Basal: KO vs WT",
  CycBasal_KOvWT = "Cycling Basal: KO vs WT",
  LP_KOvWT       = "LP: KO vs WT",
  CycLP_KOvWT    = "Cycling LP: KO vs WT",
  ML_KOvWT       = "ML: KO vs WT"
)

for (comparison in names(comparison_titles)) {
  title <- comparison_titles[comparison]
  message("Pseudobulk DE: ", title)

  tt <- topTable(efit, coef = comparison, sort.by = "P", number = Inf)
  print(head(tt, 5))
  write.csv(tt, file = file.path(GENELIST_DIR, paste0("DE_", comparison, ".csv")))

  print(plotMD(efit, coef = comparison, status = dt[, comparison], legend = "topright", main = title))

  htmlwidgets::saveWidget(
    Glimma::glimmaMA(efit, coef = comparison, status = dt[, comparison], dge = y, main = title),
    file = file.path(GENELIST_DIR, paste0("glimmaMAplot-", comparison, ".html"))
  )

  keg <- kegga(efit, coef = comparison, species = "Mm", geneid = "Entrez", FDR = 0.05)
  print(topKEGG(keg, sort = "up", truncate = 45))
  print(topKEGG(keg, sort = "down", truncate = 45))
  write.csv(topKEGG(keg, n = Inf), file = file.path(GENELIST_DIR, paste0("KEGG-", comparison, ".csv")), row.names = FALSE)

  GO <- goana(efit, coef = comparison, species = "Mm", geneid = "Entrez", FDR = 0.05)
  print(topGO(GO, sort = "up", truncate = 45))
  print(topGO(GO, sort = "down", truncate = 45))
  write.csv(topGO(GO, n = Inf), file = file.path(GENELIST_DIR, paste0("GO-", comparison, ".csv")), row.names = FALSE)
}