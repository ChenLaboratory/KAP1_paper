## 06_Pseudobulk_DE.R
##
## Pseudobulk limma-voom differential expression restricted to the cycling
## vs non-cycling epithelial subclusters (Basal/LP/ML x Cyc/parental),
## with 9 pairwise comparisons: KO vs WT within each lineage's parental and
## cycling state, plus cycling KO vs cycling WT for each lineage. Each
## comparison gets a toptable, MD plot, and interactive Glimma MA plot;
## the cycling-vs-parental comparisons also get GO/KEGG over-representation
## (the KO-vs-WT-within-cycling comparisons skip GO/KEGG).
##
## Input : RObjects/Kap1_com_rna.rds (from 02_Integration.R)
##         data/annotation/Mus_musculus.gene_info.gz (NCBI gene info table)
## Output: GeneList/DE_<comparison>.csv, GeneList/KEGG-<comparison>.csv,
##         GeneList/GO-<comparison>.csv, GeneList/glimmaMAplot-<comparison>.html

library(limma)
library(Glimma)

source("setup.R")
source("config.R")  

data_seurat <- readRDS(file.path(ROBJECT_DIR, "Kap1_com_rna.rds"))
data_seurat$cell_type <- gsub(" ", "_", data_seurat$cell_type)
data_seurat <- subset(data_seurat, subset = cell_type %in%
                         c("Basal", "Ba_Cyc", "LP", "LP_Cyc", "ML", "ML_Cyc"))

## ---- Build pseudobulk samples ---------------------------------------------
y <- Seurat2PB(data_seurat, sample = "stim", cluster = "cell_type", assay = "RNA")
rownames(y$samples) <- sub("cluster", "", rownames(y$samples))
colnames(y$counts) <- sub("cluster", "", colnames(y$counts))

cluster_order <- c("Basal", "Ba_Cyc", "LP", "LP_Cyc", "ML", "ML_Cyc")
sample_order <- unlist(lapply(c("^KO1", "^KO2", "^WT1", "^WT2"), function(cond) {
  rows_cond <- grep(cond, rownames(y$samples), value = TRUE)
  unlist(lapply(paste0(cluster_order, "$"), function(ct) grep(ct, rows_cond, value = TRUE)))
}))

y$samples <- y$samples[sample_order, ]
y$counts <- y$counts[, sample_order]
y$samples$condition <- substr(sapply(strsplit(rownames(y$samples), "_"), `[`, 1), 1, 2)
y$samples$group <- paste(y$samples$condition, y$samples$cluster, sep = "_")
before_y_samples_group <- unique(y$samples$group)

celltype_cols <- c(Basal = col.p2[1], ML = col.p2[2], LP = col.p2[3],
                    ML_Cyc = col.p2[11], LP_Cyc = col.p2[12], Ba_Cyc = col.p2[13])
celltype_of_sample <- sapply(strsplit(rownames(y$samples), "_"), function(x) paste(x[-1], collapse = "_"))
bar_cols <- celltype_cols[celltype_of_sample]

plot_lib_sizes <- function(y, title_suffix) {
  par(mar = c(12, 5, 4, 2) + 0.1)
  bp <- barplot(y$samples$lib.size / 1e6, ylim = c(0, 35), col = bar_cols,
                ylab = "Library Size (in millions)",
                main = paste("Pseudobulk Library Sizes (in millions)", title_suffix),
                cex.main = 2, xaxt = "n", yaxt = "n")
  axis(1, at = bp, labels = rownames(y$samples), las = 2, cex.axis = 0.8)
  axis(2, at = seq(0, 35, by = 1), las = 2)

  df_lib <- data.frame(Sample = rownames(y$samples), LibSize = y$samples$lib.size / 1e6,
                        Cluster = factor(y$samples$cluster, levels = cluster_order))
  print(ggplot(df_lib, aes(x = Cluster, y = LibSize)) +
          geom_violin(aes(fill = Cluster), alpha = 0.6, trim = TRUE) +
          geom_boxplot(aes(color = Cluster), width = 0.1, outlier.shape = NA, alpha = 0.8) +
          scale_fill_manual(values = celltype_cols) + scale_color_manual(values = celltype_cols) +
          theme_minimal(base_size = 14) +
          labs(title = paste("Pseudobulk library sizes (in millions)", title_suffix),
               x = "Subcluster", y = "Library size (in millions, log10 scale)") +
          scale_y_log10(breaks = c(0.01, 0.1, 1, 5, 10, 20, 35), labels = scales::comma))
}
plot_lib_sizes(y, "BEFORE filtering")

print(table(data_seurat$orig.ident, data_seurat$cell_type)[, cluster_order])
saveRDS(y, file = file.path(ROBJECT_DIR, "y_RNA_beforeFiltering.rds"))

## ---- Gene annotation + filtering ------------------------------------------
gene_info <- read.delim(file.path(annotation_dir, "Mus_musculus.gene_info.gz"), header = TRUE)[, -1]
m <- match(rownames(y), gene_info$Symbol)
genes <- cbind(y$genes$gene, gene_info[m, c(2, 6, 9)])
colnames(genes) <- c("gene", "Symbol", "Chr", "Type")

keep_annotated <- !is.na(genes$Symbol) & !genes$Chr %in% c("mt", "Un") &
  !genes$Type %in% c("pseudo", "unknown", "other")
y$genes <- genes
y <- y[keep_annotated, , keep.lib.size = FALSE]

y$genes$Entrez <- mapIds(org.Mm.eg.db, keys = y$genes$Symbol, column = "ENTREZID",
                          keytype = "SYMBOL", multiVals = "first")

y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
after_y_samples_group <- unique(y$samples$group)
message("Groups dropped by filtering: ", paste(setdiff(before_y_samples_group, after_y_samples_group), collapse = ", "))

plot_lib_sizes(y, "AFTER filtering")

y <- normLibSizes(y)
saveRDS(y, file = file.path(ROBJECT_DIR, "y_RNA_afterFiltering.rds"))

## ---- MDS plots (overall + per-lineage) ------------------------------------
group <- factor(y$samples$group, levels = paste(rep(c("KO", "WT"), each = 6), cluster_order, sep = "_"))
cell_type_f <- factor(y$samples$cluster, levels = cluster_order)
condition <- factor(y$samples$condition, levels = c("KO", "WT"))
donor <- factor(y$samples$sample, levels = c("KO1", "KO2", "WT1", "WT2"))

mds <- plotMDS(y, plot = FALSE)
par(mfrow = c(2, 2))
plotMDS(mds, col = col.p2[group], main = "MDS by group"); legend("top", legend = levels(group), pch = 16, col = col.p2)
plotMDS(mds, col = celltype_cols[cell_type_f], main = "MDS by cell type"); legend("top", legend = levels(cell_type_f), pch = 16, col = celltype_cols)
plotMDS(mds, col = col.p2[condition], main = "MDS by condition"); legend("top", legend = levels(condition), pch = 16, col = col.p2)
plotMDS(mds, col = col.p2[11:length(col.p2)][donor], main = "MDS by donor")
legend("top", legend = levels(donor), pch = 16, col = col.p2[11:length(col.p2)], ncol = 2)

cell_type_groups <- list(ML = c("ML", "ML_Cyc"), Basal = c("Basal", "Ba_Cyc"), LP = c("LP", "LP_Cyc"))
for (ct in names(cell_type_groups)) {
  y_sub <- y[, y$samples$cluster %in% cell_type_groups[[ct]], keep.lib.sizes = TRUE]
  cond_sub <- factor(y_sub$samples$condition, levels = c("KO", "WT"))
  group_sub <- factor(y_sub$samples$group, levels = paste(rep(c("KO", "WT"), each = 6), cluster_order, sep = "_"))
  mds_sub <- plotMDS(y_sub, plot = FALSE)

  par(mfrow = c(1, 2))
  plotMDS(mds_sub, col = col.p2[cond_sub], main = "MDS by condition"); legend("topleft", legend = levels(cond_sub), pch = 16, col = col.p2)
  plotMDS(mds_sub, col = col.p2[group_sub], main = "MDS by group"); legend("topleft", legend = levels(group_sub), pch = 16, col = col.p2)
}

## ---- Design, contrasts, voom fit -------------------------------------------
group <- factor(y$samples$group)
design <- model.matrix(~0 + group)
colnames(design) <- gsub("group", "", colnames(design))

contrast.matrix <- makeContrasts(
  # KO: cycling vs parental
  KO_CycBa_vs_KO_Ba = (KO_Ba_Cyc - KO_Basal),
  KO_CycLP_vs_KO_LP = (KO_LP_Cyc - KO_LP),
  KO_CycML_vs_KO_ML = (KO_ML_Cyc - KO_ML),
  # WT: cycling vs parental
  WT_CycBa_vs_WT_Ba = (WT_Ba_Cyc - WT_Basal),
  WT_CycLP_vs_WT_LP = (WT_LP_Cyc - WT_LP),
  WT_CycML_vs_WT_ML = (WT_ML_Cyc - WT_ML),
  # cycling: KO vs WT
  KO_CycBa_vs_WT_CycBa = (KO_Ba_Cyc - WT_Ba_Cyc),
  KO_CycLP_vs_WT_CycLP = (KO_LP_Cyc - WT_LP_Cyc),
  KO_CycML_vs_WT_CycML = (KO_ML_Cyc - WT_ML_Cyc),
  levels = design
)

donor <- factor(y$samples$sample)
v <- voomLmFit(y, design, block = donor, sample.weights = TRUE, plot = TRUE, save.plot = TRUE)

par(mar = c(11, 4, 4, 2))
barplot(v$targets$sample.weight, names = rownames(v$targets), main = "Sample-specific weights",
        ylab = "Weight", xlab = "", las = 2, ylim = c(0, 1.5), yaxt = "n", col = bar_cols, cex.main = 3)
axis(side = 2, at = seq(0, 1.5, by = 0.5))
abline(h = 1, col = 2, lty = 2)
mtext("Sample", side = 1, line = 9, adj = 0.5)

cfit <- contrasts.fit(v, contrast.matrix)
efit <- eBayes(cfit)
save(v, cfit, efit, file = file.path(ROBJECT_DIR, "v_fit_efit.RData"))

print(plotSA(efit, main = "Final model: Mean-variance trend"))
dt <- decideTests(efit)
print(summary(dt))

## ---- Per-comparison DE + GO/KEGG -------------------------------------------
comparisons_full <- c("KO_CycBa_vs_KO_Ba", "KO_CycLP_vs_KO_LP", "KO_CycML_vs_KO_ML",
                       "WT_CycBa_vs_WT_Ba", "WT_CycLP_vs_WT_LP", "WT_CycML_vs_WT_ML")
comparisons_de_only <- c("KO_CycBa_vs_WT_CycBa", "KO_CycLP_vs_WT_CycLP", "KO_CycML_vs_WT_CycML")

run_de <- function(comparison) {
  tt <- topTable(efit, coef = comparison, sort.by = "P", number = Inf)
  print(head(tt, 5))
  write.csv(tt, file = file.path(GENELIST_DIR, paste0("DE_", comparison, ".csv")), row.names = TRUE)

  print(plotMD(efit, coef = comparison, status = dt[, comparison], legend = "topright", main = comparison))

  htmlwidgets::saveWidget(
    Glimma::glimmaMA(efit, coef = comparison, status = dt[, comparison], dge = y, main = comparison),
    file = file.path(GENELIST_DIR, paste0("glimmaMAplot-", comparison, ".html"))
  )
  tt
}

run_ora <- function(comparison) {
  keg <- kegga(efit, coef = comparison, species = "Mm", geneid = "Entrez", FDR = 0.05)
  print(topKEGG(keg, sort = "up", truncate = 45))
  print(topKEGG(keg, sort = "down", truncate = 45))
  write.csv(topKEGG(keg, n = Inf), file = file.path(GENELIST_DIR, paste0("KEGG-", comparison, ".csv")), row.names = FALSE)

  GO <- goana(efit, coef = comparison, species = "Mm", geneid = "Entrez", FDR = 0.05)
  print(topGO(GO, sort = "up", truncate = 45))
  print(topGO(GO, sort = "down", truncate = 45))
  write.csv(topGO(GO, n = Inf), file = file.path(GENELIST_DIR, paste0("GO-", comparison, ".csv")), row.names = FALSE)
}

for (comparison in comparisons_full) {
  run_de(comparison)
  run_ora(comparison)
}
for (comparison in comparisons_de_only) {
  run_de(comparison)
}
