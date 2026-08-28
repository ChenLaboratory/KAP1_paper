## 02_Integration.R
##
## Integrate the RNA assay across the 4 QC'd samples from 01_Preprocessing_QC.R,
## cluster, annotate clusters to cell types, and generate the main RNA-side
## QC/marker figures.
##
## Input : RObjects/<sample_id>.rds (from 01_Preprocessing_QC.R)
##         data/annotation/PosSigGenes.RData (Basal/LP/ML gene signature sets:
##           objects MS2, LP2, ML2)
## Output: RObjects/Kap1_com_rna.rds
##         Markers/Kap1_KO_markers.csv, Markers/Kap1_KO_CT_markers.csv

library(stringr)
library(grid)

source("setup.R")
source("config.R")

## ---- Load QC'd samples ---------------------------------------------------
seurat_list <- lapply(samples$sample_id, function(s) {
  readRDS(file.path(ROBJECT_DIR, paste0(s, ".rds")))
})
names(seurat_list) <- samples$sample_id
list2env(seurat_list, envir = environment())

seurat_list <- lapply(seurat_list, function(x) {
  DefaultAssay(x) <- "RNA"
  x
})

## ---- Seurat CCA integration -----------------------------------------------
seurat_list <- lapply(seurat_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
  x
})

features <- SelectIntegrationFeatures(object.list = seurat_list)
anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)
Kap1_com_rna <- IntegrateData(anchorset = anchors)

DefaultAssay(Kap1_com_rna) <- "RNA"
Kap1_com_rna <- ScaleData(Kap1_com_rna, verbose = FALSE)
Kap1_com_rna <- RunPCA(Kap1_com_rna, npcs = 30, verbose = FALSE)
Kap1_com_rna <- RunUMAP(Kap1_com_rna, reduction = "pca", dims = 1:30)
print(DimPlot(Kap1_com_rna, reduction = "umap", group.by = "stim") +
        ggtitle("RNA-seq Seurat integration"))

## ---- Clustering (integrated assay, tSNE used downstream) -----------------
DefaultAssay(Kap1_com_rna) <- "integrated"
Kap1_com_rna <- RunPCA(Kap1_com_rna, npcs = 30, verbose = FALSE)
Kap1_com_rna <- RunTSNE(Kap1_com_rna, reduction = "pca", dims = 1:30)
Kap1_com_rna <- FindNeighbors(Kap1_com_rna, dims = 1:30, verbose = FALSE)
Kap1_com_rna <- FindClusters(Kap1_com_rna, graph.name = "integrated_snn",
                              verbose = FALSE, resolution = 0.4)

ncls <- length(table(Kap1_com_rna$seurat_clusters))

print(DimPlot(Kap1_com_rna, reduction = 'tsne', label = TRUE, repel = TRUE,
              cols = col.p2[1:ncls]) + NoLegend())
print(DimPlot(Kap1_com_rna, split.by = "stim", reduction = 'tsne', ncol = 2,
              cols = col.p2[1:ncls]))

Kap1_com_rna$Condition <- ifelse(Kap1_com_rna$stim %in% c("WT1", "WT2"), "WT", "KO")

## ---- Cluster composition by sample / condition ----------------------------
cluster_composition <- function(group_var) {
  df <- as.data.frame(table(Kap1_com_rna$seurat_clusters, Kap1_com_rna[[group_var]][, 1]))
  colnames(df) <- c("Seurat_cluster", "Replicate", "Frequency")
  df$Seurat_cluster <- factor(df$Seurat_cluster)
  df <- df %>%
    group_by(Replicate) %>%
    mutate(Percentage = 100 * Frequency / sum(Frequency)) %>%
    ungroup()
  df
}

df_by_sample <- cluster_composition("stim")
df_by_condition <- cluster_composition("Condition")

p_sample <- ggplot(df_by_sample, aes(x = Replicate, y = Percentage, fill = Seurat_cluster)) +
  scale_fill_manual(values = col.p2[1:ncls], guide = "none") + geom_bar(stat = "identity") +
  labs(y = "Percentage (%)", x = "", fill = "Cluster") + theme_bw(base_size = 20)

p_condition <- ggplot(df_by_condition, aes(x = Replicate, y = Percentage, fill = Seurat_cluster)) +
  scale_fill_manual(values = col.p2[1:ncls]) + geom_bar(stat = "identity") +
  labs(y = "Percentage (%)", x = "", fill = "Cluster") + theme_bw(base_size = 20)
print(p_sample + p_condition)

## ---- Cluster markers (Wilcoxon via presto) --------------------------------
Idents(Kap1_com_rna) <- "seurat_clusters"
Kap1_KO_markers <- presto:::wilcoxauc.Seurat(X = Kap1_com_rna, group_by = 'seurat_clusters',
                                              assay = 'data', seurat_assay = 'RNA')
write.csv(Kap1_KO_markers, file = file.path(MARKER_DIR, "Kap1_KO_markers.csv"))

Kap1_KO_markers_sig <- Kap1_KO_markers[Kap1_KO_markers$padj < 0.05 &
                                          abs(Kap1_KO_markers$logFC) > 0.25 &
                                          Kap1_KO_markers$auc > 0.5, ]
top5 <- Kap1_KO_markers_sig %>% group_by(group) %>% top_n(n = 5, wt = logFC)

print(DotPlot(Kap1_com_rna, assay = "RNA", features = unique(top5$feature),
              group.by = "seurat_clusters") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)))

## ---- Cell-cycle scoring ----------------------------------------------------
m.s.genes   <- orthologs(genes = cc.genes$s.genes, species = "mouse")
m.g2m.genes <- orthologs(genes = cc.genes$g2m.genes, species = "mouse")
Kap1_com_rna <- CellCycleScoring(Kap1_com_rna, g2m.features = m.g2m.genes$symbol,
                                  s.features = m.s.genes$symbol)

print(DimPlot(Kap1_com_rna, reduction = "umap", group.by = "Phase",
              cols = c("#4D96FF", "Purple", "Red"), pt.size = 1))

## ---- Basal/LP/ML signature scores -----------------------------------------
DefaultAssay(Kap1_com_rna) <- "RNA"

# MouseSig (Basal/LP/ML gene sets) comes from setup.R (PosSigGenes.RData)
layer.ls <- lapply(MouseSig, function(x) intersect(x, rownames(Kap1_com_rna)))
Kap1_com_rna <- AvgGeneSignature(Kap1_com_rna, layer.ls)
Basal <- layer.ls$Basal; LP <- layer.ls$LP; ML <- layer.ls$ML

for (sig in names(layer.ls)) {
  p_feat <- FeaturePlot(Kap1_com_rna, features = sig, reduction = 'tsne', pt.size = 0.001) +
    theme(axis.text = element_blank(), axis.ticks = element_blank())
  p_vln <- VlnPlot(Kap1_com_rna, features = sig, cols = col.p2[1:ncls], pt.size = 0) + NoLegend()
  print(p_feat + p_vln)
}

for (i in seq(1, length(known_markers), by = 4)) {
  idx <- i:min(i + 3, length(known_markers))
  # NB: identities are still numeric seurat_clusters here (cell-type annotation
  # happens further below), so use the cluster palette rather than `palette`
  # (which is keyed by cell-type name and would be misapplied at this point).
  print(VlnPlot(Kap1_com_rna, features = known_markers[idx], cols = col.p2[1:ncls], pt.size = 0, ncol = 4))
}

## ---- Ternary plot: basal/LP/ML expression composition per cell -----------
IN2 <- matrix(0L, ncol(Kap1_com_rna), 3L)
colnames(IN2) <- c("Basal", "LP", "ML")
data_exp_mat <- Kap1_com_rna@assays$RNA@counts
IN2[, "Basal"] <- colSums(data_exp_mat[Basal, ] > 0L)
IN2[, "LP"]    <- colSums(data_exp_mat[LP, ] > 0L)
IN2[, "ML"]    <- colSums(data_exp_mat[ML, ] > 0L)

Cluster <- as.integer(Kap1_com_rna$seurat_clusters)
plotOrd <- sample(ncol(Kap1_com_rna))
ternaryplot(IN2[plotOrd, c(2, 3, 1)], cex = 0.3, pch = 16, col = col.p2[Cluster][plotOrd],
            grid = TRUE, main = "")
grid_legend(0.9, 0.5, labels = paste0(0:(ncls - 1), " - ", table(Cluster)),
            pch = 16, col = col.p2[1:ncls], vgap = unit(0.2, "lines"), title = "Cluster",
            frame = FALSE, size = 0.5)

## ---- Cluster -> cell type annotation --------------------------------------
## Manual annotation based on the markers/signatures above.
cluster_to_celltype <- c(
  '0' = 'Basal', '1' = 'ML', '2' = 'LP', '3' = 'Basal', '4' = 'LP',
  '5' = 'ML', '6' = 'Lum Int', '7' = 'ML', '8' = 'Mixed ML', '9' = 'Mixed LP',
  '10' = 'ML Cyc', '11' = 'LP Cyc', '12' = 'Ba Cyc', '13' = 'Stromal', '14' = 'Stromal'
)
Kap1_com_rna$cell_type <- cluster_to_celltype[as.character(Kap1_com_rna$seurat_clusters)]
Idents(Kap1_com_rna) <- "cell_type"

print(DimPlot(Kap1_com_rna, label = TRUE, cols = palette, reduction = "tsne", label.size = 3) +
        NoLegend())

saveRDS(Kap1_com_rna, file = file.path(ROBJECT_DIR, "Kap1_com_rna.rds"))

## ---- Cell-type markers -----------------------------------------------------
Kap1_KO_CT_markers <- presto:::wilcoxauc.Seurat(X = Kap1_com_rna, group_by = 'cell_type',
                                                 assay = 'data', seurat_assay = 'RNA')
write.csv(Kap1_KO_CT_markers, file = file.path(MARKER_DIR, "Kap1_KO_CT_markers.csv"))

Kap1_KO_CT_markers_sig <- Kap1_KO_CT_markers[Kap1_KO_CT_markers$padj < 0.05 &
                                                abs(Kap1_KO_CT_markers$logFC) > 0.25 &
                                                (Kap1_KO_CT_markers$pct_in > 10 | Kap1_KO_CT_markers$pct_out > 10), ]
top10 <- Kap1_KO_CT_markers_sig %>% group_by(group) %>% top_n(n = 10, wt = logFC)

print(DotPlot(Kap1_com_rna, assay = "RNA", features = unique(top10$feature)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              axis.text = element_text(size = 20)))

Kap1_com_rna <- ScaleData(Kap1_com_rna)
print(DoHeatmap(Kap1_com_rna, assay = "RNA", features = unique(top5$feature), size = 3) + NoLegend())

## ---- Cell-type composition by sample / condition --------------------------
celltype_composition <- function(group_var) {
  df <- as.data.frame(table(Kap1_com_rna$cell_type, Kap1_com_rna[[group_var]][, 1]))
  colnames(df) <- c("Cell_type", "Replicate", "Frequency")
  df %>% group_by(Replicate) %>% mutate(Percentage = 100 * Frequency / sum(Frequency)) %>% ungroup()
}

ct_by_sample <- celltype_composition("stim")
ct_by_condition <- celltype_composition("Condition")

p_ct_sample <- ggplot(ct_by_sample, aes(x = Replicate, y = Percentage, fill = Cell_type)) +
  scale_fill_manual(values = palette, guide = "none") + geom_bar(stat = "identity") +
  labs(y = "Percentage (%)", x = "", fill = "Cluster") + theme_bw(base_size = 20)

p_ct_condition <- ggplot(ct_by_condition, aes(x = Replicate, y = Percentage, fill = Cell_type)) +
  scale_fill_manual(values = palette) + geom_bar(stat = "identity") +
  labs(y = "Percentage (%)", x = "", fill = "Cluster") + theme_bw(base_size = 20)
print(p_ct_sample + p_ct_condition)
