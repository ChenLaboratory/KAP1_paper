## 02_Integration.R
##
## Integrate the 4 QC'd samples from 01_Preprocessing_QC.R, cluster,
## annotate clusters to cell types, and generate the QC figures. 
##
## Input : RDS/data_seurat_<sample_id>.rds (from 01_Preprocessing_QC.R)
## Output: RDS/SeuratObj.rds
##         FindAllMarkers_Integrated.csv

source("setup.R")
source("config.R")

## ---- Load QC'd samples and integrate --------------------------------
seurat_list <- lapply(samples$sample_id, function(s) {
  obj <- readRDS(file.path(ROBJECT_DIR, paste0("data_seurat_", s, ".rds")))
  DefaultAssay(obj) <- "RNA"
  obj
})
names(seurat_list) <- samples$sample_id

anchors <- FindIntegrationAnchors(
  object.list = seurat_list, dims = 1:30, anchor.features = 2000,
  scale = TRUE, k.anchor = 5, k.filter = 200, k.score = 30, max.features = 200, verbose = TRUE
)
data_seurat <- IntegrateData(anchorset = anchors, dims = 1:30, k.weight = 100, new.assay.name = "integrated")

DefaultAssay(data_seurat) <- "integrated"
data_seurat <- ScaleData(data_seurat, verbose = FALSE)
data_seurat <- RunPCA(data_seurat, dims = 1:30, verbose = FALSE)
data_seurat <- RunUMAP(data_seurat, dims = 1:30, seed.use = 2025, reduction.name = "umap_seed2025")
data_seurat <- FindNeighbors(data_seurat, dims = 1:30, verbose = FALSE)
data_seurat <- FindClusters(data_seurat, resolution = 0.1, verbose = FALSE)
print(table(Idents(data_seurat)))

integration_title <- "Integrated_KAP1_scRNA"

print(DimPlot(data_seurat, reduction = "umap_seed2025", raster = FALSE, cols = col.p2) +
        ggtitle(integration_title) + theme(plot.title = element_text(hjust = 0.5)))
print(DimPlot(data_seurat, reduction = "umap_seed2025", raster = FALSE, cols = col.p2,
              split.by = "sample_name", ncol = 2) +
        ggtitle(integration_title) + theme(plot.title = element_text(hjust = 0.5)))

print(FeaturePlot(data_seurat, pt.size = 0.1, order = TRUE, reduction = "umap_seed2025",
                   features = "nCount_RNA") + ggtitle("nCount_RNA"))
print(FeaturePlot(data_seurat, pt.size = 0.1, order = TRUE, reduction = "umap_seed2025",
                   features = "nFeature_RNA") + ggtitle("Number of genes\n(nFeature_RNA)"))
print(FeaturePlot(data_seurat, pt.size = 0.1, order = TRUE, reduction = "umap_seed2025",
                   features = "percent.mt") + ggtitle("Mitochondrial percentage"))

## ---- Cluster markers -------------------------------------------------------
DefaultAssay(data_seurat) <- "RNA"
cluster_markers <- FindAllMarkers(data_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster_markers <- cluster_markers[cluster_markers$p_val_adj < 0.05, ]
write.csv(cluster_markers, file.path(".", "FindAllMarkers_Integrated.csv"))

top5 <- cluster_markers %>% dplyr::group_by(cluster) %>% dplyr::slice_head(n = 5) %>% dplyr::ungroup()
top3 <- cluster_markers %>% dplyr::group_by(cluster) %>% dplyr::slice_head(n = 3) %>% dplyr::ungroup()

genes5 <- unique(top5$gene)
print(DotPlot(data_seurat, features = genes5, dot.scale = 8) + coord_flip())

for (i in seq(1, nrow(top3), by = 12)) {
  idx <- i:min(i + 11, nrow(top3))
  print(FeaturePlot(data_seurat, features = top3[idx, ]$gene, order = TRUE, raster = FALSE,
                     pt.size = 0.4, ncol = 3, reduction = "umap_seed2025"))
}

## ---- Cell-cycle score UMAPs -------------------------------------------------
p_cc <- lapply(c("S.Score", "G2M.Score"), function(g) {
  FeaturePlot(data_seurat, pt.size = 0.1, order = TRUE, features = g, reduction = "umap_seed2025") + NoLegend()
})
print(patchwork::wrap_plots(p_cc, ncol = 2, byrow = TRUE))

## ---- Basal/LP/ML lineage signature scores ----------------------------------
DefaultAssay(data_seurat) <- "RNA"
MouseSig <- lapply(MouseSig, function(x) x[x %in% rownames(data_seurat)])
data_seurat <- AvgGeneSignature(data_seurat, MouseSig)

lineage_labels <- c("Basal", "LP", "ML")
p_sig <- lapply(seq_along(MouseSig), function(i) {
  FeaturePlot(data_seurat, pt.size = 0.1, order = TRUE, features = names(MouseSig)[i],
              reduction = "umap_seed2025") + NoLegend() + ggtitle(lineage_labels[i])
})
print(patchwork::wrap_plots(p_sig, ncol = 2, byrow = TRUE))

p_sig_vln <- lapply(seq_along(MouseSig), function(i) {
  VlnPlot(data_seurat, pt.size = 0, features = names(MouseSig)[i], col = col.p2) +
    labs(title = lineage_labels[i], x = NULL) + NoLegend()
})
print(patchwork::wrap_plots(p_sig_vln, ncol = 2, byrow = TRUE))

## ---- Curated gene panel (MmGeneAnnotation) UMAPs ---------------------------
MmGeneAnnotation_genes <- MmGeneAnnotation$Gene[MmGeneAnnotation$Gene %in% rownames(data_seurat)]
p_anno <- lapply(MmGeneAnnotation_genes, function(g) {
  FeaturePlot(data_seurat, pt.size = 0.1, order = TRUE, reduction = "umap_seed2025", features = g) + NoLegend()
})
for (i in seq(1, length(p_anno), by = 15)) {
  print(patchwork::wrap_plots(p_anno[i:min(i + 14, length(p_anno))], ncol = 3, byrow = TRUE))
}

## ---- Cell type annotation --------------------------------------
## Manual annotation based on the markers/signatures above.
cluster_annotations <- c(
  "0" = "Basal",
  "1" = "LP",
  "2" = "Cycling.LP",
  "3" = "ML",
  "4" = "Cycling.Basal",
  "5" = "LP",
  "6" = "Cluster6"
)

cluster_ids <- as.character(data_seurat$seurat_clusters)
missing_annotations <- setdiff(unique(cluster_ids), names(cluster_annotations))
if (length(missing_annotations) > 0) {
  stop("No cell-type annotation is defined for cluster(s): ", paste(missing_annotations, collapse = ", "))
}
data_seurat$cell_type <- factor(unname(cluster_annotations[cluster_ids]), levels = unique(cluster_annotations))

Idents(data_seurat) <- "cell_type"
d <- DimPlot(data_seurat, reduction = "umap_seed2025", group.by = "cell_type", raster = FALSE,
             label = FALSE, repel = TRUE, cols = celltype_cols) +
  ggtitle(integration_title) + theme(plot.title = element_text(hjust = 0.5))
print(LabelClusters(d, id = "cell_type", fontface = "bold", color = "black"))

print(table(data_seurat$cell_type))
print(DimPlot(data_seurat, reduction = "umap_seed2025", group.by = "cell_type", raster = FALSE,
              cols = col.p2, split.by = "sample_name", ncol = 2) +
        ggtitle(integration_title) + theme(plot.title = element_text(hjust = 0.5)))
print(table(data_seurat$cell_type, data_seurat$sample_name))

condition_map <- setNames(samples$condition, samples$sample_id)
data_seurat$condition <- unname(condition_map[data_seurat$sample_name])
if (anyNA(data_seurat$condition)) {
  stop("Condition is missing for one or more samples. Check the sample sheet in config.R.")
}
condition_cols <- c("KO" = "#BA094D", "WT" = "#303577")

print(DimPlot(data_seurat, reduction = "umap_seed2025", group.by = "condition", raster = FALSE,
              cols = condition_cols, label = FALSE, repel = TRUE) +
        ggtitle(integration_title) + theme(plot.title = element_text(hjust = 0.5)))

p_split <- DimPlot(data_seurat, reduction = "umap_seed2025", group.by = "condition", split.by = "condition",
                    raster = FALSE, cols = condition_cols, label = FALSE, repel = TRUE, ncol = 2) +
  ggtitle(integration_title) + theme(plot.title = element_text(hjust = 0.5))
print(patchwork::wrap_plots(p_split, nrow = 1, byrow = FALSE))

## ---- Ternary plots: Basal/LP/ML composition (all cells, WT only, KO only) --
plot_ternary <- function(obj, title) {
  TN <- matrix(0L, ncol(obj), 3L)
  colnames(TN) <- c("LP", "ML", "Basal")
  for (ct in colnames(TN)) TN[, ct] <- colSums(obj@assays$RNA$counts[MouseSig[[ct]], ] > 0L)

  cluster <- as.integer(obj$cell_type)
  ncls <- length(table(cluster))
  ternaryplot(TN, cex = 0.2, pch = 16, col = col.p2[cluster], grid = TRUE, main = title)
  grid_legend(0.9, 0.5, labels = levels(as.factor(obj$cell_type)), pch = 16, col = col.p2[1:ncls],
              vgap = unit(0.2, "lines"), title = "Cell Type", frame = FALSE, size = 0.5)
}

plot_ternary(data_seurat, "Ternary Plot")

data_seurat_WT <- subset(data_seurat, subset = sample_name %in% c("SK-K1", "SK-K2"))
data_seurat_KO <- subset(data_seurat, subset = sample_name %in% c("SK-K3", "SK-K4"))
plot_ternary(data_seurat_WT, "Ternary Plot for WT")
plot_ternary(data_seurat_KO, "Ternary Plot for KO")

## ---- Save ------------------------------------------------------------------
saveRDS(data_seurat, file = file.path(ROBJECT_DIR, "SeuratObj.rds"))
