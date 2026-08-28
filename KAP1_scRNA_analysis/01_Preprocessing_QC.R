## 01_Preprocessing_QC.R
##
## Per-sample QC for scRNA-seq: load CellRanger output, build a Seurat
## object, apply QC filters, cluster, and score cell-cycle / lineage
## signatures.
##
## Run this script once per sample: set `sample_name` below to the sample
## you're processing, then run top to bottom. The QC cutoffs and
## clustering resolution are picked automatically based on `sample_name`
##
## Input : CellRanger "outs/filtered_feature_bc_matrix/" for `sample_name`.
## Output: RDS/data_seurat_<sample_name>.rds

source("setup.R")
source("config.R")

# Set the sample to process. Run this script once per sample.
sample_name <- "SK-K1"

# Sample-specific QC cutoffs and clustering resolution.
n_mito_cutoff <- 10

if (sample_name %in% c("SK-K1", "SK-K2", "SK-K3")) {
  n_gene_cutoffs     <- c(500, 7500)
  n_gex_size_cutoffs <- c(1000, 50000)
} else if (sample_name == "SK-K4") {
  n_gene_cutoffs     <- c(500, 7000)
  n_gex_size_cutoffs <- c(1000, 40000)
} else {
  stop("Unrecognised sample_name '", sample_name, "': set n_gene_cutoffs / n_gex_size_cutoffs for it above.")
}

n_resolution <- switch(sample_name,
  "SK-K1" = 0.2,
  "SK-K2" = 0.3,
  "SK-K3" = 0.7,
  "SK-K4" = 0.7,
  stop("Unrecognised sample_name '", sample_name, "': set n_resolution for it above.")
)

n_var_genes <- 2000
n_dims <- 30

## ---- Load + QC ---------------------------------------------------------
cellranger_dir <- file.path(DATA_DIR, sample_name, "outs", "filtered_feature_bc_matrix")
if (!dir.exists(cellranger_dir)) {
  stop("CellRanger output not found: ", cellranger_dir)
}

rna_counts <- Read10X(cellranger_dir)

data_seurat_raw <- CreateSeuratObject(counts = rna_counts, min.cells = 3, assay = "RNA")
data_seurat_raw$sample_name <- sample_name
data_seurat_raw[["percent.mt"]] <- PercentageFeatureSet(data_seurat_raw, pattern = "^mt-")

print(VlnPlot(data_seurat_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              group.by = "sample_name", ncol = 3))

print(FeatureScatter(data_seurat_raw, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                      group.by = "percent.mt"))

n_genes_min    <- n_gene_cutoffs[1]
n_genes_max    <- n_gene_cutoffs[2]
n_lib_size_min <- n_gex_size_cutoffs[1]
n_lib_size_max <- n_gex_size_cutoffs[2]
n_mito_perc    <- n_mito_cutoff

data_seurat <- subset(
  x = data_seurat_raw,
  subset = nFeature_RNA > n_genes_min &
    nFeature_RNA < n_genes_max &
    nCount_RNA > n_lib_size_min &
    nCount_RNA < n_lib_size_max &
    percent.mt < n_mito_perc
)
message(sample_name, ": ", ncol(data_seurat), " / ", ncol(data_seurat_raw), " cells pass QC")

print(VlnPlot(data_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              group.by = "sample_name", ncol = 3))

## ---- Normalise, cluster --------------------------------------------------
DefaultAssay(data_seurat) <- "RNA"
data_seurat <- NormalizeData(data_seurat, verbose = FALSE)
data_seurat <- FindVariableFeatures(data_seurat, nfeatures = n_var_genes, verbose = FALSE)
data_seurat <- ScaleData(data_seurat, verbose = FALSE)
data_seurat <- RunPCA(data_seurat, verbose = FALSE)
data_seurat <- RunUMAP(data_seurat, dims = 1:n_dims, verbose = FALSE,
                        reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
data_seurat <- FindNeighbors(data_seurat, dims = 1:n_dims, verbose = FALSE)
data_seurat <- FindClusters(data_seurat, resolution = n_resolution, verbose = FALSE)

print(VariableFeaturePlot(data_seurat))
print(DimPlot(data_seurat, reduction = "umap.rna", label = TRUE))
print(FeaturePlot(data_seurat, reduction = "umap.rna", features = c("nCount_RNA", "nFeature_RNA")))
print(FeaturePlot(data_seurat, reduction = "umap.rna", features = "percent.mt"))

## ---- Cell-cycle scoring ----------------------------------------------------
DefaultAssay(data_seurat) <- "RNA"
s_genes    <- cc.genes$s.genes
g2m_genes  <- cc.genes$g2m.genes
m.s.genes   <- orthologs(genes = s_genes, species = "mouse")
m.g2m.genes <- orthologs(genes = g2m_genes, species = "mouse")

data_seurat <- CellCycleScoring(data_seurat, g2m.features = m.g2m.genes$symbol,
                                 s.features = m.s.genes$symbol)

print(DimPlot(data_seurat, reduction = "umap.rna", group.by = "Phase"))
print(FeaturePlot(data_seurat, reduction = "umap.rna", features = c("S.Score", "G2M.Score")))
print(FeaturePlot(data_seurat, reduction = "umap.rna",
                   features = kap1_hormone_genes[kap1_hormone_genes %in% rownames(data_seurat)]))
print(FeaturePlot(data_seurat, reduction = "umap.rna",
                   features = known_markers[known_markers %in% rownames(data_seurat)]))

## ---- Cluster markers -------------------------------------------------------
markers <- FindAllMarkers(data_seurat, assay = "RNA", only.pos = TRUE, logfc.threshold = 0.25)

top_markers <- dplyr::group_by(markers, cluster)
top_markers <- dplyr::top_n(top_markers, n = 3, wt = avg_log2FC)
print(as.data.frame(top_markers))

print(DotPlot(data_seurat, features = unique(top_markers$gene)) + Seurat::RotatedAxis())

## ---- Lineage signature scores ----------------------------------------------
mouse_sig_markers <- lapply(mouse_sig_genes, function(x) x[x %in% rownames(data_seurat)])
data_seurat <- AvgGeneSignature(data_seurat = data_seurat, markers_ls = mouse_sig_markers)

print(FeaturePlot(data_seurat, reduction = "umap.rna", features = names(mouse_sig_genes)))
print(VlnPlot(data_seurat, features = names(mouse_sig_genes), group.by = "seurat_clusters"))

## ---- Basal/LP/ML ternary plot ----------------------------------------------
IN2 <- matrix(0L, ncol(data_seurat), 3L)
colnames(IN2) <- c("Basal", "LP", "ML")
data_exp_mat <- data_seurat@assays$RNA@counts

Basal <- intersect(MouseSig$Basal, rownames(data_seurat))
LP    <- intersect(MouseSig$LP, rownames(data_seurat))
ML    <- intersect(MouseSig$ML, rownames(data_seurat))

IN2[, "Basal"] <- colSums(data_exp_mat[Basal, ] > 0L)
IN2[, "LP"]    <- colSums(data_exp_mat[LP, ] > 0L)
IN2[, "ML"]    <- colSums(data_exp_mat[ML, ] > 0L)

Cluster <- as.integer(data_seurat@meta.data$seurat_clusters)

plotOrd <- sample(ncol(data_seurat))
ncls <- length(table(Cluster))

vcd::ternaryplot(
  IN2[plotOrd, ],
  col = Cluster[plotOrd],
  pch = 19, cex = 0.3,
  main = paste0(sample_name, ": Basal / LP / ML ternary plot by cluster")
)

## ---- Save ------------------------------------------------------------------
saveRDS(data_seurat, file = file.path(ROBJECT_DIR, paste0("data_seurat_", sample_name, ".rds")))
