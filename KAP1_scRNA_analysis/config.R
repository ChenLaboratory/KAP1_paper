## config.R
##
## Sourced (after setup.R) at the start of every script in this
## pipeline (01-03). Defines input/output directories, the sample sheet,
## the random seed, the cell-type colour map, and the Seurat2PB()
## pseudobulk helper.
##
## This file is sourced automatically by scripts 01-03.

if (!exists("annotation_dir")) source("setup.R")

set.seed(2025)

## ---- Directories -----------------------------------------------------
## Edit these to point at your own copies of the data / outputs.
## DATA_DIR     : per-sample CellRanger output (outs/filtered_feature_bc_matrix/).
## ROBJECT_DIR  : intermediate Seurat/DGEList .rds objects.
## GENELIST_DIR : per-comparison DE/GO/KEGG tables from script 03.
## (annotation_dir is set in setup.R.)
DATA_DIR     <- "CellRanger"
ROBJECT_DIR  <- "RDS"
GENELIST_DIR <- "GeneList"

for (d in c(ROBJECT_DIR, GENELIST_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

## ---- Sample sheet ------------------------------------------------------
## One row per sample. Condition labels are used by 02_Integration.R and
## 03_Pseudobulk.R. Sample-specific QC thresholds and clustering
## resolution are set in 01_Preprocessing_QC.R.
samples <- data.frame(
  sample_id = c("SK-K1", "SK-K2", "SK-K3", "SK-K4"),
  condition = c("WT", "WT", "KO", "KO"),
  stringsAsFactors = FALSE
)

## ---- Cell-type colour palette -------------------------------------------
## Assigned once clusters have been annotated to biological cell types
## (script 02).
celltype_cols <- c(
  "Basal"          = col.p2[1],
  "LP"             = col.p2[2],
  "Cycling.LP"     = col.p2[3],
  "ML"             = col.p2[4],
  "Cycling.Basal"  = col.p2[5],
  "Cluster6"       = col.p2[7]
)

## Short labels used for pseudobulk sample names / bar-plot colours in 03.
pb_cluster_labels <- c(
  "Basal" = "Ba", "Cycling.Basal" = "CycBa",
  "LP" = "LP", "Cycling.LP" = "CycLP",
  "ML" = "ML"
)
pb_cluster_cols <- c(
  Ba = "#67BF5C", CycBa = "#ED665D", LP = "#AD8BC9", CycLP = "#A8786E", ML = "#ED97CA"
)

## ---- Helper: pseudobulk aggregation -------------------------------------
## Aggregates a Seurat object's raw counts into per-sample x per-cluster
## pseudobulk profiles, returned as a DGEList. Used by script 03.
Seurat2PB <- function(object, sample, cluster = "seurat_clusters", assay = "RNA") {
  counts <- SeuratObject::GetAssayData(object, assay = assay, slot = "counts")
  meta <- object@meta.data
  sp <- meta[, sample]
  clst <- meta[, cluster]
  genes <- data.frame(gene = rownames(object[[assay]]))
  genes <- cbind(genes, object[[assay]][[]])
  sp_clst <- factor(paste(sp, clst, sep = "_cluster"))
  group_mat <- Matrix::sparse.model.matrix(~0 + sp_clst)
  colnames(group_mat) <- gsub("^sp_clst", "", colnames(group_mat))
  counts.pb <- counts %*% group_mat
  sp.pb <- gsub("_cluster.*$", "", levels(sp_clst))
  clst.pb <- gsub("^.*_cluster", "", levels(sp_clst))
  sample.pb <- data.frame(sample = sp.pb, cluster = clst.pb)
  DGEList(counts = as.matrix(counts.pb), samples = sample.pb, genes = genes)
}

## ---- Helper: bar colour lookup for pseudobulk sample plots --------------
## Maps a pseudobulk sample name (e.g. "SK-K1_Ba") to its cluster colour.
bar_colors_for <- function(sample_names) {
  suffix <- sub("^.*_", "", sample_names)
  unname(pb_cluster_cols[suffix])
}
