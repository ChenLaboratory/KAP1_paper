## config.R
##
## Sourced (after setup.R) at the start of every script in this
## pipeline (01-06). setup.R defines packages, palettes, and marker
## gene panels; this file adds the pipeline-specific bits: input/output
## directories, the sample sheet, the random seed, the cell-type colour
## map, the Seurat2PB() pseudobulk helper, and the nearestFeature() peak
## annotation helper.
##
## Run:  source("setup.R"); source("config.R")  at the top of every script.

if (!exists("annotation_dir")) source("setup.R")

set.seed(1234)

## ---- Directories -----------------------------------------------------
## Edit these to point at your own copies of the data / outputs.
## RAW_DATA_DIR : per-sample CellRanger ARC output (filtered_feature_bc_matrix.h5,
##                 atac_fragments.tsv.gz) - one subfolder per sample.
## MACS2_DIR    : per-cell-type MACS2 narrowPeak / bed output from 03.
## ROBJECT_DIR  : intermediate Seurat/DGEList .rds objects.
## MARKER_DIR   : differential expression / accessibility result tables (.csv).
## GENELIST_DIR : per-comparison DE/GO/KEGG tables from script 06.
## FIGURE_DIR   : saved figures, if figures are written to disk rather than shown.
## (annotation_dir is set in setup.R.)
RAW_DATA_DIR <- "data/raw"
MACS2_DIR    <- "data/macs2"
ROBJECT_DIR  <- "RObjects"
MARKER_DIR   <- "Markers"
GENELIST_DIR <- "GeneList"
FIGURE_DIR   <- "figures"

for (d in c(ROBJECT_DIR, MARKER_DIR, GENELIST_DIR, FIGURE_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

## ---- Sample sheet ------------------------------------------------------
## One row per multiome sample. Update `frag_file` / `h5_file` to the
## real paths on your system (see 01_Preprocessing_QC.R).
samples <- data.frame(
  sample_id = c("WT1", "WT2", "KO1", "KO2"),
  condition = c("WT", "WT", "KO", "KO"),
  stringsAsFactors = FALSE
)

## ---- Cell-type colour palette -------------------------------------------
## Used from script 02 onwards, once clusters have been annotated to
## biological cell types. (Distinct from col.p/col.p2 in setup.R, which
## are cluster/sample-index palettes rather than named cell-type colours.)
palette <- c(
  'Basal'    = '#1F77B4',
  'Ba Cyc'   = '#FF7F0E',
  'LP'       = '#2CA02C',
  'LP Cyc'   = '#D62728',
  'Lum Int'  = '#9467BD',
  'Mixed LP' = '#8C564B',
  'Mixed ML' = '#E377C2',
  'ML'       = "#9EDAE5",
  'ML Cyc'   = '#BCBD22',
  'Stromal'  = '#7F7F7F'
)

## ---- Helper: pseudobulk aggregation -------------------------------------
## Aggregates a Seurat object's raw counts into per-sample x per-cluster
## pseudobulk profiles, returned as a DGEList. Used by scripts 05 and 06.
Seurat2PB <- function(object, sample, cluster = "seurat_clusters", assay) {
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

## ---- Helper: nearest-gene / TSS-distance annotation ---------------------
## Annotates each region in `regions` with its nearest feature in `object`
## (e.g. TSS positions), following Signac's ClosestFeature approach. Used by
## scripts 03 and 05 to classify ATAC peaks as promoter-proximal or distal.
nearestFeature <- function(regions, object) {
  nearest_feature <- distanceToNearest(x = regions, subject = object)
  feature_hits <- object[subjectHits(nearest_feature)]
  df <- as.data.frame(mcols(feature_hits))
  df$closest_region <- GRangesToString(feature_hits)
  df$query_region <- GRangesToString(regions)
  df$distance <- as.numeric(mcols(nearest_feature)$distance)
  df
}
