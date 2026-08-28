## 03_ATAC_Peak_Calling_and_Integration.R
##
## Cell-type-aware ATAC peak calling (MACS2, run externally) + Harmony
## integration of the ATAC assay across samples, then transfer of RNA-derived
## cell type labels onto the ATAC cells, and peak-to-gene linking.
##
## Note: MACS2 peak calling itself is run outside R (one call per cell type,
## using the per-cell-type barcode lists exported below). This script covers
## everything else in R.
##
## Input : RObjects/Kap1_com_rna.rds (from 02_Integration.R)
##         RObjects/<sample_id>.rds (from 01_Preprocessing_QC.R)
##         atac_fragments.tsv.gz per sample
##         data/macs2/<cell_type>_peaks.narrowPeak (from external MACS2 run)
## Output: RObjects/comb_atac.rds, RObjects/WT_links.rds, RObjects/KO_links.rds

library(EnsDb.Mmusculus.v79)
library(GenomicRanges)
library(GenomicTools.fileHandler)
library(harmony)

source("setup.R")
source("config.R")

Kap1_com_rna <- readRDS(file.path(ROBJECT_DIR, "Kap1_com_rna.rds"))
seurat_list <- lapply(samples$sample_id, function(s) readRDS(file.path(ROBJECT_DIR, paste0(s, ".rds"))))
names(seurat_list) <- samples$sample_id
list2env(seurat_list, envir = environment())

## ---- Export per-cell-type barcode lists for external MACS2 peak calling --
cell_ids <- data.frame(
  umi      = stringr::str_sub(colnames(Kap1_com_rna), 1, -3),
  cellType = stringr::str_replace(Kap1_com_rna$cell_type, " ", "_"),
  stim     = Kap1_com_rna$stim
)

cellid_dir <- file.path(MACS2_DIR, "cellIDs")
dir.create(cellid_dir, showWarnings = FALSE, recursive = TRUE)
for (s in samples$sample_id) {
  exportBed(cell_ids[cell_ids$stim == s, c("umi", "cellType")],
            file = file.path(cellid_dir, paste0(s, "_cells.bed")), header = FALSE)
}

## ---- Load MACS2 peaks per cell type and build a consensus peak set -------
## Run MACS2 externally (per cell type, using the barcode lists above) before
## this step, writing narrowPeak files to MACS2_DIR/<cell_type>_peaks.narrowPeak.
cell_types <- c("Ba_Cyc", "Basal", "LP", "LP_Cyc", "Lum_Int", "Mixed_LP", "Mixed_ML", "ML", "ML_Cyc", "Stromal")

peaks_by_celltype <- lapply(cell_types, function(ct) {
  df <- read.table(file.path(MACS2_DIR, sprintf("%s_peaks.narrowPeak", ct)),
                    col.names = c("chr", "start", "end", "name", "score", "strand",
                                  "fold_change", "neg_log10pvalue_summit",
                                  "neg_log10qvalue_summit", "relative_summit_position"))
  gr <- makeGRangesFromDataFrame(df)
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- subsetByOverlaps(gr, blacklist_mm10, invert = TRUE)
  rtracklayer::export.bed(gr, file.path(MACS2_DIR, sprintf("%s_peaks.bed", ct)))
  gr
})
names(peaks_by_celltype) <- cell_types

combined_peaks <- GenomicRanges::reduce(x = do.call(c, unname(peaks_by_celltype)))
peakwidths <- width(combined_peaks)
combined_peaks <- combined_peaks[peakwidths < 10000 & peakwidths > 20]

## ---- Re-count fragments over the consensus peak set, per sample ----------
frags <- lapply(samples$sample_id, function(s) {
  CreateFragmentObject(
    path = file.path(RAW_DATA_DIR, s, "outs", "atac_fragments.tsv.gz"),
    cells = colnames(seurat_list[[s]])
  )
})
names(frags) <- samples$sample_id

atac_me_list <- lapply(samples$sample_id, function(s) {
  counts <- FeatureMatrix(fragments = frags[[s]], features = combined_peaks,
                           cells = colnames(seurat_list[[s]]))
  assay <- CreateChromatinAssay(counts, fragments = frags[[s]], genome = 'mm10')
  obj <- CreateSeuratObject(assay, assay = "ATAC")
  obj$stim <- s
  obj
})
names(atac_me_list) <- samples$sample_id

comb_atac <- merge(x = atac_me_list[[1]], y = atac_me_list[-1],
                    add.cell.ids = samples$sample_id)

## ---- Harmony integration of the ATAC assay --------------------------------
DefaultAssay(comb_atac) <- "ATAC"
comb_atac <- FindTopFeatures(comb_atac)
comb_atac <- RunTFIDF(comb_atac)
comb_atac <- RunSVD(comb_atac)
comb_atac <- RunHarmony(object = comb_atac, group.by.vars = 'stim',
                         reduction.use = 'lsi', assay.use = 'ATAC', project.dim = FALSE)
comb_atac <- RunTSNE(comb_atac, dims = 2:30, reduction = 'harmony')
comb_atac <- FindNeighbors(object = comb_atac, reduction = 'harmony', dims = 2:30)
comb_atac <- FindClusters(object = comb_atac, verbose = FALSE, algorithm = 3)

print(DimPlot(comb_atac, pt.size = 0.1, group.by = "stim") +
        ggtitle("Harmony integration") + NoAxes())

## ---- Transfer RNA-derived cell type labels onto ATAC cells ---------------
meta_atac <- data.frame(umi = colnames(comb_atac))
meta_rna <- data.frame(
  cols = colnames(Kap1_com_rna),
  umi = substr(colnames(Kap1_com_rna), 1, nchar(colnames(Kap1_com_rna)) - 2),
  cell_type = Kap1_com_rna$cell_type,
  stim = Kap1_com_rna$stim
)
meta_rna$umi <- paste0(meta_rna$stim, "_", meta_rna$umi)
meta_atac <- dplyr::left_join(meta_atac, meta_rna, by = "umi")
comb_atac$cell_type <- meta_atac$cell_type

Idents(comb_atac) <- "cell_type"
print(DimPlot(comb_atac, label = TRUE, cols = palette, reduction = "tsne", pt.size = 0.1) + NoLegend())

saveRDS(comb_atac, file = file.path(ROBJECT_DIR, "comb_atac.rds"))

cols_cells <- levels(as.factor(comb_atac$cell_type))
p_highlight <- lapply(cols_cells, function(ct) {
  cells_int <- WhichCells(comb_atac, idents = ct)
  DimPlot(comb_atac, label = TRUE, cells.highlight = cells_int,
          cols.highlight = "darkred", cols = "grey") + ggtitle(ct)
})
for (i in seq(1, length(p_highlight), by = 2)) {
  print(patchwork::wrap_plots(p_highlight[i:min(i + 1, length(p_highlight))], ncol = 2))
}

## ---- Peak annotation: nearest gene / promoter classification -------------
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
tss <- resize(annotations, width = 1, fix = 'start')

## nearestFeature() (defined in config.R) does the nearest-gene lookup below,
## following Signac's ClosestFeature approach.
Kap1_gr <- StringToGRanges(rownames(comb_atac))
Kap1_closest_gene <- nearestFeature(regions = Kap1_gr, object = tss)
Kap1_closest_gene$classification <- ifelse(Kap1_closest_gene$distance < 1500, "promoter", "distal")

Kap1_anno <- dplyr::left_join(data.frame(query_region = rownames(comb_atac)),
                               Kap1_closest_gene[, c("query_region", "classification")],
                               by = "query_region")

comb_atac_promoters <- comb_atac[rownames(comb_atac) %in%
                                    Kap1_anno[Kap1_anno$classification == "promoter", ]$query_region, ]

DefaultAssay(comb_atac_promoters) <- "ATAC"
comb_atac_promoters <- FindTopFeatures(comb_atac_promoters)
comb_atac_promoters <- RunTFIDF(comb_atac_promoters)
comb_atac_promoters <- RunSVD(comb_atac_promoters)
comb_atac_promoters <- RunHarmony(object = comb_atac_promoters, group.by.vars = 'stim',
                                   reduction.use = 'lsi', assay.use = 'ATAC', project.dim = FALSE)
comb_atac_promoters <- RunTSNE(comb_atac_promoters, dims = 2:30, reduction = 'harmony')

p1 <- DimPlot(comb_atac_promoters, pt.size = 0.1, group.by = "stim") +
  ggtitle("Harmony integration (promoter peaks only)") + NoAxes()
Idents(comb_atac_promoters) <- "cell_type"
p2 <- DimPlot(comb_atac_promoters, label = TRUE, cols = palette, reduction = "tsne", pt.size = 0.1) + NoLegend()
print(p1 + p2)

## ---- Add ATAC assay + embeddings onto the RNA object -----------------------
comb_atac <- RenameCells(comb_atac, new.names = meta_atac$cols)

Kap1_com_rna[["ATAC"]] <- comb_atac@assays$ATAC
Kap1_com_rna[["lsi"]] <- CreateDimReducObject(
  embeddings = comb_atac[['lsi']]@cell.embeddings, key = "LSI_", assay = 'ATAC')
Kap1_com_rna[["tsne.atac"]] <- CreateDimReducObject(
  embeddings = comb_atac[['tsne']]@cell.embeddings, key = "tSNE_", assay = 'ATAC')

DefaultAssay(Kap1_com_rna) <- "ATAC"
Annotation(Kap1_com_rna@assays$ATAC) <- annotations

## ---- Peak-to-gene linking --------------------------
DefaultAssay(Kap1_com_rna) <- "ATAC"
library(BSgenome.Mmusculus.UCSC.mm10)
Kap1_com_rna <- RegionStats(Kap1_com_rna, genome = BSgenome.Mmusculus.UCSC.mm10)
Annotation(Kap1_com_rna) <- annotations

Idents(Kap1_com_rna) <- "Condition"
link_peaks <- function(obj) {
  LinkPeaks(obj, peak.assay = "ATAC", expression.assay = "RNA", expression.slot = "data",
             distance = 1e6, min.cells = 10, n_sample = 200,
             pvalue_cutoff = 0.05, score_cutoff = 0.05, gene.id = FALSE, verbose = TRUE)
}

WT_links <- link_peaks(subset(Kap1_com_rna, idents = "WT"))
KO_links <- link_peaks(subset(Kap1_com_rna, idents = "KO"))
saveRDS(WT_links, file = file.path(ROBJECT_DIR, "WT_links.rds"))
saveRDS(KO_links, file = file.path(ROBJECT_DIR, "KO_links.rds"))
