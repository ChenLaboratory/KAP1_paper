## 01_Preprocessing_QC.R
##
## Per-sample QC for a scMultiome (RNA + ATAC) dataset: load CellRanger ARC
## output, build a Seurat object with paired RNA/ATAC assays, apply QC
## filters, remove ATAC doublets (AMULET), and save one .rds per sample.
##
## Input : CellRanger ARC "outs/" folder per sample (filtered_feature_bc_matrix.h5,
##         atac_fragments.tsv.gz, per_barcode_metrics.csv), listed in `samples`
##         (config.R).
## Output: RObjects/<sample_id>.rds - one QC'd Seurat object per sample.

library(EnsDb.Mmusculus.v79)
library(GenomicRanges)

source("setup.R")
source("config.R")

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

## Sample-specific QC thresholds. These were tuned per sample based on the
## QC metric distributions (see plots below) - adjust for your own data.
qc_thresholds <- list(
  WT1 = list(atac_max = 15e4, atac_min = 1000, rna_max = 30000, rna_min = 500,
             feat_rna_min = 500, feat_rna_max = 7000, feat_atac_min = 1000, feat_atac_max = 70000),
  WT2 = list(atac_max = 15e4, atac_min = 1000, rna_max = 25000, rna_min = 500,
             feat_rna_min = 500, feat_rna_max = 7000, feat_atac_min = 1000, feat_atac_max = 70000),
  KO1 = list(atac_max = 12e4, atac_min = 1000, rna_max = 30000, rna_min = 500,
             feat_rna_min = 700, feat_rna_max = 7000, feat_atac_min = 500, feat_atac_max = 50000),
  KO2 = list(atac_max = 15e4, atac_min = 1000, rna_max = 30000, rna_min = 1000,
             feat_rna_min = 700, feat_rna_max = 7000, feat_atac_min = 500, feat_atac_max = 50000)
)

## Doublet-exclusion regions for AMULET (blacklist + sex/mito chromosomes).
otherChroms <- GRanges(c("chrM", "chrX", "chrY"), IRanges(1L, width = 10^8))
toExclude <- suppressWarnings(otherChroms)

process_sample <- function(sample_id) {

  message("Processing ", sample_id)
  th <- qc_thresholds[[sample_id]]

  h5_file      <- file.path(RAW_DATA_DIR, sample_id, "outs", "filtered_feature_bc_matrix.h5")
  frag_file    <- file.path(RAW_DATA_DIR, sample_id, "outs", "atac_fragments.tsv.gz")
  metrics_file <- file.path(RAW_DATA_DIR, sample_id, "outs", "per_barcode_metrics.csv")

  raw_data   <- Read10X_h5(h5_file)
  rna_counts <- raw_data$`Gene Expression`
  atac_counts <- raw_data$Peaks

  obj <- CreateSeuratObject(counts = rna_counts)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

  # keep ATAC peaks on standard chromosomes only
  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts <- atac_counts[as.vector(grange.use), ]

  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = 'mm10',
    fragments = frag_file,
    min.cells = 10,
    annotation = annotations
  )
  obj[["ATAC"]] <- chrom_assay
  obj$stim <- sample_id

  # Reads-in-peaks QC metrics used by the filter below, taken from CellRanger
  # ARC's own per-barcode metrics rather than recomputed in R.
  peak_fragments <- read.csv(metrics_file)
  colnames(peak_fragments)[1] <- "UMI"
  umis <- data.frame(UMI = colnames(obj))
  peak_frag <- left_join(umis, peak_fragments[, c("UMI", "atac_fragments", "atac_peak_region_fragments")],
                          by = "UMI")
  obj$peak_region_fragments <- peak_frag$atac_peak_region_fragments
  obj$atac_fragments <- peak_frag$atac_fragments
  obj$pct_reads_in_peaks <- obj$peak_region_fragments / obj$atac_fragments * 100

  print(VlnPlot(obj, features = c("nCount_ATAC", "nFeature_ATAC",
                                   "nCount_RNA", "nFeature_RNA", "percent.mt"),
                ncol = 5, log = TRUE, pt.size = 0) + NoLegend())

  qc.metrics <- obj@meta.data
  p1 <- qc.metrics %>% arrange(percent.mt) %>%
    ggplot(aes(nCount_RNA, nCount_ATAC, colour = percent.mt)) +
    geom_point(size = 1) +
    scale_color_gradientn(colors = c("black", "blue", "green2", "red", "yellow")) +
    ggtitle(paste(sample_id, "QC metrics (log scale)")) +
    theme_classic() + scale_x_log10() + scale_y_log10()

  p2 <- qc.metrics %>% arrange(percent.mt) %>%
    ggplot(aes(nFeature_RNA, nFeature_ATAC, colour = percent.mt)) +
    geom_point(size = 1) +
    scale_color_gradientn(colors = c("black", "blue", "green2", "red", "yellow")) +
    ggtitle(paste(sample_id, "QC metrics (log scale)")) +
    theme_classic() + scale_x_log10() + scale_y_log10()
  print(p1 + p2)

  DefaultAssay(obj) <- "ATAC"
  obj <- NucleosomeSignal(object = obj)
  obj <- TSSEnrichment(object = obj, fast = FALSE)
  print(VlnPlot(obj, features = c('TSS.enrichment', 'nucleosome_signal'), pt.size = 0.1, ncol = 2))

  # QC filter: fixed thresholds + distributional (mean +/- 2SD) thresholds
  meta <- obj[[]]
  obj <- subset(
    obj,
    subset = nCount_ATAC < th$atac_max & nCount_ATAC > th$atac_min &
      nCount_RNA < th$rna_max & nCount_RNA > th$rna_min &
      nFeature_RNA > th$feat_rna_min & nFeature_RNA < th$feat_rna_max &
      nFeature_ATAC > th$feat_atac_min & nFeature_ATAC < th$feat_atac_max &
      percent.mt < 10 &
      peak_region_fragments > max(0, mean(meta$peak_region_fragments) - 2 * sd(meta$peak_region_fragments)) &
      peak_region_fragments < mean(meta$peak_region_fragments) + 2 * sd(meta$peak_region_fragments) &
      pct_reads_in_peaks > max(0, mean(meta$pct_reads_in_peaks) - 2 * sd(meta$pct_reads_in_peaks)) &
      nucleosome_signal < mean(meta$nucleosome_signal) + 2 * sd(meta$nucleosome_signal) &
      TSS.enrichment > mean(meta$TSS.enrichment) - 2 * sd(meta$TSS.enrichment)
  )

  message(sample_id, ": ", ncol(obj), " cells pass fixed/distributional QC")

  # ATAC doublet detection (AMULET)
  amulet_res <- amulet(frag_file, regionsToExclude = toExclude)
  amulet_res$type <- ifelse(amulet_res$q.value < 0.01, "doublet", "singlet")
  amulet_res$umi <- rownames(amulet_res)

  cell_meta <- data.frame(umi = colnames(obj)) %>%
    left_join(amulet_res[, c("type", "umi")], by = "umi")
  obj$atac_db <- factor(cell_meta$type, levels = c("singlet", "doublet"))

  print(FeatureScatter(obj, feature1 = "nCount_ATAC", feature2 = "nFeature_ATAC",
                        group.by = "atac_db", cols = c("blue", "tan")) +
          labs(color = "db_type") + ggtitle(paste(sample_id, "ATAC-seq doublets")))
  print(table(obj$atac_db))

  obj <- subset(obj, subset = atac_db == "singlet")
  message(sample_id, ": ", ncol(obj), " singlet cells retained after doublet removal")

  saveRDS(obj, file = file.path(ROBJECT_DIR, paste0(sample_id, ".rds")))
  obj
}

## Run for every sample listed in the sample sheet.
seurat_list <- lapply(samples$sample_id, process_sample)
names(seurat_list) <- samples$sample_id
