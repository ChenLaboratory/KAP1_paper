## 05_Pseudobulk_DE_DA.R
##
## Pseudobulk (cell-type level) differential expression (edgeR, RNA) and
## differential accessibility (csaw + edgeR, ATAC) between KO and WT, plus
## GO/KEGG over-representation and a DEG/DA-promoter overlap test.
##
## Input : RObjects/Kap1_com_rna.rds (from 02_Integration.R)
##         RObjects/comb_atac.rds (from 03_ATAC_Peak_Calling_and_Integration.R)
##         data/macs2/<cell_type>_peaks.narrowPeak (from 03)
##         atac_fragments.tsv.gz per sample
## Output: Markers/*.csv (DE/DA tables, GO/KEGG results)

library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(stringr)
library(csaw)
library(GenomicFeatures)  # for genes()/promoters() on the TxDb object below
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(annotate)
library(broom)

source("setup.R")
source("config.R")  # provides Seurat2PB(), nearestFeature(), palettes, directories

Kap1_com_rna <- readRDS(file.path(ROBJECT_DIR, "Kap1_com_rna.rds"))
comb_atac <- readRDS(file.path(ROBJECT_DIR, "comb_atac.rds"))

## =============================================================================
## Pseudobulk RNA-seq DE (per cell type, excluding Stromal)
## =============================================================================
Idents(Kap1_com_rna) <- "cell_type"
y_rna <- Seurat2PB(subset(Kap1_com_rna, idents = "Stromal", invert = TRUE),
                    sample = "stim", cluster = "cell_type", assay = "RNA")

ann <- AnnotationDbi::select(org.Mm.eg.db, keys = rownames(y_rna),
                              columns = c("GENETYPE", "ENTREZID"), keytype = "SYMBOL")
colnames(ann) <- c("Symbol", "Type", "GeneID")
ann <- dplyr::distinct(ann, Symbol, .keep_all = TRUE)
y_rna$genes <- dplyr::left_join(data.frame(Symbol = rownames(y_rna)), ann, by = "Symbol")

colnames(y_rna) <- str_remove(str_replace(colnames(y_rna), " ", "_"), "cluster")
y_rna$samples$group <- paste0(substr(colnames(y_rna), 1, 2), substr(colnames(y_rna), 4, nchar(colnames(y_rna))))

y_rna <- y_rna[filterByExpr(y_rna), , keep.lib.sizes = FALSE]
y_rna <- calcNormFactors(y_rna)

group <- as.factor(y_rna$samples$group)
design <- model.matrix(~0 + group)
colnames(design) <- substr(colnames(design), 6, nchar(colnames(design)))

y_rna <- estimateDisp(y_rna, design, robust = TRUE)
fit_rna <- glmQLFit(y_rna, design, robust = TRUE)

Contrasts <- makeContrasts(
  Ba_Cyc   = (KO_Ba_Cyc - WT_Ba_Cyc),
  Basal    = (KO_Basal - WT_Basal),
  LP       = (KO_LP - WT_LP),
  LP_Cyc   = (KO_LP_Cyc - WT_LP_Cyc),
  Lum_Int  = (KO_Lum_Int - WT_Lum_Int),
  Mixed_LP = (KO_Mixed_LP - WT_Mixed_LP),
  Mixed_ML = (KO_Mixed_ML - WT_Mixed_ML),
  ML       = (KO_ML - WT_ML),
  ML_Cyc   = (KO_ML_Cyc - WT_ML_Cyc),
  levels = design
)

## Full RNA-side DE + GO/KEGG for the three main lineages (Basal/LP/ML).
cell_type <- c("Basal", "LP", "ML")
celltype_DEG <- list()

for (i in seq_along(cell_type)) {
  foo <- stringr::str_replace(cell_type[i], " ", "_")
  message("Pseudobulk RNA-seq analysis of ", cell_type[i])

  qlf_rna <- glmQLFTest(fit_rna, contrast = Contrasts[, foo])
  KO_vs_WT_DE <- as.data.frame(topTags(qlf_rna, n = Inf))
  KO_vs_WT_DE$feature <- rownames(KO_vs_WT_DE)
  write.csv(KO_vs_WT_DE, file = file.path(MARKER_DIR, sprintf("%s_KO_vs_WT_DE.csv", foo)))

  print(plotMD(qlf_rna, main = paste(cell_type[i], "MD plot")))

  KO_vs_WT_DE_sig <- KO_vs_WT_DE[KO_vs_WT_DE$FDR < 0.05, ]
  KO_vs_WT_DE_sig$trend <- ifelse(KO_vs_WT_DE_sig$logFC > 0, "Up", "Down")
  write.csv(KO_vs_WT_DE_sig, file = file.path(MARKER_DIR, sprintf("%s_KO_vs_WT_DE_sig.csv", foo)))

  celltype_DEG[[cell_type[i]]] <- KO_vs_WT_DE_sig
}

## GO and KEGG over-representation of the KO-vs-WT DEGs, per cell type.
run_ora <- function(fun, top_fun, xlab_prefix, id_col) {
  for (i in seq_along(cell_type)) {
    foo <- stringr::str_replace(cell_type[i], " ", "_")
    sig <- celltype_DEG[[cell_type[i]]]
    list_DEGs <- list(Up = subset(sig, trend == "Up")[[id_col]], Down = subset(sig, trend == "Down")[[id_col]])

    res <- fun(list_DEGs, species = "Mm")
    write.csv(res, file = file.path(MARKER_DIR, sprintf("%s_%s.csv", foo, xlab_prefix)))

    top_down <- top_fun(res, sort = "down")
    top_down$U.perc <- top_down$Up / top_down$N
    top_down$D.perc <- top_down$Down / top_down$N
    term_col <- if ("Term" %in% colnames(top_down)) "Term" else "Pathway"

    p <- ggplot(top_down, aes(x = reorder(.data[[term_col]], U.perc), y = U.perc * 100, fill = -log(P.Up))) +
      theme_bw() + geom_col(position = "dodge") + scale_fill_gradient(low = "red", high = "blue") +
      xlab(paste(xlab_prefix, "pathways of up-regulated genes in KO")) + ylab("% of Genes in Pathway") +
      coord_flip() + theme(panel.grid = element_blank(), panel.border = element_blank())

    q <- ggplot(top_down, aes(x = reorder(.data[[term_col]], D.perc), y = D.perc * 100, fill = -log(P.Down))) +
      theme_bw() + geom_col(position = "dodge") + scale_fill_gradient(low = "red", high = "blue") +
      xlab(paste(xlab_prefix, "pathways of up-regulated genes in WT")) + ylab("% of Genes in Pathway") +
      coord_flip() + theme(panel.grid = element_blank(), panel.border = element_blank())

    print(patchwork::wrap_plots(p, q, nrow = 2) + patchwork::plot_annotation(paste(xlab_prefix, "KO vs WT in", cell_type[i])))
  }
}
run_ora(goana, topGO, "GO_rna", "GeneID")
run_ora(kegga, topKEGG, "KEGG_rna", "GeneID")

## =============================================================================
## Pseudobulk ATAC-seq DA (per cell type)
## =============================================================================
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
tss <- resize(annotations, width = 1, fix = 'start')

## nearestFeature() (defined in config.R) does the nearest-gene lookup below.

## Re-count fragments over per-cell-type MACS2 peaks, split by sample.
cell_type_all <- levels(as.factor(y_rna$samples$cluster))
cell_type_all <- setdiff(cell_type_all, "Stromal")

for (ct in cell_type_all) {
  foo <- stringr::str_replace(ct, " ", "_")

  gr_macs <- makeGRangesFromDataFrame(
    read.table(file.path(MACS2_DIR, sprintf("%s_peaks.narrowPeak", foo)),
               col.names = c("chr", "start", "end", "name", "score", "strand",
                             "fold_change", "neg_log10pvalue_summit",
                             "neg_log10qvalue_summit", "relative_summit_position")),
    keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE
  )

  counts_by_sample <- lapply(samples$sample_id, function(s) {
    frag <- CreateFragmentObject(
      path = file.path(RAW_DATA_DIR, s, "outs", "atac_fragments.tsv.gz"),
      cells = stringr::str_sub(colnames(Kap1_com_rna)[Kap1_com_rna$cell_type == ct & Kap1_com_rna$stim == s], 1, -3)
    )
    FeatureMatrix(fragments = frag, features = gr_macs,
                   cells = stringr::str_sub(colnames(Kap1_com_rna)[Kap1_com_rna$cell_type == ct & Kap1_com_rna$stim == s], 1, -3))
  })
  names(counts_by_sample) <- samples$sample_id

  atac_counts <- as.matrix(data.frame(lapply(counts_by_sample, rowSums)))
  saveRDS(atac_counts, file = file.path(ROBJECT_DIR, sprintf("%s_atac_counts.rds", foo)))
}

## Loess-normalised DA testing (Basal/LP/ML) with TSS-distance annotation.
cell_type <- c("Basal", "LP", "ML")

for (ct in cell_type) {
  foo <- stringr::str_replace(ct, " ", "_")
  message("Pseudobulk ATAC-seq analysis of ", ct)

  atac_counts <- readRDS(file.path(ROBJECT_DIR, sprintf("%s_atac_counts.rds", foo)))
  group <- samples$condition[match(colnames(atac_counts), samples$sample_id)]
  y_atac <- DGEList(atac_counts, group = group)

  gr <- StringToGRanges(rownames(atac_counts))
  closest_gene <- nearestFeature(regions = gr, object = tss)
  closest_gene$classification <- ifelse(closest_gene$distance < 500, "500 bp TSS",
                                  ifelse(closest_gene$distance < 2500, "0.5-2.5 kb TSS",
                                  ifelse(closest_gene$distance < 10000, "2.5-10 kb TSS", "> 10 kb TSS")))
  y_atac$genes <- dplyr::left_join(data.frame(query_region = rownames(y_atac)), closest_gene, by = "query_region")

  y_atac <- y_atac[filterByExpr(y_atac), , keep.lib.sizes = FALSE]
  counts_sub <- y_atac$counts

  design <- model.matrix(~0 + group)
  colnames(design) <- substr(colnames(design), 6, nchar(colnames(design)))
  Contrasts_atac <- makeContrasts(KO_vs_WT = (KO - WT), levels = design)

  offset_loess <- normOffsets(counts_sub, se.out = FALSE)
  d_meta <- dplyr::left_join(data.frame(query_region = rownames(counts_sub)), y_atac$genes, by = "query_region")

  d_loess <- DGEList(counts_sub, group = group)
  d_loess$offset <- offset_loess
  d_loess$genes <- d_meta

  d_loess <- estimateDisp(d_loess, design)
  fit_loess <- glmQLFit(d_loess, design)

  qlfLoess <- glmQLFTest(fit_loess, contrast = Contrasts_atac[, "KO_vs_WT"])
  KO_vs_WT_DA <- as.data.frame(topTags(qlfLoess, n = Inf))
  write.csv(KO_vs_WT_DA, file = file.path(MARKER_DIR, sprintf("%s_KO_vs_WT_DA.csv", foo)))
  print(plotMD(qlfLoess, main = paste(ct, "MD plot")))

  # narrow the DA list with a minimum fold-change requirement
  tr <- glmTreat(fit_loess, contrast = Contrasts_atac[, "KO_vs_WT"], lfc = log2(1.2))
  print(plotMD(tr, main = paste(ct, "MD plot after glmTreat")))

  KO_vs_WT_DA_treat <- as.data.frame(topTags(tr, n = Inf))
  write.csv(KO_vs_WT_DA_treat, file = file.path(MARKER_DIR, sprintf("%s_KO_vs_WT_DA_treat.csv", foo)))

  KO_vs_WT_DA_sig <- KO_vs_WT_DA_treat[KO_vs_WT_DA_treat$FDR < 0.05, ]
  KO_vs_WT_DA_sig$Direction <- ifelse(KO_vs_WT_DA_sig$logFC > 0, "Up regulated", "Down regulated")
  write.csv(KO_vs_WT_DA_sig, file = file.path(MARKER_DIR, sprintf("%s_KO_vs_WT_DA_sig.csv", foo)))
}

## =============================================================================
## Overlap between DEGs (RNA) and DA promoters (ATAC): Fisher's exact test
## =============================================================================
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoters.gr <- promoters(genes(txdb), upstream = 1500, downstream = 500)
promoters.gr$symbol <- lookUp(promoters.gr$gene_id, 'org.Mm.eg.db', 'SYMBOL')

deg_da_overlap <- list()
for (ct in cell_type) {
  foo <- stringr::str_replace(ct, " ", "_")

  degs <- read.csv(file.path(MARKER_DIR, sprintf("%s_KO_vs_WT_DE.csv", foo)), row.names = 1)
  dapeaks <- read.csv(file.path(MARKER_DIR, sprintf("%s_KO_vs_WT_DA_treat.csv", foo)), row.names = 1)

  significant.degs <- rownames(filter(degs, FDR < 0.05))
  significant.degs.pr <- promoters.gr[promoters.gr$symbol %in% significant.degs]

  significant.peaks <- (dapeaks %>% filter(FDR < 0.05))$query_region
  significant.peaks <- as.data.frame(t(as.data.frame(strsplit(significant.peaks, '-'))))
  colnames(significant.peaks) <- c("seqnames", "start", "end")
  significant.peaks.gr <- GRanges(significant.peaks)

  deg.dap.ols <- findOverlaps(significant.degs.pr, significant.peaks.gr)
  n.deg.dap.ols <- length(unique(queryHits(deg.dap.ols)))
  n.deg.nondap.ols <- length(significant.degs) - n.deg.dap.ols

  non.degs <- rownames(filter(degs, FDR > 0.05))
  non.degs.pr <- promoters.gr[promoters.gr$symbol %in% non.degs]
  nondeg.dap.ols <- findOverlaps(non.degs.pr, significant.peaks.gr)
  n.nondeg.dap.ols <- length(unique(queryHits(nondeg.dap.ols)))
  n.nondeg.nondap.ols <- length(non.degs) - n.nondeg.dap.ols

  mtx <- data.frame("DA promoter" = c(n.deg.dap.ols, n.deg.nondap.ols),
                     "non-DA promoter" = c(n.nondeg.dap.ols, n.nondeg.nondap.ols),
                     row.names = c("DEG", "non-DEG"))
  message(ct); print(mtx)

  res <- tidy(fisher.test(mtx))
  res$celltype <- ct
  deg_da_overlap[[ct]] <- res
}
deg_da_overlap <- dplyr::bind_rows(deg_da_overlap)
write.csv(deg_da_overlap, file = file.path(MARKER_DIR, "DEG_DA_promoter_overlap_fisher.csv"), row.names = FALSE)
