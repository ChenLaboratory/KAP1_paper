## 04_scATAC.R
##
## scATAC-side UMAP, gene-activity scoring, and marker-signature exploration
## on the Harmony-integrated ATAC object from 03_ATAC_Peak_Calling_and_Integration.R.
##
## Input : RObjects/comb_atac.rds (from 03_ATAC_Peak_Calling_and_Integration.R)
##         data/annotation/GSE227750_Mouse-SigGenes.RData 
## Output: RObjects/Kap1_com_atac_wGeneActivity.rds

library(EnsDb.Mmusculus.v79)

source("setup.R")
source("config.R")

Kap1_com_atac <- readRDS(file.path(ROBJECT_DIR, "comb_atac.rds"))

## ---- UMAP on the Harmony-integrated ATAC assay ----------------------------
Kap1_com_atac <- NormalizeData(Kap1_com_atac, verbose = FALSE)
Kap1_com_atac <- RunUMAP(Kap1_com_atac, dims = 1:30, seed.use = 2024, verbose = FALSE, reduction = 'harmony')

print(DimPlot(Kap1_com_atac, pt.size = 0.1, reduction = "umap", group.by = "cell_type",
              shuffle = TRUE, label = TRUE, repel = TRUE, cols = col.p2) +
        labs(title = "Integrated_Kap1_scATAC"))

Kap1_com_atac$Condition <- ifelse(grepl("WT", Kap1_com_atac$stim), "WT",
                                   ifelse(grepl("KO", Kap1_com_atac$stim), "KO", NA))

condition_cols <- c("KO" = "#BA094D", "WT" = "#303577")

print(DimPlot(Kap1_com_atac, pt.size = 0.1, reduction = "umap", group.by = "Condition",
              cols = condition_cols) + labs(title = "Integrated_Kap1_scATAC"))

print(DimPlot(Kap1_com_atac, pt.size = 0.1, reduction = "umap", group.by = "Condition",
              split.by = "Condition", shuffle = TRUE, cols = condition_cols, ncol = 2) +
        labs(title = "Integrated_Kap1_scATAC"))

print(DimPlot(Kap1_com_atac, split.by = "stim", pt.size = 0.1, reduction = "umap",
              group.by = "cell_type", shuffle = TRUE, label = TRUE, repel = TRUE,
              cols = col.p2, ncol = 2) + labs(title = "Integrated_Kap1_scATAC"))

## ---- Gene activity scores ------
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
Annotation(Kap1_com_atac) <- annotations

gene.activities <- GeneActivity(Kap1_com_atac, extend.upstream = 5000, extend.downstream = 5000)
Kap1_com_atac[['ATAC_gene_activity']] <- CreateAssayObject(counts = gene.activities)
Kap1_com_atac <- NormalizeData(
  Kap1_com_atac, assay = 'ATAC_gene_activity', normalization.method = 'LogNormalize',
  scale.factor = median(Kap1_com_atac$nCount_ATAC)
)

DefaultAssay(Kap1_com_atac) <- "ATAC_gene_activity"
Kap1_com_atac <- RunUMAP(Kap1_com_atac, dims = 1:30, seed.use = 2024, verbose = FALSE, reduction = 'harmony')

saveRDS(Kap1_com_atac, file = file.path(ROBJECT_DIR, "Kap1_com_atac_wGeneActivity.rds"))

## ---- Basal/LP/ML gene-activity signature scores (GSE227750) --------------
load(file.path(annotation_dir, "GSE227750_Mouse-SigGenes.RData"))  
mouseSigGenes <- list(Basal_up = Basal_up, LP_up = LP_up, ML_up = ML_up)
mouseSigGenes <- lapply(mouseSigGenes, function(x) x[x %in% rownames(Kap1_com_atac)])

AvgGeneSignatureATAC <- function(data_seurat, markers_ls) {
  mean.ls <- markers_ls %>% map_dfc(~ colMeans(as.matrix(data_seurat@assays$ATAC_gene_activity@data[.x, ]), na.rm = TRUE))
  mean.ls <- as.data.frame(mean.ls)
  rownames(mean.ls) <- rownames(data_seurat@meta.data)
  AddMetaData(data_seurat, mean.ls, col.name = colnames(mean.ls))
}

Kap1_com_atac <- AvgGeneSignatureATAC(Kap1_com_atac, mouseSigGenes)

sig_labels <- c("Basal", "LP", "ML")
p_sig_score <- lapply(seq_along(mouseSigGenes), function(i) {
  FeaturePlot(Kap1_com_atac, pt.size = 0.1, order = TRUE, features = names(mouseSigGenes)[i]) +
    NoLegend() + ggtitle(sig_labels[i])
})
print(patchwork::wrap_plots(p_sig_score, ncol = 3, byrow = TRUE))

p_sig_vln <- lapply(seq_along(mouseSigGenes), function(i) {
  VlnPlot(Kap1_com_atac, pt.size = 0, features = names(mouseSigGenes)[i]) +
    labs(title = sig_labels[i], x = NULL) + NoLegend()
})
print(patchwork::wrap_plots(p_sig_vln, ncol = 3))
