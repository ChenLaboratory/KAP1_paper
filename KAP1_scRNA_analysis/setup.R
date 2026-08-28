## setup.R
##
## Sourced before every other script. Defines packages, the annotation
## directory, colour palettes, marker gene panels, the Basal/LP/ML
## signature gene sets, and the AvgGeneSignature() scoring helper used
## throughout the pipeline.
##
## Edit `annotation_dir` to point at your own copy of the annotation files.

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyverse)
library(vcd)
library(edgeR)
library(org.Mm.eg.db)
library(purrr)
library(grid)

annotation_dir <- "Annotation"  # path to annotation files

## ---- Colour palettes -------------------------------------------------
col.pMedium <- c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D", "#AD8BC9",
                  "#A8786E", "#ED97CA", "#A2A2A2", "#CDCC5D", "#6DCCDA")
col.pDark <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
               "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF")
col.pLight <- c("#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
                "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5")
col.p <- col.pMedium
col.p2 <- c(col.pDark, col.pLight, col.pMedium)

## ---- Marker gene panels ------------------------------------------------
mouse_sig_genes <- list(
  Epithelial  = c("Epcam", "Itga6"),
  Basal       = c("Krt14", "Acta2"),
  LP          = c("Elf5", "Plet1"),
  ML          = c("Prlr", "Areg"),
  Cycling     = c("Mki67", "Hmgb2"),
  Fibroblast  = c("Dcn", "Lum", "Col1a2"),
  Stromal     = c("Igfbp7"),
  Endothelial = c("Pecam1", "Vwf"),
  Macrophage  = c("C1qa", "Cd68", "Cd3d", "Cd4")
)
known_markers <- unique(unlist(mouse_sig_genes))

marker_table <- data.frame(
  Cell_type    = names(mouse_sig_genes),
  Marker_genes = vapply(mouse_sig_genes, paste, collapse = ",", FUN.VALUE = character(1))
)
MmGeneAnnotation <- read.csv(file.path(annotation_dir, "MouseSignatureGenes_CellTypes.csv"))

## ---- Basal/LP/ML signature gene sets (from PosSigGenes.RData: MS2, LP2, ML2) --
load(file.path(annotation_dir, "PosSigGenes.RData"))
MouseSig <- list(MS2, LP2, ML2)
names(MouseSig) <- c("Basal", "LP", "ML")

## Genes of interest for KAP1/hormone signalling.
kap1_hormone_genes <- c("Trim28", "Esr1", "Pgr", "Prlr")

## ---- Signature scoring helper (RNA assay) ------------------------------
## Adds one metadata column per entry of `markers_ls`, holding each cell's
## mean expression across that gene set (RNA assay).
AvgGeneSignature <- function(data_seurat, markers_ls) {
  mean.ls <- markers_ls %>%
    purrr::map_dfc(~ colMeans(x = as.matrix(data_seurat@assays$RNA[.x, ]),
                               na.rm = TRUE))
  rownames(mean.ls) <- rownames(data_seurat@meta.data)
  AddMetaData(data_seurat, mean.ls, col.name = colnames(mean.ls))
}
