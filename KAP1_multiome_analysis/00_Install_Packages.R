## 00_Install_Packages.R
##
## One-time setup: install every package used across setup.R - 06_Pseudobulk_DE.R.

# CRAN
install.packages(c(
  "Seurat",
  "SeuratObject",
  "Signac",
  "harmony",
  "GenomicTools.fileHandler",
  "ggplot2",
  "dplyr",
  "patchwork",
  "tidyverse",
  "vcd",
  "purrr",
  "babelgene",
  "broom",
  "scales",
  "htmlwidgets"
))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c(
  "edgeR",
  "limma",
  "Glimma",
  "org.Mm.eg.db",
  "scDblFinder",
  "GenomicRanges",
  "GenomicFeatures",
  "GenomeInfoDb",
  "EnsDb.Mmusculus.v79",
  "BSgenome.Mmusculus.UCSC.mm10",
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "annotate",
  "csaw",
  "rtracklayer"
))

# GitHub
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
remotes::install_github("immunogenomics/presto")  # wilcoxauc.Seurat(), used in 02_Integration.R

## Note: MACS2 (https://github.com/macs3-project/MACS2) is a separate,
## externally-run command-line tool used between the two halves of
## 03_ATAC_Peak_Calling_and_Integration.R - install it with pip/conda,
## not from R.
