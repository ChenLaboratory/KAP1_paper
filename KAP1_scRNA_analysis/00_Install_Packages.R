## 00_Install_Packages.R
##
## One-time setup: install every package used across setup.R - 03_Pseudobulk.R.

# CRAN
install.packages(c(
  "Seurat",
  "ggplot2",
  "dplyr",
  "patchwork",
  "tidyverse",
  "vcd",
  "purrr"
))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "edgeR",
  "limma",
  "Glimma",
  "org.Mm.eg.db"
))
