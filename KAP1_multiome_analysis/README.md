# Kap1 scMultiome Analysis Pipeline

Analysis workflow for a single-cell multiome dataset (paired scRNA-seq + scATAC-seq) from mouse mammary gland (WT vs. Kap1/Trim28 KO). 

## Workflow

Scripts, run in order:

| Script | Description 
|---|---|
| `00_Install_Packages.R` | One-time R package installation | — |
| `01_Preprocessing_QC.R` | Per-sample QC (RNA + ATAC), doublet removal (AMULET) 
| `02_Integration.R` | RNA-side CCA integration, clustering and cell-type annotation 
| `03_ATAC_Peak_Calling_and_Integration.R` | Cell-type-aware MACS2 peak calling (external) and Harmony ATAC integration, label transfer, peak-to-gene linking 
| `04_scATAC.R` | ATAC UMAP, gene-activity scores and signature/marker exploration |
| `05_Pseudobulk_DE_DA.R` | Cell-type pseudobulk RNA DE (edgeR) and ATAC DA (csaw + edgeR), GO/KEGG analysis 
| `06_Pseudobulk_DE.R` | Pseudobulk limma-voom DE restricted to cycling vs. parental epithelial subclusters 

`setup.R` and `config.R` are sourced in Scripts 01–06 automatically:  
`setup.R` defines packages, annotations, palettes and scoring helper functions;  
`config.R` defines paths, sample metadata, and the pseudobulk and peak-annotation helper functions.

`03_ATAC_Peak_Calling_and_Integration.R` is intentionally split by an
**external MACS2 run**: the first half exports per-cell-type barcode lists,
MACS2 is then run outside R (one call per cell type) to produce
`narrowPeak` files in `data/macs2/`, and the second half of the script
picks up from there.

## Requirements

- R installation compatible with the package versions available from CRAN/Bioconductor
- CellRanger ARC output for each sample
- Annotation files 
- [MACS2](https://github.com/macs3-project/MACS2), run externally for cell-type-aware peak calling (script 03)

Run `00_Install_Packages.R` once before the first analysis.

## Directory structure

The repository should contain the scripts and annotation directory. Large
CellRanger/MACS2 outputs and generated RObjects/results should normally be
kept outside GitHub and referenced through `config.R`.

```text
<repository>/
├── 00_Install_Packages.R
├── setup.R
├── config.R
├── 01_Preprocessing_QC.R
├── 02_Integration.R
├── 03_ATAC_Peak_Calling_and_Integration.R
├── 04_scATAC.R
├── 05_Pseudobulk_DE_DA.R
├── 06_Pseudobulk_DE.R
└── Annotation/
    ├── PosSigGenes.RData
    ├── MouseSignatureGenes_CellTypes.csv
    └── Mus_musculus.gene_info.gz
    └── GSE227750_Mouse-SigGenes.RData
```

Set `RAW_DATA_DIR` and `MACS2_DIR` in `config.R` to the directories
containing the per-sample CellRanger ARC output and the MACS2 peak calls:

```text
<RAW_DATA_DIR>/<sample_id>/outs/filtered_feature_bc_matrix.h5
<RAW_DATA_DIR>/<sample_id>/outs/atac_fragments.tsv.gz
<RAW_DATA_DIR>/<sample_id>/outs/per_barcode_metrics.csv
<MACS2_DIR>/<cell_type>_peaks.narrowPeak
```

## Running the pipeline

From the repository root:

```r
setwd("path/to/this/repository")

# Run once
source("00_Install_Packages.R")

source("01_Preprocessing_QC.R")
source("02_Integration.R")

# 03 pauses partway through for an external MACS2 run using the barcode
# lists it exports; run MACS2, then continue sourcing the rest of the script.
source("03_ATAC_Peak_Calling_and_Integration.R")

source("04_scATAC.R")
source("05_Pseudobulk_DE_DA.R")
source("06_Pseudobulk_DE.R")
```
