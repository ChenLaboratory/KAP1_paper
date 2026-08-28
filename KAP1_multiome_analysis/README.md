# Kap1 scMultiome Analysis Pipeline

Analysis workflow for the single-cell multiome dataset (paired scRNA-seq + scATAC-seq).

## Workflow

Scripts, run in order:

| Script | Description 
|---|---|
| `00_Install_Packages.R` | R package installation | — |
| `01_Preprocessing_QC.R` | Per-sample QC
| `02_Integration.R` | RNA integration, clustering and cell-type annotation 
| `03_ATAC_Peak_Calling_and_Integration.R` | MACS2 peak calling and ATAC integration, label transfer, peak-to-gene linking 
| `04_scATAC.R` | ATAC UMAP, gene-activity scores and marker gene exploration |
| `05_Pseudobulk_DE_DA.R` | Pseudobulk RNA DE and ATAC DA of KO vs WT, GO/KEGG analysis 
| `06_Pseudobulk_DE.R` | Pseudobulk RNA DE of cycling vs. parental epithelial populations

`setup.R` and `config.R` are sourced in Scripts 01–06 automatically:  
- `setup.R` defines packages, annotations, palettes and scoring helper functions;  
- `config.R` defines paths, sample metadata, and the pseudobulk and peak-annotation helper functions.

## Directory structure

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

source("00_Install_Packages.R")
source("01_Preprocessing_QC.R")
source("02_Integration.R")
source("03_ATAC_Peak_Calling_and_Integration.R")
source("04_scATAC.R")
source("05_Pseudobulk_DE_DA.R")
source("06_Pseudobulk_DE.R")
```
