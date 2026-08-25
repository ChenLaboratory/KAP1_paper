# Kap1 scRNA-seq Analysis Pipeline

Analysis workflow for the scRNA-seq dataset from mouse mammary gland (WT vs. Kap1/Trim28 KO).  

## Workflow

Scripts, run in order:

| Script | Description | Main outputs |
|---|---|---|
| `00_Install_Packages.R` | One-time R package installation 
| `01_Preprocessing_QC.R` | Per-sample QC, normalisation, clustering and scoring 
| `02_Integration.R` | CCA integration, clustering and cell-type annotation 
| `03_Pseudobulk.R` | Cell-type pseudobulk limma-voom DE and GO/KEGG analysis

`setup.R` and `config.R` are sourced in Scripts 01–03 automatically:  
`setup.R` defines packages, annotations, palettes and scoring helper functions;  
`config.R` defines paths, sample metadata and the pseudobulk helper.

## Requirements

- R installation compatible with the package versions available from CRAN/Bioconductor
- CellRanger output for each sample
- Annotation files listed below
- Internet access for the initial package installation

Run `00_Install_Packages.R` once before the first analysis.

## Directory structure

The repository should contain the scripts and annotation directory. Large CellRanger outputs and generated RDS/results should normally be kept outside GitHub and referenced through `config.R`.

```text
<repository>/
├── 00_Install_Packages.R
├── setup.R
├── config.R
├── 01_Preprocessing_QC.R
├── 02_Integration.R
├── 03_Pseudobulk.R
└── Annotation/
    ├── PosSigGenes.RData
    ├── MouseSignatureGenes_CellTypes.csv
    └── Mus_musculus.gene_info.gz
    └── GSE227750_Mouse-SigGenes.RData
```

Set `DATA_DIR` in `config.R` to the directory containing the per-sample CellRanger outputs:

```text
<DATA_DIR>/<sample>/outs/filtered_feature_bc_matrix/
```

### Annotation files

The following files are provided:

- `PosSigGenes.RData` — contains the `MS2`, `LP2` and `ML2` signature gene sets used for Basal/LP/ML scoring.
- `MouseSignatureGenes_CellTypes.csv` — mouse cell-type marker annotation.
- `Mus_musculus.gene_info.gz` — NCBI mouse gene annotation used for gene filtering and Entrez ID mapping.
- `GSE227750_Mouse-SigGenes.RData` — contains the `Basal_up`, `LP_up` and `ML_up` signature gene sets used for ATAC gene-activity scoring.

## Running the pipeline

From the repository root:

```r
setwd("path/to/this/repository")

# Run once
source("00_Install_Packages.R")

# Run once per sample. Change sample_name in 01_Preprocessing_QC.R
# before each run.
source("01_Preprocessing_QC.R")

# Run after all four samples have been processed
source("02_Integration.R")

# Run after integration
source("03_Pseudobulk.R")
```