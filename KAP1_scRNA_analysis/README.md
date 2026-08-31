# Kap1 scRNA-seq Analysis Pipeline

Analysis workflow for the scRNA-seq dataset.

## Workflow

Scripts, run in order:

| Script | Description
|---|---|
| `00_Install_Packages.R` | R package installation 
| `01_Preprocessing_QC.R` | Per-sample QC
| `02_Integration.R` | Integration, clustering and cell-type annotation 
| `03_Pseudobulk.R` | Pseudobulk DE, GO/KEGG analysis

`setup.R` and `config.R` are sourced in Scripts 01–03 automatically:  
- `setup.R` defines packages, annotations, palettes and scoring helper functions;  
- `config.R` defines paths, sample metadata and the pseudobulk helper.

## Directory structure

The repository should contain the scripts and annotation directory. 

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

## Running the pipeline

From the repository root:

```r
setwd("path/to/this/repository")

source("00_Install_Packages.R")
source("01_Preprocessing_QC.R") # Run once per sample. Change sample_name before each run.
source("02_Integration.R")
source("03_Pseudobulk.R")
```
