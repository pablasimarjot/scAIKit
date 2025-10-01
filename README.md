# scAIkit

AI-assisted toolkit for single-cell RNA-seq in R.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "scater",
    "scran",
    "scDblFinder",
    "SingleCellExperiment",
    "batchelor",
    "AUCell"
))

install.packages("remotes")
install.packages("devtools")

devtools::install_github("jinworks/CellChat")

remotes::install_github("pablasimarjot/scAIkit")
```

## Features

- Adaptive QC filtering
- Doublet detection (scDblFinder)
- Normalization & HVG selection
- Batch integration (Harmony, fastMNN)
- Clustering with resolution sweep
- Auto-annotation via ML (xgboost)
- Signature scoring (AddModuleScore / AUCell)
- Cell–cell interaction analysis (CellChat or LR scoring)
- Self-contained HTML report generator

## Example

```r
library(scAIkit)
obj <- scai_load_10x("path/to/10x_matrix")
obj <- scai_pipeline(obj, batch_var = "sample")
scai_umap_plot(obj)
scai_report(obj, file = "scAIkit_report.html")
```
