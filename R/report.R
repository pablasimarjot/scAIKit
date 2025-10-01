# HTML report generator
scai_report <- function(obj, file = "scai_report.html", include_markers = TRUE, include_signatures = TRUE,
                        include_cellchat = FALSE, cellchat_species = "Human", group.by = "seurat_clusters",
                        signatures = NULL, max_markers = 50, seed = 777) {
  .check_seurat(obj); .set_seed(seed)
  if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("rmarkdown is required for reporting")
  if (!requireNamespace("knitr", quietly = TRUE)) stop("knitr is required for reporting")
  if (!requireNamespace("glue", quietly = TRUE)) stop("glue is required for reporting")

  tmp_rmd <- tempfile(fileext = ".Rmd")
  rmd <- glue::glue(
'---
title: "scAIkit Report"
output:
  html_document:
    toc: true
    toc_float: true
    df_print: paged
    self_contained: true
---

```{{r setup, include=FALSE}}
knitr::opts_chunk$set(message=FALSE, warning=FALSE, fig.width=8, fig.height=6)
library(Seurat); library(ggplot2)
obj <- scAI_obj
```

# Overview
- Cells: `r ncol(obj)`  
- Genes: `r nrow(obj)`  
- Assay: `r DefaultAssay(obj)`  

# QC
```{{r}}
plots <- scai_qc_plot(obj); plots$qc_violin
```
```{{r}}
plots$counts_vs_mito
```

# Embeddings & Clusters
```{{r}}
scai_umap_plot(obj, group.by = "{group.by}")
```

# Signature Scores
{if (include_signatures) "```{r}
if (!is.null(signature_list)) {
  obj <<- scai_score_signatures(obj, signature_list)
  print(scai_signature_dotplot(obj))
}
```" else ""}

# Markers
{if (include_markers) "```{r}
mk <- scai_markers(obj); head(mk, {max_markers})
```" else ""}

# Cell–Cell Interactions
{if (include_cellchat) "```{r}
if (requireNamespace('CellChat', quietly = TRUE)) {
  out <- scai_cc_interactions_cellchat(obj, species = '{cellchat_species}', group.by = '{group.by}')
  print(out$plots$circle)
  print(out$plots$role_scatter)
} else {
  cat('CellChat is not installed; skipping.')
}
```" else ""}
')
  writeLines(rmd, con = tmp_rmd)
  env <- new.env(parent = globalenv())
  env$scAI_obj <- obj
  env$signature_list <- signatures
  rmarkdown::render(tmp_rmd, output_file = file, envir = env, quiet = TRUE)
  message("Report written to ", normalizePath(file))
  invisible(file)
}
