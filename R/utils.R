# Utilities
scAI_install_and_load <- function() {
  pkgs <- c(
    "Seurat", "SeuratObject", "Matrix", "dplyr", "purrr", "tibble",
    "ggplot2", "rlang", "scater", "scran", "scDblFinder", "harmony",
    "uwot", "tidymodels", "xgboost"
  )
  opt_pkgs <- c("SoupX", "SingleCellExperiment", "batchelor")
  all <- unique(c(pkgs, opt_pkgs))
  inst <- rownames(installed.packages())
  to_install <- setdiff(all, inst)
  if (length(to_install) > 0) {
    message("Installing: ", paste(to_install, collapse = ", "))
    suppressWarnings(suppressMessages(
      install.packages(setdiff(to_install, c("Seurat", "SeuratObject")), repos = "https://cloud.r-project.org")
    ))
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    bioc_needed <- intersect(to_install, c("scater", "scran", "scDblFinder", "SingleCellExperiment", "batchelor"))
    if (length(bioc_needed) > 0) BiocManager::install(bioc_needed, ask = FALSE)
  }
  invisible(lapply(pkgs, require, character.only = TRUE))
}

.check_seurat <- function(obj) {
  if (!inherits(obj, "Seurat")) stop("Expected a Seurat object.")
}

.set_seed <- function(seed = 777) {
  set.seed(seed); try(uwot::set.seed(seed), silent = TRUE)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
