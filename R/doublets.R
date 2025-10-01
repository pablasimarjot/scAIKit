# Doublet detection
scai_doublets <- function(obj, batch_var = NULL, seed = 777) {
  .check_seurat(obj); .set_seed(seed)
  if (!requireNamespace("scDblFinder", quietly = TRUE)) stop("scDblFinder required")
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment")
  sce <- Seurat::as.SingleCellExperiment(obj)
  if (!is.null(batch_var) && batch_var %in% colnames(Seurat::FetchData(obj, vars = batch_var))) {
    sce$batch <- Seurat::FetchData(obj, vars = batch_var)[,1]
  }
  sce <- scDblFinder::scDblFinder(sce, samples = if (!is.null(batch_var)) "batch" else NULL)
  obj$doublet_score <- sce$scDblFinder.score
  obj$doublet_call  <- sce$scDblFinder.class
  obj@misc$scai_doublets <- list(rate = mean(obj$doublet_call == "doublet"))
  obj
}
