# Batch integration
scai_integrate <- function(obj, batch_var, method = c("harmony", "none"), pcs = 50, seed = 777) {
  .check_seurat(obj); .set_seed(seed)
  method <- match.arg(method)
  hvgs <- Seurat::VariableFeatures(obj)
  obj <- Seurat::RunPCA(obj, features = hvgs, npcs = pcs, verbose = FALSE)
  if (method == "harmony") {
    if (!requireNamespace("harmony", quietly = TRUE)) stop("harmony required")
    obj <- harmony::RunHarmony(obj, group.by.vars = batch_var, dims.use = 1:pcs, assay.use = Seurat::DefaultAssay(obj))
    obj@reductions$harmony@misc <- list(batch_var = batch_var)
  }
  obj
}
