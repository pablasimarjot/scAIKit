# Normalization & HVGs
scai_normalize <- function(obj, method = c("sctransform", "lognorm"), n_hvgs = 3000, seed = 777) {
  .check_seurat(obj); .set_seed(seed)
  method <- match.arg(method)
  if (method == "sctransform") {
    if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat required")
    obj <- Seurat::SCTransform(obj, vars.to.regress = intersect(c("percent.mt", "nCount_RNA"), colnames(obj@meta.data)),
                               verbose = FALSE)
    hvgs <- Seurat::VariableFeatures(obj)
  } else {
    obj <- Seurat::NormalizeData(obj)
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = n_hvgs)
    hvgs <- Seurat::VariableFeatures(obj)
    obj <- Seurat::ScaleData(obj, features = hvgs)
  }
  obj@misc$scai_hvgs <- hvgs
  obj
}
