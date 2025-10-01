# Data loading & QC
scai_load_10x <- function(path, sample = basename(path)) {
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat required")
  mat <- Seurat::Read10X(path)
  obj <- Seurat::CreateSeuratObject(mat, project = sample, min.cells = 3, min.features = 200)
  obj
}

scai_qc_filter <- function(obj, mito_prefix = c("MT-", "mt-"), ribo_prefix = c("RPL", "RPS"),
                           min_features = NULL, max_features = NULL, max_mito = NULL, max_counts = NULL,
                           regress_vars = c("percent.mt")) {
  .check_seurat(obj)
  genes <- rownames(obj)
  mito_genes <- unique(unlist(lapply(mito_prefix, function(p) grep(p, genes, value = TRUE))))
  ribo_genes <- unique(unlist(lapply(ribo_prefix, function(p) grep(p, genes, value = TRUE))))
  obj[["percent.mt"]] <- Matrix::colSums(obj@assays$RNA@counts[mito_genes, , drop = FALSE]) / Matrix::colSums(obj@assays$RNA@counts) * 100
  obj[["percent.ribo"]] <- Matrix::colSums(obj@assays$RNA@counts[ribo_genes, , drop = FALSE]) / Matrix::colSums(obj@assays$RNA@counts) * 100

  df <- Seurat::FetchData(obj, vars = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
  med_feat <- stats::median(df$nFeature_RNA); mad_feat <- stats::mad(df$nFeature_RNA)
  med_cnts <- stats::median(df$nCount_RNA); mad_cnts <- stats::mad(df$nCount_RNA)
  med_mito <- stats::median(df$percent.mt); mad_mito <- stats::mad(df$percent.mt)

  min_features <- min_features %||% max(200, round(med_feat - 3 * mad_feat))
  max_features <- max_features %||% round(med_feat + 4 * mad_feat)
  max_counts   <- max_counts   %||% round(med_cnts + 4 * mad_cnts)
  max_mito     <- max_mito     %||% max(5, round(med_mito + 3 * mad_mito))

  keep <- with(df, nFeature_RNA >= min_features & nFeature_RNA <= max_features &
                    nCount_RNA <= max_counts & percent.mt <= max_mito)
  message(sprintf("QC filter: keeping %d / %d cells", sum(keep), nrow(df)))
  obj <- obj[, keep]
  obj@misc$scai_qc <- list(thresholds = list(min_features=min_features, max_features=max_features,
                                             max_counts=max_counts, max_mito=max_mito))
  obj
}

scai_qc_plot <- function(obj) {
  .check_seurat(obj)
  p1 <- Seurat::VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + ggplot2::theme_bw()
  p2 <- Seurat::FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt") + ggplot2::theme_bw()
  list(qc_violin = p1, counts_vs_mito = p2)
}
