# Dimensionality reduction & clustering
scai_dimred_cluster <- function(obj, use = c("harmony", "pca"), pcs = 50, resolutions = c(0.4, 0.6, 0.8, 1.0, 1.2),
                                neighbors_k = 30, seed = 777) {
  .check_seurat(obj); .set_seed(seed)
  use <- match.arg(use)
  dims_to_use <- 1:pcs
  if (use == "harmony" && !is.null(obj@reductions$harmony)) {
    obj <- Seurat::RunUMAP(obj, reduction = "harmony", dims = dims_to_use)
    obj <- Seurat::FindNeighbors(obj, reduction = "harmony", dims = dims_to_use, k.param = neighbors_k)
  } else {
    if (is.null(obj@reductions$pca)) obj <- Seurat::RunPCA(obj, npcs = pcs, verbose = FALSE)
    obj <- Seurat::RunUMAP(obj, reduction = "pca", dims = dims_to_use)
    obj <- Seurat::FindNeighbors(obj, reduction = "pca", dims = dims_to_use, k.param = neighbors_k)
  }
  res_scores <- purrr::map_dfr(resolutions, function(r) {
    cl <- Seurat::FindClusters(obj, resolution = r, verbose = FALSE)$seurat_clusters
    emb <- Seurat::Embeddings(obj, "umap")
    sil <- tryCatch({
      cluster::silhouette(as.integer(cl), dist(emb))
    }, error = function(e) NULL)
    mean_sil <- if (!is.null(sil)) mean(sil[, 3]) else NA_real_
    tibble::tibble(resolution = r, n_clusters = dplyr::n_distinct(cl), mean_silhouette = mean_sil)
  })
  best <- res_scores |> dplyr::arrange(dplyr::desc(mean_silhouette), n_clusters) |> dplyr::slice(1)
  best_r <- best$resolution[1]
  obj <- Seurat::FindClusters(obj, resolution = best_r, verbose = FALSE)
  obj@misc$scai_cluster_sweep <- res_scores
  obj@misc$scai_best_resolution <- best_r
  obj
}

scai_umap_plot <- function(obj, group.by = "seurat_clusters", label = TRUE) {
  .check_seurat(obj)
  Seurat::DimPlot(obj, reduction = "umap", group.by = group.by, label = label) + ggplot2::theme_bw()
}
