# Marker discovery
scai_markers <- function(obj, logfc = 0.25, min.pct = 0.1, test = "wilcox", top_n = 20) {
  .check_seurat(obj)
  markers <- Seurat::FindAllMarkers(obj, only.pos = TRUE, logfc.threshold = logfc, min.pct = min.pct, test.use = test)
  top <- markers |>
    dplyr::group_by(cluster) |>
    dplyr::arrange(dplyr::desc(avg_log2FC), p_val_adj) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::ungroup()
  obj@misc$scai_markers <- top
  top
}
