# Signature scoring
scai_score_signatures <- function(obj, gene_sets, method = c("AddModuleScore","AUCell"), assay = NULL, seed = 777) {
  .check_seurat(obj); .set_seed(seed)
  method <- match.arg(method)
  if (is.null(assay)) assay <- Seurat::DefaultAssay(obj)
  if (method == "AddModuleScore") {
    obj <- Seurat::AddModuleScore(obj, features = gene_sets, assay = assay, name = paste0(names(gene_sets),"_score"))
    score_cols <- grep("_score[0-9]+$", colnames(obj@meta.data), value = TRUE)
    obj@misc$scai_signature_scores <- obj@meta.data[, score_cols, drop = FALSE]
  } else {
    if (!requireNamespace("AUCell", quietly = TRUE)) stop("AUCell not installed; install AUCell or use method='AddModuleScore'")
    expr <- as.matrix(Seurat::GetAssayData(obj, assay = assay, slot = if (!is.null(obj@assays$SCT) && assay == "SCT") "data" else "data"))
    cells_rankings <- AUCell::AUCell_buildRankings(expr, plotStats = FALSE)
    gs <- AUCell::GeneSet(names = names(gene_sets), geneIds = unname(gene_sets))
    cells_AUC <- AUCell::AUCell_calcAUC(geneSets = gs, cells_rankings)
    auc_mat <- as.data.frame(Matrix::t(AUCell::getAUC(cells_AUC)))
    obj@meta.data <- cbind(obj@meta.data, auc_mat[colnames(obj)])
    obj@misc$scai_signature_scores <- auc_mat
  }
  obj
}

scai_signature_dotplot <- function(obj, signatures = NULL, group.by = "seurat_clusters") {
  .check_seurat(obj)
  if (is.null(signatures)) {
    signatures <- grep("(_score[0-9]+$)|(^.+AUC$)", colnames(obj@meta.data), value = TRUE)
  }
  df <- obj@meta.data[, c(group.by, signatures), drop = FALSE]
  df[[group.by]] <- factor(df[[group.by]])
  tidy <- tidyr::pivot_longer(tibble::as_tibble(df), cols = tidyselect::all_of(signatures), names_to = "signature", values_to = "score")
  stats <- tidy |>
    dplyr::group_by(!!rlang::sym(group.by), signature) |>
    dplyr::summarise(mean = mean(score, na.rm = TRUE), frac = mean(score > stats::quantile(score, 0.75, na.rm = TRUE), na.rm = TRUE), .groups = "drop")
  ggplot2::ggplot(stats, ggplot2::aes(x = !!rlang::sym(group.by), y = signature, size = frac, color = mean)) +
    ggplot2::geom_point() + ggplot2::coord_flip() + ggplot2::theme_bw() + ggplot2::labs(size = ">75% frac", color = "mean score")
}
