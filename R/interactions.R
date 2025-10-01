# Cell–cell interactions
scai_cc_interactions_cellchat <- function(obj, species = c("Human","Mouse"), group.by = "seurat_clusters", assay = NULL, ncores = 1) {
  .check_seurat(obj)
  species <- match.arg(species)
  if (!requireNamespace("CellChat", quietly = TRUE)) stop("CellChat not installed. Install it from Bioconductor/GitHub.")
  if (is.null(assay)) assay <- Seurat::DefaultAssay(obj)
  data.input <- Seurat::GetAssayData(obj, assay = assay, slot = "data")
  meta <- obj@meta.data
  meta$group <- meta[[group.by]]
  cellchat <- CellChat::createCellChat(object = data.input, meta = meta, group.by = "group")
  CellChat::DB <- switch(species,
                         Human = CellChat::CellChatDB.human,
                         Mouse = CellChat::CellChatDB.mouse)
  cellchat@DB <- CellChat::DB
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(cellchat, raw.use = FALSE)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  plt1 <- CellChat::netVisual_circle(cellchat@net$count, vertex.weight = table(meta$group), weight.scale = TRUE, label.edge = FALSE)
  plt2 <- CellChat::netAnalysis_signalingRole_scatter(cellchat)
  list(cellchat = cellchat, plots = list(circle = plt1, role_scatter = plt2), net = cellchat@net)
}

scai_cc_interactions_simple <- function(obj, lr_table, group.by = "seurat_clusters", assay = NULL) {
  .check_seurat(obj)
  if (is.null(assay)) assay <- Seurat::DefaultAssay(obj)
  mat <- Seurat::GetAssayData(obj, assay = assay, slot = "data")
  clusters <- obj@meta.data[[group.by]]
  avg <- Seurat::AverageExpression(obj, group.by = group.by, assays = assay, slot = "data")[[assay]]
  lig_in <- intersect(rownames(avg), unique(lr_table$ligand))
  rec_in <- intersect(rownames(avg), unique(lr_table$receptor))
  lr_in <- lr_table[lr_table$ligand %in% lig_in & lr_table$receptor %in% rec_in, , drop = FALSE]
  res <- purrr::map_dfr(seq_len(nrow(lr_in)), function(i) {
    lig <- lr_in$ligand[i]; rec <- lr_in$receptor[i]
    out_rows <- list()
    k <- 1
    for (s in colnames(avg)) for (r in colnames(avg)) {
      score <- avg[lig, s] * avg[rec, r]
      if (!is.na(score) && score > 0) {
        out_rows[[k]] <- tibble::tibble(sender = s, receiver = r, ligand = lig, receptor = rec, score = score); k <- k + 1
      }
    }
    if (length(out_rows) == 0) return(NULL) else dplyr::bind_rows(out_rows)
  })
  dplyr::arrange(res, dplyr::desc(score))
}
