# Ambient RNA (optional, SoupX)
scai_ambient_correct <- function(obj, soup_fraction = NULL) {
  if (!requireNamespace("SoupX", quietly = TRUE)) {
    warning("SoupX not installed; skipping ambient RNA correction.")
    return(obj)
  }
  se <- obj
  sc <- SoupX::SoupChannel(tod = se@assays$RNA@counts)
  sc <- SoupX::estimateSoup(sc)
  if (is.null(soup_fraction)) soup_fraction <- sc$metaData$autoEstCont
  out <- SoupX::adjustCounts(sc, round = TRUE, contamination = soup_fraction)
  se@assays$RNA@counts <- out
  se
}

# End-to-end convenience pipeline
scai_pipeline <- function(obj,
                          batch_var = NULL,
                          qc_args = list(),
                          normalize_method = "sctransform",
                          integration = if (!is.null(batch_var)) "harmony" else "none",
                          seed = 777) {
  .check_seurat(obj); .set_seed(seed)
  obj <- do.call(scai_qc_filter, c(list(obj = obj), qc_args))
  obj <- scai_doublets(obj, batch_var = batch_var, seed = seed)
  obj <- scai_normalize(obj, method = normalize_method, seed = seed)
  if (!is.null(batch_var)) obj <- scai_integrate(obj, batch_var = batch_var, method = integration, seed = seed)
  obj <- scai_dimred_cluster(obj, use = if (!is.null(batch_var) && integration == "harmony") "harmony" else "pca", seed = seed)
  obj
}
