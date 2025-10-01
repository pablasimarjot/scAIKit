# Auto-annotation via ML
scai_auto_annotate <- function(target_obj, ref_obj, label_col = "celltype", features = NULL, seed = 777) {
  .check_seurat(target_obj); .check_seurat(ref_obj); .set_seed(seed)
  if (!label_col %in% colnames(ref_obj@meta.data)) stop("label_col not found in ref_obj@meta.data")

  if (is.null(features)) features <- intersect(Seurat::VariableFeatures(ref_obj), rownames(target_obj))
  if (length(features) < 100) features <- head(intersect(Seurat::VariableFeatures(target_obj), rownames(ref_obj)), 2000)

  ref_mat <- t(as.matrix(Seurat::GetAssayData(ref_obj, slot = if (!is.null(ref_obj@assays$SCT)) "data" else "scale.data")[features, ]))
  tgt_mat <- t(as.matrix(Seurat::GetAssayData(target_obj, slot = if (!is.null(target_obj@assays$SCT)) "data" else "scale.data")[features, ]))

  ref_df <- as.data.frame(ref_mat)
  ref_df[[label_col]] <- factor(ref_obj@meta.data[[label_col]])

  tidymodels::tidymodels_prefer()
  set.seed(seed)
  rec <- recipes::recipe(as.formula(paste(label_col, "~ .")), data = ref_df) |>
    recipes::step_zv(all_predictors()) |>
    recipes::step_normalize(all_predictors())

  xgb_spec <- parsnip::boost_tree(trees = 1200, learn_rate = 0.05, tree_depth = 6, loss_reduction = 1e-3,
                                  sample_size = 0.8, mtry = floor(sqrt(ncol(ref_df) - 1))) |>
    parsnip::set_engine("xgboost") |>
    parsnip::set_mode("classification")

  wflow <- workflows::workflow() |>
    workflows::add_model(xgb_spec) |>
    workflows::add_recipe(rec)

  folds <- rsample::vfold_cv(ref_df, v = 5, strata = !!rlang::sym(label_col))
  fitrs <- tune::fit_resamples(wflow, resamples = folds, control = tune::control_resamples(save_pred = TRUE))

  final_fit <- parsnip::fit(wflow, data = ref_df)

  tgt_df <- as.data.frame(tgt_mat)
  preds <- predict(final_fit, new_data = tgt_df, type = "prob")
  cls   <- predict(final_fit, new_data = tgt_df, type = "class")$.pred_class
  maxp  <- do.call(pmax, as.list(preds))
  target_obj$scai_label <- cls
  target_obj$scai_label_maxprob <- maxp
  target_obj@misc$scai_autolabel <- list(features = features, n_features = length(features))
  target_obj
}
