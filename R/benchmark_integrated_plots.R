# ============================================================
# Integrated benchmark plots
# ============================================================
# Choose varying benchmark-condition dimensions from the grid, mirroring the
# single-dataset benchmark. Integrated condition grids are named lists, so the
# metric tables carry their *_name columns.
.bp_integrated_factorization_layout <- function(
    bench,
    df,
    x_var = "ng",
    facet_var = "confounding_name",
    average_over = "seed",
    split_vars = NULL
) {
  if (is.null(split_vars)) {
    grid_name_map <- c(
      seeds = "seed",
      ng_grid = "ng",
      confounding_grid = "confounding_name",
      batch_to_biology_ratio_grid = "batch_to_biology_ratio_name",
      feature_noise_shape_grid = "feature_noise_shape_name",
      batch_biology_overlap_fraction_grid = "batch_biology_overlap_fraction_name",
      time_imbalance_grid = "time_imbalance",
      family_grid = "family",
      sigma1_normalize_all_grid = "sigma1_normalize_all"
    )
    grid_keys <- intersect(names(grid_name_map), names(bench$grid))
    varying <- vapply(grid_keys, function(k) {
      x <- bench$grid[[k]]
      length(unique(x)) > 1L
    }, logical(1))
    grid_vars <- unname(grid_name_map[grid_keys[varying]])
    grid_vars <- intersect(grid_vars, names(df))
    split_vars <- setdiff(grid_vars, c(x_var, facet_var, average_over))
  } else {
    split_vars <- intersect(split_vars, names(df))
    split_vars <- split_vars[vapply(split_vars, function(v) {
      x <- df[[v]]
      length(unique(x[!is.na(x)])) > 1L
    }, logical(1))]
    split_vars <- setdiff(split_vars, c(x_var, facet_var, average_over))
  }

  list(
    x_var = x_var,
    facet_var = facet_var,
    split_vars = unique(split_vars)
  )
}

# Integrated plots preserve dataset by default and use the integrated
# metadata columns confounding_name and batch_to_biology_ratio_name.
# If dataset/module/etc. are included in average_over, they are removed
# from grouping, splitting, names, and titles.
# ============================================================

plot_integrated_factorization_metrics <- function(
    bench,
    metrics = c(
      "loading_feature_auc_signal_vs_non_signal",
      "loading_feature_auc_batch_vs_non_batch",
      "loading_feature_auc_signal_vs_batch_only",
      "reconstruction_cor_truth",
      "relative_reconstruction_error_truth"
    ),
    average_over = "seed",
    show_se = TRUE
) {
  df <- bench$factorization_metrics
  layout <- .bp_integrated_factorization_layout(bench, df, average_over = average_over)
  .bp_check_columns(
    df,
    c("method", "dataset", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance"),
    "bench$factorization_metrics"
  )

  .bp_make_metric_plots(
    df = df,
    metric_cols = metrics,
    id_cols = unique(c("method", "dataset", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance", "seed", layout$split_vars)),
    split_vars = unique(c("dataset", layout$split_vars)),
    facet_var = "confounding_name",
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Integrated factorization metrics"
  )
}

plot_integrated_module_auc <- function(
    bench,
    average_over = "seed",
    show_se = TRUE
) {
  df <- bench$factorization_module_auc
  layout <- .bp_integrated_factorization_layout(bench, df, average_over = average_over)
  .bp_check_columns(
    df,
    c("method", "dataset", "module", "module_type", "auc_module_vs_background", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance"),
    "bench$factorization_module_auc"
  )

  .bp_make_value_plots(
    df = df,
    value_col = "auc_module_vs_background",
    id_cols = unique(c("method", "dataset", "module", "module_type", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance", "seed", layout$split_vars)),
    split_vars = unique(c("dataset", "module", "module_type", layout$split_vars)),
    facet_var = "confounding_name",
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Integrated module AUC",
    ylab = "Module AUC vs background"
  )
}

plot_integrated_module_recall_summary <- function(
    bench,
    metric = "aurc",
    average_over = "seed",
    show_se = TRUE
) {
  df <- bench$factorization_module_recall_summary
  layout <- .bp_integrated_factorization_layout(bench, df, average_over = average_over)
  .bp_check_columns(
    df,
    c("method", "dataset", "module", "module_type", metric, "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance"),
    "bench$factorization_module_recall_summary"
  )

  .bp_make_value_plots(
    df = df,
    value_col = metric,
    id_cols = unique(c("method", "dataset", "module", "module_type", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance", "seed", layout$split_vars)),
    split_vars = unique(c("dataset", "module", "module_type", layout$split_vars)),
    facet_var = "confounding_name",
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Integrated module recall summary",
    ylab = metric
  )
}

plot_integrated_loading_enrichment <- function(
    bench,
    metric = "odds_ratio",
    average_over = "seed",
    show_se = TRUE
) {
  df <- bench$factorization_loading_enrichment
  layout <- .bp_integrated_factorization_layout(bench, df, average_over = average_over)
  .bp_check_columns(
    df,
    c("method", "dataset", "module", "module_type", metric, "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance"),
    "bench$factorization_loading_enrichment"
  )

  .bp_make_value_plots(
    df = df,
    value_col = metric,
    id_cols = unique(c("method", "dataset", "module", "module_type", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance", "seed", layout$split_vars)),
    split_vars = unique(c("dataset", "module", "module_type", layout$split_vars)),
    facet_var = "confounding_name",
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Integrated loading enrichment",
    ylab = metric
  )
}

plot_integrated_inference_metrics <- function(
    bench,
    metrics = c("inference_auc_signal_vs_null", "power", "type1_error", "empirical_fdr"),
    average_over = "seed",
    show_se = TRUE
) {
  df <- bench$inference_metrics
  layout <- .bp_integrated_factorization_layout(bench, df, average_over = average_over)
  .bp_check_columns(
    df,
    c("method", "dataset", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance"),
    "bench$inference_metrics"
  )

  .bp_make_metric_plots(
    df = df,
    metric_cols = metrics,
    id_cols = unique(c("method", "dataset", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance", "seed", layout$split_vars)),
    split_vars = unique(c("dataset", layout$split_vars)),
    facet_var = "confounding_name",
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Integrated inference metrics"
  )
}

plot_integrated_inference_module_enrichment <- function(
    bench,
    metric = "odds_ratio",
    average_over = "seed",
    show_se = TRUE,
    scope = NULL,
    component = NULL
) {
  df <- bench$inference_module_enrichment
  layout <- .bp_integrated_factorization_layout(bench, df, average_over = average_over)
  .bp_check_columns(
    df,
    c("method", "dataset", "module", "module_type", metric, "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance"),
    "bench$inference_module_enrichment"
  )

  if (!is.null(scope) && "scope" %in% names(df)) {
    df <- dplyr::filter(df, .data$scope %in% scope)
  }
  if (!is.null(component) && "component" %in% names(df)) {
    df <- dplyr::filter(df, .data$component %in% component)
  }

  .bp_make_value_plots(
    df = df,
    value_col = metric,
    id_cols = unique(c("method", "dataset", "module", "module_type", "scope", "component", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance", "seed", layout$split_vars)),
    split_vars = unique(c("dataset", "module", "module_type", "scope", "component", layout$split_vars)),
    facet_var = "confounding_name",
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Integrated inference module enrichment",
    ylab = metric
  )
}

plot_integrated_runtime <- function(
    bench,
    time_col = "time",
    average_over = "seed",
    show_se = TRUE
) {
  df <- bench$factorization_timing
  layout <- .bp_integrated_factorization_layout(bench, df, average_over = average_over)
  .bp_check_columns(df, c("method", time_col, "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance"), "bench$factorization_timing")

  .bp_make_value_plots(
    df = df,
    value_col = time_col,
    id_cols = unique(c("method", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance", "seed", layout$split_vars)),
    split_vars = unique(c(layout$split_vars)),
    facet_var = "confounding_name",
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Integrated runtime",
    ylab = time_col
  )
}

#' Plot integrated module-type AUC summaries.
plot_integrated_module_type_auc <- function(
    bench,
    average_over = "seed",
    show_se = TRUE
) {
  df <- bench$factorization_module_type_auc
  if (is.null(df) || nrow(df) == 0L) {
    df <- summarize_module_type_auc(bench$factorization_module_auc)
  }
  layout <- .bp_integrated_factorization_layout(bench, df, average_over = average_over)
  .bp_check_columns(
    df,
    c("method", "dataset", "model_type", "mean_auc", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance"),
    "bench$factorization_module_type_auc"
  )

  .bp_make_value_plots(
    df = df,
    value_col = "mean_auc",
    id_cols = unique(c("method", "dataset", "model_type", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance", "seed", layout$split_vars)),
    split_vars = unique(c("dataset", "model_type", layout$split_vars)),
    facet_var = "confounding_name",
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Integrated module-type recovery",
    ylab = "Mean module AUC vs background"
  )
}

#' Plot integrated module-type recall summaries.
plot_integrated_module_type_recall <- function(
    bench,
    average_over = "seed",
    show_se = TRUE
) {
  df <- bench$factorization_module_type_recall
  if (is.null(df) || nrow(df) == 0L) {
    df <- summarize_module_type_recall(bench$factorization_module_recall_summary)
  }
  layout <- .bp_integrated_factorization_layout(bench, df, average_over = average_over)
  .bp_check_columns(
    df,
    c("method", "dataset", "model_type", "mean_recall", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance"),
    "bench$factorization_module_type_recall"
  )

  .bp_make_value_plots(
    df = df,
    value_col = "mean_recall",
    id_cols = unique(c("method", "dataset", "model_type", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance", "seed", layout$split_vars)),
    split_vars = unique(c("dataset", "model_type", layout$split_vars)),
    facet_var = "confounding_name",
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Integrated module-type recall",
    ylab = "Mean module recall"
  )
}

plot_integrated_loading_module_type_enrichment <- function(
    bench,
    metric = "odds_ratio",
    average_over = "seed",
    show_se = TRUE
) {
  df <- bench$factorization_loading_module_type_enrichment
  layout <- .bp_integrated_factorization_layout(bench, df, average_over = average_over)
  .bp_check_columns(
    df,
    c("method", "dataset", "model_type", metric, "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance"),
    "bench$factorization_loading_module_type_enrichment"
  )

  .bp_make_value_plots(
    df = df,
    value_col = metric,
    id_cols = unique(c("method", "dataset", "model_type", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance", "seed", layout$split_vars)),
    split_vars = unique(c("dataset", "model_type", layout$split_vars)),
    facet_var = "confounding_name",
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Integrated loading module-type enrichment",
    ylab = metric
  )
}

plot_integrated_inference_module_auc <- function(
    bench,
    average_over = "seed",
    show_se = TRUE
) {
  df <- bench$inference_module_auc
  layout <- .bp_integrated_factorization_layout(bench, df, average_over = average_over)
  .bp_check_columns(
    df,
    c("method", "dataset", "module", "module_type", "auc_module_vs_background", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance"),
    "bench$inference_module_auc"
  )

  .bp_make_value_plots(
    df = df,
    value_col = "auc_module_vs_background",
    id_cols = unique(c("method", "dataset", "module", "module_type", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance", "seed", layout$split_vars)),
    split_vars = unique(c("dataset", "module", "module_type", layout$split_vars)),
    facet_var = "confounding_name",
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Integrated inference module AUC",
    ylab = "Inference module AUC vs background"
  )
}

plot_integrated_inference_module_recall_summary <- function(
    bench,
    metric = "aurc",
    average_over = "seed",
    show_se = TRUE
) {
  df <- bench$inference_module_recall_summary
  layout <- .bp_integrated_factorization_layout(bench, df, average_over = average_over)
  .bp_check_columns(
    df,
    c("method", "dataset", "module", "module_type", metric, "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance"),
    "bench$inference_module_recall_summary"
  )

  .bp_make_value_plots(
    df = df,
    value_col = metric,
    id_cols = unique(c("method", "dataset", "module", "module_type", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance", "seed", layout$split_vars)),
    split_vars = unique(c("dataset", "module", "module_type", layout$split_vars)),
    facet_var = "confounding_name",
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Integrated inference module recall summary",
    ylab = metric
  )
}

plot_integrated_inference_module_type_auc <- function(
    bench,
    average_over = "seed",
    show_se = TRUE
) {
  df <- bench$inference_module_type_auc
  if (is.null(df) || nrow(df) == 0L) {
    df <- summarize_module_type_auc(bench$inference_module_auc)
  }
  layout <- .bp_integrated_factorization_layout(bench, df, average_over = average_over)
  .bp_check_columns(
    df,
    c("method", "dataset", "model_type", "mean_auc", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance"),
    "bench$inference_module_type_auc"
  )

  .bp_make_value_plots(
    df = df,
    value_col = "mean_auc",
    id_cols = unique(c("method", "dataset", "model_type", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance", "seed", layout$split_vars)),
    split_vars = unique(c("dataset", "model_type", layout$split_vars)),
    facet_var = "confounding_name",
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Integrated inference module-type recovery",
    ylab = "Mean inference module AUC vs background"
  )
}

plot_integrated_inference_module_type_recall <- function(
    bench,
    average_over = "seed",
    show_se = TRUE
) {
  df <- bench$inference_module_type_recall
  if (is.null(df) || nrow(df) == 0L) {
    df <- summarize_module_type_recall(bench$inference_module_recall_summary)
  }
  layout <- .bp_integrated_factorization_layout(bench, df, average_over = average_over)
  .bp_check_columns(
    df,
    c("method", "dataset", "model_type", "mean_recall", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance"),
    "bench$inference_module_type_recall"
  )

  .bp_make_value_plots(
    df = df,
    value_col = "mean_recall",
    id_cols = unique(c("method", "dataset", "model_type", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance", "seed", layout$split_vars)),
    split_vars = unique(c("dataset", "model_type", layout$split_vars)),
    facet_var = "confounding_name",
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Integrated inference module-type recall",
    ylab = "Mean inference module recall"
  )
}

plot_integrated_inference_module_type_enrichment <- function(
    bench,
    metric = "odds_ratio",
    average_over = "seed",
    show_se = TRUE,
    scope = NULL,
    component = NULL
) {
  df <- bench$inference_module_type_enrichment
  layout <- .bp_integrated_factorization_layout(bench, df, average_over = average_over)
  .bp_check_columns(
    df,
    c("method", "dataset", "model_type", metric, "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance"),
    "bench$inference_module_type_enrichment"
  )

  if (!is.null(scope) && "scope" %in% names(df)) {
    df <- dplyr::filter(df, .data$scope %in% scope)
  }
  if (!is.null(component) && "component" %in% names(df)) {
    df <- dplyr::filter(df, .data$component %in% component)
  }

  .bp_make_value_plots(
    df = df,
    value_col = metric,
    id_cols = unique(c("method", "dataset", "model_type", "scope", "component", "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance", "seed", layout$split_vars)),
    split_vars = unique(c("dataset", "model_type", "scope", "component", layout$split_vars)),
    facet_var = "confounding_name",
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Integrated inference module-type enrichment",
    ylab = metric
  )
}

plot_integrated_inference_module_recall <- function(bench, ...) {
  plot_integrated_inference_module_recall_summary(bench, ...)
}

# Restrict an integrated plot to explicitly requested fitted method names.
# NULL preserves the complete benchmarked set.
.bp_filter_integrated_plot_methods <- function(df, methods = NULL) {
  if (is.null(methods)) return(df)
  if (!is.character(methods) || length(methods) == 0L ||
      anyNA(methods) || any(!nzchar(methods))) {
    stop("`methods` must be NULL or a nonempty character vector of method names.",
         call. = FALSE)
  }
  available <- unique(as.character(df$method))
  # Accept either an exact fitted method or a base method that expands to
  # its tuning-specific variants (e.g., concat_sPCA_keepX_*).
  matches <- lapply(methods, function(m) {
    if (m %in% available) return(m)
    available[startsWith(available, paste0(m, "_"))]
  })
  missing <- methods[lengths(matches) == 0L]
  if (length(missing)) {
    stop("Requested methods not found: ", paste(missing, collapse = ", "),
         ". Available methods: ", paste(available, collapse = ", "),
         call. = FALSE)
  }
  selected <- unique(unlist(matches, use.names = FALSE))
  df[df$method %in% selected, , drop = FALSE]
}

# ------------------------------------------------------------
# Additive metric plots (existing plot helpers only)
# ------------------------------------------------------------

plot_integrated_reconstruction_error <- function(
    bench, feature_set = "allFeatures", k_scope = "allK", module = NULL,
    average_over = "seed", show_se = TRUE, methods = NULL) {
  df <- bench$factorization_reconstruction_metrics
  layout <- .bp_integrated_factorization_layout(bench, df, average_over = average_over)
  .bp_check_columns(df, c("method", "dataset", "feature_set", "k_scope", "reconstruction_error",
                          "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance"),
                    "bench$factorization_reconstruction_metrics")
  df <- .bp_filter_integrated_plot_methods(df, methods)
  df <- dplyr::filter(df, .data$feature_set %in% .env$feature_set, .data$k_scope %in% .env$k_scope)
  if (!is.null(module)) {
    module_requested <- sub("^module([0-9]+)$", "module_\\1", as.character(module))
    available_modules <- sort(unique(as.character(df$module[!is.na(df$module)])))
    missing_modules <- setdiff(module_requested, available_modules)
    if (length(missing_modules) > 0L) {
      stop(
        "Requested module(s) not found: ", paste(module, collapse = ", "),
        ". Available modules: ", paste(available_modules, collapse = ", "),
        call. = FALSE
      )
    }
    df <- dplyr::filter(df, .data$module %in% .env$module_requested)
  }
  df <- dplyr::mutate(df, module = dplyr::coalesce(as.character(.data$module), "all"))
  .bp_make_value_plots(
    df = df, value_col = "reconstruction_error",
    id_cols = unique(c("method", "dataset", "feature_set", "k_scope", "module", "reconstruction_source",
                "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance", "seed", layout$split_vars)),
    split_vars = unique(c("dataset", "feature_set", "k_scope", "module", layout$split_vars)),
    facet_var = "confounding_name", average_over = average_over, show_se = show_se,
    title_prefix = "Integrated reconstruction error", ylab = "Squared Frobenius error / number of features"
  )
}

plot_integrated_all_signal_pooled_auc <- function(
    bench, k_scope = c("allK", "isvdxK"), average_over = "seed", show_se = TRUE, methods = NULL) {
  df <- bench$factorization_loading_auc_metrics
  layout <- .bp_integrated_factorization_layout(bench, df, average_over = average_over)
  .bp_check_columns(df, c("method", "dataset", "k_scope", "allSignalFeaturesPooled_auc",
                          "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance"),
                    "bench$factorization_loading_auc_metrics")
  df <- .bp_filter_integrated_plot_methods(df, methods)
  df <- dplyr::filter(df, .data$k_scope %in% .env$k_scope)
  .bp_make_value_plots(
    df = df, value_col = "allSignalFeaturesPooled_auc",
    id_cols = unique(c("method", "dataset", "k_scope", "ng", "confounding_name",
                "batch_to_biology_ratio_name", "time_imbalance", "seed", layout$split_vars)),
    split_vars = unique(c("dataset", "k_scope", layout$split_vars)),
    facet_var = "confounding_name", average_over = average_over, show_se = show_se,
    title_prefix = "Integrated pooled signal loading AUC", ylab = "Signal vs non-signal loading AUC"
  )
}

plot_integrated_all_batch_pooled_auc <- function(
    bench, k_scope = c("allK", "isvdxK"), average_over = "seed", show_se = TRUE, methods = NULL) {
  df <- bench$factorization_loading_auc_metrics
  layout <- .bp_integrated_factorization_layout(bench, df, average_over = average_over)
  .bp_check_columns(df, c("method", "dataset", "k_scope", "allBatchFeaturesPooled_auc",
                          "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance"),
                    "bench$factorization_loading_auc_metrics")
  df <- .bp_filter_integrated_plot_methods(df, methods)
  df <- dplyr::filter(df, .data$k_scope %in% .env$k_scope)
  .bp_make_value_plots(
    df = df, value_col = "allBatchFeaturesPooled_auc",
    id_cols = unique(c("method", "dataset", "k_scope", "ng", "confounding_name",
                "batch_to_biology_ratio_name", "time_imbalance", "seed", layout$split_vars)),
    split_vars = unique(c("dataset", "k_scope", layout$split_vars)),
    facet_var = "confounding_name", average_over = average_over, show_se = show_se,
    title_prefix = "Integrated pooled batch loading AUC", ylab = "Batch vs non-batch loading AUC"
  )
}

plot_integrated_inference_filtered_reconstruction <- function(
    bench, feature_set = "allFeatures", module = NULL, average_over = "seed", show_se = TRUE) {
  message("Inference-filtered reconstruction uses held-out inference information to filter the training reconstruction.")
  df <- bench$inference_filtered_reconstruction_metrics
  layout <- .bp_integrated_factorization_layout(bench, df, average_over = average_over)
  .bp_check_columns(df, c("method", "dataset", "feature_set", "reconstruction_error", "reconstruction_source",
                          "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance"),
                    "bench$inference_filtered_reconstruction_metrics")
  df <- dplyr::filter(df, .data$feature_set %in% .env$feature_set)
  if (!is.null(module)) {
    module_requested <- sub("^module([0-9]+)$", "module_\\1", as.character(module))
    available_modules <- sort(unique(as.character(df$module[!is.na(df$module)])))
    missing_modules <- setdiff(module_requested, available_modules)
    if (length(missing_modules) > 0L) {
      stop(
        "Requested module(s) not found: ", paste(module, collapse = ", "),
        ". Available modules: ", paste(available_modules, collapse = ", "),
        call. = FALSE
      )
    }
    df <- dplyr::filter(df, .data$module %in% .env$module_requested)
  }
  df <- dplyr::mutate(df, module = dplyr::coalesce(as.character(.data$module), "all"))
  .bp_make_value_plots(
    df = df, value_col = "reconstruction_error",
    id_cols = unique(c("method", "dataset", "feature_set", "module", "reconstruction_source",
                "ng", "confounding_name", "batch_to_biology_ratio_name", "time_imbalance", "seed", layout$split_vars)),
    split_vars = unique(c("dataset", "feature_set", "module", layout$split_vars)),
    facet_var = "confounding_name", average_over = average_over, show_se = show_se,
    title_prefix = "Integrated training reconstruction filtered by held-out inference",
    ylab = "Squared Frobenius error / number of features"
  )
}

plot_integrated_module_type_auc_with_warning <- function(bench, ...) {
  warning("Module AUC uses all features outside the target module as negatives, including features in other true signal modules; interpret this metric cautiously.", call. = FALSE)
  plot_integrated_module_type_auc(bench, ...)
}
