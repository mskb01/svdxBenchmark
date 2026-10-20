# ============================================================
# Single-dataset benchmark plots
# ============================================================

# Choose plotting dimensions for Gaussian versus generalized single benchmarks.
.bp_single_factorization_layout <- function(
    bench,
    df,
    x_var = "ng",
    facet_var = "confounding",
    average_over = "seed",
    split_vars = NULL
) {
  if (is.null(split_vars)) {
    grid_name_map <- c(
      seeds = "seed",
      ng_grid = "ng",
      confounding_grid = "confounding",
      batch_to_biology_ratio_grid = "batch_to_biology_ratio",
      feature_noise_shape_grid = "feature_noise_shape",
      batch_biology_overlap_fraction_grid = "batch_biology_overlap_fraction",
      time_imbalance_grid = "time_imbalance",
      family_grid = "family",
      count_regime_mean_grid = "count_regime_mean"
    )
    grid_keys <- intersect(names(grid_name_map), names(bench$grid))
    varying <- vapply(grid_keys, function(k) {
      x <- bench$grid[[k]]
      length(unique(x[!is.na(x)])) > 1L
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

plot_single_factorization_metrics <- function(
    bench,
    metrics = NULL,
    average_over = "seed",
    show_se = TRUE,
    x_var = "ng",
    facet_var = "confounding",
    split_vars = NULL
) {
  df <- bench$factorization_metrics
  layout <- .bp_single_factorization_layout(
    bench, df,
    x_var = x_var,
    facet_var = facet_var,
    average_over = average_over,
    split_vars = split_vars
  )

  if (is.null(metrics)) {
    if (identical(layout$x_var, "count_regime_mean")) {
      metrics <- c("deviance_explained", "mean_poisson_deviance")
    } else {
      metrics <- c(
        "loading_feature_auc_signal_vs_non_signal",
        "loading_feature_auc_batch_vs_non_batch",
        "loading_feature_auc_signal_vs_batch_only",
        "reconstruction_cor_truth",
        "reconstruction_error_truth",
        "relative_reconstruction_error_truth"
      )
    }
  }

  required <- unique(c(
    "method",
    layout$x_var,
    layout$facet_var,
    layout$split_vars
  ))
  required <- required[!is.na(required) & nzchar(required)]
  .bp_check_columns(df, required, "bench$factorization_metrics")

  id_cols <- unique(c(
    "method", layout$x_var, layout$facet_var, layout$split_vars,
    "ng", "confounding", "batch_to_biology_ratio", "time_imbalance",
    "family", "count_regime_mean", "seed"
  ))

  .bp_make_metric_plots(
    df = df,
    metric_cols = metrics,
    id_cols = id_cols,
    split_vars = layout$split_vars,
    facet_var = layout$facet_var,
    x_var = layout$x_var,
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Single factorization metrics"
  )
}

plot_single_module_auc <- function(
    bench,
    average_over = "seed",
    show_se = TRUE,
    x_var = "ng",
    facet_var = "confounding",
    split_vars = NULL) {
  df <- bench$factorization_module_auc
  layout <- .bp_single_factorization_layout(
    bench, df, x_var = x_var, facet_var = facet_var,
    average_over = average_over, split_vars = split_vars
  )
  plot_split_vars <- c("module", "module_type", layout$split_vars)

  .bp_check_columns(
    df,
    c("method", "module", "module_type", "auc_module_vs_background", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance"),
    "bench$factorization_module_auc"
  )

  .bp_make_value_plots(
    df = df,
    value_col = "auc_module_vs_background",
    id_cols = unique(c("method", "module", "module_type", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance", "seed", layout$x_var, layout$facet_var, layout$split_vars)),
    split_vars = plot_split_vars,
    facet_var = layout$facet_var,
    x_var = layout$x_var,
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Single module AUC",
    ylab = "Module AUC vs background"
  )
}

plot_single_module_recall_summary <- function(
    bench,
    metric = "aurc",
    average_over = "seed",
    show_se = TRUE,
    x_var = "ng",
    facet_var = "confounding",
    split_vars = NULL) {
  df <- bench$factorization_module_recall_summary
  layout <- .bp_single_factorization_layout(
    bench, df, x_var = x_var, facet_var = facet_var,
    average_over = average_over, split_vars = split_vars
  )
  plot_split_vars <- c("module", "module_type", layout$split_vars)

  .bp_check_columns(
    df,
    c("method", "module", "module_type", metric, "ng", "confounding", "batch_to_biology_ratio", "time_imbalance"),
    "bench$factorization_module_recall_summary"
  )

  .bp_make_value_plots(
    df = df,
    value_col = metric,
    id_cols = unique(c("method", "module", "module_type", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance", "seed", layout$x_var, layout$facet_var, layout$split_vars)),
    split_vars = plot_split_vars,
    facet_var = layout$facet_var,
    x_var = layout$x_var,
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Single module recall summary",
    ylab = metric
  )
}

plot_single_loading_enrichment <- function(
    bench,
    metric = "odds_ratio",
    average_over = "seed",
    show_se = TRUE,
    x_var = "ng",
    facet_var = "confounding",
    split_vars = NULL) {
  df <- bench$factorization_loading_enrichment
  layout <- .bp_single_factorization_layout(
    bench, df, x_var = x_var, facet_var = facet_var,
    average_over = average_over, split_vars = split_vars
  )
  plot_split_vars <- c("module", "module_type", layout$split_vars)

  .bp_check_columns(
    df,
    c("method", "module", "module_type", metric, "ng", "confounding", "batch_to_biology_ratio", "time_imbalance"),
    "bench$factorization_loading_enrichment"
  )

  .bp_make_value_plots(
    df = df,
    value_col = metric,
    id_cols = unique(c("method", "module", "module_type", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance", "seed", layout$x_var, layout$facet_var, layout$split_vars)),
    split_vars = plot_split_vars,
    facet_var = layout$facet_var,
    x_var = layout$x_var,
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Single loading enrichment",
    ylab = metric
  )
}

plot_single_inference_metrics <- function(
    bench,
    metrics = c("inference_auc_signal_vs_null", "power", "type1_error", "empirical_fdr"),
    average_over = "seed",
    show_se = TRUE,
    exclude_feature_lm = TRUE,
    x_var = "ng",
    facet_var = "confounding",
    split_vars = NULL) {
  df <- bench$inference_metrics
  layout <- .bp_single_factorization_layout(
    bench, df, x_var = x_var, facet_var = facet_var,
    average_over = average_over, split_vars = split_vars
  )
  plot_split_vars <- layout$split_vars

  .bp_check_columns(
    df,
    c("method", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance"),
    "bench$inference_metrics"
  )

  if (exclude_feature_lm) {
    df <- dplyr::filter(df, .data$method != "feature_lm")
  }

  .bp_make_metric_plots(
    df = df,
    metric_cols = metrics,
    id_cols = unique(c("method", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance", "seed", layout$x_var, layout$facet_var, layout$split_vars)),
    split_vars = plot_split_vars,
    facet_var = layout$facet_var,
    x_var = layout$x_var,
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Single inference metrics"
  )
}

plot_single_feature_lm_inference_metrics <- function(
    bench,
    metrics = c("inference_auc_signal_vs_null", "power", "type1_error", "empirical_fdr"),
    average_over = "seed",
    show_se = TRUE,
    x_var = "ng",
    facet_var = "confounding",
    split_vars = NULL) {
  df <- bench$inference_metrics |>
    dplyr::filter(.data$method == "feature_lm")

  layout <- .bp_single_factorization_layout(
    bench, df, x_var = x_var, facet_var = facet_var,
    average_over = average_over, split_vars = split_vars
  )
  plot_split_vars <- layout$split_vars

  .bp_make_metric_plots(
    df = df,
    metric_cols = metrics,
    id_cols = unique(c("method", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance", "seed", layout$x_var, layout$facet_var, layout$split_vars)),
    split_vars = plot_split_vars,
    facet_var = layout$facet_var,
    x_var = layout$x_var,
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Feature-wise LM inference metrics"
  )
}

plot_single_inference_module_enrichment <- function(
    bench,
    metric = "odds_ratio",
    average_over = "seed",
    show_se = TRUE,
    exclude_feature_lm = TRUE,
    x_var = "ng",
    facet_var = "confounding",
    split_vars = NULL) {
  df <- bench$inference_module_enrichment
  layout <- .bp_single_factorization_layout(
    bench, df, x_var = x_var, facet_var = facet_var,
    average_over = average_over, split_vars = split_vars
  )
  plot_split_vars <- c("module", "module_type", layout$split_vars)

  .bp_check_columns(
    df,
    c("method", "module", "module_type", metric, "ng", "confounding", "batch_to_biology_ratio", "time_imbalance"),
    "bench$inference_module_enrichment"
  )

  if (exclude_feature_lm) {
    df <- dplyr::filter(df, .data$method != "feature_lm")
  }

  .bp_make_value_plots(
    df = df,
    value_col = metric,
    id_cols = unique(c("method", "module", "module_type", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance", "seed", layout$x_var, layout$facet_var, layout$split_vars)),
    split_vars = plot_split_vars,
    facet_var = layout$facet_var,
    x_var = layout$x_var,
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Single inference module enrichment",
    ylab = metric
  )
}

plot_single_feature_lm_inference_module_enrichment <- function(
    bench,
    metric = "odds_ratio",
    average_over = "seed",
    show_se = TRUE,
    x_var = "ng",
    facet_var = "confounding",
    split_vars = NULL) {
  df <- bench$inference_module_enrichment |>
    dplyr::filter(.data$method == "feature_lm")

  layout <- .bp_single_factorization_layout(
    bench, df, x_var = x_var, facet_var = facet_var,
    average_over = average_over, split_vars = split_vars
  )
  plot_split_vars <- c("module", "module_type", layout$split_vars)

  .bp_make_value_plots(
    df = df,
    value_col = metric,
    id_cols = unique(c("method", "module", "module_type", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance", "seed", layout$x_var, layout$facet_var, layout$split_vars)),
    split_vars = plot_split_vars,
    facet_var = layout$facet_var,
    x_var = layout$x_var,
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Feature-wise LM module enrichment",
    ylab = metric
  )
}

plot_single_runtime <- function(
    bench,
    time_col = "time",
    average_over = "seed",
    show_se = TRUE,
    x_var = "ng",
    facet_var = "confounding",
    split_vars = NULL) {
  df <- bench$factorization_timing
  layout <- .bp_single_factorization_layout(
    bench, df, x_var = x_var, facet_var = facet_var,
    average_over = average_over, split_vars = split_vars
  )
  plot_split_vars <- layout$split_vars

  .bp_check_columns(df, c("method", time_col, "ng", "confounding", "batch_to_biology_ratio", "time_imbalance"), "bench$factorization_timing")

  .bp_make_value_plots(
    df = df,
    value_col = time_col,
    id_cols = unique(c("method", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance", "seed", layout$x_var, layout$facet_var, layout$split_vars)),
    split_vars = plot_split_vars,
    facet_var = layout$facet_var,
    x_var = layout$x_var,
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Single runtime",
    ylab = time_col
  )
}

#' Plot single-dataset module-type AUC summaries.
plot_single_module_type_auc <- function(
    bench,
    average_over = "seed",
    show_se = TRUE,
    x_var = "ng",
    facet_var = "confounding",
    split_vars = NULL
) {
  df <- bench$factorization_module_type_auc
  if (is.null(df) || nrow(df) == 0L) {
    df <- summarize_module_type_auc(bench$factorization_module_auc)
  }

  layout <- .bp_single_factorization_layout(
    bench, df,
    x_var = x_var,
    facet_var = facet_var,
    average_over = average_over,
    split_vars = split_vars
  )
  layout$split_vars <- unique(c("model_type", layout$split_vars))

  required <- unique(c(
    "method", "model_type", "mean_auc",
    layout$x_var, layout$facet_var, layout$split_vars
  ))
  required <- required[!is.na(required) & nzchar(required)]
  .bp_check_columns(df, required, "bench$factorization_module_type_auc")

  id_cols <- unique(c(
    "method", "model_type", layout$x_var, layout$facet_var,
    layout$split_vars, "ng", "confounding", "batch_to_biology_ratio",
    "time_imbalance", "family", "count_regime_mean", "seed"
  ))

  .bp_make_value_plots(
    df = df,
    value_col = "mean_auc",
    id_cols = id_cols,
    split_vars = layout$split_vars,
    facet_var = layout$facet_var,
    x_var = layout$x_var,
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Single module-type recovery",
    ylab = "Mean module AUC vs background"
  )
}

#' Plot single-dataset module-type recall summaries.
plot_single_module_type_recall <- function(
    bench,
    average_over = "seed",
    show_se = TRUE,
    x_var = "ng",
    facet_var = "confounding",
    split_vars = NULL) {
  df <- bench$factorization_module_type_recall
  if (is.null(df) || nrow(df) == 0L) {
    df <- summarize_module_type_recall(bench$factorization_module_recall_summary)
  }
  layout <- .bp_single_factorization_layout(
    bench, df, x_var = x_var, facet_var = facet_var,
    average_over = average_over, split_vars = split_vars
  )
  plot_split_vars <- c("model_type", layout$split_vars)

  .bp_check_columns(
    df,
    c("method", "model_type", "mean_recall", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance"),
    "bench$factorization_module_type_recall"
  )

  .bp_make_value_plots(
    df = df,
    value_col = "mean_recall",
    id_cols = unique(c("method", "model_type", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance", "seed", layout$x_var, layout$facet_var, layout$split_vars)),
    split_vars = plot_split_vars,
    facet_var = layout$facet_var,
    x_var = layout$x_var,
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Single module-type recall",
    ylab = "Mean module recall"
  )
}


plot_single_loading_module_type_enrichment <- function(
    bench,
    metric = "odds_ratio",
    average_over = "seed",
    show_se = TRUE,
    x_var = "ng",
    facet_var = "confounding",
    split_vars = NULL) {
  df <- bench$factorization_loading_module_type_enrichment
  layout <- .bp_single_factorization_layout(
    bench, df, x_var = x_var, facet_var = facet_var,
    average_over = average_over, split_vars = split_vars
  )
  plot_split_vars <- c("model_type", layout$split_vars)

  .bp_check_columns(
    df,
    c("method", "model_type", metric, "ng", "confounding", "batch_to_biology_ratio", "time_imbalance"),
    "bench$factorization_loading_module_type_enrichment"
  )

  .bp_make_value_plots(
    df = df,
    value_col = metric,
    id_cols = unique(c("method", "model_type", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance", "seed", layout$x_var, layout$facet_var, layout$split_vars)),
    split_vars = plot_split_vars,
    facet_var = layout$facet_var,
    x_var = layout$x_var,
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Single loading module-type enrichment",
    ylab = metric
  )
}

plot_single_inference_module_auc <- function(
    bench,
    average_over = "seed",
    show_se = TRUE,
    x_var = "ng",
    facet_var = "confounding",
    split_vars = NULL) {
  df <- bench$inference_module_auc
  layout <- .bp_single_factorization_layout(
    bench, df, x_var = x_var, facet_var = facet_var,
    average_over = average_over, split_vars = split_vars
  )
  plot_split_vars <- c("module", "module_type", layout$split_vars)

  .bp_check_columns(
    df,
    c("method", "module", "module_type", "auc_module_vs_background", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance"),
    "bench$inference_module_auc"
  )

  .bp_make_value_plots(
    df = df,
    value_col = "auc_module_vs_background",
    id_cols = unique(c("method", "module", "module_type", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance", "seed", layout$x_var, layout$facet_var, layout$split_vars)),
    split_vars = plot_split_vars,
    facet_var = layout$facet_var,
    x_var = layout$x_var,
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Single inference module AUC",
    ylab = "Inference module AUC vs background"
  )
}

plot_single_inference_module_recall_summary <- function(
    bench,
    metric = "aurc",
    average_over = "seed",
    show_se = TRUE,
    x_var = "ng",
    facet_var = "confounding",
    split_vars = NULL) {
  df <- bench$inference_module_recall_summary
  layout <- .bp_single_factorization_layout(
    bench, df, x_var = x_var, facet_var = facet_var,
    average_over = average_over, split_vars = split_vars
  )
  plot_split_vars <- c("module", "module_type", layout$split_vars)

  .bp_check_columns(
    df,
    c("method", "module", "module_type", metric, "ng", "confounding", "batch_to_biology_ratio", "time_imbalance"),
    "bench$inference_module_recall_summary"
  )

  .bp_make_value_plots(
    df = df,
    value_col = metric,
    id_cols = unique(c("method", "module", "module_type", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance", "seed", layout$x_var, layout$facet_var, layout$split_vars)),
    split_vars = plot_split_vars,
    facet_var = layout$facet_var,
    x_var = layout$x_var,
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Single inference module recall summary",
    ylab = metric
  )
}

plot_single_inference_module_type_auc <- function(
    bench,
    average_over = "seed",
    show_se = TRUE,
    x_var = "ng",
    facet_var = "confounding",
    split_vars = NULL) {
  df <- bench$inference_module_type_auc
  if (is.null(df) || nrow(df) == 0L) {
    df <- summarize_module_type_auc(bench$inference_module_auc)
  }
  layout <- .bp_single_factorization_layout(
    bench, df, x_var = x_var, facet_var = facet_var,
    average_over = average_over, split_vars = split_vars
  )
  plot_split_vars <- c("model_type", layout$split_vars)

  .bp_check_columns(
    df,
    c("method", "model_type", "mean_auc", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance"),
    "bench$inference_module_type_auc"
  )

  .bp_make_value_plots(
    df = df,
    value_col = "mean_auc",
    id_cols = unique(c("method", "model_type", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance", "seed", layout$x_var, layout$facet_var, layout$split_vars)),
    split_vars = plot_split_vars,
    facet_var = layout$facet_var,
    x_var = layout$x_var,
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Single inference module-type recovery",
    ylab = "Mean inference module AUC vs background"
  )
}

plot_single_inference_module_type_recall <- function(
    bench,
    average_over = "seed",
    show_se = TRUE,
    x_var = "ng",
    facet_var = "confounding",
    split_vars = NULL) {
  df <- bench$inference_module_type_recall
  if (is.null(df) || nrow(df) == 0L) {
    df <- summarize_module_type_recall(bench$inference_module_recall_summary)
  }
  layout <- .bp_single_factorization_layout(
    bench, df, x_var = x_var, facet_var = facet_var,
    average_over = average_over, split_vars = split_vars
  )
  plot_split_vars <- c("model_type", layout$split_vars)

  .bp_check_columns(
    df,
    c("method", "model_type", "mean_recall", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance"),
    "bench$inference_module_type_recall"
  )

  .bp_make_value_plots(
    df = df,
    value_col = "mean_recall",
    id_cols = unique(c("method", "model_type", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance", "seed", layout$x_var, layout$facet_var, layout$split_vars)),
    split_vars = plot_split_vars,
    facet_var = layout$facet_var,
    x_var = layout$x_var,
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Single inference module-type recall",
    ylab = "Mean inference module recall"
  )
}

plot_single_inference_module_type_enrichment <- function(
    bench,
    metric = "odds_ratio",
    average_over = "seed",
    show_se = TRUE,
    x_var = "ng",
    facet_var = "confounding",
    split_vars = NULL) {
  df <- bench$inference_module_type_enrichment
  layout <- .bp_single_factorization_layout(
    bench, df, x_var = x_var, facet_var = facet_var,
    average_over = average_over, split_vars = split_vars
  )
  plot_split_vars <- c("model_type", "scope", "component", layout$split_vars)

  .bp_check_columns(
    df,
    c("method", "model_type", "scope", "component", metric, "ng", "confounding", "batch_to_biology_ratio", "time_imbalance"),
    "bench$inference_module_type_enrichment"
  )

  .bp_make_value_plots(
    df = df,
    value_col = metric,
    id_cols = unique(c("method", "model_type", "scope", "component", "ng", "confounding", "batch_to_biology_ratio", "time_imbalance", "seed", layout$x_var, layout$facet_var, layout$split_vars)),
    split_vars = plot_split_vars,
    facet_var = layout$facet_var,
    x_var = layout$x_var,
    average_over = average_over,
    show_se = show_se,
    title_prefix = "Single inference module-type enrichment",
    ylab = metric
  )
}

plot_single_inference_module_recall <- function(bench, ...) {
  plot_single_inference_module_recall_summary(bench, ...)
}

# Backward-compatible aliases for module-type plot names used in reports/smoke tests.
plot_module_type_auc <- function(bench, ...) {
  if (is.list(bench) && is.data.frame(bench$factorization_module_auc)) {
    return(plot_single_module_type_auc(bench, ...))
  }
  stop("plot_module_type_auc() expects a single benchmark object with factorization_module_auc. For integrated benchmarks, use plot_integrated_module_type_auc().", call. = FALSE)
}

plot_module_type_recall <- function(bench, ...) {
  if (is.list(bench) && (is.data.frame(bench$factorization_module_recall_summary) || is.data.frame(bench$factorization_module_type_recall))) {
    return(plot_single_module_type_recall(bench, ...))
  }
  stop("plot_module_type_recall() expects a single benchmark object with factorization_module_recall_summary. For integrated benchmarks, use plot_integrated_module_type_recall().", call. = FALSE)
}

plot_single_module_recall <- function(bench, ...) {
  plot_single_module_recall_summary(bench, ...)
}

# ------------------------------------------------------------
# Additive metric plots (existing plot helpers only)
# ------------------------------------------------------------

plot_single_reconstruction_error <- function(
    bench, feature_set = "allFeatures", k_scope = "allK", module = NULL,
    average_over = "seed", show_se = TRUE, x_var = "ng",
    facet_var = "confounding", split_vars = NULL) {
  df <- bench$factorization_reconstruction_metrics
  .bp_check_columns(df, c("method", "feature_set", "k_scope", "reconstruction_error"),
                    "bench$factorization_reconstruction_metrics")
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
  layout <- .bp_single_factorization_layout(
    bench, df, x_var = x_var, facet_var = facet_var,
    average_over = average_over, split_vars = split_vars
  )
  .bp_make_value_plots(
    df = df, value_col = "reconstruction_error",
    id_cols = unique(c("method", "feature_set", "k_scope", "module", "reconstruction_source",
                       "ng", "confounding", "batch_to_biology_ratio", "time_imbalance",
                       "family", "count_regime_mean", "seed", layout$x_var, layout$facet_var, layout$split_vars)),
    split_vars = unique(c("feature_set", "k_scope", "module", layout$split_vars)),
    facet_var = layout$facet_var, x_var = layout$x_var,
    average_over = average_over, show_se = show_se,
    title_prefix = "Single reconstruction error", ylab = "Squared Frobenius error / number of features"
  )
}

plot_single_all_signal_pooled_auc <- function(
    bench, k_scope = c("allK", "svdxK"), average_over = "seed", show_se = TRUE,
    x_var = "ng", facet_var = "confounding", split_vars = NULL) {
  df <- bench$factorization_loading_auc_metrics
  .bp_check_columns(df, c("method", "k_scope", "allSignalFeaturesPooled_auc"),
                    "bench$factorization_loading_auc_metrics")
  df <- dplyr::filter(df, .data$k_scope %in% .env$k_scope)
  layout <- .bp_single_factorization_layout(
    bench, df, x_var = x_var, facet_var = facet_var,
    average_over = average_over, split_vars = split_vars
  )
  .bp_make_value_plots(
    df = df, value_col = "allSignalFeaturesPooled_auc",
    id_cols = unique(c("method", "k_scope", "ng", "confounding", "batch_to_biology_ratio",
                       "time_imbalance", "family", "count_regime_mean", "seed",
                       layout$x_var, layout$facet_var, layout$split_vars)),
    split_vars = unique(c("k_scope", layout$split_vars)), facet_var = layout$facet_var,
    x_var = layout$x_var, average_over = average_over, show_se = show_se,
    title_prefix = "Single pooled signal loading AUC", ylab = "Signal vs non-signal loading AUC"
  )
}

plot_single_inference_filtered_reconstruction <- function(
    bench, feature_set = "allFeatures", module = NULL, average_over = "seed", show_se = TRUE,
    x_var = "ng", facet_var = "confounding", split_vars = NULL) {
  message("Inference-filtered reconstruction uses held-out inference information to filter the training reconstruction.")
  df <- bench$inference_filtered_reconstruction_metrics
  .bp_check_columns(df, c("method", "feature_set", "reconstruction_error", "reconstruction_source"),
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
  layout <- .bp_single_factorization_layout(
    bench, df, x_var = x_var, facet_var = facet_var,
    average_over = average_over, split_vars = split_vars
  )
  .bp_make_value_plots(
    df = df, value_col = "reconstruction_error",
    id_cols = unique(c("method", "feature_set", "module", "reconstruction_source", "ng", "confounding",
                       "batch_to_biology_ratio", "time_imbalance", "family", "count_regime_mean", "seed",
                       layout$x_var, layout$facet_var, layout$split_vars)),
    split_vars = unique(c("feature_set", "module", layout$split_vars)), facet_var = layout$facet_var,
    x_var = layout$x_var, average_over = average_over, show_se = show_se,
    title_prefix = "Single training reconstruction filtered by held-out inference",
    ylab = "Squared Frobenius error / number of features"
  )
}

plot_single_module_type_auc_with_warning <- function(bench, ...) {
  warning("Module AUC uses all features outside the target module as negatives, including features in other true signal modules; interpret this metric cautiously.", call. = FALSE)
  plot_single_module_type_auc(bench, ...)
}
