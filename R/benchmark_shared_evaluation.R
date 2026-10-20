# .resolve_signal_features
# new_factorization_benchmark_result
# as_factorization_benchmark_result
# evaluate_factorization_result
# evaluate_method_results
# plot_method_benchmark_metrics

.resolve_signal_features <- function(Y, signal_features) {
  feature_names <- rownames(Y) %||% paste0("feature_", seq_len(nrow(Y)))
  if (is.null(signal_features)) {
    stop("signal_features must be supplied.", call. = FALSE)
  }
  if (is.character(signal_features)) signal_names <- signal_features
  else signal_names <- feature_names[as.integer(signal_features)]
  feature_names %in% signal_names
}

new_factorization_benchmark_result <- function(
    method,
    U,
    V,
    Yhat = NULL,
    feature_score = NULL,
    extra = list()
) {
  U <- as.matrix(U)
  V <- as.matrix(V)
  
  if (is.null(Yhat)) Yhat <- V %*% t(U)
  Yhat <- as.matrix(Yhat)
  
  if (is.null(feature_score)) {
    feature_score <- if (is.null(V) || is.null(U) || all(is.na(V))) {
      rep(NA_real_, nrow(Yhat))
    } else {
      .loading_union_score(U = U, V = V)
    }
  }
  
  out <- list(
    method = method,
    U = U,
    V = V,
    Yhat = Yhat,
    feature_score = as.numeric(feature_score),
    extra = extra
  )
  class(out) <- "factorizationBenchmarkResult"
  out
}

#' Wrap arbitrary factorization output into a benchmark result
#' @export
as_factorization_benchmark_result <- function(
    method,
    U,
    V,
    Yhat = NULL,
    feature_score = NULL,
    extra = list()
) {
  new_factorization_benchmark_result(
    method = method,
    U = U,
    V = V,
    Yhat = Yhat,
    feature_score = feature_score,
    extra = extra
  )
}


#' Evaluate one standardized factorization result
#' @export
evaluate_factorization_result <- function(
    result,
    Y,
    truth_signal_matrix = NULL,
    signal_features = NULL,
    feature_sets = NULL,
    signal_modules = NULL
) {
  if (!inherits(result, "factorizationBenchmarkResult")) {
    stop("result must inherit from factorizationBenchmarkResult.", call. = FALSE)
  }
  
  Y <- as.matrix(Y)
  Yhat <- as.matrix(result$Yhat)
  if (!all(dim(Yhat) == dim(Y))) {
    stop("result$Yhat must have the same dimensions as Y: feature x sample.", call. = FALSE)
  }
  
  recon_cor_observed <- .safe_cor_vec(Yhat, Y)
  rel_recon_error_observed <- .safe_rel_frob(Yhat, Y)
  
  recon_cor_truth <- NA_real_
  rel_recon_error_truth <- NA_real_
  if (!is.null(truth_signal_matrix)) {
    truth_signal_matrix <- as.matrix(truth_signal_matrix)
    if (!all(dim(truth_signal_matrix) == dim(Y))) {
      stop("truth_signal_matrix must have same dimensions as Y.", call. = FALSE)
    }
    recon_cor_truth <- .safe_cor_vec(Yhat, truth_signal_matrix)
    rel_recon_error_truth <- .safe_rel_frob(Yhat, truth_signal_matrix)
  }
  
  feature_auc <- NA_real_
  if (!is.null(signal_features)) {
    truth_signal <- .resolve_signal_features(Y, signal_features)
    feature_auc <- .safe_auc_score(result$feature_score, truth_signal)
  }
  
  metrics <- tibble::tibble(
    method = result$method,
    reconstruction_cor_observed = recon_cor_observed,
    relative_reconstruction_error_observed = rel_recon_error_observed,
    reconstruction_cor_truth = recon_cor_truth,
    relative_reconstruction_error_truth = rel_recon_error_truth,
    feature_auc_signal_vs_null = feature_auc
  )
  
  enrichment <- NULL
  if (!is.null(feature_sets) && !is.null(signal_modules)) {
    enrichment <- evaluate_module_recovery_from_scores(
      scores = result$feature_score,
      feature_names = rownames(Y) %||% paste0("feature_", seq_len(nrow(Y))),
      feature_sets = feature_sets,
      signal_modules = signal_modules,
      method = result$method
    )
  }
  
  list(metrics = metrics, enrichment = enrichment)
}

#' Evaluate module recovery from method-level feature scores
#' @export
evaluate_module_recovery_from_scores <- function(
    scores,
    feature_names,
    feature_sets,
    signal_modules,
    method = NA_character_,
    score_direction = c("larger", "smaller")
) {
  score_direction <- match.arg(score_direction)
  scores <- as.numeric(scores)
  names(scores) <- feature_names
  
  module_scores <- vapply(names(feature_sets), function(module) {
    genes <- intersect(feature_sets[[module]], feature_names)
    if (length(genes) == 0) return(NA_real_)
    stats::median(scores[genes], na.rm = TRUE)
  }, numeric(1))
  
  if (score_direction == "smaller") module_scores <- -module_scores
  
  truth <- names(module_scores) %in% signal_modules
  auc <- .safe_auc_score(module_scores, truth)
  
  tibble::tibble(
    method = method,
    module = names(module_scores),
    module_score = as.numeric(module_scores),
    is_signal_module = truth,
    module_auc_signal_vs_null = auc
  )
}

#' Evaluate multiple standardized factorization results
#' @export
evaluate_method_results <- function(
    results,
    Y,
    truth_signal_matrix = NULL,
    signal_features = NULL,
    feature_sets = NULL,
    signal_modules = NULL
) {
  evals <- lapply(
    results,
    evaluate_factorization_result,
    Y = Y,
    truth_signal_matrix = truth_signal_matrix,
    signal_features = signal_features,
    feature_sets = feature_sets,
    signal_modules = signal_modules
  )
  
  list(
    metrics = dplyr::bind_rows(lapply(evals, `[[`, "metrics")),
    enrichment = dplyr::bind_rows(lapply(evals, `[[`, "enrichment"))
  )
}

#' Plot across-method factorization benchmark metrics
#' @export
plot_method_benchmark_metrics <- function(metrics) {
  df <- metrics |>
    tidyr::pivot_longer(
      cols = c(
        reconstruction_cor_observed,
        relative_reconstruction_error_observed,
        reconstruction_cor_truth,
        relative_reconstruction_error_truth,
        feature_auc_signal_vs_null
      ),
      names_to = "metric",
      values_to = "value"
    ) |>
    dplyr::filter(is.finite(value))
  
  ggplot2::ggplot(df, ggplot2::aes(x = method, y = value, fill = method)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.4) +
    ggplot2::facet_wrap(~ metric, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    ggplot2::labs(x = NULL, y = "Value", title = "Across-method factorization benchmark")
}

#' Plot module recovery benchmarks across methods
#' @export
plot_method_module_recovery <- function(enrichment) {
  if (is.null(enrichment) || nrow(enrichment) == 0) {
    stop("No module recovery results supplied.", call. = FALSE)
  }
  
  auc_df <- enrichment |>
    dplyr::distinct(method, module_auc_signal_vs_null)
  
  ggplot2::ggplot(auc_df, ggplot2::aes(x = method, y = module_auc_signal_vs_null, fill = method)) +
    ggplot2::geom_col() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    ggplot2::labs(x = NULL, y = "Module AUC", title = "Module recovery across methods")
}
