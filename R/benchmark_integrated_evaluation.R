# ============================================================
# Integrated benchmark evaluation
# ============================================================
# Mirrors benchmark_single_evaluation.R, with dataset added to all
# dataset-specific output tables.
#
# Expected integrated simulator structure:
#   sim$Y[[dataset]]
#   sim$truth$biology[[dataset]]
#   sim$truth$signal_features[[dataset]]
#   sim$truth$module_features[[dataset]]
#   sim$truth$module_type
#
# Expected integrated fitter result structure:
#   result$method
#   result$U
#   result$V_list[[dataset]]
#   result$Yhat_list[[dataset]]
#   result$feature_score_list[[dataset]]
# ============================================================

.bench_get_dataset_item <- function(x, dataset, index, name = "object") {
  if (is.null(x)) {
    stop(name, " is NULL.", call. = FALSE)
  }

  if (!is.null(names(x)) && dataset %in% names(x)) {
    return(x[[dataset]])
  }

  if (length(x) >= index) {
    return(x[[index]])
  }

  stop(
    "Could not extract dataset '", dataset, "' from ", name, ".",
    call. = FALSE
  )
}



.benchmark_integrated_reconstruct_at_k <- function(result, sim, dataset, index, K) {
  Yd <- as.matrix(sim$Y[[index]])
  Vd <- as.matrix(.bench_get_dataset_item(result$V_list, dataset, index, "result$V_list"))
  Ud <- if (is.list(result$U)) {
    as.matrix(.bench_get_dataset_item(result$U, dataset, index, "result$U"))
  } else {
    as.matrix(result$U)
  }

  K_use <- min(as.integer(K), ncol(Ud), ncol(Vd))
  if (!is.finite(K_use) || K_use < 1L) return(matrix(NA_real_, nrow(Yd), ncol(Yd)))
  U <- Ud[, seq_len(K_use), drop = FALSE]
  V <- Vd[, seq_len(K_use), drop = FALSE]
  method <- result$method %||% ""

  # NMF/sNMF require an additive shift solely to satisfy nonnegativity.
  # All other methods return U and V on the benchmark reconstruction scale.
  if (identical(method, "concat_NMF") || grepl("^concat_sNMF", method)) {
    shift <- result$extra$preprocessing$nonnegative_shift %||% 0
    return(V %*% t(U) + as.numeric(shift))
  }

  V %*% t(U)
}

# ------------------------------------------------------------
# Loading-based factorization evaluation
# ------------------------------------------------------------


.benchmark_integrated_fixed_isvdxK_metrics <- function(result, sim, dataset, index, isvdx_K) {
  method <- result$method
  Yd <- sim$Y[[index]]
  truth <- sim$truth

  Vd <- as.matrix(.bench_get_dataset_item(result$V_list, dataset, index, "result$V_list"))
  Ud <- if (is.list(result$U)) {
    as.matrix(.bench_get_dataset_item(result$U, dataset, index, "result$U"))
  } else {
    as.matrix(result$U)
  }
  K_available <- min(ncol(Ud), ncol(Vd))
  K_target <- as.integer(isvdx_K)[1]
  K_eval <- if (is.finite(K_target) && K_target >= 1L && K_available >= 1L) min(K_target, K_available) else NA_integer_

  out <- tibble::tibble(
    method = method, dataset = dataset,
    isvdx_K = K_target,
    n_components_returned = K_available,
    n_components_evaluated_at_isvdxK = K_eval,
    has_isvdxK_components = is.finite(K_target) && K_target >= 1L && K_available >= K_target,
    reconstruction_cor_truth_at_isvdxK = NA_real_,
    relative_reconstruction_error_truth_at_isvdxK = NA_real_,
    loading_feature_auc_signal_vs_non_signal_at_isvdxK = NA_real_,
    loading_feature_auc_batch_vs_non_batch_at_isvdxK = NA_real_,
    loading_feature_auc_signal_vs_batch_only_at_isvdxK = NA_real_
  )
  if (!is.finite(K_eval) || K_eval < 1L) return(out)

  Ud_k <- Ud[, seq_len(K_eval), drop = FALSE]
  Vd_k <- Vd[, seq_len(K_eval), drop = FALSE]

  # Same loading statistic as the existing benchmark, now with a common cap.
  score_k <- .loading_union_score(U = Ud_k, V = Vd_k)
  fm <- .bench_evaluate_feature_score(
    method = method, Y = Yd, score = score_k,
    signal_features = truth$signal_features[[dataset]] %||% truth$signal_features[[index]],
    batch_features = truth$batch_features[[dataset]] %||% truth$batch_features[[index]],
    batch_only_features = truth$batch_only_features[[dataset]] %||% truth$batch_only_features[[index]],
    dataset = dataset
  )
  out$loading_feature_auc_signal_vs_non_signal_at_isvdxK <- fm$loading_feature_auc_signal_vs_non_signal %||% NA_real_
  out$loading_feature_auc_batch_vs_non_batch_at_isvdxK <- fm$loading_feature_auc_batch_vs_non_batch %||% NA_real_
  out$loading_feature_auc_signal_vs_batch_only_at_isvdxK <- fm$loading_feature_auc_signal_vs_batch_only %||% NA_real_

  Yhat_k <- .benchmark_integrated_reconstruct_at_k(
    result = result, sim = sim, dataset = dataset, index = index, K = K_eval
  )
  biology_d <- truth$biology[[dataset]] %||% truth$biology[[index]]
  out$reconstruction_cor_truth_at_isvdxK <- .bench_safe_cor_vec(Yhat_k, biology_d)
  out$relative_reconstruction_error_truth_at_isvdxK <- .bench_safe_rel_frob(Yhat_k, biology_d)
  out
}


evaluate_integrated_factorization_result <- function(
    result,
    sim,
    top_frac = seq(0.01, 0.50, by = 0.01),
    enrichment_top_n = 100,
    isvdx_K = NULL
) {
  if (!inherits(result, "integratedBenchmarkResult")) {
    stop("result must inherit from integratedBenchmarkResult.", call. = FALSE)
  }

  Y <- sim$Y
  truth <- sim$truth
  dataset_names <- names(Y) %||% paste0("dataset_", seq_along(Y))

  evals <- lapply(seq_along(dataset_names), function(d) {
    ds <- dataset_names[d]
    Yd <- Y[[d]]

    Vd <- .bench_get_dataset_item(result$V_list, ds, d, "result$V_list")
    Yhatd <- .bench_get_dataset_item(result$Yhat_list, ds, d, "result$Yhat_list")
    score_d <- .bench_get_dataset_item(result$feature_score_list, ds, d, "result$feature_score_list")

    if (!all(dim(Yhatd) == dim(Yd))) {
      stop(
        "For method '", result$method, "' and dataset '", ds,
        "', Yhat dimensions do not match Y.",
        call. = FALSE
      )
    }

    module_metadata <- .bench_module_metadata_integrated(
      Yd = Yd,
      truth = truth,
      dataset = ds
    )

    signal_metrics <- tibble::tibble(
      method = result$method,
      dataset = ds,
      reconstruction_cor_observed = .bench_safe_cor_vec(Yhatd, Yd),
      relative_reconstruction_error_observed = .bench_safe_rel_frob(Yhatd, Yd),
      reconstruction_cor_truth = .bench_safe_cor_vec(Yhatd, truth$biology[[ds]] %||% truth$biology[[d]]),
      relative_reconstruction_error_truth = .bench_safe_rel_frob(Yhatd, truth$biology[[ds]] %||% truth$biology[[d]])
    )

    feature_metrics <- .bench_evaluate_feature_score(
      method = result$method,
      Y = Yd,
      score = score_d,
      signal_features = truth$signal_features[[ds]] %||% truth$signal_features[[d]],
      batch_features = truth$batch_features[[ds]] %||% truth$batch_features[[d]],
      batch_only_features = truth$batch_only_features[[ds]] %||% truth$batch_only_features[[d]],
      dataset = ds
    )

    fixed_isvdxK_metrics <- .benchmark_integrated_fixed_isvdxK_metrics(
      result = result, sim = sim, dataset = ds, index = d, isvdx_K = isvdx_K
    )

    metrics <- dplyr::left_join(
      signal_metrics,
      dplyr::select(feature_metrics, -dataset),
      by = "method"
    )
    metrics <- dplyr::left_join(
      metrics, fixed_isvdxK_metrics, by = c("method", "dataset")
    )

    feature_scores <- .bench_feature_score_tbl(
      method = result$method,
      Y = Yd,
      score = score_d,
      dataset = ds
    )

    # module_auc <- .bench_evaluate_module_auc(
    #   method = result$method,
    #   Y = Yd,
    #   score = score_d,
    #   module_metadata = module_metadata,
    #   dataset = ds
    # )
    #
    # module_recall_curves <- .bench_evaluate_module_recall_curves(
    #   method = result$method,
    #   Y = Yd,
    #   score = score_d,
    #   module_metadata = module_metadata,
    #   top_frac = top_frac,
    #   dataset = ds
    # )
    module_auc <- .bench_evaluate_module_auc(
      method = result$method,
      Y = Yd,
      V = Vd,
      module_metadata = module_metadata,
      dataset = ds
    )

    module_recall_curves <- .bench_evaluate_module_recall_curves(
      method = result$method,
      Y = Yd,
      V = Vd,
      module_metadata = module_metadata,
      top_frac = top_frac,
      dataset = ds
    )

    module_recall_summary <- .bench_summarise_module_recall_curves(
      module_recall_curves
    )

    loading_enrichment <- .bench_evaluate_loading_enrichment(
      method = result$method,
      Y = Yd,
      score = score_d,
      module_metadata = module_metadata,
      top_n = enrichment_top_n,
      dataset = ds
    )

    loading_module_type_enrichment <- .bench_evaluate_loading_module_type_enrichment(
      method = result$method,
      Y = Yd,
      score = score_d,
      module_metadata = module_metadata,
      top_n = enrichment_top_n,
      dataset = ds
    )

    biology_d <- truth$biology[[ds]] %||% truth$biology[[d]]
    reconstruction_metrics <- .bench_reconstruction_error_rows(
      method = result$method, Yhat = Yhatd, truth_matrix = biology_d, truth = truth,
      k_scope = "allK", dataset = ds, index = d
    )
    K_target <- as.integer(isvdx_K)[1]
    Ud <- if (is.list(result$U)) as.matrix(.bench_get_dataset_item(result$U, ds, d, "result$U")) else as.matrix(result$U)
    Vd <- as.matrix(Vd)
    K_available <- min(ncol(Ud), ncol(Vd))
    if (is.finite(K_target) && K_target >= 1L && K_available >= 1L) {
      K_eval <- min(K_target, K_available)
      Yhat_k <- .benchmark_integrated_reconstruct_at_k(
        result = result, sim = sim, dataset = ds, index = d, K = K_eval
      )
      reconstruction_metrics <- dplyr::bind_rows(
        reconstruction_metrics,
        .bench_reconstruction_error_rows(result$method, Yhat_k, biology_d, truth, "isvdxK", dataset = ds, index = d)
      )
    }

    list(
      metrics = metrics,
      reconstruction_metrics = reconstruction_metrics,
      feature_scores = feature_scores,
      module_auc = module_auc,
      module_recall_curves = module_recall_curves,
      module_recall_summary = module_recall_summary,
      loading_enrichment = loading_enrichment,
      loading_module_type_enrichment = loading_module_type_enrichment
    )
  })

  list(
    metrics = dplyr::bind_rows(lapply(evals, `[[`, "metrics")),
    reconstruction_metrics = dplyr::bind_rows(lapply(evals, `[[`, "reconstruction_metrics")),
    feature_scores = dplyr::bind_rows(lapply(evals, `[[`, "feature_scores")),
    module_auc = dplyr::bind_rows(lapply(evals, `[[`, "module_auc")),
    module_recall_curves = dplyr::bind_rows(lapply(evals, `[[`, "module_recall_curves")),
    module_recall_summary = dplyr::bind_rows(lapply(evals, `[[`, "module_recall_summary")),
    loading_enrichment = dplyr::bind_rows(lapply(evals, `[[`, "loading_enrichment")),
    loading_module_type_enrichment = dplyr::bind_rows(lapply(evals, `[[`, "loading_module_type_enrichment"))
  )
}


evaluate_integrated_factorization_results <- function(
    results,
    sim,
    top_frac = seq(0.01, 0.50, by = 0.01),
    enrichment_top_n = 100
) {
  isvdx_idx <- which(vapply(results, function(x) {
    identical(tolower(x$method %||% ""), "isvdx")
  }, logical(1)))

  isvdx_K <- if (length(isvdx_idx) > 0L) {
    ref <- results[[isvdx_idx[[1]]]]
    min(ncol(as.matrix(ref$U)), min(vapply(ref$V_list, ncol, integer(1))))
  } else {
    NA_integer_
  }

  evals <- lapply(
    results,
    evaluate_integrated_factorization_result,
    sim = sim,
    top_frac = top_frac,
    enrichment_top_n = enrichment_top_n,
    isvdx_K = isvdx_K
  )

  list(
    metrics = dplyr::bind_rows(lapply(evals, `[[`, "metrics")),
    reconstruction_metrics = dplyr::bind_rows(lapply(evals, `[[`, "reconstruction_metrics")),
    feature_scores = dplyr::bind_rows(lapply(evals, `[[`, "feature_scores")),
    module_auc = dplyr::bind_rows(lapply(evals, `[[`, "module_auc")),
    module_recall_curves = dplyr::bind_rows(lapply(evals, `[[`, "module_recall_curves")),
    module_recall_summary = dplyr::bind_rows(lapply(evals, `[[`, "module_recall_summary")),
    loading_enrichment = dplyr::bind_rows(lapply(evals, `[[`, "loading_enrichment")),
    loading_module_type_enrichment = dplyr::bind_rows(lapply(evals, `[[`, "loading_module_type_enrichment"))
  )
}


# ------------------------------------------------------------
# Inference-based evaluation for integrated SVDX objects
# ------------------------------------------------------------

evaluate_integrated_svdx_inference <- function(
    method,
    fit_obj,
    sim,
    alpha = 0.1,
    top_frac = seq(0.01, 0.50, by = 0.01)
) {
  dataset_names <- names(sim$Y) %||% paste0("dataset_", seq_along(sim$Y))

  evals <- lapply(dataset_names, function(ds) {
    if (method %in% c("isvdx", "isvdx_adjusted")) {
      inf_d <- fit_obj$inf$dataset_results[[ds]]
    } else if (identical(method, "separate_svdx")) {
      inf_d <- fit_obj$inf[[ds]]
    } else {
      stop("Unsupported integrated inference method: ", method, call. = FALSE)
    }

    Yd <- sim$Y[[ds]]
    module_metadata <- .bench_module_metadata_integrated(
      Yd = Yd,
      truth = sim$truth,
      dataset = ds
    )

    union_tbl <- .bench_extract_union_inference(inf_d, alpha = alpha)

    module_eval <- .bench_evaluate_inference_modules(
      method = method,
      Y = Yd,
      union_tbl = union_tbl,
      module_metadata = module_metadata,
      top_frac = top_frac,
      dataset = ds
    )

    list(
      feature_inference = .bench_union_feature_table(
        method = method,
        Y = Yd,
        union_tbl = union_tbl,
        dataset = ds
      ),
      metrics = .bench_evaluate_union_inference(
        method = method,
        Y = Yd,
        union_tbl = union_tbl,
        signal_features = sim$truth$signal_features[[ds]],
        alpha = alpha,
        dataset = ds
      ),
      module_auc = module_eval$module_auc,
      module_recall_curves = module_eval$module_recall_curves,
      module_recall_summary = module_eval$module_recall_summary,
      module_enrichment = .bench_evaluate_inference_enrichment(
        method = method,
        Y = Yd,
        union_tbl = union_tbl,
        module_metadata = module_metadata,
        alpha = alpha,
        dataset = ds
      ),
      module_type_enrichment = .bench_evaluate_inference_module_type_enrichment(
        method = method,
        Y = Yd,
        union_tbl = union_tbl,
        module_metadata = module_metadata,
        alpha = alpha,
        dataset = ds
      )
    )
  })

  list(
    feature_inference = dplyr::bind_rows(lapply(evals, `[[`, "feature_inference")),
    metrics = dplyr::bind_rows(lapply(evals, `[[`, "metrics")),
    module_auc = dplyr::bind_rows(lapply(evals, `[[`, "module_auc")),
    module_recall_curves = dplyr::bind_rows(lapply(evals, `[[`, "module_recall_curves")),
    module_recall_summary = dplyr::bind_rows(lapply(evals, `[[`, "module_recall_summary")),
    module_enrichment = dplyr::bind_rows(lapply(evals, `[[`, "module_enrichment")),
    module_type_enrichment = dplyr::bind_rows(lapply(evals, `[[`, "module_type_enrichment"))
  )
}


# ------------------------------------------------------------
# Timing/sparsity extraction
# ------------------------------------------------------------

extract_integrated_timing <- function(results) {
  dplyr::bind_rows(lapply(results, function(res) {
    if (is.null(res)) return(NULL)

    tibble::tibble(
      method = res$method %||% NA_character_,
      time = res$time %||% NA_real_
    )
  }))
}


extract_integrated_sparsity <- function(results) {
  dplyr::bind_rows(lapply(results, function(res) {
    if (is.null(res) || is.null(res$extra$sparsity)) return(NULL)

    res$extra$sparsity |>
      dplyr::mutate(
        method = res$method,
        hyperparameter_name = res$extra$hyperparameter_name %||% NA_character_,
        hyperparameter = res$extra$hyperparameter %||% NA_real_,
        .before = 1
      )
  }))
}
