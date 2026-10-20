# ============================================================
# Single-dataset benchmark evaluation
# ============================================================
# Uses simulator structure from simulate_svdx_interactions():
#   sim$Y
#   sim$meta
#   sim$truth$biology
#   sim$truth$signal_features
#   sim$truth$module_features
#   sim$truth$module_type
# ============================================================


.benchmark_poisson_deviance <- function(Y, mu, mu_lower = 1e-8, mu_upper = 1e10) {
  Y <- as.matrix(Y)
  mu <- as.matrix(mu)
  if (!all(dim(Y) == dim(mu))) {
    stop("Y and mu must have identical dimensions for Poisson deviance.", call. = FALSE)
  }
  if (any(!is.finite(Y)) || any(Y < 0)) {
    stop("Y must contain finite non-negative values for Poisson deviance.", call. = FALSE)
  }
  mu <- pmin(pmax(mu, mu_lower), mu_upper)
  term <- ifelse(Y == 0, 0, Y * log(Y / mu))
  2 * sum(term - (Y - mu))
}

.benchmark_single_poisson_mu_hat <- function(result, sim) {
  Yhat <- as.matrix(result$Yhat)
  scale <- result$extra$fitter_extra$reconstruction_scale %||% NULL

  mu_lower <- 1e-8
  mu_upper <- 1e10

  if (identical(scale, "response_mean")) {
    mu_hat <- Yhat
  } else {
    nuisance <- sim$Y_baseline_nuisance %||% matrix(0, nrow(Yhat), ncol(Yhat))
    log1p_mu_hat <- Yhat + nuisance
    log1p_mu_hat <- pmin(
      pmax(log1p_mu_hat, log1p(mu_lower)),
      log1p(mu_upper)
    )
    mu_hat <- expm1(log1p_mu_hat)
  }

  mu_hat <- pmin(pmax(mu_hat, mu_lower), mu_upper)
  if (any(!is.finite(mu_hat))) {
    stop("Non-finite expected counts produced for method '", result$method, "'.", call. = FALSE)
  }
  mu_hat
}

.benchmark_single_poisson_mu_hat_at_k <- function(result, sim, K_eval) {
  U <- as.matrix(result$U)
  V <- as.matrix(result$V)
  Yhat <- as.matrix(result$Yhat)
  K_eval <- as.integer(K_eval)[1]

  if (!is.finite(K_eval) || K_eval < 1L) {
    stop("K_eval must be a positive finite integer.", call. = FALSE)
  }
  K_available <- min(ncol(U), ncol(V))
  if (K_available < 1L) {
    stop("No factorization components are available for count reconstruction.", call. = FALSE)
  }
  K_eval <- min(K_eval, K_available)

  U_k <- U[, seq_len(K_eval), drop = FALSE]
  V_k <- V[, seq_len(K_eval), drop = FALSE]
  full_component <- V %*% t(U)
  truncated_component <- V_k %*% t(U_k)

  scale <- result$extra$fitter_extra$reconstruction_scale %||% NULL
  mu_lower <- 1e-8
  mu_upper <- 1e10

  if (identical(scale, "response_mean")) {
    mu_full <- pmin(pmax(Yhat, mu_lower), mu_upper)
    eta_full <- log(mu_full)
    eta_offset <- eta_full - full_component
    eta_k <- truncated_component + eta_offset
    eta_k <- pmin(pmax(eta_k, log(mu_lower)), log(mu_upper))
    mu_hat_k <- exp(eta_k)
  } else {
    model_offset <- Yhat - full_component
    model_k <- truncated_component + model_offset
    nuisance <- sim$Y_baseline_nuisance %||% matrix(0, nrow(Yhat), ncol(Yhat))
    log1p_mu_hat_k <- model_k + nuisance
    log1p_mu_hat_k <- pmin(
      pmax(log1p_mu_hat_k, log1p(mu_lower)),
      log1p(mu_upper)
    )
    mu_hat_k <- expm1(log1p_mu_hat_k)
  }

  mu_hat_k <- pmin(pmax(mu_hat_k, mu_lower), mu_upper)
  if (any(!is.finite(mu_hat_k))) {
    stop("Non-finite truncated expected counts produced for method '",
         result$method, "'.", call. = FALSE)
  }
  mu_hat_k
}

# ------------------------------------------------------------
# Loading-based factorization evaluation
# ------------------------------------------------------------

.benchmark_single_fixed_svdxK_metrics <- function(result, sim, svdx_K) {
  Y <- sim$Y
  truth <- sim$truth
  method <- result$method

  U <- as.matrix(result$U)
  V <- as.matrix(result$V)
  Yhat <- as.matrix(result$Yhat)

  K_available <- min(ncol(U), ncol(V))
  K_target <- as.integer(svdx_K)[1]

  empty <- tibble::tibble(
    method = method,
    svdx_K = K_target,
    n_components_returned = K_available,
    n_components_evaluated = if (is.finite(K_target) && K_target >= 1L) min(K_target, K_available) else NA_integer_,
    has_svdxK_components = is.finite(K_target) && K_target >= 1L && K_available >= K_target,
    reconstruction_cor_truth_at_svdxK = NA_real_,
    reconstruction_error_truth_at_svdxK = NA_real_,
    relative_reconstruction_error_truth_at_svdxK = NA_real_,
    loading_feature_auc_signal_vs_non_signal_at_svdxK = NA_real_,
    loading_feature_auc_batch_vs_non_batch_at_svdxK = NA_real_,
    loading_feature_auc_signal_vs_batch_only_at_svdxK = NA_real_
  )

  if (!is.finite(K_target) || K_target < 1L || K_available < 1L) {
    return(empty)
  }

  # Cap at SVDX K, but if a method returned fewer components, use all it has.
  # Preserve natural order: no refitting, matching, or truth-based reordering.
  K_eval <- min(K_target, K_available)
  U_k <- U[, seq_len(K_eval), drop = FALSE]
  V_k <- V[, seq_len(K_eval), drop = FALSE]

  score_k <- .loading_union_score(U = U_k, V = V_k)
  feature_metrics <- .bench_evaluate_feature_score(
    method = method,
    Y = Y,
    score = score_k,
    signal_features = truth$signal_features,
    batch_features = truth$batch_features,
    batch_only_features = truth$batch_only_features
  )

  empty$loading_feature_auc_signal_vs_non_signal_at_svdxK <-
    feature_metrics$loading_feature_auc_signal_vs_non_signal %||% NA_real_
  empty$loading_feature_auc_batch_vs_non_batch_at_svdxK <-
    feature_metrics$loading_feature_auc_batch_vs_non_batch %||% NA_real_
  empty$loading_feature_auc_signal_vs_batch_only_at_svdxK <-
    feature_metrics$loading_feature_auc_signal_vs_batch_only %||% NA_real_

  # For Gaussian fits, reconstruct from exactly the first K_target returned
  # components. Preserve any additive non-factor term already present in the
  # fitter's full reconstruction (e.g. PCA centering) so the only quantity
  # being truncated is the component contribution.
  if (!.benchmark_family_has_count(sim$family %||% NULL)) {
    full_component_reconstruction <- V %*% t(U)
    reconstruction_offset <- Yhat - full_component_reconstruction
    Yhat_k <- V_k %*% t(U_k) + reconstruction_offset

    empty$reconstruction_cor_truth_at_svdxK <-
      .bench_safe_cor_vec(Yhat_k, truth$biology)
    empty$reconstruction_error_truth_at_svdxK <-
      .bench_safe_frob(Yhat_k, truth$biology)
    empty$relative_reconstruction_error_truth_at_svdxK <-
      .bench_safe_rel_frob(Yhat_k, truth$biology)
  }

  empty
}

evaluate_single_factorization_result <- function(
    result,
    sim,
    top_frac = seq(0.01, 0.50, by = 0.01),
    enrichment_top_n = 100,
    svdx_K = NULL
) {
  Y <- sim$Y
  truth <- sim$truth

  method <- result$method
  Yhat <- as.matrix(result$Yhat)

  if (!all(dim(Yhat) == dim(Y))) {
    stop("result$Yhat must have same dimensions as sim$Y.", call. = FALSE)
  }

  score <- result$feature_score
  if (is.null(score)) {
    stop("result must contain feature_score from loading union score.", call. = FALSE)
  }

  module_metadata <- .bench_module_metadata_single(Y, truth)

  if (.benchmark_family_has_count(sim$family %||% NULL)) {
    mu_hat <- .benchmark_single_poisson_mu_hat(result, sim)
    model_deviance <- .benchmark_poisson_deviance(Y, mu_hat)
    null_mu <- matrix(rowMeans(Y), nrow = nrow(Y), ncol = ncol(Y))
    null_deviance <- .benchmark_poisson_deviance(Y, null_mu)
    deviance_explained <- if (is.finite(null_deviance) && null_deviance > 0) {
      1 - model_deviance / null_deviance
    } else {
      NA_real_
    }

    signal_metrics <- tibble::tibble(
      method = method,
      poisson_deviance = model_deviance,
      mean_poisson_deviance = model_deviance / length(Y),
      null_poisson_deviance = null_deviance,
      deviance_explained = deviance_explained,
      deviance_explained_percent = 100 * deviance_explained,
      reconstruction_cor_observed = .bench_safe_cor_vec(mu_hat, Y),
      relative_reconstruction_error_observed = NA_real_,
      reconstruction_cor_truth = if (!is.null(truth$mu)) .bench_safe_cor_vec(mu_hat, truth$mu) else NA_real_,
      reconstruction_error_truth = if (!is.null(truth$mu)) .bench_safe_frob(mu_hat, truth$mu) else NA_real_,
      relative_reconstruction_error_truth = if (!is.null(truth$mu)) .bench_safe_rel_frob(mu_hat, truth$mu) else NA_real_
    )
  } else {
    signal_metrics <- tibble::tibble(
      method = method,
      poisson_deviance = NA_real_,
      mean_poisson_deviance = NA_real_,
      null_poisson_deviance = NA_real_,
      deviance_explained = NA_real_,
      deviance_explained_percent = NA_real_,
      reconstruction_cor_observed = .bench_safe_cor_vec(Yhat, Y),
      relative_reconstruction_error_observed = .bench_safe_rel_frob(Yhat, Y),
      reconstruction_cor_truth = .bench_safe_cor_vec(Yhat, truth$biology),
      reconstruction_error_truth = .bench_safe_frob(Yhat, truth$biology),
      relative_reconstruction_error_truth = .bench_safe_rel_frob(Yhat, truth$biology)
    )
  }

  feature_metrics <- .bench_evaluate_feature_score(
    method = method,
    Y = Y,
    score = score,
    signal_features = truth$signal_features,
    batch_features = truth$batch_features,
    batch_only_features = truth$batch_only_features
  )

  fixed_svdxK_metrics <- .benchmark_single_fixed_svdxK_metrics(
    result = result,
    sim = sim,
    svdx_K = svdx_K
  )

  metrics <- dplyr::bind_cols(
    signal_metrics,
    dplyr::select(feature_metrics, -method),
    dplyr::select(fixed_svdxK_metrics, -method)
  )

  feature_scores <- .bench_feature_score_tbl(method, Y, score)

  # module_auc <- .bench_evaluate_module_auc(
  #   method = method,
  #   Y = Y,
  #   score = score,
  #   module_metadata = module_metadata
  # )
  #
  # module_recall_curves <- .bench_evaluate_module_recall_curves(
  #   method = method,
  #   Y = Y,
  #   score = score,
  #   module_metadata = module_metadata,
  #   top_frac = top_frac
  # )
  module_auc <- .bench_evaluate_module_auc(
    method = method,
    Y = Y,
    V = result$V,
    module_metadata = module_metadata
  )

  module_recall_curves <- .bench_evaluate_module_recall_curves(
    method = method,
    Y = Y,
    V = result$V,
    module_metadata = module_metadata,
    top_frac = top_frac
  )


  module_recall_summary <- .bench_summarise_module_recall_curves(module_recall_curves)

  loading_enrichment <- .bench_evaluate_loading_enrichment(
    method = method,
    Y = Y,
    score = score,
    module_metadata = module_metadata,
    top_n = enrichment_top_n
  )

  loading_module_type_enrichment <- .bench_evaluate_loading_module_type_enrichment(
    method = method,
    Y = Y,
    score = score,
    module_metadata = module_metadata,
    top_n = enrichment_top_n
  )

  truth_target <- if (.benchmark_family_has_count(sim$family %||% NULL)) {
    .benchmark_single_poisson_mu_hat(result, sim)
  } else {
    truth$biology
  }
  reconstruction_metrics <- .bench_reconstruction_error_rows(
    method = method, Yhat = if (.benchmark_family_has_count(sim$family %||% NULL)) truth_target else Yhat,
    truth_matrix = if (.benchmark_family_has_count(sim$family %||% NULL)) truth$mu else truth$biology,
    truth = truth, k_scope = "allK"
  )

  K_target <- as.integer(svdx_K)[1]
  U <- as.matrix(result$U); V <- as.matrix(result$V)
  K_available <- min(ncol(U), ncol(V))
  if (is.finite(K_target) && K_target >= 1L && K_available >= 1L) {
    K_eval <- min(K_target, K_available)
    if (.benchmark_family_has_count(sim$family %||% NULL)) {
      Yhat_k <- .benchmark_single_poisson_mu_hat_at_k(result, sim, K_eval)
      truth_k <- truth$mu
    } else {
      U_k <- U[, seq_len(K_eval), drop = FALSE]
      V_k <- V[, seq_len(K_eval), drop = FALSE]
      Yhat_k <- V_k %*% t(U_k) + (Yhat - V %*% t(U))
      truth_k <- truth$biology
    }
    reconstruction_metrics <- dplyr::bind_rows(
      reconstruction_metrics,
      .bench_reconstruction_error_rows(method, Yhat_k, truth_k, truth, "svdxK")
    )
  } else {
    na_rows <- .bench_reconstruction_error_rows(method, matrix(NA_real_, nrow(Y), ncol(Y)),
                                                 matrix(NA_real_, nrow(Y), ncol(Y)), truth, "svdxK")
    reconstruction_metrics <- dplyr::bind_rows(reconstruction_metrics, na_rows)
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
}


evaluate_single_factorization_results <- function(
    results,
    sim,
    top_frac = seq(0.01, 0.50, by = 0.01),
    enrichment_top_n = 100
) {
  svdx_idx <- which(vapply(results, function(x) {
    identical(tolower(x$method %||% ""), "svdx")
  }, logical(1)))

  svdx_K <- if (length(svdx_idx) > 0L) {
    svdx_result <- results[[svdx_idx[[1]]]]
    min(ncol(as.matrix(svdx_result$U)), ncol(as.matrix(svdx_result$V)))
  } else {
    NA_integer_
  }

  evals <- lapply(results, evaluate_single_factorization_result,
                  sim = sim, top_frac = top_frac,
                  enrichment_top_n = enrichment_top_n,
                  svdx_K = svdx_K)

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
# Inference-based evaluation for svdx/inferV objects
# ------------------------------------------------------------

evaluate_single_svdx_inference <- function(
    inf,
    sim,
    method = "svdx",
    alpha = inf$alpha %||% 0.1,
    inference_scope = c("global", "union"),
    top_frac = seq(0.01, 0.50, by = 0.01)
) {
  inference_scope <- match.arg(inference_scope)
  Y <- sim$Y
  truth <- sim$truth
  module_metadata <- .bench_module_metadata_single(Y, truth)

  union_tbl <- .bench_extract_union_inference(
    inf,
    alpha = alpha,
    inference_scope = inference_scope
  )

  feature_inference <- .bench_union_feature_table(
    method = method,
    Y = Y,
    union_tbl = union_tbl
  )

  metrics <- .bench_evaluate_union_inference(
    method = method,
    Y = Y,
    union_tbl = union_tbl,
    signal_features = truth$signal_features,
    alpha = alpha
  )

  module_eval <- .bench_evaluate_inference_modules(
    method = method,
    Y = Y,
    union_tbl = union_tbl,
    module_metadata = module_metadata,
    top_frac = top_frac
  )

  module_enrichment <- .bench_evaluate_inference_enrichment(
    method = method,
    Y = Y,
    union_tbl = union_tbl,
    module_metadata = module_metadata,
    alpha = alpha
  )

  module_type_enrichment <- .bench_evaluate_inference_module_type_enrichment(
    method = method,
    Y = Y,
    union_tbl = union_tbl,
    module_metadata = module_metadata,
    alpha = alpha
  )

  list(
    feature_inference = feature_inference,
    metrics = metrics,
    module_auc = module_eval$module_auc,
    module_recall_curves = module_eval$module_recall_curves,
    module_recall_summary = module_eval$module_recall_summary,
    module_enrichment = module_enrichment,
    module_type_enrichment = module_type_enrichment
  )
}

# ------------------------------------------------------------
# Inference evaluation for feature-wise LM
# ------------------------------------------------------------

evaluate_single_feature_lm_inference <- function(
    coef_mat,
    p_mat,
    sim,
    method = "feature_lm",
    terms = NULL,
    alpha = 0.1,
    inference_scope = c("global", "union"),
    fit = NULL,
    top_frac = seq(0.01, 0.50, by = 0.01)
) {
  inference_scope <- match.arg(inference_scope)

  # coef_mat and p_mat are expected to be p x q matrices with matching
  # columns. Intercept terms are excluded by default and should not be
  # evaluated for signal-feature recovery.
  Y <- sim$Y
  truth <- sim$truth
  feature_names <- .bench_feature_names(Y)

  coef_mat <- as.matrix(coef_mat)
  p_mat <- as.matrix(p_mat)

  if (!all(dim(coef_mat) == dim(p_mat))) {
    stop("coef_mat and p_mat must have identical dimensions.", call. = FALSE)
  }

  if (is.null(rownames(p_mat))) rownames(p_mat) <- feature_names
  if (is.null(colnames(p_mat))) colnames(p_mat) <- paste0("term_", seq_len(ncol(p_mat)))

  eval_terms <- colnames(p_mat)
  eval_terms <- setdiff(eval_terms, c("(Intercept)", "Intercept", "intercept"))
  if (!is.null(terms)) eval_terms <- intersect(eval_terms, terms)

  if (length(eval_terms) == 0) {
    stop("No non-intercept feature-wise LM terms available for evaluation.", call. = FALSE)
  }

  p_sub <- p_mat[, eval_terms, drop = FALSE]
  stat_sub <- coef_mat[, eval_terms, drop = FALSE]

  union_p <- apply(p_sub, 1, function(x) {
    ok <- is.finite(x) & !is.na(x)
    if (!any(ok)) return(NA_real_)
    pmin(1, sum(ok) * min(x[ok], na.rm = TRUE))
  })

  max_abs_stat <- apply(stat_sub, 1, function(x) {
    if (all(!is.finite(x))) return(NA_real_)
    max(abs(x), na.rm = TRUE)
  })

  global_p <- union_p
  if (!is.null(fit) && requireNamespace("limma", quietly = TRUE)) {
    tt <- tryCatch(
      limma::topTable(
        fit,
        coef = eval_terms,
        number = Inf,
        sort.by = "none"
      ),
      error = function(e) NULL
    )
    if (!is.null(tt) && "P.Value" %in% names(tt)) {
      global_p <- tt$P.Value[match(rownames(p_mat), rownames(tt))]
    }
  }

  union_q <- .bench_safe_bh(union_p)
  global_q <- .bench_safe_bh(global_p)

  scope_p <- if (identical(inference_scope, "global")) global_p else union_p
  scope_q <- if (identical(inference_scope, "global")) global_q else union_q
  if (!any(is.finite(scope_p) & !is.na(scope_p)) ||
      !any(is.finite(scope_q) & !is.na(scope_q))) {
    stop(
      "Requested inference scope '", inference_scope,
      "' is unavailable for method '", method, "'.",
      call. = FALSE
    )
  }

  union_tbl <- tibble::tibble(
    feature = rownames(p_mat),
    K_tested = length(eval_terms),
    min_p = apply(p_sub, 1, min, na.rm = TRUE),
    union_p = union_p,
    union_q = union_q,
    global_p = global_p,
    global_q = global_q,
    max_abs_stat = max_abs_stat,
    score = -log10(pmax(scope_p, .Machine$double.xmin)),
    discovered = !is.na(scope_q) & is.finite(scope_q) & scope_q <= alpha,
    scope = inference_scope,
    component = inference_scope
  )

  module_metadata <- .bench_module_metadata_single(Y, truth)

  module_eval <- .bench_evaluate_inference_modules(
    method = method,
    Y = Y,
    union_tbl = union_tbl,
    module_metadata = module_metadata,
    top_frac = top_frac
  )

  list(
    feature_inference = .bench_union_feature_table(method, Y, union_tbl),
    metrics = .bench_evaluate_union_inference(method, Y, union_tbl, truth$signal_features, alpha = alpha),
    module_auc = module_eval$module_auc,
    module_recall_curves = module_eval$module_recall_curves,
    module_recall_summary = module_eval$module_recall_summary,
    module_enrichment = .bench_evaluate_inference_enrichment(method, Y, union_tbl, module_metadata, alpha = alpha),
    module_type_enrichment = .bench_evaluate_inference_module_type_enrichment(method, Y, union_tbl, module_metadata, alpha = alpha),
    evaluated_terms = eval_terms
  )
}
