if (!exists(".benchmark_family_has_count", mode = "function")) {
  .benchmark_family_has_count <- function(family) {
    if (is.null(family)) return(FALSE)
    any(unlist(family, use.names = FALSE) %in% c("poisson", "quasipoisson"), na.rm = TRUE)
  }
}

# ============================================================
# Single-dataset benchmark grid runner
# ============================================================
# Public API:
#   run_single_benchmark_once()
#   run_single_benchmark_grid()
#
# This file intentionally exposes one single-dataset grid runner.
# Factorization and inference evaluations are both orchestrated here,
# but returned as separate top-level tables.
# ============================================================


# ------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------


.benchmark_collect_once_args <- function(...) {
  args <- list(...)

  if ("methods" %in% names(args)) {
    stop(
      "Argument 'methods' is not used by benchmark runners. ",
      "Use 'factorization_methods' and/or 'inference_methods' instead.",
      call. = FALSE
    )
  }

  args
}


.benchmark_single_overlap_feasibility <- function(
    p,
    n_modules,
    pi_signal,
    batch_feature_fraction,
    batch_biology_overlap_fraction
) {
  module_size <- floor(p * pi_signal / n_modules)
  n_signal <- n_modules * module_size
  n_null <- p - n_signal
  n_batch_features <- as.integer(round(batch_feature_fraction * p))
  n_overlap <- n_modules * as.integer(round(batch_biology_overlap_fraction * module_size))
  n_batch_only <- n_batch_features - n_overlap

  if (n_overlap > n_batch_features) {
    return(list(
      feasible = FALSE,
      reason = paste0(
        "requested batch-biology overlap requires ", n_overlap,
        " biological batch-contaminated features, but batch_feature_fraction allows only ",
        n_batch_features, " total batch features"
      )
    ))
  }

  if (n_batch_only < 0L) {
    return(list(
      feasible = FALSE,
      reason = "requested batch-biology overlap leaves a negative batch-only feature count"
    ))
  }

  if (n_batch_only > n_null) {
    return(list(
      feasible = FALSE,
      reason = paste0(
        "requested batch_feature_fraction requires ", n_batch_only,
        " batch-only features, but only ", n_null,
        " non-biological features are available"
      )
    ))
  }

  list(feasible = TRUE, reason = NA_character_)
}

.benchmark_single_grid_overlap_feasibility <- function(
    batch_biology_overlap_fraction,
    once_args,
    simulator_args
) {
  get_arg <- function(name, default) {
    if (!is.null(once_args[[name]])) return(once_args[[name]])
    if (!is.null(simulator_args[[name]])) return(simulator_args[[name]])
    default
  }

  .benchmark_single_overlap_feasibility(
    p = get_arg("p", 2000),
    n_modules = get_arg("n_modules", 8),
    pi_signal = get_arg("pi_signal", 0.40),
    batch_feature_fraction = get_arg("batch_feature_fraction", 0.10),
    batch_biology_overlap_fraction = batch_biology_overlap_fraction
  )
}

.benchmark_add_single_grid_metadata <- function(
    x,
    seed,
    confounding,
    batch_to_biology_ratio,
    feature_noise_shape,
    batch_biology_overlap_fraction,
    time_imbalance,
    G,
    ng,
    K,
    n_modules,
    family = NULL,
    count_regime_mean = NULL,
    xformula = NULL,
    zformula = NULL
) {
  if (is.null(x) || nrow(x) == 0) return(x)

  x |>
    dplyr::mutate(
      seed = seed,
      confounding = confounding,
      batch_to_biology_ratio = batch_to_biology_ratio,
      feature_noise_shape = feature_noise_shape,
      batch_biology_overlap_fraction = batch_biology_overlap_fraction,
      time_imbalance = time_imbalance,
      G = G,
      ng = ng,
      n_samples = G * ng,
      K = K,
      n_modules = n_modules,
      family = if (is.null(family)) NA_character_ else as.character(family)[1L],
      count_regime_mean = if (is.null(count_regime_mean)) NA_real_ else as.numeric(count_regime_mean)[1L],
      xformula = if (is.null(xformula)) NA_character_ else paste(deparse(xformula), collapse = ""),
      zformula = if (is.null(zformula)) NA_character_ else paste(deparse(zformula), collapse = "")
    )
}


.benchmark_preprocess_single_sim <- function(
    sim,
    standardize = TRUE,
    zformula = NULL,
    family = NULL
) {
  if (isTRUE(standardize) && .benchmark_family_has_count(family)) {
    message("standardize=TRUE ignored for poisson/quasipoisson family; using raw count data.")
    standardize <- FALSE
  }

  if (!standardize) {
    sim$Y_fit <- sim$Y

    Y_baseline_input <- if (.benchmark_family_has_count(family)) {
      log1p(sim$Y)
    } else {
      sim$Y
    }

    sim$Y_for_baselines <- if (!is.null(zformula)) {
      adjust_Y_by_formula(Y_baseline_input, sim$meta, zformula)
    } else {
      Y_baseline_input
    }
    sim$Y_baseline_nuisance <- Y_baseline_input - sim$Y_for_baselines
    sim$standardization <- NULL
    return(sim)
  }

  std <- standardize_features(sim$Y)
  Y <- std$Y
  truth <- sim$truth

  # Put additive truth components on the same standardized scale as Y.
  # Each component must be centered by its own feature-wise mean; using the
  # observed-Y center would inject other components (e.g. batch) into the
  # biological truth.
  if (!is.null(truth$biology)) {
    biology_center <- rowMeans(truth$biology, na.rm = TRUE)
    truth$biology <- sweep(truth$biology, 1, biology_center, "-")
    truth$biology <- sweep(truth$biology, 1, std$scale, "/")
  }

  if (!is.null(truth$batch_effect)) {
    batch_center <- rowMeans(truth$batch_effect, na.rm = TRUE)
    truth$batch_effect <- sweep(truth$batch_effect, 1, batch_center, "-")
    truth$batch_effect <- sweep(truth$batch_effect, 1, std$scale, "/")
  }

  sim$Y <- Y
  sim$truth <- truth
  sim$Y_fit <- Y
  sim$Y_for_baselines <- if (!is.null(zformula)) {
    adjust_Y_by_formula(Y, sim$meta, zformula)
  } else {
    Y
  }
  sim$standardization <- std

  sim
}


.benchmark_subset_single_sim <- function(sim, idx) {
  truth <- sim$truth
  for (nm in c("biology", "batch_effect", "linear_predictor", "mu")) {
    if (!is.null(truth[[nm]]) && is.matrix(truth[[nm]]) && ncol(truth[[nm]]) >= max(idx)) {
      truth[[nm]] <- truth[[nm]][, idx, drop = FALSE]
    }
  }
  if (!is.null(truth$module_signal) && is.matrix(truth$module_signal) && ncol(truth$module_signal) >= max(idx)) {
    truth$module_signal <- truth$module_signal[, idx, drop = FALSE]
  }

  out <- list(
    Y = sim$Y[, idx, drop = FALSE],
    meta = sim$meta[idx, , drop = FALSE],
    data = sim$meta[idx, , drop = FALSE],
    truth = truth,
    params = sim$params,
    family = sim$family %||% NULL,
    rowOffsets = sim$rowOffsets %||% NULL,
    colOffsets = if (!is.null(sim$colOffsets)) sim$colOffsets[idx] else NULL
  )

  if (!is.null(sim$Y_fit)) {
    out$Y_fit <- sim$Y_fit[, idx, drop = FALSE]
  }
  if (!is.null(sim$Y_for_baselines)) {
    out$Y_for_baselines <- sim$Y_for_baselines[, idx, drop = FALSE]
  }
  if (!is.null(sim$Y_baseline_nuisance)) {
    out$Y_baseline_nuisance <- sim$Y_baseline_nuisance[, idx, drop = FALSE]
  }
  if (!is.null(sim$standardization)) {
    out$standardization <- sim$standardization
  }

  out
}


.benchmark_as_single_factorization_result <- function(x) {
  if (is.null(x)) return(NULL)

  # feature_lm is intentionally excluded from loading-based
  # factorization benchmarks. It should be handled separately in
  # inference benchmarking, excluding intercept terms.
  if (identical(x$method, "feature_lm")) return(NULL)

  U <- x$scores
  V <- x$loadings

  feature_score <- if (!is.null(U) && !is.null(V) && exists(".loading_union_score", mode = "function")) {
    .loading_union_score(U = U, V = V)
  } else {
    rep(NA_real_, nrow(x$Yhat))
  }

  as_factorization_benchmark_result(
    method = x$method,
    U = U,
    V = V,
    Yhat = x$Yhat,
    feature_score = feature_score,
    extra = list(
      fit = x$fit %||% NULL,
      time = x$time %||% NA_real_,
      hyperparameter = x$hyperparameter %||% NULL,
      hyperparameter_name = x$hyperparameter_name %||% NULL,
      sparsity = x$sparsity %||% NULL,
      fitter_extra = x$extra %||% list()
    )
  )
}


.benchmark_extract_single_timing <- function(fits) {
  dplyr::bind_rows(lapply(fits, function(x) {
    if (is.null(x)) return(NULL)
    tibble::tibble(
      method = x$method %||% NA_character_,
      time = x$time %||% NA_real_
    )
  }))
}


.benchmark_extract_single_sparsity <- function(fits) {
  dplyr::bind_rows(lapply(fits, function(x) {
    if (is.null(x) || is.null(x$sparsity)) return(NULL)

    x$sparsity |>
      dplyr::mutate(
        method = x$method,
        hyperparameter_name = x$hyperparameter_name %||% NA_character_,
        hyperparameter = x$hyperparameter %||% NA_real_,
        .before = 1
      )
  }))
}



.benchmark_fit_single_factorization_methods <- function(
    Y,
    Y_base,
    meta,
    K,
    xformula,
    zformula,
    methods,
    family = NULL,
    fit_family = NULL,
    rowOffsets = NULL,
    colOffsets = NULL,
    graphData = NULL,
    graphDistFn = NULL,
    spca_keepX_grid = NULL,
    sparse_nmf_lambda_grid = 10^c(-3, -2, -1),
    mefisto_covariates = NULL,
    numThreads = 4,
    seed = 1,
    modelType = NULL,
    svdx_args = list(),
    svdx_adjusted_use_z = FALSE
) {
  raw_fits <- list()

  if (any(c("svdx", "svdx_adjusted") %in% methods)) {
    if (!is.list(svdx_args)) {
      stop("svdx_args must be a list.", call. = FALSE)
    }
    if (length(svdx_args) > 0L &&
        (is.null(names(svdx_args)) || any(names(svdx_args) == ""))) {
      stop("All elements of svdx_args must be named.", call. = FALSE)
    }

    protected_svdx_args <- c(
      "Y", "K", "xformula", "zformula", "data", "zdata",
      "numThreads", "modelType", "family", "rowOffsets", "colOffsets",
      "graphData", "graphDistFn"
    )
    conflicting_svdx_args <- intersect(names(svdx_args), protected_svdx_args)
    if (length(conflicting_svdx_args) > 0L) {
      stop(
        "svdx_args cannot override benchmark-managed arguments: ",
        paste(conflicting_svdx_args, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if ("svdx" %in% methods) {
    raw_fits$svdx <- do.call(
      fit_svdx_wrapper,
      c(
        list(
          Y = Y,
          K = K,
          xformula = xformula,
          zformula = zformula,
          data = meta,
          zdata = meta,
          numThreads = numThreads,
          modelType = modelType,
          family = fit_family,
          rowOffsets = rowOffsets,
          colOffsets = colOffsets,
          graphData = graphData,
          graphDistFn = graphDistFn
        ),
        svdx_args
      )
    )
  }

  if ("svdx_adjusted" %in% methods) {
    if (.benchmark_family_has_count(family)) {
      stop(
        "svdx_adjusted is currently defined only for Gaussian benchmarks: ",
        "Y_for_baselines is transformed/residualized for count-family benchmarks.",
        call. = FALSE
      )
    }

    raw_fits$svdx_adjusted <- do.call(
      fit_svdx_wrapper,
      c(
        list(
          Y = Y_base,
          K = K,
          xformula = xformula,
          zformula = if (isTRUE(svdx_adjusted_use_z)) zformula else ~1,
          data = meta,
          zdata = meta,
          numThreads = numThreads,
          modelType = modelType,
          family = fit_family,
          rowOffsets = rowOffsets,
          colOffsets = colOffsets,
          graphData = graphData,
          graphDistFn = graphDistFn
        ),
        svdx_args
      )
    )
    raw_fits$svdx_adjusted$method <- "svdx_adjusted"
  }

  if ("glm_pca" %in% methods) {
    if (!.benchmark_family_has_count(family)) {
      stop("glm_pca is only available for count-family benchmarks.", call. = FALSE)
    }
    glm_pca_X <- NULL
    if (!is.null(zformula)) {
      glm_pca_X <- stats::model.matrix(zformula, data = meta)
      if ("(Intercept)" %in% colnames(glm_pca_X)) {
        glm_pca_X <- glm_pca_X[, colnames(glm_pca_X) != "(Intercept)", drop = FALSE]
      }
      if (ncol(glm_pca_X) == 0L) glm_pca_X <- NULL
    }

    raw_fits$glm_pca <- fit_glmpca_wrapper(
      Y = Y,
      K = K,
      X = glm_pca_X,
      colOffsets = colOffsets,
      seed = seed
    )
  }

  if ("svd" %in% methods) raw_fits$svd <- fit_svd(Y_base, K = K)
  if ("pca" %in% methods) raw_fits$pca <- fit_pca(Y_base, K = K)
  if ("nmf" %in% methods) raw_fits$nmf <- fit_nmf(Y_base, K = K, seed = seed)

  if ("mixomics_pls" %in% methods) {
    raw_fits$mixomics_pls <- fit_mixomics_pls_formula(
      Y = Y_base,
      meta = meta,
      K = K,
      xformula = xformula
    )
  }

  if ("spca" %in% methods) {
    if (is.null(spca_keepX_grid)) {
      spca_keepX_grid <- round(c(0.01, 0.25, 0.75) * nrow(Y_base))
    }

    raw_fits <- c(
      raw_fits,
      fit_mixomics_spca_grid(
        Y = Y_base,
        K = K,
        keepX_grid = spca_keepX_grid
      )
    )
  }

  if ("sparse_nmf" %in% methods) {
    raw_fits <- c(
      raw_fits,
      fit_sparse_nmf_sklearn_grid(
        Y = Y_base,
        K = K,
        lambda_grid = sparse_nmf_lambda_grid
      )
    )
  }

  if ("mefisto_duplicate_view" %in% methods) {
    raw_fits$mefisto_duplicate_view <- fit_mefisto_duplicate_view_wrapper(
      Y = Y_base,
      data = meta,
      K = K,
      covariates = mefisto_covariates,
      xformula = xformula,
      seed = seed
    )
  }

  Filter(Negate(is.null), raw_fits)
}


.benchmark_fit_single_svdx_inference <- function(
    sim_train,
    sim_test,
    K,
    xformula,
    zformula,
    trained_fit = NULL,
    family = NULL,
    fit_family = NULL,
    inference_family = "gaussian",
    graphData = NULL,
    graphDistFn = NULL,
    alpha = 0.1,
    z_adjustment = c("fitted", "reestimate"),
    numThreads = 4,
    seed = 1,
    modelType = NULL,
    use_baseline = FALSE,
    method = "svdx"
) {
  z_adjustment <- match.arg(z_adjustment)

  Ytrain <- if (isTRUE(use_baseline)) sim_train$Y_for_baselines else sim_train$Y
  Ytest <- if (isTRUE(use_baseline)) sim_test$Y_for_baselines else sim_test$Y
  meta_train <- sim_train$meta
  meta_test <- sim_test$meta

  if (is.null(trained_fit)) {
    tm_fit <- system.time({
      fit <- svdx(
        Y = Ytrain,
        formula = xformula,
        data = meta_train,
        zformula = zformula,
        zdata = meta_train,
        K = K,
        numThreads = numThreads,
        modelType = modelType,
        family = fit_family,
        rowOffsets = sim_train$rowOffsets,
        colOffsets = sim_train$colOffsets,
        graphData = graphData,
        graphDistFn = graphDistFn
      )
    })
    time_fit <- unname(tm_fit["elapsed"])
  } else {
    fit <- trained_fit
    time_fit <- 0
  }

  tm_inf <- system.time({
    inf <- inferV(
      fit = fit,
      Ytest = Ytest,
      newdata = meta_test,
      newzdata = meta_test,
      zformula = zformula,
      family = inference_family,
      rowOffsets = sim_test$rowOffsets,
      colOffsets = sim_test$colOffsets,
      alpha = alpha,
      z_adjustment = z_adjustment
    )
  })

  list(
    fit = fit,
    inf = inf,
    timing = tibble::tibble(
      method = method,
      time_fit = time_fit,
      time_inference = unname(tm_inf["elapsed"]),
      time_total = time_fit + unname(tm_inf["elapsed"])
    )
  )
}


# ------------------------------------------------------------
# Public single-dataset runner
# ------------------------------------------------------------

run_single_benchmark_once <- function(
    seed = 1,
    p = 2000,
    G = 4,
    ng = 40,
    B = 6,
    K = 8,
    n_modules = 8,
    xformula = ~ group * t + severity * t,
    zformula = ~ batch,
    confounding = 0.4,
    time_imbalance = 0.4,
    batch_to_biology_ratio = 1,
    feature_noise_shape = 1,
    batch_biology_overlap_fraction = 0.25,
    factorization_methods = c("svdx", "svdx_adjusted", "svd", "pca", "nmf", "mixomics_pls", "spca", "sparse_nmf", "mefisto_duplicate_view"),
    inference_methods = c("svdx", "svdx_adjusted"),
    run_factorization = TRUE,
    run_inference = TRUE,
    standardize = TRUE,
    family = NULL,
    fit_family = NULL,
    inference_family = "gaussian",
    graphData = NULL,
    graphDistFn = NULL,
    spca_keepX_grid = NULL,
    sparse_nmf_lambda_grid = 10^c(-3, -2, -1),
    mefisto_covariates = NULL,
    test_frac = 0.30,
    alpha = 0.1,
    inference_scope = c("global", "union"),
    z_adjustment = c("fitted", "reestimate"),
    numThreads = 4,
    top_frac = seq(0.01, 0.50, by = 0.01),
    loading_top_n = 100,
    keep_fits = FALSE,
    simulator_args = list(),
    modelType = NULL,
    svdx_args = list(),
    svdx_adjusted_use_z = FALSE,
    ...
) {
  z_adjustment <- match.arg(z_adjustment)
  inference_scope <- match.arg(inference_scope)

  valid_factorization_methods <- c("svdx", "svdx_adjusted", "glm_pca", "svd", "pca", "nmf", "mixomics_pls", "spca", "sparse_nmf", "mefisto_duplicate_view")
  invalid_factorization_methods <- setdiff(factorization_methods, valid_factorization_methods)
  if (length(invalid_factorization_methods) > 0) {
    stop(
      "Unknown factorization_methods: ", paste(invalid_factorization_methods, collapse = ", "),
      ". Valid options are: ", paste(valid_factorization_methods, collapse = ", "),
      call. = FALSE
    )
  }

  valid_inference_methods <- c("svdx", "svdx_adjusted", "feature_lm")
  invalid_inference_methods <- setdiff(inference_methods, valid_inference_methods)
  if (length(invalid_inference_methods) > 0) {
    stop(
      "Unknown inference_methods: ", paste(invalid_inference_methods, collapse = ", "),
      ". Valid options are: ", paste(valid_inference_methods, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.null(graphDistFn) && !is.function(graphDistFn)) {
    stop("graphDistFn must be a function when supplied.", call. = FALSE)
  }

  if ("family" %in% names(simulator_args)) {
    if (!is.null(family) && !identical(simulator_args$family, family)) {
      stop(
        "Conflicting family values supplied through 'family' and 'simulator_args$family'.",
        call. = FALSE
      )
    }
  } else {
    simulator_args$family <- family
  }

  simulator_args <- c(
    list(
      p = p,
      G = G,
      ng = ng,
      B = B,
      n_modules = n_modules,
      confounding = confounding,
      time_imbalance = time_imbalance,
      batch_to_biology_ratio = batch_to_biology_ratio,
      feature_noise_shape = feature_noise_shape,
      batch_biology_overlap_fraction = batch_biology_overlap_fraction,
      seed = seed
    ),
    simulator_args
  )

  sim <- do.call(simulate_svdx_interactions, simulator_args)

  resolved_family <- sim$family
  if (is.null(resolved_family)) {
    resolved_family <- if ("family" %in% names(simulator_args)) simulator_args$family else family
  }

  sim <- .benchmark_preprocess_single_sim(
    sim = sim,
    standardize = standardize,
    zformula = zformula,
    family = resolved_family
  )

  svdx_model_type <- .benchmark_complete_model_type(
    modelType = modelType,
    formula = xformula,
    data = sim$meta
  )

  if (!is.null(graphData)) {
    if (!is.data.frame(graphData)) {
      stop("graphData must be a data.frame when supplied.", call. = FALSE)
    }
    if (nrow(graphData) < 1L) {
      stop("graphData must contain at least one row.", call. = FALSE)
    }
    if (ncol(graphData) != ncol(sim$Y)) {
      stop("ncol(graphData) must equal ncol(Y).", call. = FALSE)
    }
  }

  split <- make_train_test_split(sim$meta, test_frac = test_frac, seed = seed)
  sim_train <- .benchmark_subset_single_sim(sim, split$train)
  sim_test <- .benchmark_subset_single_sim(sim, split$test)

  add_meta <- function(x) {
    .benchmark_add_single_grid_metadata(
      x,
      seed = seed,
      confounding = confounding,
      batch_to_biology_ratio = batch_to_biology_ratio,
      feature_noise_shape = feature_noise_shape,
      batch_biology_overlap_fraction = batch_biology_overlap_fraction,
      time_imbalance = time_imbalance,
      G = G,
      ng = ng,
      K = K,
      n_modules = n_modules,
      family = resolved_family,
      count_regime_mean = sim$params$count_regime_mean %||% NULL,
      xformula = xformula,
      zformula = zformula
    )
  }

  factorization <- list()
  raw_fits <- NULL
  benchmark_results <- NULL

  if (run_factorization) {
    raw_fits <- .benchmark_fit_single_factorization_methods(
      Y = sim_train$Y_fit,
      Y_base = sim_train$Y_for_baselines,
      meta = sim_train$meta,
      K = K,
      xformula = xformula,
      zformula = zformula,
      methods = factorization_methods,
      family = resolved_family,
      fit_family = fit_family,
      rowOffsets = sim_train$rowOffsets,
      colOffsets = sim_train$colOffsets,
      graphData = graphData,
      graphDistFn = graphDistFn,
      spca_keepX_grid = spca_keepX_grid,
      sparse_nmf_lambda_grid = sparse_nmf_lambda_grid,
      mefisto_covariates = mefisto_covariates,
      numThreads = numThreads,
      seed = seed,
      modelType = svdx_model_type,
      svdx_args = svdx_args,
      svdx_adjusted_use_z = svdx_adjusted_use_z
    )

    benchmark_results <- raw_fits |>
      lapply(.benchmark_as_single_factorization_result) |>
      Filter(Negate(is.null), x = _)

    eval <- evaluate_single_factorization_results(
      results = benchmark_results,
      sim = sim_train,
      top_frac = top_frac,
      enrichment_top_n = loading_top_n
    )

    factorization <- list(
      metrics = add_meta(eval$metrics),
      reconstruction_metrics = add_meta(dplyr::select(eval$reconstruction_metrics, -dataset)),
      loading_auc_metrics = add_meta(dplyr::select(.bench_loading_auc_alias_rows(eval$metrics, "svdxK"), -dataset)),
      feature_scores = add_meta(eval$feature_scores),
      module_auc = add_meta(eval$module_auc),
      module_recall_curves = add_meta(eval$module_recall_curves),
      module_recall_summary = add_meta(eval$module_recall_summary),
      loading_enrichment = add_meta(eval$loading_enrichment),
      loading_module_type_enrichment = add_meta(eval$loading_module_type_enrichment),
      timing = add_meta(.benchmark_extract_single_timing(raw_fits)),
      sparsity_metrics = add_meta(.benchmark_extract_single_sparsity(raw_fits))
    )
  }

  inference <- list()
  inf_fits <- list()
  inf_evals <- list()
  inf_timing <- list()

  if (run_inference) {
    if ("svdx" %in% inference_methods) {
      trained_svdx_fit <- NULL
      if (!is.null(raw_fits) && !is.null(raw_fits$svdx) && !is.null(raw_fits$svdx$fit)) {
        trained_svdx_fit <- raw_fits$svdx$fit
      }

      inf_fits$svdx <- .benchmark_fit_single_svdx_inference(
        sim_train = sim_train,
        sim_test = sim_test,
        K = K,
        xformula = xformula,
        zformula = zformula,
        trained_fit = trained_svdx_fit,
        family = resolved_family,
        fit_family = fit_family,
        inference_family = inference_family,
        graphData = graphData,
        graphDistFn = graphDistFn,
        alpha = alpha,
        z_adjustment = z_adjustment,
        numThreads = numThreads,
        seed = seed,
        modelType = svdx_model_type
      )

      inf_evals$svdx <- evaluate_single_svdx_inference(
        inf = inf_fits$svdx$inf,
        sim = sim_test,
        method = "svdx",
        alpha = alpha,
        inference_scope = inference_scope,
        top_frac = top_frac
      )

      inf_timing$svdx <- inf_fits$svdx$timing
    }

    if ("svdx_adjusted" %in% inference_methods) {
      if (.benchmark_family_has_count(resolved_family)) {
        stop(
          "svdx_adjusted inference is currently defined only for Gaussian benchmarks.",
          call. = FALSE
        )
      }

      trained_svdx_adjusted_fit <- NULL
      if (!is.null(raw_fits) &&
          !is.null(raw_fits$svdx_adjusted) &&
          !is.null(raw_fits$svdx_adjusted$fit)) {
        trained_svdx_adjusted_fit <- raw_fits$svdx_adjusted$fit
      }

      adjusted_zformula <- if (isTRUE(svdx_adjusted_use_z)) zformula else ~1

      inf_fits$svdx_adjusted <- .benchmark_fit_single_svdx_inference(
        sim_train = sim_train,
        sim_test = sim_test,
        K = K,
        xformula = xformula,
        zformula = adjusted_zformula,
        trained_fit = trained_svdx_adjusted_fit,
        family = resolved_family,
        fit_family = fit_family,
        inference_family = inference_family,
        graphData = graphData,
        graphDistFn = graphDistFn,
        alpha = alpha,
        z_adjustment = z_adjustment,
        numThreads = numThreads,
        seed = seed,
        modelType = svdx_model_type,
        use_baseline = TRUE,
        method = "svdx_adjusted"
      )

      inf_evals$svdx_adjusted <- evaluate_single_svdx_inference(
        inf = inf_fits$svdx_adjusted$inf,
        sim = sim_test,
        method = "svdx_adjusted",
        alpha = alpha,
        inference_scope = inference_scope,
        top_frac = top_frac
      )

      inf_timing$svdx_adjusted <- inf_fits$svdx_adjusted$timing
    }

    if ("feature_lm" %in% inference_methods) {
      lm_fit <- fit_feature_lm_formula(
        Y = sim_train$Y_for_baselines,
        meta = sim_train$meta,
        xformula = xformula
      )

      if (is.null(lm_fit$p_mat)) {
        stop(
          "fit_feature_lm_formula() must return p_mat for feature_lm inference benchmarking.",
          call. = FALSE
        )
      }

      inf_fits$feature_lm <- lm_fit

      inf_evals$feature_lm <- evaluate_single_feature_lm_inference(
        coef_mat = lm_fit$coef_mat,
        p_mat = lm_fit$p_mat,
        sim = sim_train,
        method = "feature_lm",
        alpha = alpha,
        inference_scope = inference_scope,
        fit = lm_fit$fit,
        top_frac = top_frac
      )

      inf_timing$feature_lm <- tibble::tibble(
        method = "feature_lm",
        time_fit = lm_fit$time,
        time_inference = NA_real_,
        time_total = lm_fit$time
      )
    }

    inference <- list(
      feature_inference = add_meta(dplyr::bind_rows(lapply(inf_evals, `[[`, "feature_inference"))),
      metrics = add_meta(dplyr::bind_rows(lapply(inf_evals, `[[`, "metrics"))),
      module_auc = add_meta(dplyr::bind_rows(lapply(inf_evals, `[[`, "module_auc"))),
      module_recall_curves = add_meta(dplyr::bind_rows(lapply(inf_evals, `[[`, "module_recall_curves"))),
      module_recall_summary = add_meta(dplyr::bind_rows(lapply(inf_evals, `[[`, "module_recall_summary"))),
      module_enrichment = add_meta(dplyr::bind_rows(lapply(inf_evals, `[[`, "module_enrichment"))),
      module_type_enrichment = add_meta(dplyr::bind_rows(lapply(inf_evals, `[[`, "module_type_enrichment"))),
      timing = add_meta(dplyr::bind_rows(inf_timing))
    )
  }

  if (run_factorization && run_inference && !is.null(benchmark_results)) {
    filtered_rows <- list()

    for (method_name in c("svdx", "svdx_adjusted")) {
      if (is.null(inf_evals[[method_name]])) next

      method_res <- Filter(
        function(x) identical(tolower(x$method %||% ""), method_name),
        benchmark_results
      )
      if (length(method_res) == 0L) next
      method_res <- method_res[[1]]

      fi <- inf_evals[[method_name]]$feature_inference
      keep <- rep(FALSE, nrow(method_res$V))
      if (!is.null(fi) && all(c("feature", "global_q") %in% names(fi))) {
        idx <- match(rownames(method_res$V) %||% rownames(sim_train$Y), fi$feature)
        keep <- !is.na(idx) & !is.na(fi$global_q[idx]) & fi$global_q[idx] < alpha
      }

      Vf <- as.matrix(method_res$V)
      Vf[!keep, ] <- 0
      U <- as.matrix(method_res$U)
      V <- as.matrix(method_res$V)
      Yhat <- as.matrix(method_res$Yhat)

      if (.benchmark_family_has_count(resolved_family)) {
        eta_f <- Vf %*% t(U)
        fit0 <- method_res$extra$fit %||% NULL
        Zeta <- if (!is.null(fit0)) fit0$fitAtOptCV$Zeta %||% NULL else NULL
        if (!is.null(zformula) && !is.null(Zeta)) {
          Z <- stats::model.matrix(zformula, data = sim_train$zdata %||% sim_train$meta)
          Zeta <- as.matrix(Zeta)
          if (!is.null(colnames(Zeta)) && !is.null(colnames(Z))) {
            Z <- Z[, colnames(Zeta), drop = FALSE]
          }
          eta_f <- eta_f + Zeta %*% t(Z)
        }
        if (!is.null(sim_train$rowOffsets)) {
          eta_f <- sweep(eta_f, 1, sim_train$rowOffsets, "+")
        }
        if (!is.null(sim_train$colOffsets)) {
          eta_f <- sweep(eta_f, 2, sim_train$colOffsets, "+")
        }
        eta_f <- pmin(
          pmax(eta_f, log(.Machine$double.xmin)),
          log(.Machine$double.xmax) - 1
        )
        Yhat_f <- exp(eta_f)
        truth_f <- sim_train$truth$mu
      } else {
        Yhat_f <- Vf %*% t(U) + (Yhat - V %*% t(U))
        truth_f <- sim_train$truth$biology
      }

      filtered_rows[[method_name]] <- .bench_reconstruction_error_rows(
        method_name,
        Yhat_f,
        truth_f,
        sim_train$truth,
        "svdxK",
        reconstruction_source = "training_filtered_by_heldout_inference"
      )
    }

    if (length(filtered_rows) > 0L) {
      filtered <- dplyr::bind_rows(filtered_rows)
      factorization$inference_filtered_reconstruction_metrics <-
        add_meta(dplyr::select(filtered, -dataset))
    }
  }

  prediction_metrics <- NULL
  svdx_pred_fit <- NULL
  if (!is.null(raw_fits) && !is.null(raw_fits$svdx) && !is.null(raw_fits$svdx$fit)) {
    svdx_pred_fit <- raw_fits$svdx$fit
  } else if (!is.null(inf_fits$svdx) && !is.null(inf_fits$svdx$fit)) {
    svdx_pred_fit <- inf_fits$svdx$fit
  }
  if (!is.null(svdx_pred_fit)) {
    prediction_metrics <- tryCatch(
      add_meta(.benchmark_svdx_single_prediction_metrics(
        fit = svdx_pred_fit,
        Ytrain = sim_train$Y_fit,
        Ytest = sim_test$Y_fit,
        xdata_test = sim_test$meta,
        zdata_test = sim_test$zdata %||% sim_test$meta,
        zformula = zformula,
        family = resolved_family,
        rowOffsets = sim_test$rowOffsets,
        colOffsets = sim_test$colOffsets,
        method = "svdx",
        dataset = NA_character_,
        numThreads = numThreads
      )),
      error = function(e) NULL
    )
  }

  out <- list(
    sim = sim,
    split = split,
    factorization = factorization,
    inference = inference,
    prediction_metrics = prediction_metrics,
    fits = if (keep_fits) list(
      factorization = raw_fits,
      factorization_benchmark_results = benchmark_results,
      inference = inf_fits
    ) else NULL,
    xformula = xformula,
    zformula = zformula,
    factorization_methods = factorization_methods,
    inference_methods = inference_methods,
    run_factorization = run_factorization,
    run_inference = run_inference
  )

  class(out) <- "singleBenchmarkRun"
  out
}


run_single_benchmark_grid <- function(
    seeds = 1:3,
    ng_grid = c(5, 10, 25, 50),
    confounding_grid = c(0, 0.25, 0.5, 0.75, 1),
    batch_to_biology_ratio_grid = c(1, 10),
    feature_noise_shape_grid = c(1, 5),
    batch_biology_overlap_fraction_grid = c(0.10, 0.25, 0.50, 0.75),
    time_imbalance_grid = c(0, 0.5),
    family_grid = "gaussian",
    fit_family = NULL,
    inference_family = "gaussian",
    mefisto_covariates = NULL,
    spca_keepX_grid = NULL,
    sparse_nmf_lambda_grid = 10^c(-3, -2, -1),
    count_regime_mean_grid = NULL,
    keep_runs = FALSE,
    simulator_args = list(),
    modelType = NULL,
    svdx_args = list(),
    svdx_adjusted_use_z = FALSE,
    ...
) {
  once_args <- .benchmark_collect_once_args(...)

  if ("svdx_args" %in% names(once_args)) {
    stop("Specify svdx_args through the dedicated grid argument, not through '...'.", call. = FALSE)
  }
  if ("svdx_adjusted_use_z" %in% names(once_args)) {
    stop("Specify svdx_adjusted_use_z through the dedicated grid argument, not through '...'.", call. = FALSE)
  }
  if (!is.list(svdx_args)) {
    stop("svdx_args must be a list.", call. = FALSE)
  }

  if ("family" %in% names(once_args)) {
    stop("Specify family through family_grid, not through '...'.", call. = FALSE)
  }
  if ("fit_family" %in% names(once_args) || "inference_family" %in% names(once_args) ||
      "mefisto_covariates" %in% names(once_args)) {
    stop(
      "Specify fit_family, inference_family, and mefisto_covariates through their dedicated grid arguments, not through '...'.",
      call. = FALSE
    )
  }
  if (!is.null(simulator_args$family)) {
    stop("Specify family through family_grid, not through simulator_args.", call. = FALSE)
  }

  family_grid <- unique(as.character(family_grid))
  if (length(family_grid) == 0L || anyNA(family_grid) ||
      any(!family_grid %in% c("gaussian", "poisson"))) {
    stop("family_grid must contain only 'gaussian' and/or 'poisson'.", call. = FALSE)
  }

  if (!is.null(count_regime_mean_grid)) {
    count_regime_mean_grid <- as.numeric(count_regime_mean_grid)
    if (length(count_regime_mean_grid) == 0L ||
        any(!is.finite(count_regime_mean_grid)) ||
        any(count_regime_mean_grid <= 0)) {
      stop("count_regime_mean_grid must contain positive finite values.", call. = FALSE)
    }
    if (!is.null(simulator_args$count_regime_mean)) {
      stop(
        "Specify count_regime_mean through either count_regime_mean_grid or simulator_args, not both.",
        call. = FALSE
      )
    }
  }

  runs <- if (isTRUE(keep_runs)) list() else NULL
  run_tables <- list(
    factorization_metrics = list(),
    factorization_reconstruction_metrics = list(),
    factorization_loading_auc_metrics = list(),
    inference_filtered_reconstruction_metrics = list(),
    factorization_feature_scores = list(),
    factorization_module_auc = list(),
    factorization_module_recall_curves = list(),
    factorization_module_recall_summary = list(),
    factorization_loading_enrichment = list(),
    factorization_loading_module_type_enrichment = list(),
    factorization_timing = list(),
    factorization_sparsity_metrics = list(),
    inference_feature_inference = list(),
    inference_metrics = list(),
    inference_module_auc = list(),
    inference_module_recall_curves = list(),
    inference_module_recall_summary = list(),
    inference_module_enrichment = list(),
    inference_module_type_enrichment = list(),
    inference_timing = list(),
    prediction_metrics = list()
  )
  errors <- list()
  idx <- 1L
  err_idx <- 1L

  for (family_value in family_grid) {
    count_values <- if (identical(family_value, "poisson")) {
      if (!is.null(count_regime_mean_grid)) {
        count_regime_mean_grid
      } else if (!is.null(simulator_args$count_regime_mean)) {
        as.numeric(simulator_args$count_regime_mean)
      } else {
        stop(
          "count_regime_mean_grid must be supplied when family_grid includes 'poisson'.",
          call. = FALSE
        )
      }
    } else {
      NA_real_
    }

    for (count_regime_mean in count_values) {
      run_simulator_args <- simulator_args
      if (identical(family_value, "poisson")) {
        run_simulator_args$count_regime_mean <- count_regime_mean
      } else {
        run_simulator_args$count_regime_mean <- NULL
      }

      feature_noise_values <- if (identical(family_value, "poisson")) {
        1
      } else {
        feature_noise_shape_grid
      }

      for (ng in ng_grid) {
        for (batch_to_biology_ratio in batch_to_biology_ratio_grid) {
          for (feature_noise_shape in feature_noise_values) {
            for (batch_biology_overlap_fraction in batch_biology_overlap_fraction_grid) {
              feasibility <- .benchmark_single_grid_overlap_feasibility(
                batch_biology_overlap_fraction = batch_biology_overlap_fraction,
                once_args = once_args,
                simulator_args = run_simulator_args
              )
              if (!isTRUE(feasibility$feasible)) next

              for (confounding in confounding_grid) {
                for (time_imbalance in time_imbalance_grid) {
                  for (seed in seeds) {
                    run_args <- c(
                      list(
                        seed = seed,
                        family = family_value,
                        fit_family = fit_family,
                        inference_family = inference_family,
                        mefisto_covariates = mefisto_covariates,
                        spca_keepX_grid = spca_keepX_grid,
                        sparse_nmf_lambda_grid = sparse_nmf_lambda_grid,
                        ng = ng,
                        confounding = confounding,
                        batch_to_biology_ratio = batch_to_biology_ratio,
                        feature_noise_shape = feature_noise_shape,
                        batch_biology_overlap_fraction = batch_biology_overlap_fraction,
                        time_imbalance = time_imbalance,
                        simulator_args = run_simulator_args,
                        modelType = modelType,
                        svdx_args = svdx_args,
                        svdx_adjusted_use_z = svdx_adjusted_use_z
                      ),
                      once_args
                    )
                    res <- tryCatch(
                      do.call(run_single_benchmark_once, run_args),
                      error = function(e) e
                    )
                    if (inherits(res, "error")) {
                      errors[[err_idx]] <- tibble::tibble(
                        seed = seed,
                        family = family_value,
                        count_regime_mean = if (is.null(count_regime_mean)) NA_real_ else count_regime_mean,
                        ng = ng,
                        confounding = confounding,
                        batch_to_biology_ratio = batch_to_biology_ratio,
                        feature_noise_shape = feature_noise_shape,
                        batch_biology_overlap_fraction = batch_biology_overlap_fraction,
                        time_imbalance = time_imbalance,
                        error_class = paste(class(res), collapse = ","),
                        error_message = conditionMessage(res)
                      )
                      err_idx <- err_idx + 1L
                    } else {
                      run_tables$factorization_metrics[idx] <- list(res$factorization$metrics)
                      run_tables$factorization_reconstruction_metrics[idx] <- list(res$factorization$reconstruction_metrics)
                      run_tables$factorization_loading_auc_metrics[idx] <- list(res$factorization$loading_auc_metrics)
                      run_tables$inference_filtered_reconstruction_metrics[idx] <- list(res$factorization$inference_filtered_reconstruction_metrics)
                      run_tables$factorization_feature_scores[idx] <- list(res$factorization$feature_scores)
                      run_tables$factorization_module_auc[idx] <- list(res$factorization$module_auc)
                      run_tables$factorization_module_recall_curves[idx] <- list(res$factorization$module_recall_curves)
                      run_tables$factorization_module_recall_summary[idx] <- list(res$factorization$module_recall_summary)
                      run_tables$factorization_loading_enrichment[idx] <- list(res$factorization$loading_enrichment)
                      run_tables$factorization_loading_module_type_enrichment[idx] <- list(res$factorization$loading_module_type_enrichment)
                      run_tables$factorization_timing[idx] <- list(res$factorization$timing)
                      run_tables$factorization_sparsity_metrics[idx] <- list(res$factorization$sparsity_metrics)

                      run_tables$inference_feature_inference[idx] <- list(res$inference$feature_inference)
                      run_tables$inference_metrics[idx] <- list(res$inference$metrics)
                      run_tables$inference_module_auc[idx] <- list(res$inference$module_auc)
                      run_tables$inference_module_recall_curves[idx] <- list(res$inference$module_recall_curves)
                      run_tables$inference_module_recall_summary[idx] <- list(res$inference$module_recall_summary)
                      run_tables$inference_module_enrichment[idx] <- list(res$inference$module_enrichment)
                      run_tables$inference_module_type_enrichment[idx] <- list(res$inference$module_type_enrichment)
                      run_tables$inference_timing[idx] <- list(res$inference$timing)
                      run_tables$prediction_metrics[idx] <- list(res$prediction_metrics)

                      if (isTRUE(keep_runs)) runs[[idx]] <- res
                      rm(res)
                      idx <- idx + 1L
                    }
                  }
                }
              }
            }
          }
        }
      }

    }
  }
  errors <- dplyr::bind_rows(errors)

  factorization_module_auc <- dplyr::bind_rows(run_tables$factorization_module_auc)
  factorization_module_recall_summary <- dplyr::bind_rows(run_tables$factorization_module_recall_summary)
  inference_module_auc <- dplyr::bind_rows(run_tables$inference_module_auc)
  inference_module_recall_summary <- dplyr::bind_rows(run_tables$inference_module_recall_summary)

  out <- list(
    runs = if (keep_runs) runs else NULL,

    factorization_metrics = dplyr::bind_rows(run_tables$factorization_metrics),
    factorization_reconstruction_metrics = dplyr::bind_rows(run_tables$factorization_reconstruction_metrics),
    factorization_loading_auc_metrics = dplyr::bind_rows(run_tables$factorization_loading_auc_metrics),
    inference_filtered_reconstruction_metrics = dplyr::bind_rows(run_tables$inference_filtered_reconstruction_metrics),
    factorization_feature_scores = dplyr::bind_rows(run_tables$factorization_feature_scores),
    factorization_module_auc = factorization_module_auc,
    factorization_module_type_auc = summarize_module_type_auc(factorization_module_auc),
    factorization_module_recall_curves = dplyr::bind_rows(run_tables$factorization_module_recall_curves),
    factorization_module_recall_summary = factorization_module_recall_summary,
    factorization_module_type_recall = summarize_module_type_recall(factorization_module_recall_summary),
    factorization_loading_enrichment = dplyr::bind_rows(run_tables$factorization_loading_enrichment),
    factorization_loading_module_type_enrichment = dplyr::bind_rows(run_tables$factorization_loading_module_type_enrichment),
    factorization_timing = dplyr::bind_rows(run_tables$factorization_timing),
    factorization_sparsity_metrics = dplyr::bind_rows(run_tables$factorization_sparsity_metrics),

    inference_feature_inference = dplyr::bind_rows(run_tables$inference_feature_inference),
    inference_metrics = dplyr::bind_rows(run_tables$inference_metrics),
    inference_module_auc = inference_module_auc,
    inference_module_type_auc = summarize_module_type_auc(inference_module_auc),
    inference_module_recall_curves = dplyr::bind_rows(run_tables$inference_module_recall_curves),
    inference_module_recall_summary = inference_module_recall_summary,
    inference_module_type_recall = summarize_module_type_recall(inference_module_recall_summary),
    inference_module_enrichment = dplyr::bind_rows(run_tables$inference_module_enrichment),
    inference_module_type_enrichment = dplyr::bind_rows(run_tables$inference_module_type_enrichment),
    inference_timing = dplyr::bind_rows(run_tables$inference_timing),
    prediction_metrics = dplyr::bind_rows(run_tables$prediction_metrics),
    grid_errors = errors,
    grid = list(
      seeds = seeds,
      ng_grid = ng_grid,
      confounding_grid = confounding_grid,
      batch_to_biology_ratio_grid = batch_to_biology_ratio_grid,
      feature_noise_shape_grid = feature_noise_shape_grid,
      batch_biology_overlap_fraction_grid = batch_biology_overlap_fraction_grid,
      time_imbalance_grid = time_imbalance_grid,
      family_grid = family_grid,
      count_regime_mean_grid = count_regime_mean_grid,
      svdx_args = svdx_args,
      svdx_adjusted_use_z = svdx_adjusted_use_z
    )
  )

  class(out) <- "singleBenchmarkGrid"
  out
}
