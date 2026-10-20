# ============================================================
# Integrated benchmark grid runner
# ============================================================
# Public API:
#   run_integrated_benchmark_once()
#   run_integrated_benchmark_grid()
#
# This file intentionally exposes one integrated grid runner.
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


.benchmark_integrated_overlap_feasibility <- function(
    p,
    n_modules,
    pi_signal,
    batch_feature_fraction,
    batch_biology_overlap_fraction
) {
  D <- length(p)
  recycle_arg <- function(x, name) {
    if (length(x) == 1L) return(rep(x, D))
    if (length(x) != D) {
      stop(name, " must have length 1 or length D for integrated overlap feasibility checks.", call. = FALSE)
    }
    x
  }

  pi_signal <- recycle_arg(pi_signal, "pi_signal")
  batch_feature_fraction <- recycle_arg(batch_feature_fraction, "batch_feature_fraction")
  batch_biology_overlap_fraction <- recycle_arg(
    batch_biology_overlap_fraction,
    "batch_biology_overlap_fraction"
  )

  for (d in seq_len(D)) {
    module_size <- floor(p[[d]] * pi_signal[[d]] / n_modules)
    n_signal <- n_modules * module_size
    n_null <- p[[d]] - n_signal
    n_batch_features <- as.integer(round(batch_feature_fraction[[d]] * p[[d]]))
    n_overlap <- n_modules * as.integer(round(batch_biology_overlap_fraction[[d]] * module_size))
    n_batch_only <- n_batch_features - n_overlap

    if (n_overlap > n_batch_features) {
      return(list(
        feasible = FALSE,
        reason = paste0(
          "dataset ", d, ": requested batch-biology overlap requires ", n_overlap,
          " biological batch-contaminated features, but batch_feature_fraction allows only ",
          n_batch_features, " total batch features"
        )
      ))
    }

    if (n_batch_only < 0L) {
      return(list(
        feasible = FALSE,
        reason = paste0(
          "dataset ", d,
          ": requested batch-biology overlap leaves a negative batch-only feature count"
        )
      ))
    }

    if (n_batch_only > n_null) {
      return(list(
        feasible = FALSE,
        reason = paste0(
          "dataset ", d, ": requested batch_feature_fraction requires ", n_batch_only,
          " batch-only features, but only ", n_null,
          " non-biological features are available"
        )
      ))
    }
  }

  list(feasible = TRUE, reason = NA_character_)
}

.benchmark_integrated_grid_overlap_feasibility <- function(
    batch_biology_overlap_fraction,
    once_args,
    simulator_args
) {
  get_arg <- function(name, default) {
    if (!is.null(once_args[[name]])) return(once_args[[name]])
    if (!is.null(simulator_args[[name]])) return(simulator_args[[name]])
    default
  }

  D <- get_arg("D", 3)
  p <- get_arg("p", c(2000, 1500, 1000))
  if (length(p) == 1L) p <- rep(p, D)

  .benchmark_integrated_overlap_feasibility(
    p = p,
    n_modules = get_arg("n_modules", 8),
    pi_signal = get_arg("pi_signal", rep(0.40, D)),
    batch_feature_fraction = get_arg("batch_feature_fraction", rep(0.10, D)),
    batch_biology_overlap_fraction = batch_biology_overlap_fraction
  )
}

.benchmark_add_integrated_grid_metadata <- function(
    x,
    seed,
    confounding_name,
    confounding,
    batch_to_biology_ratio_name,
    batch_to_biology_ratio,
    feature_noise_shape_name,
    feature_noise_shape,
    batch_biology_overlap_fraction_name,
    batch_biology_overlap_fraction,
    time_imbalance,
    sigma1_normalize_all,
    G,
    ng,
    K,
    n_modules,
    family = NULL,
    xformula = NULL
) {
  if (is.null(x) || nrow(x) == 0) return(x)

  x |>
    dplyr::mutate(
      seed = seed,
      confounding_name = confounding_name,
      confounding_pattern = paste(confounding, collapse = ","),
      batch_to_biology_ratio_name = batch_to_biology_ratio_name,
      batch_to_biology_ratio_pattern = paste(batch_to_biology_ratio, collapse = ","),
      feature_noise_shape_name = feature_noise_shape_name,
      feature_noise_shape_pattern = paste(feature_noise_shape, collapse = ","),
      batch_biology_overlap_fraction_name = batch_biology_overlap_fraction_name,
      batch_biology_overlap_fraction_pattern = paste(batch_biology_overlap_fraction, collapse = ","),
      time_imbalance = time_imbalance,
      sigma1_normalize_all = isTRUE(sigma1_normalize_all),
      G = G,
      ng = ng,
      n_samples = G * ng,
      K = K,
      n_modules = n_modules,
      family = if (is.null(family)) NA_character_ else paste(unlist(family, use.names = FALSE), collapse = ","),
      xformula = if (is.null(xformula)) NA_character_ else paste(deparse(xformula), collapse = "")
    )
}


.benchmark_make_integrated_train_test_split <- function(
    Y,
    train_fraction = 0.70,
    seed = NULL
) {
  n <- ncol(Y[[1]])
  all_cols <- seq_len(n)

  if (!is.null(seed)) set.seed(seed)

  train <- sort(sample(all_cols, size = round(train_fraction * n)))
  test <- setdiff(all_cols, train)

  list(
    train = train,
    test = test,
    train_names = colnames(Y[[1]])[train],
    test_names = colnames(Y[[1]])[test]
  )
}


.benchmark_subset_integrated_sim <- function(sim, idx) {
  truth <- sim$truth

  subset_matrix_cols <- function(x) {
    if (is.matrix(x) && ncol(x) >= max(idx)) x[, idx, drop = FALSE] else x
  }
  subset_list_matrix_cols <- function(x) {
    if (!is.list(x)) return(x)
    lapply(x, subset_matrix_cols)
  }

  for (nm in c(
    "biology", "batch_effect", "module_signal", "shared_biology",
    "dataset_specific_biology", "dataset_specific_module_signal", "total_biology"
  )) {
    if (!is.null(truth[[nm]])) truth[[nm]] <- subset_list_matrix_cols(truth[[nm]])
  }
  if (!is.null(truth$shared_module_signal) && is.matrix(truth$shared_module_signal) &&
      ncol(truth$shared_module_signal) >= max(idx)) {
    truth$shared_module_signal <- truth$shared_module_signal[, idx, drop = FALSE]
  }

  list(
    Y = lapply(sim$Y, function(Yd) Yd[, idx, drop = FALSE]),
    xdata = sim$xdata[idx, , drop = FALSE],
    data = sim$xdata[idx, , drop = FALSE],
    zdata = lapply(sim$zdata, function(zd) zd[idx, , drop = FALSE]),
    zformula = sim$zformula,
    family = sim$family,
    truth = truth
  )
}


.benchmark_fit_integrated_factorization_methods <- function(
    Y,
    xdata,
    zdata,
    zformula,
    family,
    fit_family = NULL,
    K,
    Y_base = NULL,
    xformula,
    methods,
    graphData = NULL,
    graphDistFn = NULL,
    spca_keepX_grid = NULL,
    sparse_nmf_lambda_grid = 10^c(-3, -2, -1),
    block_spls_keepX_grid = c(0.01, 0.10, 0.25, 0.50, 0.75),
    mefisto_covariates = c("t", "severity"),
    mefisto_maxiter = 1000,
    numThreads = 4,
    seed = 1,
    modelType = NULL,
    isvdx_adjusted_use_z = FALSE,
    isvdx_wd = NULL
) {
  results <- list()
  Y_for_baselines <- Y_base %||% if (.benchmark_family_has_count(family)) lapply(Y, log1p) else Y

  if ("isvdx" %in% methods) {
    results$isvdx <- fit_isvdx_wrapper(
      Y = Y,
      K = K,
      xformula = xformula,
      zformula = zformula,
      xdata = xdata,
      zdata = zdata,
      family = fit_family,
      wd = isvdx_wd,
      graphData = graphData,
      graphDistFn = graphDistFn,
      numThreads = numThreads,
      modelType = modelType
    )
  }

  if ("isvdx_adjusted" %in% methods) {
    if (any(family %in% c("poisson", "quasipoisson"), na.rm = TRUE)) {
      stop(
        "isvdx_adjusted is currently defined only for Gaussian benchmarks because ",
        "the adjusted baseline matrices are not count-scale observations.",
        call. = FALSE
      )
    }

    results$isvdx_adjusted <- fit_isvdx_wrapper(
      Y = Y_for_baselines,
      K = K,
      xformula = xformula,
      zformula = if (isTRUE(isvdx_adjusted_use_z)) zformula else rep(list(~1), length(Y_for_baselines)),
      xdata = xdata,
      zdata = zdata,
      family = fit_family,
      wd = isvdx_wd,
      graphData = graphData,
      graphDistFn = graphDistFn,
      numThreads = numThreads,
      modelType = modelType
    )
    results$isvdx_adjusted$method <- "isvdx_adjusted"
  }

  if ("separate_svdx" %in% methods) {
    results$separate_svdx <- fit_separate_svdx_wrapper(
      Y = Y,
      K = K,
      xformula = xformula,
      zformula = zformula,
      xdata = xdata,
      zdata = zdata,
      family = fit_family,
      graphData = graphData,
      graphDistFn = graphDistFn,
      numThreads = numThreads,
      modelType = modelType
    )
  }

  if ("MEFISTO" %in% methods) {
    results$MEFISTO <- fit_mefisto_wrapper(
      Y = Y_for_baselines,
      xdata = xdata,
      covariates = mefisto_covariates,
      K = K,
      maxiter = mefisto_maxiter,
      seed = seed
    )
  }

  if ("concat_PCA" %in% methods) results$concat_PCA <- fit_concat_pca(Y_for_baselines, K = K)
  if ("MFA" %in% methods) results$MFA <- fit_mfa(Y_for_baselines, K = K)
  if ("concat_NMF" %in% methods) results$concat_NMF <- fit_concat_nmf(Y_for_baselines, K = K, seed = seed)

  if ("concat_sPCA" %in% methods) {
    if (is.null(spca_keepX_grid)) {
      spca_keepX_grid <- round(c(0.01, 0.10, 0.25, 0.50, 0.75) * sum(vapply(Y_for_baselines, nrow, integer(1))))
    }

    results <- c(
      results,
      fit_concat_spca_grid(
        Y = Y_for_baselines,
        K = K,
        keepX_grid = spca_keepX_grid
      )
    )
  }

  if ("concat_sNMF" %in% methods) {
    results <- c(
      results,
      fit_concat_sparse_nmf_grid(
        Y = Y_for_baselines,
        K = K,
        lambda_grid = sparse_nmf_lambda_grid,
        seed = seed
      )
    )
  }

  if ("block_PLS" %in% methods) {
    results$block_PLS <- fit_block_pls_wrapper(
      Y = Y_for_baselines, K = K, xformula = xformula, xdata = xdata
    )
  }

  if ("block_sPLS" %in% methods) {
    results <- c(results, fit_block_spls_grid(
      Y = Y_for_baselines, K = K, xformula = xformula, xdata = xdata,
      keepX_grid = block_spls_keepX_grid
    ))
  }

  Filter(Negate(is.null), results)
}


.benchmark_fit_integrated_inference_isvdx <- function(
    Ytrain,
    Ytest,
    xdata_train,
    xdata_test,
    zdata_train,
    zdata_test,
    xformula,
    zformula,
    fit_family = NULL,
    inference_family = "gaussian",
    graphData = NULL,
    graphDistFn = NULL,
    trained_fit = NULL,
    K = 5,
    numThreads = 4,
    modelType = NULL,
    alpha = 0.1,
    z_adjustment = c("fitted", "reestimate"),
    isvdx_wd = NULL,
    return_design = TRUE,
    verbose = FALSE,
    ...
) {
  z_adjustment <- match.arg(z_adjustment)

  if (is.null(trained_fit)) {
    tm_fit <- system.time({
      fit <- isvdx(
        Y = Ytrain,
        K = K,
        formula = xformula,
        data = xdata_train,
        zformula = zformula,
        zdata = zdata_train,
        family = fit_family,
        wd = isvdx_wd,
        graphData = graphData,
        graphDistFn = graphDistFn,
        numThreads = numThreads,
        modelType = modelType,
        ...
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
      newdata = xdata_test,
      newzdata = zdata_test,
      zformula = zformula,
      family = inference_family,
      alpha = alpha,
      z_adjustment = z_adjustment,
      return_design = return_design,
      verbose = verbose
    )
  })

  list(
    fit = fit,
    inf = inf,
    time_fit = time_fit,
    time_inference = unname(tm_inf["elapsed"]),
    time_total = time_fit + unname(tm_inf["elapsed"])
  )
}


.benchmark_fit_integrated_inference_separate_svdx <- function(
    Ytrain,
    Ytest,
    xdata_train,
    xdata_test,
    zdata_train,
    zdata_test,
    xformula,
    zformula,
    fit_family = NULL,
    inference_family = "gaussian",
    graphData = NULL,
    graphDistFn = NULL,
    trained_fit = NULL,
    K = 5,
    numThreads = 4,
    modelType = NULL,
    alpha = 0.1,
    z_adjustment = c("fitted", "reestimate"),
    return_design = TRUE,
    verbose = FALSE,
    ...
) {
  z_adjustment <- match.arg(z_adjustment)

  dataset_names <- names(Ytrain) %||% paste0("dataset_", seq_along(Ytrain))

  fits <- list()
  infs <- list()
  time_fit <- numeric(length(dataset_names))
  time_inference <- numeric(length(dataset_names))
  names(time_fit) <- names(time_inference) <- dataset_names

  for (ds in dataset_names) {
    fit_fam_d <- if (is.null(fit_family)) NULL else fit_family[[ds]]
    inference_fam_d <- if (is.null(inference_family)) NULL else inference_family[[ds]]

    if (!is.null(trained_fit) && !is.null(trained_fit[[ds]])) {
      fits[[ds]] <- trained_fit[[ds]]
      time_fit[ds] <- 0
    } else {
      tm_fit <- system.time({
        fits[[ds]] <- svdx(
          Y = Ytrain[[ds]],
          K = K,
          formula = xformula,
          data = xdata_train,
          zformula = zformula[[ds]],
          zdata = zdata_train[[ds]],
          family = fit_fam_d,
          graphData = graphData,
          graphDistFn = graphDistFn,
          numThreads = numThreads,
          modelType = modelType,
          ...
        )
      })
      time_fit[ds] <- unname(tm_fit["elapsed"])
    }

    tm_inf <- system.time({
      infs[[ds]] <- inferV(
        fit = fits[[ds]],
        Ytest = Ytest[[ds]],
        newdata = xdata_test,
        newzdata = zdata_test[[ds]],
        zformula = zformula[[ds]],
        family = inference_fam_d,
        alpha = alpha,
        z_adjustment = z_adjustment,
        return_design = return_design,
        verbose = verbose
      )
    })

    time_inference[ds] <- unname(tm_inf["elapsed"])
  }

  list(
    fit = fits,
    inf = infs,
    time_fit = sum(time_fit, na.rm = TRUE),
    time_inference = sum(time_inference, na.rm = TRUE),
    time_total = sum(time_fit + time_inference, na.rm = TRUE),
    timing_by_dataset = tibble::tibble(
      dataset = dataset_names,
      time_fit = as.numeric(time_fit),
      time_inference = as.numeric(time_inference),
      time_total = as.numeric(time_fit + time_inference)
    )
  )
}


# ------------------------------------------------------------
# Public integrated runner
# ------------------------------------------------------------

run_integrated_benchmark_once <- function(
    seed = 1,
    D = 3,
    p = c(2000, 1500, 1000),
    G = 4,
    ng = 40,
    B = c(6, 5, 4),
    K = 8,
    n_modules = 8,
    xformula = ~ group * t + severity * t,
    confounding = rep(0.4, D),
    confounding_name = paste(confounding, collapse = "_"),
    time_imbalance = 0.4,
    batch_to_biology_ratio = rep(1, D),
    batch_to_biology_ratio_name = paste(batch_to_biology_ratio, collapse = "_"),
    feature_noise_shape = rep(1, D),
    feature_noise_shape_name = paste(feature_noise_shape, collapse = "_"),
    batch_biology_overlap_fraction = rep(0.25, D),
    batch_biology_overlap_fraction_name = paste(batch_biology_overlap_fraction, collapse = "_"),
    module_presence = NULL,
    pi_signal = NULL,
    include_dataset_specific_processes = FALSE,
    n_dataset_specific_modules = 0,
    pi_dataset_specific_signal = NULL,
    dataset_specific_module_presence = NULL,
    dataset_specific_module_sd = NULL,
    dataset_specific_strength = 0.5,
    shared_process_sd = NULL,
    severity_group_shift = 0.7,
    sigma = NULL,
    module_sd = NULL,
    phase_shift_sd = pi / 4,
    batch_feature_fraction = NULL,
    family = NULL,
    fit_family = NULL,
    inference_family = "gaussian",
    graphData = NULL,
    graphDistFn = NULL,
    factorization_methods = c("isvdx", "isvdx_adjusted", "separate_svdx", "MEFISTO", "concat_PCA", "MFA", "concat_NMF", "concat_sPCA", "concat_sNMF", "block_PLS", "block_sPLS"),
    inference_methods = c("isvdx", "isvdx_adjusted", "separate_svdx"),
    run_factorization = TRUE,
    run_inference = TRUE,
    standardize_features = TRUE,
    block_scale = c("none", "leading_sv"),
    sigma1_normalize_all = FALSE,
    spca_keepX_grid = NULL,
    sparse_nmf_lambda_grid = 10^c(-3, -2, -1),
    block_spls_keepX_grid = c(0.01, 0.10, 0.25, 0.50, 0.75),
    mefisto_covariates = c("t", "severity"),
    mefisto_maxiter = 1000,
    train_fraction = 0.70,
    alpha = 0.1,
    z_adjustment = c("fitted", "reestimate"),
    numThreads = 4,
    top_frac = seq(0.01, 0.50, by = 0.01),
    loading_top_n = 100,
    return_design = TRUE,
    verbose = FALSE,
    keep_fits = FALSE,
    simulator_args = list(),
    modelType = NULL,
    isvdx_adjusted_use_z = FALSE,
    ...
) {
  block_scale <- match.arg(block_scale)
  dataset_names <- paste0("dataset_", seq_len(D))
  family <- normalize_isvdx_family(family, D, dataset_names = dataset_names)
  fit_family <- normalize_isvdx_family(fit_family, D, dataset_names = dataset_names)
  inference_family <- normalize_isvdx_family(inference_family, D, dataset_names = dataset_names)

  z_adjustment <- match.arg(z_adjustment)

  valid_factorization_methods <- c("isvdx", "isvdx_adjusted", "separate_svdx", "MEFISTO", "concat_PCA", "MFA", "concat_NMF", "concat_sPCA", "concat_sNMF", "block_PLS", "block_sPLS")
  invalid_factorization_methods <- setdiff(factorization_methods, valid_factorization_methods)
  if (length(invalid_factorization_methods) > 0) {
    stop(
      "Unknown factorization_methods: ", paste(invalid_factorization_methods, collapse = ", "),
      ". Valid options are: ", paste(valid_factorization_methods, collapse = ", "),
      call. = FALSE
    )
  }

  valid_inference_methods <- c("isvdx", "isvdx_adjusted", "separate_svdx")
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

  simulator_args <- c(
    list(
      D = D,
      p = p,
      G = G,
      ng = ng,
      B = B,
      n_modules = n_modules,
      module_presence = module_presence,
      include_dataset_specific_processes = include_dataset_specific_processes,
      n_dataset_specific_modules = n_dataset_specific_modules,
      dataset_specific_module_presence = dataset_specific_module_presence,
      dataset_specific_strength = dataset_specific_strength,
      confounding = confounding,
      time_imbalance = time_imbalance,
      severity_group_shift = severity_group_shift,
      batch_to_biology_ratio = batch_to_biology_ratio,
      feature_noise_shape = feature_noise_shape,
      batch_biology_overlap_fraction = batch_biology_overlap_fraction,
      phase_shift_sd = phase_shift_sd,
      family = family,
      seed = seed
    ),
    if (!is.null(pi_signal)) list(pi_signal = pi_signal) else list(),
    if (!is.null(sigma)) list(sigma = sigma) else list(),
    if (!is.null(module_sd)) list(module_sd = module_sd) else list(),
    if (!is.null(shared_process_sd)) list(shared_process_sd = shared_process_sd) else list(),
    if (!is.null(pi_dataset_specific_signal)) list(pi_dataset_specific_signal = pi_dataset_specific_signal) else list(),
    if (!is.null(dataset_specific_module_sd)) list(dataset_specific_module_sd = dataset_specific_module_sd) else list(),
    if (!is.null(batch_feature_fraction)) list(batch_feature_fraction = batch_feature_fraction) else list(),
    simulator_args
  )

  sim <- do.call(simulate_isvdx_interactions, simulator_args)

  sim <- preprocess_isvdx_sim(
    sim,
    standardize_features = standardize_features,
    block_scale = block_scale,
    sigma1_normalize_all = sigma1_normalize_all,
    family = family
  )

  if (!is.null(graphData)) {
    if (!is.data.frame(graphData)) {
      stop("graphData must be a data.frame when supplied.", call. = FALSE)
    }
    if (nrow(graphData) < 1L) {
      stop("graphData must contain at least one row.", call. = FALSE)
    }
    if (ncol(graphData) != ncol(sim$Y[[1]])) {
      stop("ncol(graphData) must equal ncol(Y_list[[1]]).", call. = FALSE)
    }
  }

  svdx_model_type <- .benchmark_complete_model_type(
    modelType = modelType,
    formula = xformula,
    data = sim$xdata
  )

  split <- .benchmark_make_integrated_train_test_split(
    Y = sim$Y,
    train_fraction = train_fraction,
    seed = seed
  )
  sim_train <- .benchmark_subset_integrated_sim(sim, split$train)
  sim_test <- .benchmark_subset_integrated_sim(sim, split$test)
  sim_train_baseline_input <- if (.benchmark_family_has_count(sim$family)) lapply(sim_train$Y, log1p) else sim_train$Y
  sim_train_Y_base <- .benchmark_adjust_Y_list_by_formula(
    Y = sim_train_baseline_input,
    zdata = sim_train$zdata,
    zformula = sim$zformula
  )

  add_meta <- function(x) {
    .benchmark_add_integrated_grid_metadata(
      x,
      seed = seed,
      confounding_name = confounding_name,
      confounding = confounding,
      batch_to_biology_ratio_name = batch_to_biology_ratio_name,
      batch_to_biology_ratio = batch_to_biology_ratio,
      feature_noise_shape_name = feature_noise_shape_name,
      feature_noise_shape = feature_noise_shape,
      batch_biology_overlap_fraction_name = batch_biology_overlap_fraction_name,
      batch_biology_overlap_fraction = batch_biology_overlap_fraction,
      time_imbalance = time_imbalance,
      sigma1_normalize_all = sigma1_normalize_all,
      G = G,
      ng = ng,
      K = K,
      n_modules = n_modules,
      family = sim$family,
      xformula = xformula
    )
  }

  isvdx_wd <- if (isTRUE(sigma1_normalize_all)) rep(1, length(sim$Y)) else NULL

  factorization <- list()
  fact_fits <- NULL

  if (run_factorization) {
    fact_fits <- .benchmark_fit_integrated_factorization_methods(
      Y = sim_train$Y,
      xdata = sim_train$xdata,
      zdata = sim_train$zdata,
      zformula = sim$zformula,
      family = sim$family,
      fit_family = fit_family,
      K = K,
      Y_base = sim_train_Y_base,
      xformula = xformula,
      methods = factorization_methods,
      graphData = graphData,
      graphDistFn = graphDistFn,
      spca_keepX_grid = spca_keepX_grid,
      sparse_nmf_lambda_grid = sparse_nmf_lambda_grid,
      block_spls_keepX_grid = block_spls_keepX_grid,
      mefisto_covariates = mefisto_covariates,
      mefisto_maxiter = mefisto_maxiter,
      numThreads = numThreads,
      seed = seed,
      modelType = svdx_model_type,
      isvdx_adjusted_use_z = isvdx_adjusted_use_z,
      isvdx_wd = isvdx_wd
    )

    eval <- evaluate_integrated_factorization_results(
      results = fact_fits,
      sim = sim_train,
      top_frac = top_frac,
      enrichment_top_n = loading_top_n
    )

    factorization <- list(
      metrics = add_meta(eval$metrics),
      reconstruction_metrics = add_meta(eval$reconstruction_metrics),
      loading_auc_metrics = add_meta(.bench_loading_auc_alias_rows(eval$metrics, "isvdxK")),
      feature_scores = add_meta(eval$feature_scores),
      module_auc = add_meta(eval$module_auc),
      module_recall_curves = add_meta(eval$module_recall_curves),
      module_recall_summary = add_meta(eval$module_recall_summary),
      loading_enrichment = add_meta(eval$loading_enrichment),
      loading_module_type_enrichment = add_meta(eval$loading_module_type_enrichment),
      timing = add_meta(extract_integrated_timing(fact_fits)),
      sparsity_metrics = add_meta(extract_integrated_sparsity(fact_fits))
    )
  }

  inference <- list()
  inf_fits <- NULL

  sim_test_Y_base <- NULL
  if (run_inference && "isvdx_adjusted" %in% inference_methods) {
    sim_test_Y_base <- .benchmark_adjust_Y_list_by_formula(
      Y = sim_test$Y,
      zdata = sim_test$zdata,
      zformula = sim$zformula
    )
  }

  if (run_inference) {
    inf_fits <- list()

    if ("isvdx" %in% inference_methods) {
      inf_fits$isvdx <- .benchmark_fit_integrated_inference_isvdx(
        Ytrain = sim_train$Y,
        Ytest = sim_test$Y,
        xdata_train = sim_train$xdata,
        xdata_test = sim_test$xdata,
        zdata_train = sim_train$zdata,
        zdata_test = sim_test$zdata,
        xformula = xformula,
        zformula = sim$zformula,
        fit_family = fit_family,
        inference_family = inference_family,
        graphData = graphData,
        graphDistFn = graphDistFn,
        trained_fit = if (!is.null(fact_fits) && !is.null(fact_fits$isvdx)) fact_fits$isvdx$fit else NULL,
        K = K,
        numThreads = numThreads,
        modelType = svdx_model_type,
        alpha = alpha,
        z_adjustment = z_adjustment,
        isvdx_wd = isvdx_wd,
        return_design = return_design,
        verbose = verbose
      )
    }

    if ("isvdx_adjusted" %in% inference_methods) {
      if (.benchmark_family_has_count(sim$family)) {
        stop(
          "isvdx_adjusted inference is currently defined only for Gaussian benchmarks.",
          call. = FALSE
        )
      }

      adjusted_zformula <- if (isTRUE(isvdx_adjusted_use_z)) {
        sim$zformula
      } else {
        rep(list(~1), length(sim_train_Y_base))
      }

      inf_fits$isvdx_adjusted <- .benchmark_fit_integrated_inference_isvdx(
        Ytrain = sim_train_Y_base,
        Ytest = sim_test_Y_base,
        xdata_train = sim_train$xdata,
        xdata_test = sim_test$xdata,
        zdata_train = sim_train$zdata,
        zdata_test = sim_test$zdata,
        xformula = xformula,
        zformula = adjusted_zformula,
        fit_family = fit_family,
        inference_family = inference_family,
        graphData = graphData,
        graphDistFn = graphDistFn,
        trained_fit = if (!is.null(fact_fits) && !is.null(fact_fits$isvdx_adjusted)) {
          fact_fits$isvdx_adjusted$fit
        } else {
          NULL
        },
        K = K,
        numThreads = numThreads,
        modelType = svdx_model_type,
        alpha = alpha,
        z_adjustment = z_adjustment,
        isvdx_wd = isvdx_wd,
        return_design = return_design,
        verbose = verbose
      )
    }

    if ("separate_svdx" %in% inference_methods) {
      inf_fits$separate_svdx <- .benchmark_fit_integrated_inference_separate_svdx(
        Ytrain = sim_train$Y,
        Ytest = sim_test$Y,
        xdata_train = sim_train$xdata,
        xdata_test = sim_test$xdata,
        zdata_train = sim_train$zdata,
        zdata_test = sim_test$zdata,
        xformula = xformula,
        zformula = sim$zformula,
        fit_family = fit_family,
        inference_family = inference_family,
        graphData = graphData,
        graphDistFn = graphDistFn,
        trained_fit = if (!is.null(fact_fits) && !is.null(fact_fits$separate_svdx)) fact_fits$separate_svdx$fit else NULL,
        K = K,
        numThreads = numThreads,
        modelType = svdx_model_type,
        alpha = alpha,
        z_adjustment = z_adjustment,
        return_design = return_design,
        verbose = verbose
      )
    }

    evals <- lapply(names(inf_fits), function(method) {
      evaluate_integrated_svdx_inference(
        method = method,
        fit_obj = inf_fits[[method]],
        sim = sim_test,
        alpha = alpha,
        top_frac = top_frac
      )
    })
    names(evals) <- names(inf_fits)

    timing <- dplyr::bind_rows(lapply(names(inf_fits), function(method) {
      tibble::tibble(
        method = method,
        time_fit = inf_fits[[method]]$time_fit,
        time_inference = inf_fits[[method]]$time_inference,
        time_total = inf_fits[[method]]$time_total
      )
    }))

    inference <- list(
      feature_inference = add_meta(dplyr::bind_rows(lapply(evals, `[[`, "feature_inference"))),
      metrics = add_meta(dplyr::bind_rows(lapply(evals, `[[`, "metrics"))),
      module_auc = add_meta(dplyr::bind_rows(lapply(evals, `[[`, "module_auc"))),
      module_recall_curves = add_meta(dplyr::bind_rows(lapply(evals, `[[`, "module_recall_curves"))),
      module_recall_summary = add_meta(dplyr::bind_rows(lapply(evals, `[[`, "module_recall_summary"))),
      module_enrichment = add_meta(dplyr::bind_rows(lapply(evals, `[[`, "module_enrichment"))),
      module_type_enrichment = add_meta(dplyr::bind_rows(lapply(evals, `[[`, "module_type_enrichment"))),
      timing = add_meta(timing)
    )
  }

  if (run_factorization && run_inference && !is.null(fact_fits)) {
    filtered_all <- list()

    for (method_name in c("isvdx", "isvdx_adjusted")) {
      if (is.null(inf_fits[[method_name]]) || is.null(fact_fits[[method_name]])) next

      method_res <- fact_fits[[method_name]]
      filtered_rows <- lapply(seq_along(sim_train$Y), function(d) {
        ds <- names(sim_train$Y)[d]
        ev <- evals[[method_name]]$feature_inference
        fi <- if (!is.null(ev) && "dataset" %in% names(ev)) {
          dplyr::filter(ev, .data$dataset == ds)
        } else {
          ev
        }

        Vd <- as.matrix(.bench_get_dataset_item(method_res$V_list, ds, d, "result$V_list"))
        Ud <- if (is.list(method_res$U)) {
          as.matrix(.bench_get_dataset_item(method_res$U, ds, d, "result$U"))
        } else {
          as.matrix(method_res$U)
        }
        keep <- rep(FALSE, nrow(Vd))
        if (!is.null(fi) && all(c("feature", "global_q") %in% names(fi))) {
          idx <- match(rownames(Vd) %||% rownames(sim_train$Y[[d]]), fi$feature)
          keep <- !is.na(idx) & !is.na(fi$global_q[idx]) & fi$global_q[idx] < alpha
        }

        Vf <- Vd
        Vf[!keep, ] <- 0
        Yhat_f <- Vf %*% t(Ud)
        biology_d <- sim_train$truth$biology[[ds]] %||% sim_train$truth$biology[[d]]

        .bench_reconstruction_error_rows(
          method_name,
          Yhat_f,
          biology_d,
          sim_train$truth,
          "isvdxK",
          dataset = ds,
          index = d,
          reconstruction_source = "training_filtered_by_heldout_inference"
        )
      })

      filtered_all[[method_name]] <- dplyr::bind_rows(filtered_rows)
    }

    if (length(filtered_all) > 0L) {
      factorization$inference_filtered_reconstruction_metrics <-
        add_meta(dplyr::bind_rows(filtered_all))
    }
  }

  prediction_metrics <- list()
  if (!is.null(fact_fits) && !is.null(fact_fits$isvdx) && !is.null(fact_fits$isvdx$fit)) {
    prediction_metrics$isvdx <- tryCatch(
      .benchmark_isvdx_prediction_metrics(
        fit = fact_fits$isvdx$fit,
        Ytrain = sim_train$Y,
        Ytest = sim_test$Y,
        xdata_test = sim_test$xdata,
        zdata_test = sim_test$zdata,
        zformula = sim$zformula,
        method = "isvdx",
        numThreads = numThreads
      ),
      error = function(e) NULL
    )
  } else if (!is.null(inf_fits$isvdx) && !is.null(inf_fits$isvdx$fit)) {
    prediction_metrics$isvdx <- tryCatch(
      .benchmark_isvdx_prediction_metrics(
        fit = inf_fits$isvdx$fit,
        Ytrain = sim_train$Y,
        Ytest = sim_test$Y,
        xdata_test = sim_test$xdata,
        zdata_test = sim_test$zdata,
        zformula = sim$zformula,
        method = "isvdx",
        numThreads = numThreads
      ),
      error = function(e) NULL
    )
  }
  sep_fit <- NULL
  if (!is.null(fact_fits) && !is.null(fact_fits$separate_svdx) && !is.null(fact_fits$separate_svdx$fit)) sep_fit <- fact_fits$separate_svdx$fit
  else if (!is.null(inf_fits$separate_svdx) && !is.null(inf_fits$separate_svdx$fit)) sep_fit <- inf_fits$separate_svdx$fit
  if (!is.null(sep_fit)) {
    prediction_metrics$separate_svdx <- tryCatch(
      .benchmark_separate_svdx_prediction_metrics(
        fits = sep_fit,
        Ytrain = sim_train$Y,
        Ytest = sim_test$Y,
        xdata_test = sim_test$xdata,
        zdata_test = sim_test$zdata,
        zformula = sim$zformula,
        method = "separate_svdx",
        numThreads = numThreads
      ),
      error = function(e) NULL
    )
  }
  prediction_metrics <- dplyr::bind_rows(prediction_metrics)
  if (nrow(prediction_metrics) > 0L) prediction_metrics <- add_meta(prediction_metrics) else prediction_metrics <- NULL

  out <- list(
    sim = sim,
    split = split,
    factorization = factorization,
    inference = inference,
    prediction_metrics = prediction_metrics,
    fits = if (keep_fits) list(
      factorization = fact_fits,
      inference = inf_fits
    ) else NULL,
    xformula = xformula,
    zformula = sim$zformula,
    factorization_methods = factorization_methods,
    inference_methods = inference_methods,
    run_factorization = run_factorization,
    run_inference = run_inference
  )

  class(out) <- "integratedBenchmarkRun"
  out
}


run_integrated_benchmark_grid <- function(
    seeds = 1:3,
    ng_grid = c(5, 10, 25, 50),
    confounding_grid = list(
      all_clean = c(0, 0, 0),
      one_confounded = c(1, 0, 0),
      two_confounded = c(1, 1, 0),
      all_confounded = c(1, 1, 1)
    ),
    batch_to_biology_ratio_grid = list(
      balanced = c(1, 1, 1),
      strong_batch = c(10, 10, 10)
    ),
    feature_noise_shape_grid = list(
      shape_1 = c(1, 1, 1),
      shape_5 = c(5, 5, 5)
    ),
    batch_biology_overlap_fraction_grid = list(
      overlap_10 = c(0.10, 0.10, 0.10),
      overlap_25 = c(0.25, 0.25, 0.25),
      overlap_50 = c(0.50, 0.50, 0.50),
      overlap_75 = c(0.75, 0.75, 0.75)
    ),
    time_imbalance_grid = c(0, 0.5),
    family_grid = "gaussian",
    sigma1_normalize_all_grid = FALSE,
    fit_family = NULL,
    inference_family = "gaussian",
    mefisto_covariates = c("t", "severity"),
    spca_keepX_grid = NULL,
    sparse_nmf_lambda_grid = 10^c(-3, -2, -1),
    block_spls_keepX_grid = c(0.01, 0.10, 0.25, 0.50, 0.75),
    keep_runs = FALSE,
    simulator_args = list(),
    modelType = NULL,
    isvdx_adjusted_use_z = FALSE,
    ...
) {
  once_args <- .benchmark_collect_once_args(...)

  if ("isvdx_adjusted_use_z" %in% names(once_args)) {
    stop("Specify isvdx_adjusted_use_z through the dedicated grid argument, not through '...'.", call. = FALSE)
  }

  dedicated <- c("family", "fit_family", "inference_family", "sigma1_normalize_all", "mefisto_covariates", "spca_keepX_grid", "sparse_nmf_lambda_grid", "block_spls_keepX_grid")
  dup <- intersect(dedicated, names(once_args))
  if (length(dup)) stop("Specify ", paste(dup, collapse = ", "), " through dedicated grid arguments, not through '...'.", call. = FALSE)
  if (!is.null(simulator_args$family)) stop("Specify family through family_grid, not simulator_args.", call. = FALSE)
  family_grid <- unique(as.character(family_grid))
  if (!length(family_grid) || anyNA(family_grid) || any(!family_grid %in% c("gaussian", "poisson", "quasipoisson"))) {
    stop("family_grid must contain only gaussian, poisson, and/or quasipoisson.", call. = FALSE)
  }
  sigma1_normalize_all_grid <- unique(as.logical(sigma1_normalize_all_grid))
  if (!length(sigma1_normalize_all_grid) || anyNA(sigma1_normalize_all_grid)) {
    stop("sigma1_normalize_all_grid must contain only TRUE and/or FALSE.", call. = FALSE)
  }
  if (any(sigma1_normalize_all_grid) && any(family_grid %in% c("poisson", "quasipoisson"))) {
    stop("sigma1_normalize_all_grid=TRUE is supported only for gaussian benchmarks.", call. = FALSE)
  }
  once_args$mefisto_covariates <- mefisto_covariates
  once_args$spca_keepX_grid <- spca_keepX_grid
  once_args$sparse_nmf_lambda_grid <- sparse_nmf_lambda_grid
  once_args$block_spls_keepX_grid <- block_spls_keepX_grid

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
    for (sigma1_normalize_all in sigma1_normalize_all_grid) {
    for (ng in ng_grid) {
    for (confounding_name in names(confounding_grid)) {
      for (batch_to_biology_ratio_name in names(batch_to_biology_ratio_grid)) {
        for (feature_noise_shape_name in names(feature_noise_shape_grid)) {
          for (batch_biology_overlap_fraction_name in names(batch_biology_overlap_fraction_grid)) {
            for (time_imbalance in time_imbalance_grid) {
              for (seed in seeds) {
                confounding <- confounding_grid[[confounding_name]]
                batch_to_biology_ratio <- batch_to_biology_ratio_grid[[batch_to_biology_ratio_name]]
                feature_noise_shape <- feature_noise_shape_grid[[feature_noise_shape_name]]
                batch_biology_overlap_fraction <- batch_biology_overlap_fraction_grid[[batch_biology_overlap_fraction_name]]
                feasibility <- .benchmark_integrated_grid_overlap_feasibility(
                  batch_biology_overlap_fraction = batch_biology_overlap_fraction,
                  once_args = once_args,
                  simulator_args = simulator_args
                )
                if (!isTRUE(feasibility$feasible)) next

                run_args <- c(
                  list(
                    seed = seed,
                    ng = ng,
                    family = family_value,
                    sigma1_normalize_all = sigma1_normalize_all,
                    fit_family = fit_family,
                    inference_family = inference_family,
                    confounding = confounding,
                    confounding_name = confounding_name,
                    batch_to_biology_ratio = batch_to_biology_ratio,
                    batch_to_biology_ratio_name = batch_to_biology_ratio_name,
                    feature_noise_shape = feature_noise_shape,
                    feature_noise_shape_name = feature_noise_shape_name,
                    batch_biology_overlap_fraction = batch_biology_overlap_fraction,
                    batch_biology_overlap_fraction_name = batch_biology_overlap_fraction_name,
                    time_imbalance = time_imbalance,
                    simulator_args = simulator_args,
                    modelType = modelType,
                    isvdx_adjusted_use_z = isvdx_adjusted_use_z
                  ),
                  once_args
                )
                res <- tryCatch(
                  do.call(run_integrated_benchmark_once, run_args),
                  error = function(e) e
                )
                if (inherits(res, "error")) {
                  errors[[err_idx]] <- tibble::tibble(
                    seed = seed,
                    ng = ng,
                    family = family_value,
                    sigma1_normalize_all = sigma1_normalize_all,
                    confounding_name = confounding_name,
                    confounding_pattern = paste(confounding, collapse = ","),
                    batch_to_biology_ratio_name = batch_to_biology_ratio_name,
                    batch_to_biology_ratio_pattern = paste(batch_to_biology_ratio, collapse = ","),
                    feature_noise_shape_name = feature_noise_shape_name,
                    feature_noise_shape_pattern = paste(feature_noise_shape, collapse = ","),
                    batch_biology_overlap_fraction_name = batch_biology_overlap_fraction_name,
                    batch_biology_overlap_fraction_pattern = paste(batch_biology_overlap_fraction, collapse = ","),
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
      sigma1_normalize_all_grid = sigma1_normalize_all_grid,
      mefisto_covariates = mefisto_covariates,
      spca_keepX_grid = spca_keepX_grid,
      sparse_nmf_lambda_grid = sparse_nmf_lambda_grid,
      block_spls_keepX_grid = block_spls_keepX_grid,
      isvdx_adjusted_use_z = isvdx_adjusted_use_z
    )
  )

  class(out) <- "integratedBenchmarkGrid"
  out
}
