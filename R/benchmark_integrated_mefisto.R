fit_mefisto_wrapper <- function(
    Y,
    xdata,
    covariates = NULL,
    K = 8,
    likelihoods = NULL,
    scale_views = FALSE,
    seed = 1,
    maxiter = 1000,
    outfile = NULL,
    verbose = FALSE,
    ...
) {
  if (!requireNamespace("MOFA2", quietly = TRUE)) {
    return(NULL)
  }

  if (is.null(covariates) || length(covariates) == 0L) {
    stop(
      "MEFISTO requires at least one explicit continuous covariate.",
      call. = FALSE
    )
  }
  covariates <- unique(as.character(covariates))
  missing_covariates <- setdiff(covariates, colnames(xdata))
  if (length(missing_covariates) > 0L) {
    stop(
      "MEFISTO covariates are missing from xdata: ",
      paste(missing_covariates, collapse = ", "),
      call. = FALSE
    )
  }
  non_numeric_covariates <- covariates[!vapply(xdata[covariates], is.numeric, logical(1))]
  if (length(non_numeric_covariates) > 0L) {
    stop(
      "MEFISTO covariates must be numeric/continuous. Non-numeric covariates: ",
      paste(non_numeric_covariates, collapse = ", "),
      call. = FALSE
    )
  }
  
  D <- length(Y)
  view_names <- names(Y) %||% paste0("dataset_", seq_len(D))
  
  if (is.null(likelihoods)) {
    likelihoods <- rep("gaussian", D)
  }
  if (length(likelihoods) == 1) {
    likelihoods <- rep(likelihoods, D)
  }
  
  sample_names <- colnames(Y[[1]])
  for (d in seq_along(Y)) {
    if (is.null(colnames(Y[[d]]))) {
      colnames(Y[[d]]) <- sample_names
    }
    if (!identical(colnames(Y[[d]]), sample_names)) {
      stop(
        "All Y datasets must have identical sample columns in the same order.",
        call. = FALSE
      )
    }
  }
  
  data_list <- Y
  names(data_list) <- view_names
  
  covariate_df <- xdata |>
    dplyr::mutate(sample = sample_names) |>
    dplyr::select(sample, dplyr::all_of(covariates)) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(covariates),
      names_to = "covariate",
      values_to = "value"
    )
  
  tm <- system.time({
    mofa <- MOFA2::create_mofa(data_list)
    
    mofa <- MOFA2::set_covariates(mofa, covariates = covariate_df)
    
    # mofa <- MOFA2::set_smooth_options(
    #   object = mofa,
    #   scale_cov = TRUE,
    #   start_opt = 10,
    #   opt_freq = 10,
    #   sparseGP = FALSE
    # )
    
    data_opts <- MOFA2::get_default_data_options(mofa)
    data_opts$scale_views <- scale_views
    
    model_opts <- MOFA2::get_default_model_options(mofa)
    model_opts$num_factors <- K
    
    train_opts <- MOFA2::get_default_training_options(mofa)
    train_opts$seed <- seed
    train_opts$maxiter <- maxiter
    train_opts$verbose <- verbose
    
    mefisto_opts <- MOFA2::get_default_mefisto_options(mofa)
    # mefisto_opts$covariates <- covariates
    # mefisto_opts$warping <- FALSE #does not work with two covariates?
    # mefisto_opts$sparseGP <- FALSE
    # mefisto_opts$frac_inducing <- 0.75
    # mefisto_opts$model_groups <- FALSE
    # mefisto_opts$new_values <- NULL
    
    mofa <- MOFA2::prepare_mofa(
      object = mofa,
      data_options = data_opts,
      model_options = model_opts,
      training_options = train_opts,
      mefisto_options = mefisto_opts
    )
    
    fit <- MOFA2::run_mofa(
      mofa,
      outfile = outfile,
      use_basilisk = FALSE
    )
  })
  
  factors <- MOFA2::get_factors(fit, factors = "all", as.data.frame = FALSE)
  weights <- MOFA2::get_weights(fit, views = "all", factors = "all", as.data.frame = FALSE)
  
  U <- as.matrix(factors[[1]])
  
  V_list <- lapply(view_names, function(vn) {
    W <- as.matrix(weights[[vn]])

    # MOFA/MEFISTO may internally modify feature names when the same
    # feature identifiers occur in more than one view (for example,
    # appending '_view1'/'_view2' in the duplicate-view benchmark).
    # Orient the loading matrix first, then restore the exact feature
    # identifiers supplied to the benchmark so that all downstream
    # feature-level evaluation uses the original feature namespace.
    if (nrow(W) != nrow(Y[[vn]]) && ncol(W) == nrow(Y[[vn]])) {
      W <- t(W)
    }

    if (nrow(W) != nrow(Y[[vn]])) {
      stop(
        "MEFISTO weight matrix for view '", vn,
        "' has ", nrow(W), " rows after orientation, but the input view has ",
        nrow(Y[[vn]]), " features.",
        call. = FALSE
      )
    }

    original_feature_ids <- rownames(Y[[vn]])
    if (!is.null(original_feature_ids)) {
      rownames(W) <- original_feature_ids
    }

    W
  })
  names(V_list) <- view_names
  
  # ------------------------------------------------------------
  # Benchmark normalization convention
  # ------------------------------------------------------------
  # MEFISTO provides shared factors and view-specific weights, but
  # the scale of factors/weights is not guaranteed to match the
  # integrated benchmark convention.
  #
  # For each component k, we use the same convention as the rest of
  # the integrated benchmark:
  #
  #   || c(v_1k, v_2k, ..., v_Dk) ||_2 = 1
  #
  # and transfer the concatenated loading norm into the shared U:
  #
  #   v_dk <- v_dk / t_k
  #   u_k  <- u_k  * t_k
  #
  # This preserves V_d %*% t(U) for each dataset while making
  # loading-based feature scores comparable across integrated methods.
  if (!exists(".normalize_integrated_shared_uv", mode = "function")) {
    stop(
      "Function '.normalize_integrated_shared_uv()' must be available before calling fit_mefisto_wrapper().",
      call. = FALSE
    )
  }

  normed <- .normalize_integrated_shared_uv(
    U = U,
    V_list = V_list
  )

  U <- normed$U
  V_list <- normed$V_list

  Yhat_list <- lapply(view_names, function(vn) {
    V_list[[vn]] %*% t(U)
  })
  names(Yhat_list) <- view_names
  
  as_integrated_benchmark_result(
    method = "MEFISTO",
    U = U,
    V_list = V_list,
    Yhat_list = Yhat_list,
    fit = fit,
    time = unname(tm["elapsed"]),
    extra = list(
      covariates = covariates,
      likelihoods = likelihoods,
      scale_views = scale_views
    )
  )
}