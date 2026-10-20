# fit_pca
# fit_svd
# fit_nmf
# fit_mixomics_pls_formula
# fit_feature_lm_formula
# fit_svdx_wrapper
# fit_glmpca_wrapper
# as_benchmark_from_fitter


# ------------------------------------------------------------
# Method fitters in common old-style format:
# list(method, scores, loadings, Yhat, time, fit)
# ------------------------------------------------------------
.normalize_loading_scale <- function(U, V, eps = 1e-12) {
  # For scale-ambiguous methods, enforce ||v_k||_2 = 1 and
  # transfer the loading scale into u_k. This keeps V %*% t(U)
  # unchanged while making loading-based feature scores comparable.
  U <- as.matrix(U)
  V <- as.matrix(V)

  K_use <- min(ncol(U), ncol(V))
  if (!is.finite(K_use) || K_use < 1) {
    return(list(U = U[, FALSE, drop = FALSE], V = V[, FALSE, drop = FALSE]))
  }

  U <- U[, seq_len(K_use), drop = FALSE]
  V <- V[, seq_len(K_use), drop = FALSE]

  for (k in seq_len(K_use)) {
    tk <- sqrt(sum(V[, k]^2, na.rm = TRUE))

    if (is.finite(tk) && tk > eps) {
      V[, k] <- V[, k] / tk
      U[, k] <- U[, k] * tk
    }
  }

  list(U = U, V = V)
}


.loading_union_score <- function(U, V, eps = 1e-12, weighted=FALSE) {
  # Feature-level loading score:
  #   score_j = sqrt((1 / K_use) * sum_k sqrt(||u_k||_2) * v_jk^2)
  #
  # Assumes V columns have already been normalized to unit L2 norm.
  U <- as.matrix(U)
  V <- as.matrix(V)

  score <- NA
  if(weighted){
    K_use <- min(ncol(U), ncol(V))
    if (!is.finite(K_use) || K_use < 1) {
      return(rep(NA_real_, nrow(V)))
    }

    U <- U[, seq_len(K_use), drop = FALSE]
    V <- V[, seq_len(K_use), drop = FALSE]

    d <- vapply(seq_len(K_use), function(k) {
      uk_norm <- sqrt(sum(U[, k]^2, na.rm = TRUE))
      if (!is.finite(uk_norm) || uk_norm <= eps) 0 else sqrt(uk_norm)
    }, numeric(1))

    score <- sqrt(rowSums(sweep(V^2, 2, d, `*`), na.rm = TRUE) / K_use)
  } else {
    score <- apply(abs(V), 1, max )
  }

  as.numeric(score)

}


fit_pca <- function(Y, K, center = TRUE) {
  tm <- system.time({
    pc <- stats::prcomp(t(Y), center = center, scale. = FALSE, rank. = K)
    k_use <- min(K, ncol(pc$x), ncol(pc$rotation))
    scores <- pc$x[, seq_len(k_use), drop = FALSE]
    loadings <- pc$rotation[, seq_len(k_use), drop = FALSE]
    Yhat <- t(scores %*% t(loadings))
    if (is.numeric(pc$center)) {
      Yhat <- sweep(Yhat, 1, pc$center, "+")
    }
  })

  list(
    method = "PCA",
    fit = pc,
    scores = scores,
    loadings = loadings,
    Yhat = Yhat,
    time = unname(tm["elapsed"])
  )
}

fit_svd <- function(Y, K) {
  tm <- system.time({
    s <- base::svd(Y, nu = K, nv = K)
    k_use <- min(K, length(s$d), ncol(s$u), ncol(s$v))
    loadings <- s$u[, seq_len(k_use), drop = FALSE]
    scores <- sweep(s$v[, seq_len(k_use), drop = FALSE], 2, s$d[seq_len(k_use)], "*")
    Yhat <- loadings %*% t(scores)
  })

  list(
    method = "SVD",
    fit = s,
    scores = scores,
    loadings = loadings,
    Yhat = Yhat,
    time = unname(tm["elapsed"])
  )
}

fit_nmf <- function(Y, K, seed = 1) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(NULL)
  if (!reticulate::py_module_available("sklearn")) return(NULL)

  sklearn_decomposition <- reticulate::import("sklearn.decomposition", delay_load = FALSE)

  ymin <- min(Y, na.rm = TRUE)
  Y_nonneg <- Y - ymin + 1e-6
  X <- t(Y_nonneg)

  tm <- system.time({
    fit <- sklearn_decomposition$NMF(
      n_components = as.integer(K)
    )
    W <- fit$fit_transform(X)
    H <- fit$components_

    normed <- .normalize_loading_scale(
      U = W,
      V = t(H)
    )
    scores <- normed$U
    loadings <- normed$V

    rownames(scores) <- colnames(Y)
    rownames(loadings) <- rownames(Y)
    component_names <- paste0("component_", seq_len(ncol(loadings)))
    colnames(scores) <- component_names
    colnames(loadings) <- component_names

    Yhat <- t(scores %*% t(loadings)) + ymin - 1e-6
  })

  list(
    method = "NMF",
    fit = fit,
    scores = scores,
    loadings = loadings,
    Yhat = Yhat,
    time = unname(tm["elapsed"]),
    preprocessing = list(nonnegative_shift = ymin - 1e-6)
  )
}

fit_mixomics_pls_formula <- function(Y, meta, K, xformula, scale = TRUE) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) return(NULL)

  X_design <- make_design(meta, xformula)

  #to avoid the zero variance columns for mixomics pls
  keep <- apply(X_design, 2, function(x) {
    stats::sd(x, na.rm = TRUE) > 0
  })
  X_design <- X_design[, keep, drop = FALSE]
  if (ncol(X_design) == 0) return(NULL)

  tm <- system.time({
    fit <- mixOmics::pls(
      X = t(Y),
      Y = X_design,
      ncomp = min(K, ncol(X_design), ncol(Y) - 1),
      mode = "regression",
      scale = scale
    )

    k_use <- ncol(fit$variates$X)
    scores <- fit$variates$X[, seq_len(k_use), drop = FALSE]
    loadings <- fit$loadings$X[, seq_len(k_use), drop = FALSE]

    normed <- .normalize_loading_scale(
      U = scores,
      V = loadings
    )
    scores <- normed$U
    loadings <- normed$V

    # Approximate latent reconstruction in transformed X-space.
    Yhat <- t(scores %*% t(loadings))
  })

  list(
    method = "mixOmics PLS",
    fit = fit,
    scores = scores,
    loadings = loadings,
    Yhat = Yhat,
    time = unname(tm["elapsed"])
  )
}

fit_feature_lm_formula <- function(Y, meta, xformula) {
  .require_namespace("limma")
  X <- make_design(meta, xformula)

  tm <- system.time({
    fit <- limma::lmFit(Y, X)
    fit <- limma::eBayes(fit)
    coef_mat <- fit$coefficients         # p x q
    p_mat <- fit$p.value                 # p x q
    Yhat <- coef_mat %*% t(X)            # p x n
  })

  list(
    method = "feature_lm",
    fit = fit,
    scores = X,
    loadings = t(coef_mat),
    coef_mat = coef_mat,
    p_mat = p_mat,
    Yhat = Yhat,
    time = unname(tm["elapsed"])
  )
}

# ------------------------------------------------------------
# Sparse PCA / sparse NMF fitters
# ------------------------------------------------------------

.effective_sparsity <- function(V, tol = 1e-8) {
  V <- as.matrix(V)

  tibble::tibble(
    mean_nonzero_fraction = mean(abs(V) > tol, na.rm = TRUE),
    mean_zero_fraction = mean(abs(V) <= tol, na.rm = TRUE),
    min_component_nonzero_fraction = min(colMeans(abs(V) > tol, na.rm = TRUE), na.rm = TRUE),
    max_component_nonzero_fraction = max(colMeans(abs(V) > tol, na.rm = TRUE), na.rm = TRUE)
  )
}

fit_mixomics_spca <- function(
    Y,
    K,
    keepX,
    scale = FALSE
) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) return(NULL)

  p <- nrow(Y)
  keepX <- pmax(1, pmin(as.integer(keepX), p))

  # mixOmics requires unique feature identifiers (columns of transposed Y).
  # Keep the numeric data unchanged; repair only missing/duplicate names.
  X <- t(Y)
  feature_ids <- colnames(X)
  if (is.null(feature_ids)) {
    feature_ids <- paste0("feature_", seq_len(ncol(X)))
  } else {
    missing_ids <- is.na(feature_ids) | !nzchar(feature_ids)
    feature_ids[missing_ids] <- paste0("feature_", which(missing_ids))
  }
  colnames(X) <- make.unique(feature_ids)

  tm <- system.time({
    fit <- mixOmics::spca(
      X = X,
      ncomp = K,
      keepX = rep(keepX, K),
      scale = scale
    )

    k_use <- min(K, ncol(fit$variates$X), ncol(fit$loadings$X))
    scores <- fit$variates$X[, seq_len(k_use), drop = FALSE]
    loadings <- fit$loadings$X[, seq_len(k_use), drop = FALSE]
    normed <- .normalize_loading_scale(
      U = scores,
      V = loadings
    )
    scores <- normed$U
    loadings <- normed$V

    Yhat <- t(scores %*% t(loadings))
  })

  sparsity <- .effective_sparsity(loadings)

  list(
    method = paste0("sPCA_keepX_", keepX),
    fit = fit,
    scores = scores,
    loadings = loadings,
    Yhat = Yhat,
    time = unname(tm["elapsed"]),
    hyperparameter = keepX,
    hyperparameter_name = "keepX",
    sparsity = sparsity
  )
}

fit_mixomics_spca_grid <- function(
    Y,
    K,
    keepX_grid = round(c(0.01, 0.25, 0.75) * nrow(Y)),
    scale = FALSE
) {
  keepX_grid <- sort(unique(pmax(1, pmin(as.integer(keepX_grid), nrow(Y)))))

  fits <- lapply(keepX_grid, function(kx) {
    fit_mixomics_spca(
      Y = Y,
      K = K,
      keepX = kx,
      scale = scale
    )
  })

  fits <- Filter(Negate(is.null), fits)
  names(fits) <- vapply(fits, `[[`, character(1), "method")
  fits
}

fit_sparse_nmf <- function(
    Y,
    K,
    lambda,
    seed = 1,
    method = "snmf/r"
) {
  fit_sparse_nmf_sklearn(
    Y = Y,
    K = K,
    lambda = lambda
  )
}

fit_sparse_nmf_grid <- function(
    Y,
    K,
    lambda_grid = 10^c(-3, -2, -1),
    seed = 1,
    method = "snmf/r"
) {
  fit_sparse_nmf_sklearn_grid(
    Y = Y,
    K = K,
    lambda_grid = lambda_grid
  )
}


# Fast sparse NMF for the single simulation benchmark.
# Uses sklearn defaults except for rank and feature-loading L1 regularization.
fit_sparse_nmf_sklearn <- function(
    Y,
    K,
    lambda
) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(NULL)
  if (!reticulate::py_module_available("sklearn")) return(NULL)

  sklearn_decomposition <- reticulate::import("sklearn.decomposition", delay_load = FALSE)

  ymin <- min(Y, na.rm = TRUE)
  Y_nonneg <- Y - ymin + 1e-6
  X <- t(Y_nonneg)

  tm <- system.time({
    fit <- sklearn_decomposition$NMF(
      n_components = as.integer(K),
      alpha_W = 0,
      alpha_H = as.numeric(lambda),
      l1_ratio = 1
    )

    W <- fit$fit_transform(X)
    H <- fit$components_

    normed <- .normalize_loading_scale(
      U = W,
      V = t(H)
    )
    scores <- normed$U
    loadings <- normed$V

    rownames(scores) <- colnames(Y)
    rownames(loadings) <- rownames(Y)
    component_names <- paste0("component_", seq_len(ncol(loadings)))
    colnames(scores) <- component_names
    colnames(loadings) <- component_names

    Yhat <- t(scores %*% t(loadings)) + ymin - 1e-6
  })

  sparsity <- .effective_sparsity(loadings)
  lambda_label <- gsub("\\.", "p", format(lambda, scientific = FALSE))

  list(
    method = paste0("sNMF_lambda_", lambda_label),
    fit = fit,
    scores = scores,
    loadings = loadings,
    Yhat = Yhat,
    time = unname(tm["elapsed"]),
    preprocessing = list(nonnegative_shift = ymin - 1e-6),
    hyperparameter = lambda,
    hyperparameter_name = "lambda",
    sparsity = sparsity
  )
}

fit_sparse_nmf_sklearn_grid <- function(
    Y,
    K,
    lambda_grid = 10^c(-3, -2, -1)
) {
  lambda_grid <- sort(unique(as.numeric(lambda_grid)))

  fits <- lapply(lambda_grid, function(lambda) {
    fit_sparse_nmf_sklearn(
      Y = Y,
      K = K,
      lambda = lambda
    )
  })

  fits <- Filter(Negate(is.null), fits)
  names(fits) <- vapply(fits, `[[`, character(1), "method")
  fits
}

fit_glmpca_wrapper <- function(
    Y,
    K,
    X = NULL,
    colOffsets = NULL,
    seed = 1,
    ...
) {
  if (!requireNamespace("glmpca", quietly = TRUE)) return(NULL)

  Y <- as.matrix(Y)
  if (any(!is.finite(Y)) || any(Y < 0) || any(Y != floor(Y))) {
    stop("GLM-PCA requires a finite non-negative integer count matrix.", call. = FALSE)
  }

  if (!is.null(X)) {
    X <- as.matrix(X)
    if (nrow(X) != ncol(Y)) {
      stop("X must have ncol(Y) rows, one per sample.", call. = FALSE)
    }
    if (any(!is.finite(X))) {
      stop("X must contain only finite values.", call. = FALSE)
    }
    if (ncol(X) == 0L) X <- NULL
  }

  if (!is.null(colOffsets)) {
    if (!is.numeric(colOffsets) || length(colOffsets) != ncol(Y) ||
        any(!is.finite(colOffsets))) {
      stop("colOffsets must be NULL or a finite numeric vector of length ncol(Y).", call. = FALSE)
    }
    size_factors <- exp(colOffsets)
    if (any(!is.finite(size_factors)) || any(size_factors <= 0)) {
      stop("exp(colOffsets) must contain finite positive GLM-PCA size factors.", call. = FALSE)
    }
  } else {
    size_factors <- rep(1, ncol(Y))
  }

  set.seed(seed)
  tm <- system.time({
    fit <- glmpca::glmpca(
      Y = Y,
      L = K,
      fam = "poi",
      X = X,
      sz = size_factors,
      ...
    )

    scores <- as.matrix(fit$factors)
    loadings <- as.matrix(fit$loadings)
    Yhat <- as.matrix(stats::predict(fit))
  })

  if (!all(dim(Yhat) == dim(Y))) {
    stop("predict.glmpca() did not return a matrix with the same dimensions as Y.", call. = FALSE)
  }

  normed <- .normalize_loading_scale(
    U = scores,
    V = loadings
  )

  list(
    method = "glm_pca",
    fit = fit,
    scores = normed$U,
    loadings = normed$V,
    Yhat = Yhat,
    time = unname(tm["elapsed"]),
    extra = list(
      family = "poisson",
      reconstruction_scale = "response_mean",
      size_factors = size_factors,
      nuisance_covariates = colnames(X) %||% character(0)
    )
  )
}


fit_svdx_wrapper <- function(
    Y,
    K,
    xformula,
    zformula = NULL,
    data,
    zdata = data,
    numThreads = 4,
    family = NULL,
    rowOffsets = NULL,
    colOffsets = NULL,
    graphData = NULL,
    graphDistFn = NULL,
    ...
) {
  tm <- system.time({
    # fit_args <- list(
    #   Y = Y,
    #   formula = xformula,
    #   data = data,
    #   K = K,
    #   numThreads = numThreads,
    #   family = family,
    #   rowOffsets = rowOffsets,
    #   colOffsets = colOffsets,
    #   graphData = graphData,
    #   graphDistFn = graphDistFn
    # )
    #
    # if (!is.null(zformula)) {
    #   fit_args$zformula <- zformula
    #   fit_args$zdata <- zdata
    # }

    fit_args <- list(
      Y = Y,
      formula = xformula,
      data = data,
      K = K,
      numThreads = numThreads,
      family = family,
      graphData = graphData,
      graphDistFn = graphDistFn
    )

    if (!is.null(rowOffsets)) {
      fit_args$rowOffsets <- rowOffsets
    }

    if (!is.null(colOffsets)) {
      fit_args$colOffsets <- colOffsets
    }

    if (!is.null(zformula)) {
      fit_args$zformula <- zformula
      fit_args$zdata <- zdata
    }

    fit <- do.call(svdx, c(fit_args, list(...)))
  })

  U <- fit$fitAtOptCV$umat
  V <- fit$fitAtOptCV$vmat
  Yhat <- V %*% t(U)

  reconstruction_scale <- "model"
  if (!is.null(family) && tolower(as.character(family)[1]) %in% c("poisson", "quasipoisson")) {
    eta_hat <- Yhat

    Zeta <- fit$fitAtOptCV$Zeta %||% NULL
    if (!is.null(zformula) && !is.null(Zeta)) {
      Zeta <- as.matrix(Zeta)
      Z <- stats::model.matrix(zformula, data = zdata)
      if (!is.null(colnames(Zeta)) && !is.null(colnames(Z))) {
        missing_cols <- setdiff(colnames(Zeta), colnames(Z))
        if (length(missing_cols) > 0) {
          stop("SVDX nuisance design is missing columns required by fit$fitAtOptCV$Zeta: ",
               paste(missing_cols, collapse = ", "), call. = FALSE)
        }
        Z <- Z[, colnames(Zeta), drop = FALSE]
      }
      eta_hat <- eta_hat + Zeta %*% t(Z)
    }

    if (!is.null(rowOffsets)) {
      eta_hat <- sweep(eta_hat, 1, rowOffsets, "+")
    }
    if (!is.null(colOffsets)) {
      eta_hat <- sweep(eta_hat, 2, colOffsets, "+")
    }

    eta_hat <- pmin(pmax(eta_hat, log(.Machine$double.xmin)), log(.Machine$double.xmax) - 1)
    Yhat <- exp(eta_hat)
    reconstruction_scale <- "response_mean"
  }

  list(
    method = "svdx",
    fit = fit,
    scores = U,
    loadings = V,
    Yhat = Yhat,
    time = unname(tm["elapsed"]),
    extra = list(reconstruction_scale = reconstruction_scale)
  )
}



.benchmark_resolve_mefisto_covariates <- function(covariates = NULL, data, xformula = NULL) {
  if (!is.null(covariates)) {
    covariates <- unique(as.character(covariates))
  } else {
    if (is.null(xformula)) {
      stop(
        "MEFISTO covariates were not supplied and cannot be inferred without xformula.",
        call. = FALSE
      )
    }
    candidates <- intersect(all.vars(xformula), colnames(data))
    covariates <- candidates[vapply(data[candidates], is.numeric, logical(1))]
  }

  if (length(covariates) == 0L) {
    stop(
      "No numeric MEFISTO covariates are available. Supply mefisto_covariates explicitly.",
      call. = FALSE
    )
  }

  missing_covariates <- setdiff(covariates, colnames(data))
  if (length(missing_covariates) > 0L) {
    stop(
      "MEFISTO covariates are missing from data: ",
      paste(missing_covariates, collapse = ", "),
      call. = FALSE
    )
  }

  non_numeric_covariates <- covariates[!vapply(data[covariates], is.numeric, logical(1))]
  if (length(non_numeric_covariates) > 0L) {
    stop(
      "MEFISTO covariates must be numeric/continuous. Non-numeric covariates: ",
      paste(non_numeric_covariates, collapse = ", "),
      call. = FALSE
    )
  }

  covariates
}

fit_mefisto_duplicate_view_wrapper <- function(
    Y,
    data,
    K,
    covariates = NULL,
    xformula = NULL,
    maxiter = 1000,
    seed = 1,
    ...
) {
  if (!exists("fit_mefisto_wrapper", mode = "function")) {
    stop(
      "Function 'fit_mefisto_wrapper()' must be available before calling ",
      "fit_mefisto_duplicate_view_wrapper().",
      call. = FALSE
    )
  }

  covariates <- .benchmark_resolve_mefisto_covariates(
    covariates = covariates,
    data = data,
    xformula = xformula
  )

  Y_list <- list(view1 = Y, view2 = Y)

  res <- fit_mefisto_wrapper(
    Y = Y_list,
    xdata = data,
    covariates = covariates,
    K = K,
    maxiter = maxiter,
    seed = seed,
    ...
  )

  if (is.null(res)) return(NULL)

  V <- res$V_list[[1]]
  Yhat <- res$Yhat_list[[1]]

  list(
    method = "mefisto_duplicate_view",
    fit = res$fit %||% NULL,
    scores = res$U,
    loadings = V,
    Yhat = Yhat,
    time = res$time %||% NA_real_,
    extra = list(
      source_method = res$method %||% "MEFISTO",
      duplicate_views = c("view1", "view2"),
      covariates = covariates,
      reconstruction_scale = "model"
    )
  )
}

# ------------------------------------------------------------
# Convert old-style fitter output to benchmark result
# ------------------------------------------------------------

as_benchmark_from_fitter <- function(x) {
  if (is.null(x)) return(NULL)

  if (identical(x$method, "feature_lm") || !is.null(x$coef_mat)) {
    stop(
      "feature_lm is an inference baseline and cannot be converted ",
      "to a factorization benchmark result.",
      call. = FALSE
    )
  }

  V <- x$loadings
  U <- x$scores

  feature_score <- if (is.null(V) || is.null(U) || all(is.na(V))) {
    rep(NA_real_, nrow(x$Yhat))
  } else {
    .loading_union_score(U = U, V = V)
  }

  as_factorization_benchmark_result(
    method = x$method,
    U = U,
    V = V,
    Yhat = x$Yhat,
    feature_score = feature_score,
    extra = list(fit = x$fit %||% NULL,
                 time = x$time,
                 hyperparameter = x$hyperparameter %||% NULL,
                 hyperparameter_name = x$hyperparameter_name %||% NULL,
                 sparsity = x$sparsity %||% NULL)
  )
}
