.normalize_integrated_shared_uv <- function(U, V_list, eps = 1e-12) {
  # Integrated shared-U convention:
  # for each component k, the concatenated loading vector
  # c(v_1k, ..., v_Dk) has unit L2 norm, and the corresponding
  # scale is transferred to the shared u_k.
  U <- as.matrix(U)
  V_list <- lapply(V_list, as.matrix)

  K_use <- min(
    ncol(U),
    min(vapply(V_list, ncol, integer(1)))
  )

  if (!is.finite(K_use) || K_use < 1) {
    return(list(
      U = U[, FALSE, drop = FALSE],
      V_list = lapply(V_list, function(V) V[, FALSE, drop = FALSE])
    ))
  }

  U <- U[, seq_len(K_use), drop = FALSE]
  V_list <- lapply(V_list, function(V) {
    V[, seq_len(K_use), drop = FALSE]
  })

  for (k in seq_len(K_use)) {
    tk2 <- sum(vapply(V_list, function(V) {
      sum(V[, k]^2, na.rm = TRUE)
    }, numeric(1)))

    tk <- sqrt(tk2)

    if (is.finite(tk) && tk > eps) {
      V_list <- lapply(V_list, function(V) {
        V[, k] <- V[, k] / tk
        V
      })

      U[, k] <- U[, k] * tk
    }
  }

  list(U = U, V_list = V_list)
}


.normalize_integrated_list_uv <- function(U_list, V_list, eps = 1e-12) {
  # Dataset-specific-U analogue used for separate_svdx. For each component
  # index k, the concatenated loading vector across datasets is normalized,
  # and the same scale is transferred to every dataset-specific u_dk.
  U_list <- lapply(U_list, as.matrix)
  V_list <- lapply(V_list, as.matrix)

  K_use <- min(
    min(vapply(U_list, ncol, integer(1))),
    min(vapply(V_list, ncol, integer(1)))
  )

  if (!is.finite(K_use) || K_use < 1) {
    return(list(
      U_list = lapply(U_list, function(U) U[, FALSE, drop = FALSE]),
      V_list = lapply(V_list, function(V) V[, FALSE, drop = FALSE])
    ))
  }

  U_list <- lapply(U_list, function(U) U[, seq_len(K_use), drop = FALSE])
  V_list <- lapply(V_list, function(V) V[, seq_len(K_use), drop = FALSE])

  for (k in seq_len(K_use)) {
    tk2 <- sum(vapply(V_list, function(V) {
      sum(V[, k]^2, na.rm = TRUE)
    }, numeric(1)))

    tk <- sqrt(tk2)

    if (is.finite(tk) && tk > eps) {
      V_list <- lapply(V_list, function(V) {
        V[, k] <- V[, k] / tk
        V
      })

      U_list <- lapply(U_list, function(U) {
        U[, k] <- U[, k] * tk
        U
      })
    }
  }

  list(U_list = U_list, V_list = V_list)
}


# .loading_union_score <- function(U, V, eps = 1e-12) {
#   U <- as.matrix(U)
#   V <- as.matrix(V)
#
#   K_use <- min(ncol(U), ncol(V))
#   if (!is.finite(K_use) || K_use < 1) {
#     return(rep(NA_real_, nrow(V)))
#   }
#
#   U <- U[, seq_len(K_use), drop = FALSE]
#   V <- V[, seq_len(K_use), drop = FALSE]
#
#   d <- vapply(seq_len(K_use), function(k) {
#     uk_norm <- sqrt(sum(U[, k]^2, na.rm = TRUE))
#     if (!is.finite(uk_norm) || uk_norm <= eps) 0 else sqrt(uk_norm)
#   }, numeric(1))
#
#   sqrt(rowSums(sweep(V^2, 2, d, `*`), na.rm = TRUE) / K_use)
# }


.integrated_feature_scores <- function(U, V_list) {
  if (is.list(U)) {
    out <- lapply(seq_along(V_list), function(d) {
      .loading_union_score(U = U[[d]], V = V_list[[d]])
    })
    names(out) <- names(V_list)
    out
  } else {
    lapply(V_list, function(V) {
      .loading_union_score(U = U, V = V)
    })
  }
}

# ------------------------------------------------------------
# Integrated benchmark result object
# ------------------------------------------------------------

as_integrated_benchmark_result <- function(
    method,
    U,
    V_list,
    Yhat_list,
    feature_score_list = NULL,
    fit = NULL,
    time = NA_real_,
    extra = list()
) {
  if (is.null(feature_score_list)) {
    feature_score_list <- .integrated_feature_scores(U = U, V_list = V_list)
  }

  out <- list(
    method = method,
    U = U,
    V_list = V_list,
    Yhat_list = Yhat_list,
    feature_score_list = feature_score_list,
    fit = fit,
    time = time,
    shared_U = !is.list(U),
    extra = extra
  )

  class(out) <- c(
    "integratedBenchmarkResult",
    "list"
  )
  out
}


# ------------------------------------------------------------
# isvdx wrapper
# ------------------------------------------------------------

normalize_isvdx_family <- function(family, D, dataset_names = NULL) {
  if (is.null(family)) return(NULL)

  allowed <- c("gaussian", "poisson", "quasipoisson")
  dataset_names <- dataset_names %||% paste0("dataset_", seq_len(D))

  if (is.list(family)) {
    fam <- family
  } else if (is.character(family)) {
    fam <- as.list(family)
  } else {
    stop("family must be NULL, a character vector, or a list.", call. = FALSE)
  }

  if (length(fam) == 1L && D > 1L && (is.null(names(fam)) || !nzchar(names(fam)[[1L]]))) {
    fam <- rep(fam, D)
  } else if (!is.null(names(fam)) && all(nzchar(names(fam)))) {
    missing <- setdiff(dataset_names, names(fam))
    if (length(missing)) {
      stop("family is missing entries for dataset(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
    fam <- fam[dataset_names]
  } else if (length(fam) != D) {
    stop("family must be NULL, length 1, or length equal to the number of datasets.", call. = FALSE)
  }

  bad_null <- vapply(fam, is.null, logical(1))
  if (any(bad_null)) {
    stop("family entries cannot be NULL; use family = NULL for all default Gaussian behavior.", call. = FALSE)
  }

  ok <- vapply(fam, function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && x %in% allowed
  }, logical(1))

  if (!all(ok)) {
    stop("Each family entry must be one of: gaussian, poisson, quasipoisson.", call. = FALSE)
  }

  names(fam) <- dataset_names
  fam
}

fit_isvdx_wrapper <- function(
    Y,
    K,
    xformula,
    zformula,
    xdata,
    zdata,
    family = NULL,
    wd = NULL,
    graphData = NULL,
    graphDistFn = NULL,
    numThreads = 4,
    ...
) {
  D <- length(Y)
  dataset_names <- names(Y) %||% paste0("dataset_", seq_len(D))
  names(Y) <- dataset_names

  if (!is.null(zformula)) {
    names(zformula) <- dataset_names
    if (is.null(zdata)) {
      stop("zdata must be supplied when zformula is non-NULL.", call. = FALSE)
    }
    names(zdata) <- dataset_names
  }

  family <- normalize_isvdx_family(family, D, dataset_names = dataset_names)

  tm <- system.time({
    fit_args <- list(
      Y = Y,
      formula = xformula,
      data = xdata,
      K = K,
      numThreads = numThreads
    )

    if (!is.null(zformula)) {
      fit_args$zformula <- zformula
      fit_args$zdata <- zdata
    }

    if (!is.null(family)) fit_args$family <- family
    if (!is.null(wd)) fit_args$wd <- wd
    if (!is.null(graphData)) fit_args$graphData <- graphData
    if (!is.null(graphDistFn)) fit_args$graphDistFn <- graphDistFn

    fit <- do.call(isvdx, c(fit_args, list(...)))
  })

  U <- fit$fitAtOptCV$umat %||% fit$umat
  V_list <- fit$fitAtOptCV$vmat %||% fit$vmat

  if (!is.list(V_list)) {
    stop("Could not extract list-valued V matrices from isvdx fit.", call. = FALSE)
  }

  names(V_list) <- names(V_list) %||% dataset_names

  Yhat_list <- lapply(V_list, function(V) {
    V %*% t(U)
  })

  names(Yhat_list) <- names(V_list)

  as_integrated_benchmark_result(
    method = "isvdx",
    U = U,
    V_list = V_list,
    Yhat_list = Yhat_list,
    fit = fit,
    time = unname(tm["elapsed"]),
    extra = list(
      fit_summary = tryCatch(summary(fit, Y = Y), error = function(e) NULL)
    )
  )
}

# ------------------------------------------------------------
# Separate svdx per dataset
# ------------------------------------------------------------

fit_separate_svdx_wrapper <- function(
    Y,
    K,
    xformula,
    zformula,
    xdata,
    zdata,
    family = NULL,
    graphData = NULL,
    graphDistFn = NULL,
    numThreads = 4,
    ...
) {
  D <- length(Y)

  dataset_names <- names(Y) %||% paste0("dataset_", seq_len(D))
  names(Y) <- dataset_names
  family <- normalize_isvdx_family(family, D, dataset_names = dataset_names)

  tm <- system.time({
    fits <- lapply(seq_len(D), function(d) {

      family_d <- if (is.null(family)) NULL else family[[d]]

      fit_svdx_wrapper(
        Y = Y[[d]],
        K = K,
        xformula = xformula,
        zformula = zformula[[d]],
        data = xdata,
        zdata = zdata[[d]],
        family = family_d,
        graphData = graphData,
        graphDistFn = graphDistFn,
        numThreads = numThreads,
        ...
      )
    })
  })

  U_list <- lapply(fits, `[[`, "scores")
  V_list <- lapply(fits, `[[`, "loadings")

  names(V_list) <- names(Y)
  names(U_list) <- names(Y)

  normed <- .normalize_integrated_list_uv(
    U_list = U_list,
    V_list = V_list
  )
  U_list <- normed$U_list
  V_list <- normed$V_list

  Yhat_list <- Map(function(V, U) {
    V %*% t(U)
  }, V_list, U_list)
  names(Yhat_list) <- names(Y)

  as_integrated_benchmark_result(
    method = "separate_svdx",
    U = U_list,
    V_list = V_list,
    Yhat_list = Yhat_list,
    fit = fits,
    time = unname(tm["elapsed"]),
    extra = list(U_list = U_list)
  )
}


# ------------------------------------------------------------
# Concatenated baseline helpers
# ------------------------------------------------------------

.concat_Y <- function(Y) {
  do.call(rbind, Y)
}

.split_by_feature_blocks <- function(M, Y) {
  sizes <- vapply(Y, nrow, integer(1))
  ends <- cumsum(sizes)
  starts <- c(1L, head(ends, -1L) + 1L)

  out <- lapply(seq_along(Y), function(d) {
    M[seq.int(starts[d], ends[d]), , drop = FALSE]
  })

  names(out) <- names(Y)
  out
}


fit_concat_pca <- function(Y, K) {
  Y_concat <- .concat_Y(Y)
  fit <- fit_pca(Y_concat, K = K, center = FALSE)

  V_list <- .split_by_feature_blocks(fit$loadings, Y)
  Yhat_list <- .split_by_feature_blocks(fit$Yhat, Y)

  as_integrated_benchmark_result(
    method = "concat_PCA",
    U = fit$scores,
    V_list = V_list,
    Yhat_list = Yhat_list,
    fit = fit$fit,
    time = fit$time,
    extra = list(preprocessing = list(block_normalization = "none"))
  )
}


# Multiple Factor Analysis baseline: conventional MFA block weighting divides
# each incoming benchmark block by its leading singular value before a common
# PCA.  Unlike sigma1_normalize_all, this normalization is part of the method.
fit_mfa <- function(Y, K) {
  pp <- .leading_sv_scale_Y_list(Y, target = "one")
  Y_fit <- pp$Y
  Y_concat <- .concat_Y(Y_fit)
  fit <- fit_pca(Y_concat, K = K, center = FALSE)

  V_fit <- .split_by_feature_blocks(fit$loadings, Y_fit)
  # Return MFA loadings on the incoming benchmark scale while retaining the
  # single shared score matrix U.  If Y_d / sigma1_d ~= V_d^* U^T, then
  # Y_d ~= (sigma1_d V_d^*) U^T.
  V_list <- Map(function(Vd, s) Vd * s, V_fit, pp$block_scale)
  names(V_list) <- names(Y)
  Yhat_list <- lapply(V_list, function(Vd) Vd %*% t(fit$scores))

  as_integrated_benchmark_result(
    method = "MFA",
    U = fit$scores,
    V_list = V_list,
    Yhat_list = Yhat_list,
    fit = fit$fit,
    time = fit$time,
    extra = list(
      preprocessing = list(
        block_normalization = "mfa_leading_sv",
        block_scale = pp$block_scale,
        leading_sv = pp$leading_sv
      )
    )
  )
}


fit_concat_nmf <- function(Y, K, seed = 1) {
  Y_concat <- .concat_Y(Y)
  fit <- fit_nmf(Y_concat, K = K, seed = seed)
  if (is.null(fit)) return(NULL)

  V_list <- .split_by_feature_blocks(fit$loadings, Y)
  Yhat_list <- .split_by_feature_blocks(fit$Yhat, Y)

  as_integrated_benchmark_result(
    method = "concat_NMF",
    U = fit$scores,
    V_list = V_list,
    Yhat_list = Yhat_list,
    fit = fit$fit,
    time = fit$time,
    extra = list(
      preprocessing = list(
        block_normalization = "none",
        nonnegative_shift = fit$preprocessing$nonnegative_shift
      )
    )
  )
}


fit_concat_spca_grid <- function(
    Y,
    K,
    keepX_grid = round(c(0.01, 0.10, 0.25, 0.50, 0.75) * sum(vapply(Y, nrow, integer(1))))
) {
  Y_concat <- .concat_Y(Y)

  fits <- fit_mixomics_spca_grid(
    Y = Y_concat,
    K = K,
    keepX_grid = keepX_grid,
    scale = FALSE
  )

  out <- lapply(fits, function(fit) {
    V_list <- .split_by_feature_blocks(fit$loadings, Y)
    # Use the factorization returned by the method directly.
    Yhat_all <- as.matrix(fit$loadings) %*% t(as.matrix(fit$scores))
    Yhat_list <- .split_by_feature_blocks(Yhat_all, Y)

    as_integrated_benchmark_result(
      method = paste0("concat_", fit$method),
      U = fit$scores,
      V_list = V_list,
      Yhat_list = Yhat_list,
      fit = fit$fit,
      time = fit$time,
      extra = list(
        hyperparameter = fit$hyperparameter,
        hyperparameter_name = fit$hyperparameter_name,
        sparsity = fit$sparsity,
        preprocessing = list(block_normalization = "none"),
        reconstruction = "returned_factorization"
      )
    )
  })

  names(out) <- vapply(out, `[[`, character(1), "method")
  out
}


fit_concat_sparse_nmf_grid <- function(
    Y,
    K,
    lambda_grid = 10^c(-3, -2, -1),
    seed = 1
) {
  Y_concat <- .concat_Y(Y)

  fits <- fit_sparse_nmf_grid(
    Y = Y_concat,
    K = K,
    lambda_grid = lambda_grid,
    seed = seed
  )

  out <- lapply(fits, function(fit) {
    V_list <- .split_by_feature_blocks(fit$loadings, Y)
    Yhat_list <- .split_by_feature_blocks(fit$Yhat, Y)

    as_integrated_benchmark_result(
      method = paste0("concat_", fit$method),
      U = fit$scores,
      V_list = V_list,
      Yhat_list = Yhat_list,
      fit = fit$fit,
      time = fit$time,
      extra = list(
        hyperparameter = fit$hyperparameter,
        hyperparameter_name = fit$hyperparameter_name,
        sparsity = fit$sparsity,
        preprocessing = list(
          block_normalization = "none",
          nonnegative_shift = fit$preprocessing$nonnegative_shift
        )
      )
    )
  })

  names(out) <- vapply(out, `[[`, character(1), "method")
  out
}


# ------------------------------------------------------------
# Multiblock PLS / sparse PLS with covariate response block
# ------------------------------------------------------------

# .benchmark_block_pls_xmatrix <- function(xformula, xdata) {
#   X <- stats::model.matrix(xformula, data = as.data.frame(xdata))
#   storage.mode(X) <- "double"
#   if (ncol(X) < 1L) {
#     stop("xformula must produce at least one design column for block PLS.", call. = FALSE)
#   }
#   X
# }

.benchmark_block_pls_xmatrix <- function(xformula, xdata) {
  X <- stats::model.matrix(
    xformula,
    data = as.data.frame(xdata)
  )

  storage.mode(X) <- "double"

  # Remove constant design columns, including the intercept.
  keep <- apply(X, 2L, function(x) {
    all(is.finite(x)) && stats::sd(x) > 0
  })

  X <- X[, keep, drop = FALSE]

  if (ncol(X) < 1L) {
    stop(
      "xformula produces no nonconstant design columns for block PLS.",
      call. = FALSE
    )
  }

  X
}

# Reproduce mixOmics block scale=TRUE preprocessing on the design block only.
# Constant columns (including an intercept) become zero, matching mixOmics.
.benchmark_mixomics_scale_design <- function(X) {
  X <- as.matrix(X)
  mu <- colMeans(X, na.rm = TRUE)
  s <- vapply(seq_len(ncol(X)), function(j) stats::sd(X[, j], na.rm = TRUE), numeric(1))
  Xs <- sweep(X, 2L, mu, "-")
  nz <- is.finite(s) & s > 0
  if (any(nz)) Xs[, nz] <- sweep(Xs[, nz, drop = FALSE], 2L, s[nz], "/")
  if (any(!nz)) Xs[, !nz] <- 0
  dimnames(Xs) <- dimnames(X)
  Xs
}

.benchmark_block_pls_repackage <- function(fit, Y, dataset_names, method, elapsed, extra = list()) {
  U_list <- vector("list", length(dataset_names))
  V_list <- vector("list", length(dataset_names))
  Yhat_list <- vector("list", length(dataset_names))
  names(U_list) <- names(V_list) <- names(Yhat_list) <- dataset_names

  for (d in seq_along(dataset_names)) {
    ds <- dataset_names[[d]]
    U <- as.matrix(fit$variates[[ds]])
    V <- as.matrix(fit$loadings[[ds]])

    normed <- .normalize_loading_scale(U = U, V = V)
    U <- normed$U
    V <- normed$V

    if (nrow(U) != ncol(Y[[d]]) || nrow(V) != nrow(Y[[d]])) {
      stop(
        "Unexpected block PLS dimensions for dataset '", ds, "': expected U = samples x K and V = features x K.",
        call. = FALSE
      )
    }

    U_list[[ds]] <- U
    V_list[[ds]] <- V

    # Use the factorization returned by the method directly.
    Yhat_list[[ds]] <- V %*% t(U)
  }

  as_integrated_benchmark_result(
    method = method,
    U = U_list,
    V_list = V_list,
    Yhat_list = Yhat_list,
    fit = fit,
    time = elapsed,
    extra = c(extra, list(reconstruction = "returned_factorization"))
  )
}

fit_block_pls_wrapper <- function(Y, K, xformula, xdata, scale = FALSE) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package 'mixOmics' is required for block_PLS.", call. = FALSE)
  }

  dataset_names <- names(Y) %||% paste0("dataset_", seq_along(Y))
  X_response_raw <- .benchmark_block_pls_xmatrix(xformula, xdata)
  X_response <- .benchmark_mixomics_scale_design(X_response_raw)
  blocks <- lapply(Y, function(Yd) t(as.matrix(Yd)))
  names(blocks) <- dataset_names
  blocks$X_design <- X_response
  indY <- length(blocks)

  K_use <- min(
    as.integer(K),
    nrow(X_response) - 1L,
    min(vapply(blocks, ncol, integer(1)))
  )
  if (!is.finite(K_use) || K_use < 1L) return(NULL)

  tm <- system.time({
    fit <- mixOmics::block.pls(
      X = blocks,
      indY = indY,
      ncomp = K_use,
      design = "full",
      mode = "regression",
      scale = scale
    )
  })

  .benchmark_block_pls_repackage(
    fit = fit,
    Y = Y,
    dataset_names = dataset_names,
    method = "block_PLS",
    elapsed = unname(tm["elapsed"]),
    extra = list(
      design = "full", response_block = "X_design", scale = scale,
      design_preprocessing = "mixomics_center_scale_only_X_design"
    )
  )
}

fit_block_spls_grid <- function(
    Y,
    K,
    xformula,
    xdata,
    keepX_grid = c(0.01, 0.10, 0.25, 0.50, 0.75),
    scale = FALSE
) {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package 'mixOmics' is required for block_sPLS.", call. = FALSE)
  }

  dataset_names <- names(Y) %||% paste0("dataset_", seq_along(Y))
  X_response_raw <- .benchmark_block_pls_xmatrix(xformula, xdata)
  X_response <- .benchmark_mixomics_scale_design(X_response_raw)
  blocks <- lapply(Y, function(Yd) t(as.matrix(Yd)))
  names(blocks) <- dataset_names
  blocks$X_design <- X_response
  indY <- length(blocks)
  p_blocks <- vapply(Y, nrow, integer(1))

  K_use <- min(
    as.integer(K),
    nrow(X_response) - 1L,
    min(vapply(blocks, ncol, integer(1)))
  )
  if (!is.finite(K_use) || K_use < 1L) return(list())

  make_keep <- function(val) {
    if (length(val) == 1L && is.finite(val) && val > 0 && val <= 1) {
      pmax(1L, pmin(p_blocks, round(val * p_blocks)))
    } else if (length(val) == 1L && is.finite(val) && val > 1) {
      pmax(1L, pmin(p_blocks, as.integer(round(val))))
    } else if (length(val) == length(Y) && all(is.finite(val))) {
      pmax(1L, pmin(p_blocks, as.integer(round(val))))
    } else {
      stop(
        "Each block_sPLS keepX setting must be a fraction, a positive scalar count, or one count per biological block.",
        call. = FALSE
      )
    }
  }

  out <- lapply(seq_along(keepX_grid), function(i) {
    setting <- keepX_grid[[i]]
    keep <- make_keep(setting)

    # Sparse selection is applied only to the biological blocks.  The complete
    # xformula design is retained in the response block on every component.
    keep_list <- lapply(seq_along(keep), function(d) rep(keep[d], K_use))
    names(keep_list) <- dataset_names
    keep_list$X_design <- rep(ncol(X_response), K_use)

    tm <- system.time({
      fit <- mixOmics::block.spls(
        X = blocks,
        indY = indY,
        ncomp = K_use,
        keepX = keep_list,
        design = "full",
        mode = "regression",
        scale = scale
      )
    })

    label <- if (length(setting) == 1L) {
      format(setting, trim = TRUE, scientific = FALSE)
    } else {
      paste(setting, collapse = "_")
    }
    method <- paste0("block_sPLS_keepX_", label)

    sparsity <- dplyr::bind_rows(lapply(seq_along(dataset_names), function(d) {
      ds <- dataset_names[[d]]
      L <- as.matrix(fit$loadings[[ds]])
      tibble::tibble(
        dataset = ds,
        n_features = nrow(L),
        n_nonzero = sum(rowSums(abs(L)) > 0),
        sparsity = 1 - n_nonzero / nrow(L)
      )
    }))

    .benchmark_block_pls_repackage(
      fit = fit,
      Y = Y,
      dataset_names = dataset_names,
      method = method,
      elapsed = unname(tm["elapsed"]),
      extra = list(
        design = "full",
        response_block = "X_design",
        scale = scale,
        design_preprocessing = "mixomics_center_scale_only_X_design",
        hyperparameter = if (length(setting) == 1L) as.numeric(setting) else NA_real_,
        hyperparameter_name = "keepX",
        keepX = keep,
        sparsity = sparsity
      )
    )
  })

  names(out) <- vapply(out, `[[`, character(1), "method")
  out
}
