# ============================================================
# SVDX held-out prediction metrics
# ============================================================
# Computes SVDX-only held-out prediction MSE and percent test-data
# variation explained relative to a training-mean null model.

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

if (!exists(".benchmark_family_has_count", mode = "function")) {
  .benchmark_family_has_count <- function(family) {
    if (is.null(family)) return(FALSE)
    family_name <- if (is.character(family)) family[[1L]] else family$family %||% family$name %||% ""
    tolower(as.character(family_name)) %in% c("poisson", "quasipoisson")
  }
}

.benchmark_get_svdx_mse_fn <- function(required = TRUE) {
  pkg_opt <- getOption("svdxBenchmark.svdx_package", NULL)
  candidates <- unique(stats::na.omit(c(pkg_opt, "mySvdx", "svdx", "mSvdx", "SVDX")))
  for (pkg in candidates) {
    ns <- tryCatch(asNamespace(pkg), error = function(e) NULL)
    if (!is.null(ns) && exists("mse", envir = ns, inherits = FALSE)) {
      return(get("mse", envir = ns, inherits = FALSE))
    }
  }
  if (isTRUE(required)) {
    stop("Could not find mse() in svdx, SVDX, mySvdx, or mSvdx.", call. = FALSE)
  }
  NULL
}

.benchmark_get_svdx_predict_fn <- function(required = TRUE) {
  # Prefer normal S3 dispatch via stats::predict().  This helper exists only
  # to make failures explicit and future package renames harmless.
  if (exists("predict", envir = asNamespace("stats"), inherits = FALSE)) {
    return(stats::predict)
  }
  if (isTRUE(required)) stop("Could not locate stats::predict().", call. = FALSE)
  NULL
}

.benchmark_svdx_fit_slot <- function(fit, name, dataset = NULL) {
  val <- fit$fitAtOptCV[[name]] %||% fit[[name]]
  if (!is.null(dataset) && is.list(val)) {
    if (!dataset %in% names(val)) {
      stop("Could not locate `", name, "` for dataset '", dataset, "'.", call. = FALSE)
    }
    val <- val[[dataset]]
  }
  val
}

.benchmark_svdx_build_Z <- function(zformula = NULL, zdata = NULL, zeta = NULL, n_samples = NULL) {
  if (is.null(zformula) || is.null(zeta)) {
    q <- if (is.null(zeta)) 0L else ncol(as.matrix(zeta))
    n <- n_samples %||% if (!is.null(zdata)) nrow(zdata) else 0L
    Z <- matrix(0, nrow = n, ncol = q)
    if (!is.null(zeta) && !is.null(colnames(zeta))) colnames(Z) <- colnames(zeta)
    return(Z)
  }
  if (is.null(zdata)) stop("zdata is required to construct test Z from zformula.", call. = FALSE)
  Z <- stats::model.matrix(zformula, data = as.data.frame(zdata))
  zeta <- as.matrix(zeta)
  if (!is.null(colnames(zeta)) && !is.null(colnames(Z)) && !identical(colnames(zeta), colnames(Z))) {
    missing_cols <- setdiff(colnames(zeta), colnames(Z))
    extra_cols <- setdiff(colnames(Z), colnames(zeta))
    if (length(missing_cols) || length(extra_cols)) {
      stop(
        "Test Z columns do not match fitted Zeta columns. Missing: ",
        paste(missing_cols, collapse = ", "), "; extra: ", paste(extra_cols, collapse = ", "),
        call. = FALSE
      )
    }
    Z <- Z[, colnames(zeta), drop = FALSE]
  }
  Z
}

.benchmark_null_mse_train_mean <- function(Ytrain, Ytest) {
  Ytrain <- as.matrix(Ytrain)
  Ytest <- as.matrix(Ytest)
  mu <- rowMeans(Ytrain, na.rm = TRUE)
  resid <- sweep(Ytest, 1L, mu, "-")
  mean(resid^2, na.rm = TRUE)
}

.benchmark_svdx_single_prediction_metrics <- function(
    fit,
    Ytrain,
    Ytest,
    xdata_test,
    zdata_test = NULL,
    zformula = NULL,
    family = NULL,
    rowOffsets = NULL,
    colOffsets = NULL,
    method = "svdx",
    dataset = NA_character_,
    numThreads = NULL,
    ...
) {
  if (is.null(fit)) return(NULL)
  Ytrain <- as.matrix(Ytrain)
  Ytest <- as.matrix(Ytest)
  predict_fn <- .benchmark_get_svdx_predict_fn(required = TRUE)
  mse_fn <- .benchmark_get_svdx_mse_fn(required = TRUE)

  Uhat <- tryCatch(
    as.matrix(predict_fn(fit, newdata = as.data.frame(xdata_test), ...)),
    error = function(e) stop("SVDX prediction failed for method '", method, "': ", e$message, call. = FALSE)
  )
  V <- as.matrix(.benchmark_svdx_fit_slot(fit, "vmat"))
  Zeta <- .benchmark_svdx_fit_slot(fit, "Zeta")
  if (is.null(Zeta)) {
    Zeta <- matrix(0, nrow = nrow(Ytest), ncol = 0L)
  } else {
    Zeta <- as.matrix(Zeta)
  }
  Ztest <- .benchmark_svdx_build_Z(
    zformula = zformula,
    zdata = zdata_test,
    zeta = Zeta,
    n_samples = ncol(Ytest)
  )

  if (.benchmark_family_has_count(family)) {
    if (nrow(Uhat) != ncol(Ytest)) {
      stop("Predicted SVDX factors do not match the number of test samples.", call. = FALSE)
    }
    if (nrow(V) != nrow(Ytest)) {
      stop("Fitted SVDX loadings do not match the number of test features.", call. = FALSE)
    }

    eta_hat <- V %*% t(Uhat)
    if (ncol(Ztest) > 0L) {
      eta_hat <- eta_hat + Zeta %*% t(Ztest)
    }

    if (!is.null(rowOffsets)) {
      rowOffsets <- as.numeric(rowOffsets)
      if (length(rowOffsets) != nrow(Ytest) || any(!is.finite(rowOffsets))) {
        stop("rowOffsets must be finite and have length nrow(Ytest).", call. = FALSE)
      }
      eta_hat <- sweep(eta_hat, 1L, rowOffsets, "+")
    }
    if (!is.null(colOffsets)) {
      colOffsets <- as.numeric(colOffsets)
      if (length(colOffsets) != ncol(Ytest) || any(!is.finite(colOffsets))) {
        stop("colOffsets must be finite and have length ncol(Ytest).", call. = FALSE)
      }
      eta_hat <- sweep(eta_hat, 2L, colOffsets, "+")
    }

    mu_lower <- 1e-8
    mu_upper <- 1e10
    eta_hat <- pmin(pmax(eta_hat, log(mu_lower)), log(mu_upper))
    mu_hat <- exp(eta_hat)

    poisson_deviance <- function(Y, mu) {
      Y <- as.matrix(Y)
      mu <- pmin(pmax(as.matrix(mu), mu_lower), mu_upper)
      if (!identical(dim(Y), dim(mu))) {
        stop("Poisson prediction and test-data dimensions do not match.", call. = FALSE)
      }
      term <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
      positive <- Y > 0
      term[positive] <- Y[positive] * log(Y[positive] / mu[positive])
      2 * sum(term - (Y - mu), na.rm = TRUE)
    }

    test_deviance <- poisson_deviance(Ytest, mu_hat)
    null_mu <- matrix(rowMeans(Ytrain, na.rm = TRUE), nrow = nrow(Ytest), ncol = ncol(Ytest))
    null_deviance <- poisson_deviance(Ytest, null_mu)
    explained <- if (is.finite(null_deviance) && null_deviance > 0) {
      1 - (test_deviance / null_deviance)
    } else {
      NA_real_
    }

    return(tibble::tibble(
      method = method,
      dataset = dataset,
      test_mse = NA_real_,
      null_mse = NA_real_,
      explained_variation = NA_real_,
      explained_variation_percent = NA_real_,
      test_poisson_deviance = test_deviance,
      mean_test_poisson_deviance = test_deviance / length(Ytest),
      null_poisson_deviance = null_deviance,
      deviance_explained = explained,
      deviance_explained_percent = 100 * explained,
      n_train = ncol(Ytrain),
      n_test = ncol(Ytest),
      n_features = nrow(Ytest)
    ))
  }

  test_mse <- do.call(
    mse_fn,
    c(
      list(
        Y = Ytest,
        Z = Ztest,
        Zeta = Zeta,
        umat = Uhat,
        vmat = V
      ),
      if (!is.null(numThreads)) list(numThreads = numThreads) else list()
    )
  )
  test_mse <- as.numeric(test_mse)[1]
  null_mse <- .benchmark_null_mse_train_mean(Ytrain, Ytest)
  explained <- if (is.finite(null_mse) && null_mse > 0) 1 - (test_mse / null_mse) else NA_real_

  tibble::tibble(
    method = method,
    dataset = dataset,
    test_mse = test_mse,
    null_mse = null_mse,
    explained_variation = explained,
    explained_variation_percent = 100 * explained,
    test_poisson_deviance = NA_real_,
    mean_test_poisson_deviance = NA_real_,
    null_poisson_deviance = NA_real_,
    deviance_explained = NA_real_,
    deviance_explained_percent = NA_real_,
    n_train = ncol(Ytrain),
    n_test = ncol(Ytest),
    n_features = nrow(Ytest)
  )

}

.benchmark_isvdx_prediction_metrics <- function(
    fit,
    Ytrain,
    Ytest,
    xdata_test,
    zdata_test,
    zformula,
    method = "isvdx",
    numThreads = NULL,
    ...
) {
  if (is.null(fit)) return(NULL)
  dataset_names <- names(Ytest) %||% paste0("dataset_", seq_along(Ytest))
  names(Ytrain) <- names(Ytrain) %||% dataset_names
  names(Ytest) <- names(Ytest) %||% dataset_names
  names(zdata_test) <- names(zdata_test) %||% dataset_names
  names(zformula) <- names(zformula) %||% dataset_names
  mse_fn <- .benchmark_get_svdx_mse_fn(required = TRUE)
  predict_fn <- .benchmark_get_svdx_predict_fn(required = TRUE)

  out <- lapply(seq_along(dataset_names), function(i) {
    ds <- dataset_names[[i]]
    Ytr <- as.matrix(Ytrain[[ds]])
    Yte <- as.matrix(Ytest[[ds]])
    Uhat <- as.matrix(predict_fn(fit, newdata = as.data.frame(xdata_test), d = i, ...))
    V_list <- .benchmark_svdx_fit_slot(fit, "vmat")
    V <- as.matrix(V_list[[ds]] %||% V_list[[i]])
    Zeta_list <- .benchmark_svdx_fit_slot(fit, "Zeta")
    Zeta <- if (is.list(Zeta_list)) Zeta_list[[ds]] %||% Zeta_list[[i]] else Zeta_list
    if (is.null(Zeta)) Zeta <- matrix(0, nrow = nrow(Yte), ncol = 0L) else Zeta <- as.matrix(Zeta)
    Ztest <- .benchmark_svdx_build_Z(zformula = zformula[[ds]], zdata = zdata_test[[ds]], zeta = Zeta, n_samples = ncol(Yte))
    test_mse <- do.call(
      mse_fn,
      c(list(Y = Yte, Z = Ztest, Zeta = Zeta, umat = Uhat, vmat = V), if (!is.null(numThreads)) list(numThreads = numThreads) else list())
    )
    test_mse <- as.numeric(test_mse)[1]
    null_mse <- .benchmark_null_mse_train_mean(Ytr, Yte)
    explained <- if (is.finite(null_mse) && null_mse > 0) 1 - (test_mse / null_mse) else NA_real_
    tibble::tibble(
      method = method,
      dataset = ds,
      test_mse = test_mse,
      null_mse = null_mse,
      explained_variation = explained,
      explained_variation_percent = 100 * explained,
      n_train = ncol(Ytr),
      n_test = ncol(Yte),
      n_features = nrow(Yte)
    )
  })
  dplyr::bind_rows(out)
}

.benchmark_separate_svdx_prediction_metrics <- function(
    fits,
    Ytrain,
    Ytest,
    xdata_test,
    zdata_test,
    zformula,
    method = "separate_svdx",
    numThreads = NULL,
    ...
) {
  if (is.null(fits)) return(NULL)
  dataset_names <- names(Ytest) %||% paste0("dataset_", seq_along(Ytest))
  names(Ytrain) <- names(Ytrain) %||% dataset_names
  names(Ytest) <- names(Ytest) %||% dataset_names
  names(zdata_test) <- names(zdata_test) %||% dataset_names
  names(zformula) <- names(zformula) %||% dataset_names
  out <- lapply(dataset_names, function(ds) {
    fit_obj <- fits[[ds]] %||% fits[[match(ds, dataset_names)]]
    if (is.list(fit_obj) && !is.null(fit_obj$fit)) fit_obj <- fit_obj$fit
    .benchmark_svdx_single_prediction_metrics(
      fit = fit_obj,
      Ytrain = Ytrain[[ds]],
      Ytest = Ytest[[ds]],
      xdata_test = xdata_test,
      zdata_test = zdata_test[[ds]],
      zformula = zformula[[ds]],
      method = method,
      dataset = ds,
      numThreads = numThreads,
      ...
    )
  })
  dplyr::bind_rows(out)
}
