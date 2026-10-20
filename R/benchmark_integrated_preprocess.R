if (!exists(".benchmark_family_has_count", mode = "function")) {
  .benchmark_family_has_count <- function(family) {
    if (is.null(family)) return(FALSE)
    any(unlist(family, use.names = FALSE) %in% c("poisson", "quasipoisson"), na.rm = TRUE)
  }
}

standardize_Y_list <- function(Y, truth = NULL) {
  std <- lapply(Y, standardize_features)
  
  Y_scaled <- lapply(std, `[[`, "Y")
  centers <- lapply(std, `[[`, "center")
  scales <- lapply(std, `[[`, "scale")
  
  out <- list(
    Y = Y_scaled,
    center = centers,
    scale = scales
  )
  
  if (!is.null(truth)) {
    truth2 <- truth
    
    # Put additive truth components on the same standardized scale as Y.
    # Center each truth component by its own feature-wise mean rather than by
    # the observed-Y center, which would mix other additive components into it.
    if (!is.null(truth$biology)) {
      truth2$biology <- Map(function(Bd, sdv) {
        component_center <- rowMeans(Bd, na.rm = TRUE)
        Bd <- sweep(Bd, 1, component_center, "-")
        sweep(Bd, 1, sdv, "/")
      }, truth$biology, scales)
    }
    
    if (!is.null(truth$batch_effect)) {
      truth2$batch_effect <- Map(function(Bd, sdv) {
        component_center <- rowMeans(Bd, na.rm = TRUE)
        Bd <- sweep(Bd, 1, component_center, "-")
        sweep(Bd, 1, sdv, "/")
      }, truth$batch_effect, scales)
    }
    
    out$truth <- truth2
  }
  
  out
}

# Leading-SV block normalization helper used by common preprocessing and MFA.
.leading_sv_scale_Y_list <- function(Y, target = c("geometric_mean", "one")) {
  target <- match.arg(target)
  if (!requireNamespace("irlba", quietly = TRUE)) {
    stop("Package 'irlba' is required for leading-SV block normalization.", call. = FALSE)
  }

  leading_sv <- vapply(Y, function(Yd) {
    sv <- irlba::irlba(as.matrix(Yd), nv = 1, nu = 1)$d[1]
    if (!is.finite(sv) || sv <= 0) 1 else sv
  }, numeric(1))

  # Equalize leading singular values without forcing them to 1.  The common
  # target is the geometric mean of the original leading singular values.
  # `block_scale` remains the divisor used by the existing fitter/inverse-
  # transform code, so Y_scaled = Y / block_scale and Y = Y_scaled * block_scale.
  target_sv <- if (identical(target, "one")) 1 else exp(mean(log(leading_sv)))
  block_scale <- leading_sv / target_sv

  list(
    Y = Map(function(Yd, s) Yd / s, Y, block_scale),
    block_scale = block_scale,
    leading_sv = leading_sv,
    target_sv = target_sv
  )
}


preprocess_isvdx_sim <- function(
    sim,
    standardize_features = TRUE,
    block_scale = c("none", "leading_sv"),
    sigma1_normalize_all = FALSE,
    family = NULL
) {
  # `block_scale` is retained for backwards compatibility only. Common block
  # normalization is controlled by sigma1_normalize_all; MFA performs its own
  # defining leading-SV normalization inside the MFA fitter.
  block_scale <- match.arg(block_scale)
  if (block_scale != "none") {
    stop(
      "Global block_scale is no longer supported; use sigma1_normalize_all for common normalization.",
      call. = FALSE
    )
  }

  Y <- sim$Y
  truth <- sim$truth
  if (is.null(family) && !is.null(sim$family)) family <- sim$family

  if (isTRUE(standardize_features) && .benchmark_family_has_count(family)) {
    message("standardize_features=TRUE ignored for poisson/quasipoisson family; using raw count data.")
    standardize_features <- FALSE
  }

  std <- NULL
  if (standardize_features) {
    std <- standardize_Y_list(Y, truth = truth)
    Y <- std$Y
    truth <- std$truth
  }

  sigma1_scale <- rep(1, length(Y))
  if (isTRUE(sigma1_normalize_all)) {
    if (!isTRUE(standardize_features)) {
      stop(
        "sigma1_normalize_all=TRUE requires standardize_features=TRUE.",
        call. = FALSE
      )
    }
    if (.benchmark_family_has_count(family)) {
      stop(
        "sigma1_normalize_all=TRUE is not supported for poisson/quasipoisson benchmarks.",
        call. = FALSE
      )
    }
    scaled <- .leading_sv_scale_Y_list(Y)
    Y <- scaled$Y
    sigma1_scale <- scaled$block_scale

    # Keep additive truth components on the same block-normalized scale as Y
    # so all reconstruction metrics remain on a common target scale within
    # this benchmark condition.
    if (!is.null(truth$biology)) {
      truth$biology <- Map(function(Bd, s) Bd / s, truth$biology, sigma1_scale)
    }
    if (!is.null(truth$batch_effect)) {
      truth$batch_effect <- Map(function(Bd, s) Bd / s, truth$batch_effect, sigma1_scale)
    }
  }

  sim$Y <- Y
  sim$truth <- truth
  sim$preprocess <- list(
    standardize_features = standardize_features,
    center = std$center %||% NULL,
    scale = std$scale %||% NULL,
    block_scale = sigma1_scale,
    block_scale_method = if (isTRUE(sigma1_normalize_all)) "leading_sv_all_methods" else "none",
    sigma1_normalize_all = isTRUE(sigma1_normalize_all)
  )

  sim
}
