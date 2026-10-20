# .benchmark_subset_simulation_features
# simulate_svdx_interactions
# make_train_test_split
# .module_feature_sets_from_sim


.benchmark_subset_simulation_features <- function(sim, keep_features) {
  p_old <- nrow(sim$Y)

  if (!is.logical(keep_features) || length(keep_features) != p_old ||
      anyNA(keep_features)) {
    stop(
      "keep_features must be a non-missing logical vector of length nrow(sim$Y).",
      call. = FALSE
    )
  }
  keep <- keep_features

  if (!any(keep)) {
    stop("Feature filtering removed every simulated feature.", call. = FALSE)
  }
  if (all(keep)) return(sim)

  index_fields <- c(
    "signal_features",
    "null_features",
    "batch_features",
    "non_batch_features",
    "batch_biology_overlap_features",
    "batch_only_features",
    "biology_only_features",
    "clean_null_features"
  )

  memberships <- lapply(index_fields, function(field) {
    out <- rep(FALSE, p_old)
    idx <- sim$truth[[field]]
    if (length(idx) > 0L) out[idx] <- TRUE
    out[keep]
  })
  names(memberships) <- index_fields

  sim$Y <- sim$Y[keep, , drop = FALSE]
  sim$rowOffsets <- sim$rowOffsets[keep]

  matrix_fields <- c(
    "biology", "batch_effect", "linear_predictor", "mu"
  )
  for (field in matrix_fields) {
    if (!is.null(sim$truth[[field]])) {
      sim$truth[[field]] <- sim$truth[[field]][keep, , drop = FALSE]
    }
  }

  vector_fields <- c(
    "module_id", "feature_noise_var", "feature_noise_sd"
  )
  for (field in vector_fields) {
    if (!is.null(sim$truth[[field]])) {
      sim$truth[[field]] <- sim$truth[[field]][keep]
    }
  }

  for (field in index_fields) {
    sim$truth[[field]] <- which(memberships[[field]])
  }

  sim$truth$module_features <- lapply(
    seq_along(sim$truth$module_type),
    function(m) which(sim$truth$module_id == m)
  )

  sim$params$p <- nrow(sim$Y)
  sim
}


simulate_svdx_interactions <- function(
    p = 2000,
    G = 4,
    ng = 40,
    B = 6,
    n_modules = 8,
    pi_signal = 0.40,
    confounding = 0.4,
    time_imbalance = 0.4,
    severity_group_shift = 0.7,
    sigma = 1,
    batch_to_biology_ratio = 1,
    c_bio = 1,
    c_batch = c_bio,
    feature_noise_shape = 1,
    batch_biology_overlap_fraction = 0.25,
    module_sd = 1.5,
    phase_shift_sd = pi / 4,
    batch_feature_fraction = 0.10,
    family = NULL,
    count_regime_mean = NULL,
    rowOffsets = NULL,
    colOffsets = NULL,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  if (is.null(family)) family <- "gaussian"
  if (!is.character(family) || length(family) != 1L ||
      is.na(family) || !family %in% c("gaussian", "poisson")) {
    stop("family must be one of: gaussian, poisson.", call. = FALSE)
  }

  if (!is.null(count_regime_mean) &&
      (!is.numeric(count_regime_mean) || length(count_regime_mean) != 1L ||
       !is.finite(count_regime_mean) || count_regime_mean <= 0)) {
    stop(
      "count_regime_mean must be NULL or a single positive finite number.",
      call. = FALSE
    )
  }
  if (family == "poisson" && is.null(count_regime_mean)) {
    stop("count_regime_mean must be supplied for Poisson simulation.", call. = FALSE)
  }

  if (!is.numeric(feature_noise_shape) || length(feature_noise_shape) != 1L ||
      !is.finite(feature_noise_shape) || feature_noise_shape <= 0) {
    stop("feature_noise_shape must be a single positive finite number.", call. = FALSE)
  }
  if (!is.numeric(batch_biology_overlap_fraction) || length(batch_biology_overlap_fraction) != 1L ||
      !is.finite(batch_biology_overlap_fraction) ||
      batch_biology_overlap_fraction < 0 || batch_biology_overlap_fraction > 1) {
    stop("batch_biology_overlap_fraction must be a single number in [0, 1].", call. = FALSE)
  }

  batch_sd <- batch_to_biology_ratio * module_sd
  
  n <- G * ng
  group <- rep(seq_len(G), each = ng)
  
  base_t <- seq(0, 2 * pi, length.out = ng)
  t <- numeric(n)
  
  for (g in seq_len(G)) {
    idx <- which(group == g)
    jitter <- rnorm(ng, 0, time_imbalance)
    t[idx] <- sort((base_t + jitter) %% (2 * pi))
  }
  
  severity <- as.numeric(scale(rnorm(n, 0, 1) + severity_group_shift * as.numeric(scale(group))))
  
  batch_random <- sample(seq_len(B), n, replace = TRUE)
  batch_confounded <- ((group - 1) %% B) + 1
  use_confounded <- rbinom(n, 1, confounding)
  batch <- ifelse(use_confounded == 1, batch_confounded, batch_random)
  
  meta <- data.frame(
    sample = paste0("s", seq_len(n)),
    group = factor(group),
    batch = factor(batch),
    t = t,
    cos_t = cos(t),
    sin_t = sin(t),
    severity = severity,
    sev_cos_t = severity * cos(t),
    sev_sin_t = severity * sin(t)
  )
  
  if (family == "gaussian") {
    feature_noise_var <- sigma^2 * stats::rgamma(
      p,
      shape = feature_noise_shape,
      rate = feature_noise_shape
    )
    feature_noise_sd <- sqrt(feature_noise_var)
    gaussian_noise <- matrix(
      stats::rnorm(p * n, mean = 0, sd = rep(feature_noise_sd, times = n)),
      p, n
    )
  } else {
    feature_noise_var <- NULL
    feature_noise_sd <- NULL
    gaussian_noise <- NULL
  }
  
  module_id <- rep(NA_integer_, p)
  module_features <- vector("list", n_modules)
  available <- seq_len(p)
  module_size <- floor(p * pi_signal / n_modules)
  
  for (m in seq_len(n_modules)) {
    fm <- sample(available, module_size)
    module_features[[m]] <- fm
    module_id[fm] <- m
    available <- setdiff(available, fm)
  }
  
  signal_features <- which(!is.na(module_id))
  null_features <- which(is.na(module_id))
  
  module_type <- rep(
    c("group_mean", "shared_oscillation", "group_phase_oscillation", "severity_modulated_oscillation"),
    length.out = n_modules
  )
  module_type <- sample(module_type)
  
  biology <- matrix(0, p, n)
  module_signal <- matrix(0, n_modules, n)
  phase_by_group <- matrix(0, n_modules, G)
  
  for (m in seq_len(n_modules)) {
    fm <- module_features[[m]]
    loading <- rnorm(length(fm), mean = 1, sd = 0.15)
    loading <- loading * sample(c(-1, 1), length(fm), replace = TRUE)
    
    if (module_type[m] == "group_mean") {
      group_effect <- rnorm(G, 0, module_sd)
      sample_signal <- group_effect[group]
    }
    
    if (module_type[m] == "shared_oscillation") {
      amp <- abs(rnorm(1, module_sd, 0.25))
      phase <- runif(1, -pi, pi)
      sample_signal <- amp * cos(t + phase)
    }
    
    if (module_type[m] == "group_phase_oscillation") {
      amp <- abs(rnorm(1, module_sd, 0.25))
      phase_g <- rnorm(G, 0, phase_shift_sd)
      phase_by_group[m, ] <- phase_g
      sample_signal <- amp * cos(t + phase_g[group])
    }
    
    if (module_type[m] == "severity_modulated_oscillation") {
      amp0 <- abs(rnorm(1, module_sd, 0.25))
      amp1 <- rnorm(1, module_sd, 0.25)
      phase <- runif(1, -pi, pi)
      sample_signal <- (amp0 + amp1 * severity) * cos(t + phase)
    }
    
    module_signal[m, ] <- sample_signal
    biology[fm, ] <- loading %*% t(sample_signal)
  }

  biology <- c_bio * biology
  
  batch_effect <- matrix(0, p, n)

  n_batch_features <- as.integer(round(batch_feature_fraction * p))
  overlap_by_module <- vapply(
    module_features,
    function(fm) as.integer(round(batch_biology_overlap_fraction * length(fm))),
    integer(1)
  )
  n_overlap <- sum(overlap_by_module)

  if (n_overlap > n_batch_features) {
    stop(
      "Requested batch_biology_overlap_fraction cannot be achieved because ",
      "the implied number of batch-contaminated biological features exceeds ",
      "round(batch_feature_fraction * p).",
      call. = FALSE
    )
  }

  batch_biology_overlap_features <- unlist(
    lapply(seq_along(module_features), function(m) {
      fm <- module_features[[m]]
      nm <- overlap_by_module[[m]]
      if (length(fm) == 0L || nm == 0L) return(integer(0))
      sample(fm, size = nm, replace = FALSE)
    }),
    use.names = FALSE
  )
  batch_biology_overlap_features <- sort(unique(as.integer(batch_biology_overlap_features)))

  n_batch_only <- n_batch_features - length(batch_biology_overlap_features)
  if (n_batch_only > length(null_features)) {
    stop(
      "Requested batch_feature_fraction cannot be achieved after controlled ",
      "batch-biology overlap because there are not enough non-biological features ",
      "for batch-only effects.",
      call. = FALSE
    )
  }

  batch_only_features <- if (n_batch_only > 0L) {
    sample(null_features, size = n_batch_only, replace = FALSE)
  } else {
    integer(0)
  }

  batch_features <- sort(unique(c(batch_biology_overlap_features, batch_only_features)))
  non_batch_features <- setdiff(seq_len(p), batch_features)
  biology_only_features <- setdiff(signal_features, batch_biology_overlap_features)
  clean_null_features <- setdiff(null_features, batch_only_features)

  if (length(batch_features) > 0L) {
    batch_loading <- matrix(
      stats::rnorm(length(batch_features) * B, 0, batch_sd),
      nrow = length(batch_features),
      ncol = B
    )
    for (b in seq_len(B)) {
      idx <- which(batch == b)
      if (length(idx) > 0L) {
        batch_effect[batch_features, idx] <- batch_loading[, b, drop = FALSE] %*%
          t(rep(1, length(idx)))
      }
    }
  }

  batch_effect <- c_batch * batch_effect
  
  if (is.null(rowOffsets)) rowOffsets <- rep(0, p)
  if (!is.numeric(rowOffsets) || length(rowOffsets) != p ||
      any(!is.finite(rowOffsets))) {
    stop("rowOffsets must be NULL or a finite numeric vector of length p.", call. = FALSE)
  }

  if (is.null(colOffsets)) colOffsets <- rep(0, n)
  if (!is.numeric(colOffsets) || length(colOffsets) != n ||
      any(!is.finite(colOffsets))) {
    stop("colOffsets must be NULL or a finite numeric vector of length n.", call. = FALSE)
  }

  offset_matrix <- outer(rowOffsets, rep(1, n)) +
    outer(rep(1, p), colOffsets)

  if (family == "poisson") {
    linear_predictor <-
      log(count_regime_mean) + biology + batch_effect + offset_matrix

    if (any(!is.finite(linear_predictor))) {
      stop("The Poisson linear predictor contains non-finite values.", call. = FALSE)
    }

    mu <- exp(linear_predictor)
    if (any(!is.finite(mu)) || any(mu <= 0)) {
      stop(
        "The Poisson mean matrix is numerically invalid. Reduce ",
        "count_regime_mean, module_sd, batch_to_biology_ratio, or the supplied offsets.",
        call. = FALSE
      )
    }

    Y <- matrix(
      stats::rpois(length(mu), lambda = as.numeric(mu)),
      nrow = p,
      ncol = n
    )
    if (any(!is.finite(Y))) {
      stop(
        "Poisson sampling produced non-finite counts. Reduce the simulated mean scale.",
        call. = FALSE
      )
    }
  } else {
    linear_predictor <- biology + batch_effect + offset_matrix
    mu <- linear_predictor
    Y <- gaussian_noise + linear_predictor
  }
  
  rownames(Y) <- paste0("feature_", seq_len(p))
  colnames(Y) <- meta$sample
  rownames(linear_predictor) <- rownames(Y)
  colnames(linear_predictor) <- colnames(Y)
  rownames(mu) <- rownames(Y)
  colnames(mu) <- colnames(Y)
  names(rowOffsets) <- rownames(Y)
  names(colOffsets) <- colnames(Y)
  
  truth <- list(
    module_id = module_id,
    module_features = module_features,
    module_type = module_type,
    module_signal = module_signal,
    signal_features = signal_features,
    null_features = null_features,
    batch_features = batch_features,
    non_batch_features = non_batch_features,
    batch_biology_overlap_features = batch_biology_overlap_features,
    batch_only_features = batch_only_features,
    biology_only_features = biology_only_features,
    clean_null_features = clean_null_features,
    feature_noise_var = feature_noise_var,
    feature_noise_sd = feature_noise_sd,
    phase_by_group = phase_by_group,
    biology = biology,
    batch_effect = batch_effect,
    linear_predictor = linear_predictor,
    mu = mu
  )
  
  sim <- list(
    Y = Y,
    meta = meta,
    family = family,
    rowOffsets = rowOffsets,
    colOffsets = colOffsets,
    truth = truth,
    params = list(
      p_requested = p,
      p = p, n = n, G = G, ng = ng, B = B, n_modules = n_modules,
      pi_signal = pi_signal, confounding = confounding,
      time_imbalance = time_imbalance,
      severity_group_shift = severity_group_shift,
      sigma = sigma, batch_to_biology_ratio = batch_to_biology_ratio,
      c_bio = c_bio, c_batch = c_batch,
      feature_noise_shape = feature_noise_shape,
      batch_biology_overlap_fraction = batch_biology_overlap_fraction,
      batch_sd = batch_sd, module_sd = module_sd,
      batch_feature_fraction = batch_feature_fraction,
      family = family,
      count_regime_mean = if (family == "poisson") count_regime_mean else NULL
    )
  )

  if (family == "poisson") {
    sim <- .benchmark_subset_simulation_features(
      sim,
      keep_features = rowSums(sim$Y) > 0
    )
  }

  sim
}

make_train_test_split <- function(meta, test_frac = 0.30, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(meta)
  test_idx <- sample(seq_len(n), size = round(test_frac * n))
  list(train = setdiff(seq_len(n), test_idx), test = test_idx)
}

.module_feature_sets_from_sim <- function(sim) {
  feature_names <- rownames(sim$Y)
  sets <- lapply(seq_along(sim$truth$module_features), function(m) {
    feature_names[sim$truth$module_features[[m]]]
  })
  names(sets) <- paste0("M", seq_along(sets), "_", sim$truth$module_type)
  sets
}