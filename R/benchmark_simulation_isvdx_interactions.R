# ------------------------------------------------------------
# Integrative SVDX simulator
# ------------------------------------------------------------
# Design:
#   - Matched samples across datasets.
#   - Dataset-specific feature spaces and batches.
#   - Shared latent biological processes are the primary target.
#   - Optional dataset-specific biological processes are nuisance structure.
#   - Evaluation should use truth$shared_* or the backward-compatible aliases
#     truth$biology, truth$signal_features, truth$module_features, which point
#     to the shared target by default.
# ------------------------------------------------------------

.recycle_to_length <- function(x, n, name = "argument") {
  if (length(x) == n) return(x)
  if (length(x) == 1) return(rep(x, n))
  stop(name, " must have length 1 or length D.", call. = FALSE)
}

.make_shared_isvdx_meta <- function(
    G = 4,
    ng = 40,
    time_imbalance = 0.4,
    severity_group_shift = 0.7
) {
  n <- G * ng
  group <- rep(seq_len(G), each = ng)

  base_t <- seq(0, 2 * pi, length.out = ng)
  t <- numeric(n)

  for (g in seq_len(G)) {
    idx <- which(group == g)
    jitter <- stats::rnorm(ng, 0, time_imbalance)
    t[idx] <- sort((base_t + jitter) %% (2 * pi))
  }

  severity <- as.numeric(scale(
    stats::rnorm(n, 0, 1) +
      severity_group_shift * as.numeric(scale(group))
  ))

  meta <- data.frame(
    sample = paste0("s", seq_len(n)),
    group = factor(group),
    t = t,
    cos_t = cos(t),
    sin_t = sin(t),
    severity = severity,
    sev_cos_t = severity * cos(t),
    sev_sin_t = severity * sin(t)
  )
  rownames(meta) <- meta$sample
  meta
}

.simulate_module_signal <- function(
    module_type,
    meta,
    G,
    module_sd = 1.5,
    phase_shift_sd = pi / 4
) {
  group <- as.integer(meta$group)
  t <- meta$t
  severity <- meta$severity

  phase_g <- rep(0, G)

  if (module_type == "group_mean") {
    group_effect <- stats::rnorm(G, 0, module_sd)
    sample_signal <- group_effect[group]
  }

  if (module_type == "shared_oscillation") {
    amp <- abs(stats::rnorm(1, module_sd, 0.25))
    phase <- stats::runif(1, -pi, pi)
    sample_signal <- amp * cos(t + phase)
  }

  if (module_type == "group_phase_oscillation") {
    amp <- abs(stats::rnorm(1, module_sd, 0.25))
    phase_g <- stats::rnorm(G, 0, phase_shift_sd)
    sample_signal <- amp * cos(t + phase_g[group])
  }

  if (module_type == "severity_modulated_oscillation") {
    amp0 <- abs(stats::rnorm(1, module_sd, 0.25))
    amp1 <- stats::rnorm(1, module_sd, 0.25)
    phase <- stats::runif(1, -pi, pi)
    sample_signal <- (amp0 + amp1 * severity) * cos(t + phase)
  }

  list(
    sample_signal = sample_signal,
    phase_by_group = phase_g
  )
}

.simulate_dataset_batches <- function(
    meta,
    B = 6,
    confounding = 0.4
) {
  n <- nrow(meta)
  group <- as.integer(meta$group)

  batch_random <- sample(seq_len(B), n, replace = TRUE)
  batch_confounded <- ((group - 1) %% B) + 1
  use_confounded <- stats::rbinom(n, 1, confounding)

  ifelse(use_confounded == 1, batch_confounded, batch_random)
}

.sample_module_features <- function(
    available,
    n_modules,
    module_size,
    module_presence
) {
  module_features <- vector("list", n_modules)
  module_id <- integer(0)

  if (n_modules == 0) {
    return(list(
      module_features = module_features,
      module_id = module_id,
      available = available
    ))
  }

  pd <- if (length(available) > 0) max(available) else 0
  module_id <- rep(NA_integer_, pd)

  for (m in seq_len(n_modules)) {
    if (!module_presence[m] || module_size <= 0 || length(available) == 0) {
      module_features[[m]] <- integer(0)
      next
    }

    size_m <- min(module_size, length(available))
    fm <- sample(available, size_m)
    module_features[[m]] <- fm
    module_id[fm] <- m
    available <- setdiff(available, fm)
  }

  list(
    module_features = module_features,
    module_id = module_id,
    available = available
  )
}

#' Simulate integrative SSANOVA-SVD benchmark data
#'
#' The primary integrated benchmark target is recovery of shared biological
#' processes across matched samples. Optional dataset-specific processes can be
#' added as nuisance biological structure. Batches are always dataset-specific.
#'
#' @export
simulate_isvdx_interactions <- function(
    D = 3,
    p = c(2000, 1500, 1000),
    G = 4,
    ng = 40,
    B = c(6, 5, 4),
    n_modules = 8,
    pi_signal = c(0.40, 0.35, 0.30),
    module_presence = NULL,
    include_dataset_specific_processes = FALSE,
    n_dataset_specific_modules = NULL,
    pi_dataset_specific_signal = 0.10,
    dataset_specific_module_presence = NULL,
    dataset_specific_module_sd = NULL,
    dataset_specific_strength = 1,
    confounding = c(0.4, 0.4, 0.4),
    time_imbalance = 0.4,
    severity_group_shift = 0.7,
    sigma = c(1, 1, 1),
    batch_to_biology_ratio = c(1, 1, 1),
    c_bio = 1,
    c_batch = c_bio,
    feature_noise_shape = 1,
    batch_biology_overlap_fraction = 0.25,
    module_sd = c(1.5, 1.5, 1.5),
    shared_process_sd = NULL,
    phase_shift_sd = pi / 4,
    batch_feature_fraction = c(0.10, 0.10, 0.10),
    family = NULL,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  p <- .recycle_to_length(p, D, "p")
  B <- .recycle_to_length(B, D, "B")
  pi_signal <- .recycle_to_length(pi_signal, D, "pi_signal")
  pi_dataset_specific_signal <- .recycle_to_length(
    pi_dataset_specific_signal, D, "pi_dataset_specific_signal"
  )
  confounding <- .recycle_to_length(confounding, D, "confounding")
  sigma <- .recycle_to_length(sigma, D, "sigma")
  batch_to_biology_ratio <- .recycle_to_length(
    batch_to_biology_ratio, D, "batch_to_biology_ratio"
  )
  c_bio <- .recycle_to_length(c_bio, D, "c_bio")
  c_batch <- .recycle_to_length(c_batch, D, "c_batch")
  feature_noise_shape <- .recycle_to_length(
    feature_noise_shape, D, "feature_noise_shape"
  )
  if (any(!is.finite(feature_noise_shape)) || any(feature_noise_shape <= 0)) {
    stop("feature_noise_shape values must be positive and finite.", call. = FALSE)
  }
  batch_biology_overlap_fraction <- .recycle_to_length(
    batch_biology_overlap_fraction, D, "batch_biology_overlap_fraction"
  )
  if (any(!is.finite(batch_biology_overlap_fraction)) ||
      any(batch_biology_overlap_fraction < 0) ||
      any(batch_biology_overlap_fraction > 1)) {
    stop("batch_biology_overlap_fraction values must be in [0, 1].", call. = FALSE)
  }
  module_sd <- .recycle_to_length(module_sd, D, "module_sd")
  if (is.null(dataset_specific_module_sd)) {
    dataset_specific_module_sd <- module_sd
  }
  dataset_specific_module_sd <- .recycle_to_length(
    dataset_specific_module_sd, D, "dataset_specific_module_sd"
  )
  batch_feature_fraction <- .recycle_to_length(
    batch_feature_fraction, D, "batch_feature_fraction"
  )

  family <- normalize_isvdx_family(family, D, dataset_names =paste0("dataset_", seq_len(D)))

  if (is.null(shared_process_sd)) {
    shared_process_sd <- mean(module_sd)
  }

  if (is.null(n_dataset_specific_modules)) {
    n_dataset_specific_modules <- if (include_dataset_specific_processes) n_modules else 0L
  }
  n_dataset_specific_modules <- as.integer(n_dataset_specific_modules)
  if (!include_dataset_specific_processes) {
    n_dataset_specific_modules <- 0L
  }

  if (is.null(module_presence)) {
    module_presence <- matrix(TRUE, nrow = D, ncol = n_modules)
  } else {
    module_presence <- as.matrix(module_presence)
    if (!all(dim(module_presence) == c(D, n_modules))) {
      stop("module_presence must be a D x n_modules logical matrix.", call. = FALSE)
    }
    module_presence <- module_presence > 0
  }

  if (n_dataset_specific_modules > 0) {
    if (is.null(dataset_specific_module_presence)) {
      dataset_specific_module_presence <- matrix(
        TRUE, nrow = D, ncol = n_dataset_specific_modules
      )
    } else {
      dataset_specific_module_presence <- as.matrix(dataset_specific_module_presence)
      if (!all(dim(dataset_specific_module_presence) == c(D, n_dataset_specific_modules))) {
        stop(
          "dataset_specific_module_presence must be D x n_dataset_specific_modules.",
          call. = FALSE
        )
      }
      dataset_specific_module_presence <- dataset_specific_module_presence > 0
    }
  } else {
    dataset_specific_module_presence <- matrix(FALSE, nrow = D, ncol = 0)
  }

  meta <- .make_shared_isvdx_meta(
    G = G,
    ng = ng,
    time_imbalance = time_imbalance,
    severity_group_shift = severity_group_shift
  )

  n <- nrow(meta)

  process_types <- c(
    "group_mean",
    "shared_oscillation",
    "group_phase_oscillation",
    "severity_modulated_oscillation"
  )

  module_type <- sample(rep(process_types, length.out = n_modules))
  dataset_specific_module_type <- if (n_dataset_specific_modules > 0) {
    sample(rep(process_types, length.out = n_dataset_specific_modules))
  } else {
    character(0)
  }

  # ----------------------------------------------------------
  # Shared sample-level processes: generated once and reused in
  # every dataset. Dataset-specific features/loadings can differ.
  # ----------------------------------------------------------
  shared_module_signal <- matrix(0, n_modules, n)
  shared_phase_by_group <- matrix(0, n_modules, G)

  for (m in seq_len(n_modules)) {
    sig <- .simulate_module_signal(
      module_type = module_type[m],
      meta = meta,
      G = G,
      module_sd = shared_process_sd,
      phase_shift_sd = phase_shift_sd
    )
    shared_module_signal[m, ] <- sig$sample_signal
    shared_phase_by_group[m, ] <- sig$phase_by_group
  }

  Y <- vector("list", D)
  zdata <- vector("list", D)
  zformula <- vector("list", D)

  truth <- list(
    module_type = module_type,
    module_presence = module_presence,
    dataset_specific_module_type = dataset_specific_module_type,
    dataset_specific_module_presence = dataset_specific_module_presence,

    # Backward-compatible aliases used by current integrated evaluation.
    # These intentionally point to the shared target.
    module_id = vector("list", D),
    module_features = vector("list", D),
    signal_features = vector("list", D),
    null_features = vector("list", D),
    module_signal = vector("list", D),
    phase_by_group = vector("list", D),
    biology = vector("list", D),

    # Explicit shared target.
    shared_module_id = vector("list", D),
    shared_module_features = vector("list", D),
    shared_signal_features = vector("list", D),
    shared_null_features = vector("list", D),
    shared_module_signal = shared_module_signal,
    shared_phase_by_group = shared_phase_by_group,
    shared_biology = vector("list", D),

    # Optional dataset-specific nuisance biology.
    dataset_specific_module_id = vector("list", D),
    dataset_specific_module_features = vector("list", D),
    dataset_specific_signal_features = vector("list", D),
    dataset_specific_module_signal = vector("list", D),
    dataset_specific_phase_by_group = vector("list", D),
    dataset_specific_biology = vector("list", D),

    total_biology = vector("list", D),
    batch_effect = vector("list", D),
    batch_features = vector("list", D),
    non_batch_features = vector("list", D),
    batch_biology_overlap_features = vector("list", D),
    batch_only_features = vector("list", D),
    biology_only_features = vector("list", D),
    clean_null_features = vector("list", D),
    feature_noise_var = vector("list", D),
    feature_noise_sd = vector("list", D)
  )

  for (d in seq_len(D)) {
    pd <- as.integer(p[d])
    Bd <- as.integer(B[d])
    batch_sd <- batch_to_biology_ratio[d] * module_sd[d]

    batch <- .simulate_dataset_batches(
      meta = meta,
      B = Bd,
      confounding = confounding[d]
    )

    zdata[[d]] <- data.frame(
      sample = meta$sample,
      batch = factor(batch),
      row.names = meta$sample
    )
    zformula[[d]] <- ~ batch

    feature_noise_var <- sigma[d]^2 * stats::rgamma(
      pd,
      shape = feature_noise_shape[d],
      rate = feature_noise_shape[d]
    )
    feature_noise_sd <- sqrt(feature_noise_var)
    noise <- matrix(
      stats::rnorm(pd * n, mean = 0, sd = rep(feature_noise_sd, times = n)),
      pd, n
    )

    available <- seq_len(pd)

    n_present_shared <- sum(module_presence[d, ])
    shared_module_size <- if (n_present_shared > 0) {
      floor(pd * pi_signal[d] / n_present_shared)
    } else {
      0L
    }

    shared_alloc <- .sample_module_features(
      available = available,
      n_modules = n_modules,
      module_size = shared_module_size,
      module_presence = module_presence[d, ]
    )

    shared_module_features <- shared_alloc$module_features
    shared_module_id <- shared_alloc$module_id
    available <- shared_alloc$available

    shared_biology <- matrix(0, pd, n)

    for (m in seq_len(n_modules)) {
      fm <- shared_module_features[[m]]
      if (length(fm) == 0) next

      loading <- stats::rnorm(length(fm), mean = 1, sd = 0.15)
      loading <- loading * sample(c(-1, 1), length(fm), replace = TRUE)

      shared_biology[fm, ] <- loading %*% t(shared_module_signal[m, ])
    }

    n_present_ds <- if (n_dataset_specific_modules > 0) {
      sum(dataset_specific_module_presence[d, ])
    } else {
      0L
    }

    ds_module_size <- if (n_present_ds > 0) {
      floor(pd * pi_dataset_specific_signal[d] / n_present_ds)
    } else {
      0L
    }

    ds_alloc <- .sample_module_features(
      available = available,
      n_modules = n_dataset_specific_modules,
      module_size = ds_module_size,
      module_presence = if (n_dataset_specific_modules > 0) {
        dataset_specific_module_presence[d, ]
      } else {
        logical(0)
      }
    )

    ds_module_features <- ds_alloc$module_features
    ds_module_id <- ds_alloc$module_id

    dataset_specific_biology <- matrix(0, pd, n)
    ds_module_signal <- matrix(0, n_dataset_specific_modules, n)
    ds_phase_by_group <- matrix(0, n_dataset_specific_modules, G)

    if (n_dataset_specific_modules > 0) {
      for (m in seq_len(n_dataset_specific_modules)) {
        fm <- ds_module_features[[m]]
        if (length(fm) == 0) next

        sig <- .simulate_module_signal(
          module_type = dataset_specific_module_type[m],
          meta = meta,
          G = G,
          module_sd = dataset_specific_module_sd[d],
          phase_shift_sd = phase_shift_sd
        )

        sample_signal <- dataset_specific_strength * sig$sample_signal
        loading <- stats::rnorm(length(fm), mean = 1, sd = 0.15)
        loading <- loading * sample(c(-1, 1), length(fm), replace = TRUE)

        ds_module_signal[m, ] <- sample_signal
        ds_phase_by_group[m, ] <- sig$phase_by_group
        dataset_specific_biology[fm, ] <- loading %*% t(sample_signal)
      }
    }

    shared_biology <- c_bio[d] * shared_biology
    dataset_specific_biology <- c_bio[d] * dataset_specific_biology
    total_biology <- shared_biology + dataset_specific_biology

    batch_effect <- matrix(0, pd, n)

    shared_signal_features_pre <- which(!is.na(shared_module_id))
    shared_null_features_pre <- setdiff(seq_len(pd), shared_signal_features_pre)

    n_batch_features <- as.integer(round(batch_feature_fraction[d] * pd))
    overlap_by_module <- vapply(
      shared_module_features,
      function(fm) as.integer(round(batch_biology_overlap_fraction[d] * length(fm))),
      integer(1)
    )
    n_overlap <- sum(overlap_by_module)

    if (n_overlap > n_batch_features) {
      stop(
        "For dataset ", d, ", requested batch_biology_overlap_fraction cannot be achieved ",
        "because the implied number of batch-contaminated biological features exceeds ",
        "round(batch_feature_fraction * p).",
        call. = FALSE
      )
    }

    batch_biology_overlap_features <- unlist(
      lapply(seq_along(shared_module_features), function(m) {
        fm <- shared_module_features[[m]]
        nm <- overlap_by_module[[m]]
        if (length(fm) == 0L || nm == 0L) return(integer(0))
        sample(fm, size = nm, replace = FALSE)
      }),
      use.names = FALSE
    )
    batch_biology_overlap_features <- sort(unique(as.integer(batch_biology_overlap_features)))

    n_batch_only <- n_batch_features - length(batch_biology_overlap_features)
    if (n_batch_only > length(shared_null_features_pre)) {
      stop(
        "For dataset ", d, ", requested batch_feature_fraction cannot be achieved after ",
        "controlled batch-biology overlap because there are not enough non-shared-signal ",
        "features for batch-only effects.",
        call. = FALSE
      )
    }

    batch_only_features <- if (n_batch_only > 0L) {
      sample(shared_null_features_pre, size = n_batch_only, replace = FALSE)
    } else {
      integer(0)
    }

    batch_features <- sort(unique(c(batch_biology_overlap_features, batch_only_features)))
    non_batch_features <- setdiff(seq_len(pd), batch_features)
    biology_only_features <- setdiff(shared_signal_features_pre, batch_biology_overlap_features)
    clean_null_features <- setdiff(shared_null_features_pre, batch_only_features)

    if (length(batch_features) > 0L) {
      batch_loading <- matrix(
        stats::rnorm(length(batch_features) * Bd, 0, batch_sd),
        nrow = length(batch_features),
        ncol = Bd
      )
      for (b in seq_len(Bd)) {
        idx <- which(batch == b)
        if (length(idx) > 0L) {
          batch_effect[batch_features, idx] <- batch_loading[, b, drop = FALSE] %*%
            t(rep(1, length(idx)))
        }
      }
    }

    batch_effect <- c_batch[d] * batch_effect

    Yd <- noise + total_biology + batch_effect

    fam_d <- family[[d]]
    if (!is.null(fam_d) && fam_d %in% c("poisson", "quasipoisson")) {
      mu <- exp(pmin(Yd, 8))
      Yd <- matrix(stats::rpois(length(mu), lambda = as.numeric(mu)), pd, n)
    }

    rownames(Yd) <- paste0("dataset_", d, "_feature_", seq_len(pd))
    colnames(Yd) <- meta$sample

    Y[[d]] <- Yd

    shared_signal_features <- which(!is.na(shared_module_id))
    # For integrated evaluation, null means "not shared signal". Dataset-specific
    # nuisance features are therefore non-target features.
    shared_null_features <- setdiff(seq_len(pd), shared_signal_features)
    ds_signal_features <- if (length(ds_module_id) > 0) {
      which(!is.na(ds_module_id))
    } else {
      integer(0)
    }

    # Backward-compatible shared-target aliases.
    truth$module_id[[d]] <- shared_module_id
    truth$module_features[[d]] <- shared_module_features
    truth$signal_features[[d]] <- shared_signal_features
    truth$null_features[[d]] <- shared_null_features
    truth$module_signal[[d]] <- shared_module_signal
    truth$phase_by_group[[d]] <- shared_phase_by_group
    truth$biology[[d]] <- shared_biology

    truth$shared_module_id[[d]] <- shared_module_id
    truth$shared_module_features[[d]] <- shared_module_features
    truth$shared_signal_features[[d]] <- shared_signal_features
    truth$shared_null_features[[d]] <- shared_null_features
    truth$shared_biology[[d]] <- shared_biology

    truth$dataset_specific_module_id[[d]] <- ds_module_id
    truth$dataset_specific_module_features[[d]] <- ds_module_features
    truth$dataset_specific_signal_features[[d]] <- ds_signal_features
    truth$dataset_specific_module_signal[[d]] <- ds_module_signal
    truth$dataset_specific_phase_by_group[[d]] <- ds_phase_by_group
    truth$dataset_specific_biology[[d]] <- dataset_specific_biology

    truth$total_biology[[d]] <- total_biology
    truth$batch_effect[[d]] <- batch_effect
    truth$batch_features[[d]] <- batch_features
    truth$non_batch_features[[d]] <- non_batch_features
    truth$batch_biology_overlap_features[[d]] <- batch_biology_overlap_features
    truth$batch_only_features[[d]] <- batch_only_features
    truth$biology_only_features[[d]] <- biology_only_features
    truth$clean_null_features[[d]] <- clean_null_features
    truth$feature_noise_var[[d]] <- feature_noise_var
    truth$feature_noise_sd[[d]] <- feature_noise_sd
  }

  names(Y) <- paste0("dataset_", seq_len(D))
  names(zdata) <- names(Y)
  names(zformula) <- names(Y)

  list_fields <- c(
    "module_id",
    "module_features",
    "signal_features",
    "null_features",
    "module_signal",
    "phase_by_group",
    "biology",
    "shared_module_id",
    "shared_module_features",
    "shared_signal_features",
    "shared_null_features",
    "shared_biology",
    "dataset_specific_module_id",
    "dataset_specific_module_features",
    "dataset_specific_signal_features",
    "dataset_specific_module_signal",
    "dataset_specific_phase_by_group",
    "dataset_specific_biology",
    "total_biology",
    "batch_effect",
    "batch_features",
    "non_batch_features",
    "batch_biology_overlap_features",
    "batch_only_features",
    "biology_only_features",
    "clean_null_features",
    "feature_noise_var",
    "feature_noise_sd"
  )

  for (nm in list_fields) {
    names(truth[[nm]]) <- names(Y)
  }

  list(
    Y = Y,
    xdata = meta,
    data = meta,
    xformula = ~ group * t + severity * t,
    formula = ~ group * t + severity * t,
    zdata = zdata,
    zformula = zformula,
    family = family,
    integrated = TRUE,
    truth = truth,
    params = list(
      D = D,
      p = p,
      n = n,
      G = G,
      ng = ng,
      B = B,
      n_modules = n_modules,
      pi_signal = pi_signal,
      module_presence = module_presence,
      include_dataset_specific_processes = include_dataset_specific_processes,
      n_dataset_specific_modules = n_dataset_specific_modules,
      pi_dataset_specific_signal = pi_dataset_specific_signal,
      dataset_specific_module_presence = dataset_specific_module_presence,
      dataset_specific_strength = dataset_specific_strength,
      confounding = confounding,
      time_imbalance = time_imbalance,
      severity_group_shift = severity_group_shift,
      sigma = sigma,
      batch_to_biology_ratio = batch_to_biology_ratio,
      c_bio = c_bio,
      c_batch = c_batch,
      feature_noise_shape = feature_noise_shape,
      batch_biology_overlap_fraction = batch_biology_overlap_fraction,
      module_sd = module_sd,
      shared_process_sd = shared_process_sd,
      dataset_specific_module_sd = dataset_specific_module_sd,
      batch_feature_fraction = batch_feature_fraction,
      family = family
    )
  )
}

.module_feature_sets_from_isvdx_sim <- function(sim, dataset) {
  d <- if (is.numeric(dataset)) dataset else match(dataset, names(sim$Y))
  if (is.na(d)) stop("Unknown dataset.", call. = FALSE)

  feature_names <- rownames(sim$Y[[d]])

  # Current integrated benchmarks evaluate shared modules only.
  sets <- sim$truth$shared_module_features[[d]]
  names(sets) <- paste0("module_", seq_along(sets))

  lapply(sets, function(idx) feature_names[idx])
}

make_isvdx_train_test_split <- function(xdata, test_frac = 0.30, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(xdata)
  test_idx <- sample(seq_len(n), size = round(test_frac * n), replace = FALSE)

  list(
    train = setdiff(seq_len(n), test_idx),
    test = test_idx
  )
}

subset_isvdx_inputs <- function(
    Y,
    xdata,
    zdata,
    graphData = NULL,
    si
) {
  list(
    Y = lapply(Y, function(Yd) Yd[, si, drop = FALSE]),
    xdata = xdata[si, , drop = FALSE],
    zdata = lapply(zdata, function(zd) zd[si, , drop = FALSE]),
    graphData = if (is.null(graphData)) NULL else graphData[, si, drop = FALSE]
  )
}
