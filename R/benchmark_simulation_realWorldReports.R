# ============================================================
# Simulation-as-real-world adapters and illustrative reports
# ============================================================
# These helpers convert a single simulated dataset into the same prepared-object
# structure consumed by run_real_world_single_benchmark().  This is intended for
# illustrative reports, not for quantitative grid benchmarking.
# ============================================================

`%||%` <- function(x, y) if (!is.null(x)) x else y

.benchmark_sim_find_function <- function(name) {
  if (exists(name, mode = "function", inherits = TRUE)) {
    return(get(name, mode = "function", inherits = TRUE))
  }
  pkg_opt <- getOption("svdxBenchmark.svdx_package", NULL)
  candidates <- unique(stats::na.omit(c(pkg_opt, "mySvdx", "svdx", "mSvdx", "SVDX")))
  for (pkg in candidates) {
    if (requireNamespace(pkg, quietly = TRUE) &&
        exists(name, envir = asNamespace(pkg), inherits = FALSE)) {
      return(getFromNamespace(name, pkg))
    }
  }
  NULL
}

.benchmark_sim_feature_ids <- function(sim, feature_prefix = "feature_") {
  if (!is.null(rownames(sim$Y))) return(rownames(sim$Y))
  paste0(feature_prefix, seq_len(nrow(sim$Y)))
}

.benchmark_sim_sample_ids <- function(sim, sample_prefix = "sample_") {
  if (!is.null(colnames(sim$Y))) return(colnames(sim$Y))
  paste0(sample_prefix, seq_len(ncol(sim$Y)))
}

#' Build target feature sets from a single SVDX interaction simulation.
#'
#' @param sim Object returned by simulate_svdx_interactions().
#' @param include_modules Include one set per simulated module.
#' @param include_module_types Include one set per module type.
#' @param feature_prefix Prefix used if sim$Y has no row names.
#' @return Named list of feature IDs/rownames.
simulation_module_target_sets <- function(
    sim,
    include_modules = TRUE,
    include_module_types = TRUE,
    feature_prefix = "feature_"
) {
  if (is.null(sim$Y) || is.null(sim$truth)) {
    stop("sim must be an object returned by simulate_svdx_interactions().", call. = FALSE)
  }
  feature_ids <- .benchmark_sim_feature_ids(sim, feature_prefix = feature_prefix)
  truth <- sim$truth
  if (is.null(truth$module_features) || is.null(truth$module_type)) {
    stop("sim$truth must contain module_features and module_type.", call. = FALSE)
  }

  out <- list()

  if (isTRUE(include_modules)) {
    for (m in seq_along(truth$module_features)) {
      idx <- truth$module_features[[m]]
      idx <- idx[!is.na(idx) & idx >= 1L & idx <= length(feature_ids)]
      nm <- paste0("module_", m)
      if (!is.null(truth$module_type[[m]])) {
        nm <- paste0(nm, "__", truth$module_type[[m]])
      }
      out[[nm]] <- feature_ids[idx]
    }
  }

  if (isTRUE(include_module_types)) {
    types <- unique(as.character(truth$module_type))
    types <- types[!is.na(types) & nzchar(types)]
    for (tp in types) {
      mods <- which(as.character(truth$module_type) == tp)
      idx <- unique(unlist(truth$module_features[mods], use.names = FALSE))
      idx <- idx[!is.na(idx) & idx >= 1L & idx <= length(feature_ids)]
      out[[paste0("module_type__", tp)]] <- feature_ids[idx]
    }
  }

  out[lengths(out) > 0L]
}

#' Create fdata for a single simulated dataset.
#'
#' The returned object is passed through as_svdx_fdata() when that function is
#' available. If the package helper is unavailable or rejects the constructed
#' data frame, the plain data frame is returned and the caller may set
#' create_fdata = FALSE.
#'
#' @param sim Object returned by simulate_svdx_interactions().
#' @param feature_prefix Prefix used if sim$Y has no row names.
#' @return fdata object or data frame.
prepare_simulation_fdata <- function(
    sim,
    feature_prefix = "feature_"
) {
  if (is.null(sim$Y) || is.null(sim$truth)) {
    stop("sim must be an object returned by simulate_svdx_interactions().", call. = FALSE)
  }
  feature_ids <- .benchmark_sim_feature_ids(sim, feature_prefix = feature_prefix)
  p <- length(feature_ids)
  truth <- sim$truth

  module_id <- truth$module_id %||% rep(NA_integer_, p)
  if (length(module_id) != p) {
    stop("sim$truth$module_id must have length nrow(sim$Y).", call. = FALSE)
  }

  module_type <- rep(NA_character_, p)
  valid <- !is.na(module_id) & module_id >= 1L & module_id <= length(truth$module_type)
  module_type[valid] <- as.character(truth$module_type[module_id[valid]])

  fdata <- data.frame(
    feature_id = feature_ids,
    feature_name = feature_ids,
    enrichment_id = feature_ids,
    module_id = module_id,
    module_type = module_type,
    is_signal_feature = seq_len(p) %in% (truth$signal_features %||% integer()),
    is_null_feature = seq_len(p) %in% (truth$null_features %||% integer()),
    is_batch_feature = seq_len(p) %in% unique(unlist(truth$batch_features %||% list(), use.names = FALSE)),
    stringsAsFactors = FALSE
  )
  rownames(fdata) <- feature_ids

  as_fdata <- .benchmark_sim_find_function("as_svdx_fdata")
  if (!is.null(as_fdata)) {
    # The package helper requires the column containing feature IDs to be
    # specified explicitly.  Simulation fdata always uses feature_id as the
    # canonical row-feature identifier and feature_name as the display label.
    fmls <- tryCatch(names(formals(as_fdata)), error = function(e) character())
    args <- list(fdata = fdata)
    if ("feature_col" %in% fmls) args$feature_col <- "feature_id"
    if ("feature_name_col" %in% fmls) args$feature_name_col <- "feature_name"
    if ("label_col" %in% fmls) args$label_col <- "feature_name"
    if ("name_col" %in% fmls) args$name_col <- "feature_name"
    if ("enrichment_id_col" %in% fmls) args$enrichment_id_col <- "enrichment_id"
    if ("enrichment_col" %in% fmls) args$enrichment_col <- "enrichment_id"
    if ("enrichment_feature_col" %in% fmls) args$enrichment_feature_col <- "enrichment_id"
    if ("id_col" %in% fmls) args$id_col <- "enrichment_id"

    converted <- tryCatch(do.call(as_fdata, args), error = function(e) e)
    if (!inherits(converted, "error")) return(converted)

    # Backward-compatible fallback for helper variants that expect the data
    # frame as the first unnamed argument rather than an argument named fdata.
    args2 <- args
    args2$fdata <- NULL
    args2 <- c(list(fdata), args2)
    converted <- tryCatch(do.call(as_fdata, args2), error = function(e) e)
    if (!inherits(converted, "error")) return(converted)

    warning(
      "as_svdx_fdata() was found but could not convert simulation fdata: ",
      conditionMessage(converted),
      call. = FALSE
    )
  }

  fdata
}

#' Convert a single simulated dataset to a prepared real-world benchmark object.
#'
#' This adapter lets a chosen simulation be treated like a real-world dataset for
#' illustrative reports. It does not run benchmark-grid truth metrics.
#'
#' @param sim Object returned by simulate_svdx_interactions().
#' @param xformula Formula for biological/design covariates.
#' @param zformula Formula for nuisance covariates.
#' @param xdata Optional sample metadata; defaults to sim$meta.
#' @param zdata Optional nuisance metadata; defaults to xdata.
#' @param family Optional family passed through to real-world benchmark runner.
#' @param modelType Optional modelType passed through to real-world benchmark runner.
#' @param create_fdata Whether to construct simulation fdata.
#' @param create_target_sets Whether to construct simulation module target sets.
#' @return Prepared-object list consumable by run_real_world_single_benchmark().
prepare_simulation_as_real_world <- function(
    sim,
    xformula = ~ group * t + severity * t,
    zformula = ~ batch,
    xdata = NULL,
    zdata = NULL,
    family = NULL,
    modelType = NULL,
    create_fdata = TRUE,
    create_target_sets = TRUE,
    feature_prefix = "feature_",
    sample_prefix = "sample_"
) {
  if (is.null(sim$Y) || is.null(sim$meta) || is.null(sim$truth)) {
    stop("sim must contain Y, meta, and truth.", call. = FALSE)
  }

  Y <- as.matrix(sim$Y)
  feature_ids <- .benchmark_sim_feature_ids(sim, feature_prefix = feature_prefix)
  sample_ids <- .benchmark_sim_sample_ids(sim, sample_prefix = sample_prefix)
  rownames(Y) <- feature_ids
  colnames(Y) <- sample_ids

  xdata <- as.data.frame(xdata %||% sim$meta)
  if (nrow(xdata) != ncol(Y)) {
    stop("xdata must have one row per sample/column of sim$Y.", call. = FALSE)
  }
  rownames(xdata) <- sample_ids

  zdata <- zdata %||% xdata
  if (!is.null(zdata)) {
    zdata <- as.data.frame(zdata)
    if (nrow(zdata) != ncol(Y)) {
      stop("zdata must have one row per sample/column of sim$Y.", call. = FALSE)
    }
    rownames(zdata) <- sample_ids
  }

  fdata <- if (isTRUE(create_fdata)) {
    prepare_simulation_fdata(sim, feature_prefix = feature_prefix)
  } else {
    NULL
  }

  target_sets <- if (isTRUE(create_target_sets)) {
    simulation_module_target_sets(sim, feature_prefix = feature_prefix)
  } else {
    NULL
  }

  # Package enrichment ultimately requires a feature -> enrichment identifier map.
  # For simulation illustrations, the enrichment identifier system is the
  # simulation feature ID itself, so module target_sets and enrichment IDs are
  # guaranteed to be in the same namespace.  This mirrors the package report
  # behavior where the map is either derived from fdata or supplied explicitly.
  enrichment_id_map <- data.frame(
    feature = feature_ids,
    enrichment_id = feature_ids,
    stringsAsFactors = FALSE
  )

  out <- list(
    Y = Y,
    xdata = xdata,
    zdata = zdata,
    xformula = xformula,
    zformula = zformula,
    family = family,
    modelType = modelType,
    fdata = fdata,
    enrichment_id_map = enrichment_id_map,
    target_sets = target_sets,
    feature_sets = target_sets,
    truth = sim$truth,
    params = sim$params,
    source = "simulation_as_real_world"
  )
  class(out) <- c("preparedSimulationRealWorldBenchmark", "preparedSVDXBenchmark", class(out))
  out
}

#' Run an illustrative real-world-style report for one simulated dataset.
#'
#' This is a convenience wrapper around prepare_simulation_as_real_world(),
#' run_real_world_single_benchmark(), and runSingleBenchmarkReports().  It is
#' intended for one chosen simulation used
#' in a paper or presentation, not for benchmark grids.
#'
#' @return Invisible list containing prepared, rw_single, and report.
runSimulationAsRealWorldReports <- function(
    sim,
    methods = c("svdx", "pca", "nmf", "mixomics_pls"),
    inference_methods = c("svdx"),
    K = 6,
    xformula = ~ group * t + severity * t,
    zformula = ~ batch,
    family = NULL,
    modelType = NULL,
    create_fdata = TRUE,
    create_target_sets = TRUE,
    target_sets = NULL,
    evaluate_target_sets = TRUE,
    combine_sets = "union",
    a_grid = c(0, 0.5, 1, 2, 4),
    outdir = "benchmark_reports/simulation_as_real_world",
    figure_format = "pdf",
    overwrite = TRUE,
    verbose = TRUE,
    ...
) {
  prepared <- prepare_simulation_as_real_world(
    sim,
    xformula = xformula,
    zformula = zformula,
    family = family,
    modelType = modelType,
    create_fdata = create_fdata,
    create_target_sets = create_target_sets
  )

  target_sets <- target_sets %||% prepared$target_sets

  rw_single <- run_real_world_single_benchmark(
    prepared,
    K = K,
    methods = methods,
    inference_methods = inference_methods,
    feature_sets = target_sets,
    enrichment_id_map = prepared$enrichment_id_map,
    target_sets = if (isTRUE(evaluate_target_sets)) target_sets else NULL,
    target_set_auc_args = list(
      combine_sets = combine_sets,
      a_grid = a_grid
    ),
    ...
  )

  report <- runSingleBenchmarkReports(
    rw_single,
    outdir = outdir,
    figure_format = figure_format,
    overwrite = overwrite,
    verbose = verbose
  )

  invisible(list(
    prepared = prepared,
    rw_single = rw_single,
    report = report
  ))
}
