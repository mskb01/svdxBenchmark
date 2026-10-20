# ============================================================
# Real-world integrated benchmarking wrappers
# ============================================================

# ------------------------------------------------------------
# Internal validation helpers for integrated real-world wrappers
# ------------------------------------------------------------


# ------------------------------------------------------------
# Feature label helpers for real-world plotting
# ------------------------------------------------------------
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

if (!exists(".benchmark_family_has_count", mode = "function")) {
  .benchmark_family_has_count <- function(family) {
    if (is.null(family)) return(FALSE)
    any(unlist(family, use.names = FALSE) %in% c("poisson", "quasipoisson"), na.rm = TRUE)
  }
}


if (!exists(".safe_fisher_enrichment", mode = "function")) {
  .safe_fisher_enrichment <- function(discovered, in_set) {
    if (exists(".bench_safe_fisher_enrichment", mode = "function", inherits = TRUE)) {
      return(.bench_safe_fisher_enrichment(discovered, in_set))
    }
    discovered <- as.logical(discovered)
    in_set <- as.logical(in_set)
    ok <- !is.na(discovered) & !is.na(in_set)
    discovered <- discovered[ok]
    in_set <- in_set[ok]
    if (length(discovered) == 0L || length(in_set) == 0L) {
      return(list(odds_ratio = NA_real_, p_value = NA_real_))
    }
    tab <- table(factor(discovered, levels = c(FALSE, TRUE)),
                 factor(in_set, levels = c(FALSE, TRUE)))
    ft <- tryCatch(stats::fisher.test(tab), error = function(e) NULL)
    if (is.null(ft)) return(list(odds_ratio = NA_real_, p_value = NA_real_))
    list(odds_ratio = unname(ft$estimate %||% NA_real_), p_value = unname(ft$p.value %||% NA_real_))
  }
}

if (!exists(".benchmark_get_svdx_internal", mode = "function")) {
  .benchmark_get_svdx_internal <- function(
      name,
      packages = unique(c(
        getOption("svdxBenchmark.svdx_package", NULL),
        "mySvdx", "svdx", "mSvdx", "SVDX"
      ))) {
    if (exists(name, envir = .GlobalEnv, inherits = TRUE)) {
      return(get(name, envir = .GlobalEnv, inherits = TRUE))
    }
    packages <- packages[!is.na(packages) & nzchar(packages)]
    for (pkg in packages) {
      if (requireNamespace(pkg, quietly = TRUE) &&
          exists(name, envir = asNamespace(pkg), inherits = FALSE)) {
        return(getFromNamespace(name, pkg))
      }
    }
    stop("Could not find SVDX package helper `", name, "`. Please load or install the SVDX package.", call. = FALSE)
  }
}

if (!exists(".benchmark_resolve_fdata_for_dataset", mode = "function")) {
  .benchmark_resolve_fdata_for_dataset <- function(fdata, dataset = NULL, required = FALSE) {
    if (is.null(fdata)) {
      if (isTRUE(required)) stop("fdata is required.", call. = FALSE)
      return(NULL)
    }
    fn <- .benchmark_get_svdx_internal(".svdx_resolve_fdata_for_dataset")
    fn(fdata, dataset = dataset, required = required)
  }
}

if (!exists(".benchmark_fdata_to_row_labels", mode = "function")) {
  .benchmark_fdata_to_row_labels <- function(fdata, dataset = NULL) {
    if (is.null(fdata)) return(NULL)
    fdata_use <- .benchmark_resolve_fdata_for_dataset(fdata, dataset = dataset, required = FALSE)
    if (is.null(fdata_use)) return(NULL)
    fn <- .benchmark_get_svdx_internal(".svdx_fdata_to_row_labels")
    fn(fdata_use)
  }
}

if (!exists(".benchmark_normalize_row_labels", mode = "function")) {
  .benchmark_normalize_row_labels <- function(row_labels = NULL, fdata = NULL, dataset = NULL) {
    if (!is.null(row_labels)) return(row_labels)
    .benchmark_fdata_to_row_labels(fdata, dataset = dataset)
  }
}

if (!exists(".benchmark_feature_label_map", mode = "function")) {
  .benchmark_feature_label_map <- function(row_labels = NULL, features = NULL) {
    if (is.null(row_labels)) return(NULL)
    if (is.vector(row_labels) && !is.list(row_labels)) {
      if (is.null(names(row_labels))) {
        if (is.null(features) || length(row_labels) != length(features)) {
          stop("Unnamed row_labels must have the same length as features.", call. = FALSE)
        }
        return(data.frame(feature = as.character(features), label = as.character(row_labels), stringsAsFactors = FALSE))
      }
      return(data.frame(feature = as.character(names(row_labels)), label = as.character(row_labels), stringsAsFactors = FALSE))
    }
    if (is.data.frame(row_labels)) {
      if (ncol(row_labels) < 2L) {
        stop("row_labels as a data frame must have at least two columns: feature and display label.", call. = FALSE)
      }
      out <- row_labels[, 1:2, drop = FALSE]
      names(out) <- c("feature", "label")
      out$feature <- as.character(out$feature)
      out$label <- as.character(out$label)
      return(out)
    }
    stop("row_labels must be NULL, a named vector, an unnamed vector matching features, or a two-column data frame.", call. = FALSE)
  }
}

if (!exists(".benchmark_resolve_feature_ids_and_labels", mode = "function")) {
  .benchmark_resolve_feature_ids_and_labels <- function(features, feature_names, row_labels = NULL) {
    feature_names <- as.character(feature_names)
    if (is.null(features)) features <- feature_names
    requested <- as.character(features)
    label_map <- .benchmark_feature_label_map(row_labels, features = feature_names)
    resolved <- requested
    matched_feature <- requested %in% feature_names
    unresolved <- !matched_feature
    if (any(unresolved) && !is.null(label_map)) {
      dup_labels <- unique(label_map$label[duplicated(label_map$label)])
      requested_labels <- intersect(requested[unresolved], label_map$label)
      ambiguous <- intersect(requested_labels, dup_labels)
      if (length(ambiguous)) {
        stop("Feature label(s) map to multiple features: ", paste(ambiguous, collapse = ", "), ". Please specify original feature IDs or unique labels.", call. = FALSE)
      }
      hits <- match(requested[unresolved], label_map$label)
      unresolved_idx <- which(unresolved)
      has_hit <- !is.na(hits)
      resolved[unresolved_idx[has_hit]] <- label_map$feature[hits[has_hit]]
    }
    missing <- setdiff(resolved, feature_names)
    if (length(missing)) {
      stop("Some requested features/labels were not found: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    labels <- resolved
    if (!is.null(label_map)) {
      hits <- match(resolved, label_map$feature)
      use <- !is.na(hits)
      labels[use] <- label_map$label[hits[use]]
    }
    list(feature = resolved, label = labels)
  }
}

.validate_real_world_integrated_benchmark <- function(rw) {
  if (!inherits(rw, "realWorldIntegratedBenchmark")) {
    stop("rw must inherit from 'realWorldIntegratedBenchmark'.", call. = FALSE)
  }
  invisible(TRUE)
}

.validate_integrated_real_world_graph_inputs <- function(graphData, graphDistFn, n_samples) {
  if (!is.null(graphData)) {
    if (!is.data.frame(graphData)) {
      stop("graphData must be a data.frame when supplied.", call. = FALSE)
    }
    if (nrow(graphData) < 1L) {
      stop("graphData must contain at least one row.", call. = FALSE)
    }
    if (ncol(graphData) != n_samples) {
      stop("ncol(graphData) must equal ncol(Y[[1]]).", call. = FALSE)
    }
  }
  if (!is.null(graphDistFn) && !is.function(graphDistFn)) {
    stop("graphDistFn must be a function when supplied.", call. = FALSE)
  }
  invisible(TRUE)
}

.validate_integrated_feature_sets_basic <- function(feature_sets, Y = NULL) {
  if (is.null(feature_sets)) return(invisible(TRUE))
  if (length(feature_sets) == 0L) {
    stop("feature_sets must contain at least one dataset-specific feature set list.", call. = FALSE)
  }
  if (!is.null(Y)) {
    feature_names <- unique(unlist(lapply(Y, rownames), use.names = FALSE))
    fs_features <- unique(unlist(feature_sets, recursive = TRUE, use.names = FALSE))
    if (length(feature_names) > 0L && length(fs_features) > 0L &&
        length(intersect(fs_features, feature_names)) == 0L) {
      warning(
        "No overlap detected between supplied feature sets and benchmark feature names.",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

.validate_integrated_real_world_output <- function(x) {
  required <- c("Y", "results", "feature_scores", "loading_enrichment", "timing", "methods")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop(
      "Internal error: integrated real-world benchmark output is missing required fields: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}


# ------------------------------------------------------------
# Canonical integrated real-world plotting helpers
# ------------------------------------------------------------
# Defined unconditionally with rwbench-specific names to avoid stale helper
# definitions from older sourced files.

.rwbench_get_integrated_method_result <- function(rw, method) {
  .validate_real_world_integrated_benchmark(rw)
  if (is.null(rw$results) || !is.list(rw$results)) {
    stop("rw$results must be a named/list object.", call. = FALSE)
  }

  if (!is.null(names(rw$results)) && method %in% names(rw$results)) {
    res <- rw$results[[method]]
    if (is.null(res)) stop("Method not found in rw$results: ", method, call. = FALSE)
    return(res)
  }

  idx <- which(vapply(rw$results, function(z) identical(z$method, method), logical(1)))
  if (length(idx) == 0L) {
    available <- names(rw$results)
    available <- available[!is.na(available) & nzchar(available)]
    msg <- paste0("Method not found in rw$results: ", method)
    if (length(available)) msg <- paste0(msg, ". Available methods: ", paste(available, collapse = ", "))
    stop(msg, call. = FALSE)
  }
  if (length(idx) > 1L) stop("Multiple results found for method: ", method, call. = FALSE)
  rw$results[[idx[1]]]
}

.rwbench_fdata_row_labels <- function(feature_names, fdata = NULL, row_labels = NULL, dataset = NULL) {
  feature_names <- as.character(feature_names)
  if (!is.null(row_labels)) {
    labs <- row_labels
  } else if (!is.null(fdata)) {
    labs <- tryCatch({
      fdata_use <- fdata
      if (!is.null(dataset)) {
        resolver <- .benchmark_get_svdx_internal(".svdx_resolve_fdata_for_dataset")
        fdata_use <- resolver(fdata, dataset = dataset, required = FALSE)
      }
      if (is.null(fdata_use)) NULL else {
        labeller <- .benchmark_get_svdx_internal(".svdx_fdata_to_row_labels")
        labeller(fdata_use)
      }
    }, error = function(e) NULL)
  } else {
    labs <- NULL
  }
  if (is.null(labs)) return(data.frame(feature = feature_names, label = feature_names, stringsAsFactors = FALSE))
  if (is.vector(labs) && !is.list(labs)) {
    if (is.null(names(labs))) {
      if (length(labs) != length(feature_names)) stop("Unnamed row_labels/fdata labels must have the same length as feature names.", call. = FALSE)
      return(data.frame(feature = feature_names, label = as.character(labs), stringsAsFactors = FALSE))
    }
    return(data.frame(feature = as.character(names(labs)), label = as.character(labs), stringsAsFactors = FALSE))
  }
  if (is.data.frame(labs)) {
    if (ncol(labs) < 2L) stop("row_labels as a data frame must have at least two columns: feature and display label.", call. = FALSE)
    out <- labs[, 1:2, drop = FALSE]
    names(out) <- c("feature", "label")
    out$feature <- as.character(out$feature)
    out$label <- as.character(out$label)
    return(out)
  }
  stop("row_labels must be NULL, a vector, or a two-column data frame.", call. = FALSE)
}

.rwbench_resolve_features_and_labels <- function(features, feature_names, fdata = NULL, row_labels = NULL, dataset = NULL) {
  feature_names <- as.character(feature_names)
  if (is.null(features)) features <- feature_names
  requested <- as.character(features)
  label_map <- .rwbench_fdata_row_labels(feature_names = feature_names, fdata = fdata, row_labels = row_labels, dataset = dataset)
  resolved <- requested
  unresolved <- !(requested %in% feature_names)
  if (any(unresolved)) {
    dup_labels <- unique(label_map$label[duplicated(label_map$label)])
    requested_labels <- intersect(requested[unresolved], label_map$label)
    ambiguous <- intersect(requested_labels, dup_labels)
    if (length(ambiguous)) {
      stop("Feature label(s) map to multiple features: ", paste(ambiguous, collapse = ", "), ". Please specify original feature IDs or unique labels.", call. = FALSE)
    }
    hits <- match(requested[unresolved], label_map$label)
    unresolved_idx <- which(unresolved)
    has_hit <- !is.na(hits)
    resolved[unresolved_idx[has_hit]] <- label_map$feature[hits[has_hit]]
  }
  missing <- setdiff(resolved, feature_names)
  if (length(missing)) stop("Some requested features/labels were not found: ", paste(missing, collapse = ", "), call. = FALSE)
  labels <- resolved
  hits <- match(resolved, label_map$feature)
  use <- !is.na(hits)
  labels[use] <- label_map$label[hits[use]]
  list(feature = resolved, label = labels)
}


extract_integrated_feature_scores <- function(integrated_results) {
  dplyr::bind_rows(lapply(integrated_results, function(res) {
    dplyr::bind_rows(lapply(seq_along(res$V_list), function(d) {
      V <- res$V_list[[d]]
      if (is.null(V) || is.null(res$feature_score_list[[d]]) || length(res$feature_score_list[[d]]) == 0L) {
        return(NULL)
      }
      feature_names <- rownames(V) %||% paste0("feature_", seq_len(nrow(V)))
      dataset <- names(res$V_list)[d] %||% paste0("dataset_", d)

      tibble::tibble(
        method = res$method,
        dataset = dataset,
        feature = feature_names,
        feature_score = res$feature_score_list[[d]]
      )
    }))
  }))
}


enrich_integrated_loadings <- function(
    integrated_results,
    feature_sets,
    top_n = 100
) {
  dplyr::bind_rows(lapply(integrated_results, function(res) {
    dplyr::bind_rows(lapply(seq_along(res$V_list), function(d) {
      dataset <- names(res$V_list)[d] %||% paste0("dataset_", d)

      if (is.null(feature_sets[[dataset]])) return(NULL)

      V <- res$V_list[[d]]
      scores <- res$feature_score_list[[d]]
      if (is.null(V) || is.null(scores) || length(scores) == 0L) return(NULL)
      feature_names <- rownames(V) %||% paste0("feature_", seq_len(nrow(V)))
      names(scores) <- feature_names

      ranked <- sort(scores, decreasing = TRUE)
      selected <- names(ranked)[seq_len(min(top_n, length(ranked)))]
      universe <- feature_names

      out <- lapply(names(feature_sets[[dataset]]), function(set_name) {
        set_features <- intersect(feature_sets[[dataset]][[set_name]], universe)

        discovered <- universe %in% selected
        in_set <- universe %in% set_features

        sf <- .safe_fisher_enrichment(discovered, in_set)

        tibble::tibble(
          method = res$method,
          dataset = dataset,
          pathway = set_name,
          set_size = length(set_features),
          n_selected = length(selected),
          overlap = sum(discovered & in_set),
          odds_ratio = sf$odds_ratio,
          p_value = sf$p_value
        )
      })

      dplyr::bind_rows(out) |>
        dplyr::mutate(q_value = .safe_bh(p_value)) |>
        dplyr::arrange(q_value, p_value)
    }))
  }))
}


.make_real_world_integrated_train_test_split <- function(
    n_samples,
    train_fraction = 0.70,
    train_cols = NULL,
    test_cols = NULL,
    seed = NULL,
    sample_names = NULL
) {
  all_cols <- seq_len(n_samples)

  resolve_idx <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.character(x)) {
      if (is.null(sample_names)) {
        stop("Character train/test columns require colnames(Y[[1]]).", call. = FALSE)
      }
      idx <- match(x, sample_names)
      if (any(is.na(idx))) stop("Unknown sample name in train/test columns.", call. = FALSE)
      return(idx)
    }
    as.integer(x)
  }

  train <- resolve_idx(train_cols)
  test <- resolve_idx(test_cols)

  if (is.null(train) && is.null(test)) {
    if (!is.null(seed)) set.seed(seed)
    train <- sort(sample(all_cols, size = round(train_fraction * n_samples)))
    test <- setdiff(all_cols, train)
  } else if (is.null(train)) {
    test <- sort(unique(test))
    train <- setdiff(all_cols, test)
  } else if (is.null(test)) {
    train <- sort(unique(train))
    test <- setdiff(all_cols, train)
  } else {
    train <- sort(unique(train))
    test <- sort(unique(test))
  }

  if (length(train) == 0L) {
    stop("Training split must contain at least one sample.", call. = FALSE)
  }
  if (any(train < 1L | train > n_samples) || any(test < 1L | test > n_samples)) {
    stop("Train/test indices are out of bounds.", call. = FALSE)
  }
  if (length(intersect(train, test)) > 0L) {
    stop("Train and test splits must be disjoint.", call. = FALSE)
  }

  list(
    train = train,
    test = test,
    train_names = sample_names[train] %||% NULL,
    test_names = sample_names[test] %||% NULL
  )
}

.subset_integrated_real_world_inputs <- function(Y, xdata, zdata, idx) {
  list(
    Y = lapply(Y, function(y) y[, idx, drop = FALSE]),
    xdata = xdata[idx, , drop = FALSE],
    zdata = lapply(zdata, function(zd) zd[idx, , drop = FALSE])
  )
}



.benchmark_validate_prepared_integrated <- function(prepared) {
  if (!is.list(prepared) || is.null(prepared$Y) || is.null(prepared$xdata) || is.null(prepared$xformula)) {
    stop("Prepared integrated object must contain Y, xdata, and xformula.", call. = FALSE)
  }
  if (!isTRUE(prepared$integrated)) {
    stop("Prepared single-dataset objects should be passed to run_real_world_single_benchmark().", call. = FALSE)
  }
  if (!is.list(prepared$Y) || length(prepared$Y) == 0L) {
    stop("prepared$Y must be a non-empty named list of feature-by-sample matrices.", call. = FALSE)
  }
  xdata <- as.data.frame(prepared$xdata)
  for (nm in names(prepared$Y) %||% seq_along(prepared$Y)) {
    Yd <- as.matrix(prepared$Y[[nm]])
    if (ncol(Yd) != nrow(xdata)) {
      stop("All prepared$Y datasets must have ncol equal to nrow(prepared$xdata).", call. = FALSE)
    }
  }
  invisible(TRUE)
}


run_real_world_integrated_benchmark <- function(
    Y,
    xdata = NULL,
    xformula = NULL,
    zformula = NULL,
    zdata = NULL,
    K = 10,
    methods = c(
      "isvdx",
      "isvdx_adjusted",
      "separate_svdx",
      "MEFISTO",
      "concat_PCA",
      "MFA",
      "concat_NMF",
      "concat_sPCA",
      "concat_sNMF",
      "block_PLS",
      "block_sPLS"
    ),
    inference_methods = c("isvdx", "isvdx_adjusted", "separate_svdx"),
    run_factorization = TRUE,
    run_inference = TRUE,
    train_fraction = 0.70,
    train_cols = NULL,
    test_cols = NULL,
    seed = 1,
    alpha = 0.1,
    modelType = NULL,
    family = NULL,
    graphData = NULL,
    graphDistFn = NULL,
    standardize_features = TRUE,
    block_scale = c("none", "leading_sv"),
    adjust_for_z = TRUE,
    isvdx_adjusted_use_z = FALSE,
    feature_sets = NULL,
    target_sets = NULL,
    target_set_auc_args = list(),
    fdata = NULL,
    enrichment_id_map = NULL,
    loading_top_n = 100,
    loading_scope = c("all", "component", "union", "global"),
    loading_enrichment = NULL,
    inference_scope = c("all", "component", "union", "global"),
    inference_enrichment = TRUE,
    enrichment_backend = c("auto", "internal", "clusterProfiler", "fgsea"),
    msigdb_db_species = NULL,
    msigdb_species = NULL,
    msigdb_collection = NULL,
    msigdb_subcollection = NULL,
    enrichment_id_type = c("symbol", "entrez", "ensembl"),
    enrichment_score = "auto",
    enrichment_term2gene = NULL,
    enrichment_universe = NULL,
    min_feature_overlap_warn = 5L,
    min_feature_overlap_fraction = 0.01,
    enrichment_extra_args = list(),
    spca_keepX_grid = NULL,
    sparse_nmf_lambda_grid = 10^c(-3, -2, -1),
    block_spls_keepX_grid = c(0.01, 0.10, 0.25, 0.50, 0.75),
    mefisto_covariates = c("t", "severity"),
    mefisto_maxiter = 1000,
    numThreads = 4,
    return_design = TRUE,
    verbose = FALSE,
    z_adjustment = c("fitted", "reestimate"),
    ...
) {
  block_scale <- match.arg(block_scale)
  inference_scope <- match.arg(inference_scope)
  loading_scope <- .benchmark_normalize_loading_scope(loading_scope)
  z_adjustment <- match.arg(z_adjustment)
  do_inference_enrichment <- isTRUE(inference_enrichment)
  if (!is.list(target_set_auc_args)) {
    stop("target_set_auc_args must be a list.", call. = FALSE)
  }

  # Support prepared integrated workflow objects, mirroring
  # runISVDXWorkflow()/runISVDXReport() semantics.  A prepared integrated
  # object should contain Y/xdata/xformula/zdata/zformula and may contain
  # family/fdata/modelType.
  if (is.list(Y) && !is.null(Y$Y) && !is.null(Y$xdata)) {
    prepared <- Y
    .benchmark_validate_prepared_integrated(prepared)
    Y <- prepared$Y
    if (is.null(xdata)) xdata <- prepared$xdata
    if (is.null(xformula)) xformula <- prepared$xformula
    if (is.null(zdata)) zdata <- prepared$zdata
    if (is.null(zformula)) zformula <- prepared$zformula
    if (is.null(family)) family <- prepared$family
    if (is.null(fdata)) fdata <- prepared$fdata
    if (is.null(enrichment_id_map)) enrichment_id_map <- prepared$enrichment_id_map
    if (is.null(target_sets)) target_sets <- prepared$target_sets
    if (is.null(feature_sets)) feature_sets <- prepared$feature_sets %||% prepared$target_sets
    if (is.null(modelType)) modelType <- prepared$modelType
  }

  do_loading_enrichment <- if (is.null(loading_enrichment)) {
    !is.null(feature_sets) || !is.null(enrichment_term2gene) || !is.null(msigdb_collection)
  } else {
    isTRUE(loading_enrichment)
  }

  if (is.null(xdata)) stop("xdata must be supplied or available as prepared$xdata.", call. = FALSE)
  if (is.null(xformula)) stop("xformula must be supplied or available as prepared$xformula.", call. = FALSE)
  if (is.null(zdata)) stop("zdata must be supplied or available as prepared$zdata.", call. = FALSE)
  if (is.null(zformula)) stop("zformula must be supplied or available as prepared$zformula.", call. = FALSE)

  if (!is.list(Y) || length(Y) == 0L) {
    stop("Y must be a non-empty list of matrices.", call. = FALSE)
  }
  Y <- lapply(Y, as.matrix)
  dataset_names <- names(Y) %||% paste0("dataset_", seq_along(Y))
  names(Y) <- dataset_names
  n_samples <- vapply(Y, ncol, integer(1))
  if (length(unique(n_samples)) != 1L) {
    stop("All entries of Y must have the same number of columns/samples.", call. = FALSE)
  }
  xdata <- as.data.frame(xdata)
  if (nrow(xdata) != n_samples[[1]]) {
    stop("nrow(xdata) must equal ncol(Y[[1]]).", call. = FALSE)
  }
  if (length(zformula) != length(Y)) {
    stop("length(zformula) must equal length(Y).", call. = FALSE)
  }
  if (length(zdata) != length(Y)) {
    stop("length(zdata) must equal length(Y).", call. = FALSE)
  }
  zdata <- lapply(zdata, as.data.frame)
  for (d in seq_along(Y)) {
    if (nrow(zdata[[d]]) != ncol(Y[[d]])) {
      stop("Each zdata[[d]] must have nrow equal to ncol(Y[[d]]).", call. = FALSE)
    }
  }
  if (any(vapply(Y, function(y) is.null(colnames(y)), logical(1)))) {
    warning("One or more Y entries have no column names; downstream sample labels may be unavailable.", call. = FALSE)
  }
  if (any(vapply(Y, function(y) is.null(rownames(y)), logical(1)))) {
    warning("One or more Y entries have no row names; feature-level enrichment may be unavailable or less informative.", call. = FALSE)
  }
  valid_methods <- c("isvdx", "isvdx_adjusted", "separate_svdx", "MEFISTO", "concat_PCA", "MFA", "concat_NMF", "concat_sPCA", "concat_sNMF", "block_PLS", "block_sPLS")
  invalid_methods <- setdiff(methods, valid_methods)
  if (length(invalid_methods) > 0L) {
    stop(
      "Unknown integrated real-world methods: ", paste(invalid_methods, collapse = ", "),
      ". Valid options are: ", paste(valid_methods, collapse = ", "),
      call. = FALSE
    )
  }

  valid_inference_methods <- c("isvdx", "isvdx_adjusted", "separate_svdx")
  invalid_inference_methods <- setdiff(inference_methods, valid_inference_methods)
  if (length(invalid_inference_methods) > 0L) {
    stop(
      "Unknown inference_methods: ", paste(invalid_inference_methods, collapse = ", "),
      ". Valid options are: ", paste(valid_inference_methods, collapse = ", "),
      call. = FALSE
    )
  }

  family <- normalize_isvdx_family(family, length(Y), dataset_names = dataset_names)
  svdx_model_type <- .benchmark_complete_model_type(
    modelType = modelType,
    formula = xformula,
    data = xdata
  )

  .validate_integrated_real_world_graph_inputs(graphData, graphDistFn, n_samples[[1]])
  .validate_integrated_feature_sets_basic(feature_sets, Y)
  if (is.null(enrichment_id_map)) enrichment_id_map <- NULL
  if (!is.null(fdata)) {
    if (!is.list(fdata) || is.data.frame(fdata)) {
      stop("For integrated real-world benchmarking, fdata must be a named list of svdx_fdata objects, one per dataset.", call. = FALSE)
    }
    if (length(fdata) != length(Y)) {
      stop("For integrated real-world benchmarking, fdata must have one entry per dataset in Y.", call. = FALSE)
    }
    if (!is.null(names(fdata)) && all(nzchar(names(fdata)))) {
      missing_fdata <- setdiff(dataset_names, names(fdata))
      if (length(missing_fdata)) {
        stop("fdata is missing entries for dataset(s): ", paste(missing_fdata, collapse = ", "), call. = FALSE)
      }
      fdata <- fdata[dataset_names]
    } else {
      names(fdata) <- dataset_names
    }
    if (is.null(enrichment_id_map)) {
      enrichment_id_map <- vector("list", length(dataset_names))
      names(enrichment_id_map) <- dataset_names
    }
    for (ds in dataset_names) {
      .benchmark_validate_fdata_matches_features(fdata[[ds]], rownames(Y[[ds]]), context = paste0("Y[[", ds, "]]"))
      if (is.null(enrichment_id_map[[ds]])) {
        enrichment_id_map[[ds]] <- .benchmark_fdata_to_enrichment_id_map(fdata[[ds]], required = FALSE)
      }
    }
  }

  enrichment_backend <- match.arg(enrichment_backend)
  enrichment_id_type <- match.arg(enrichment_id_type)

  split <- .make_real_world_integrated_train_test_split(
    n_samples = n_samples[[1]],
    train_fraction = train_fraction,
    train_cols = train_cols,
    test_cols = test_cols,
    seed = seed,
    sample_names = colnames(Y[[1]])
  )
  has_test <- length(split$test) > 0L
  run_inference_effective <- isTRUE(run_inference) && has_test
  if (isTRUE(run_inference) && !has_test && isTRUE(verbose)) {
    message("Skipping integrated real-world inference because the test split is empty (e.g., train_fraction = 1).")
  }

  if (isTRUE(standardize_features) && .benchmark_family_has_count(family)) {
    message("standardize_features=TRUE ignored for poisson/quasipoisson family; using raw count data.")
    standardize_features <- FALSE
  }
  if (block_scale != "none") {
    stop(
      "Global block_scale is no longer supported; block normalization is method-specific.",
      call. = FALSE
    )
  }

  if (standardize_features) {
    std <- standardize_Y_list(Y, truth = NULL)
    Y <- std$Y
  } else {
    std <- NULL
  }

  train_inputs <- .subset_integrated_real_world_inputs(Y, xdata, zdata, split$train)
  test_inputs <- .subset_integrated_real_world_inputs(Y, xdata, zdata, split$test)
  Y_baseline_train <- if (.benchmark_family_has_count(family)) lapply(train_inputs$Y, log1p) else train_inputs$Y
  Y_method_train <- if (adjust_for_z && !is.null(zformula)) {
    .benchmark_adjust_Y_list_by_formula(
      Y = Y_baseline_train,
      zdata = train_inputs$zdata,
      zformula = zformula
    )
  } else {
    Y_baseline_train
  }
  Y_method_test <- NULL
  if (has_test) {
    Y_baseline_test <- if (.benchmark_family_has_count(family)) lapply(test_inputs$Y, log1p) else test_inputs$Y
    Y_method_test <- if (adjust_for_z && !is.null(zformula)) {
      .benchmark_adjust_Y_list_by_formula(
        Y = Y_baseline_test,
        zdata = test_inputs$zdata,
        zformula = zformula
      )
    } else {
      Y_baseline_test
    }
  }

  results <- list()

  if (run_factorization) {
    if ("isvdx" %in% methods) {
      results$isvdx <- fit_isvdx_wrapper(
        Y = train_inputs$Y,
        K = K,
        xformula = xformula,
        zformula = zformula,
        xdata = train_inputs$xdata,
        zdata = train_inputs$zdata,
        family = family,
        graphData = graphData,
        graphDistFn = graphDistFn,
        numThreads = numThreads,
        modelType = svdx_model_type,
        ...
      )
    }

    if ("isvdx_adjusted" %in% methods) {
      if (.benchmark_family_has_count(family)) {
        stop("isvdx_adjusted is currently defined only for Gaussian integrated real-world runs.", call. = FALSE)
      }
      adjusted_zformula <- if (isTRUE(isvdx_adjusted_use_z)) zformula else rep(list(~1), length(Y_method_train))
      results$isvdx_adjusted <- fit_isvdx_wrapper(
        Y = Y_method_train, K = K, xformula = xformula, zformula = adjusted_zformula,
        xdata = train_inputs$xdata, zdata = train_inputs$zdata, family = family,
        graphData = graphData, graphDistFn = graphDistFn, numThreads = numThreads,
        modelType = svdx_model_type, ...
      )
      if (!is.null(results$isvdx_adjusted)) results$isvdx_adjusted$method <- "isvdx_adjusted"
    }

    if ("separate_svdx" %in% methods) {
      results$separate_svdx <- fit_separate_svdx_wrapper(
        Y = train_inputs$Y,
        K = K,
        xformula = xformula,
        zformula = zformula,
        xdata = train_inputs$xdata,
        zdata = train_inputs$zdata,
        family = family,
        graphData = graphData,
        graphDistFn = graphDistFn,
        numThreads = numThreads,
        modelType = svdx_model_type,
        ...
      )
    }

    if ("MEFISTO" %in% methods) {
      results$MEFISTO <- fit_mefisto_wrapper(
        Y = Y_method_train,
        xdata = train_inputs$xdata,
        covariates = mefisto_covariates,
        K = K,
        maxiter = mefisto_maxiter,
        seed = seed
      )
    }

    if ("concat_PCA" %in% methods) results$concat_PCA <- fit_concat_pca(Y_method_train, K = K)
    if ("MFA" %in% methods) results$MFA <- fit_mfa(Y_method_train, K = K)
    if ("concat_NMF" %in% methods) results$concat_NMF <- fit_concat_nmf(Y_method_train, K = K, seed = seed)

    if ("concat_sPCA" %in% methods) {
      if (is.null(spca_keepX_grid)) {
        spca_keepX_grid <- round(c(0.01, 0.10, 0.25, 0.50, 0.75) * sum(vapply(Y_method_train, nrow, integer(1))))
      }
      results <- c(results, fit_concat_spca_grid(Y = Y_method_train, K = K, keepX_grid = spca_keepX_grid))
    }

    if ("concat_sNMF" %in% methods) {
      results <- c(
        results,
        fit_concat_sparse_nmf_grid(
          Y = Y_method_train,
          K = K,
          lambda_grid = sparse_nmf_lambda_grid,
          seed = seed
        )
      )
    }

    if ("block_PLS" %in% methods) {
      results$block_PLS <- fit_block_pls_wrapper(
        Y = Y_method_train, K = K, xformula = xformula, xdata = train_inputs$xdata
      )
    }

    if ("block_sPLS" %in% methods) {
      results <- c(results, fit_block_spls_grid(
        Y = Y_method_train, K = K, xformula = xformula, xdata = train_inputs$xdata,
        keepX_grid = block_spls_keepX_grid
      ))
    }
  }

  results <- Filter(Negate(is.null), results)
  feature_scores <- extract_integrated_feature_scores(results)

  loading_enrichment <- NULL
  if (isTRUE(do_loading_enrichment)) {
    loading_enrichment <- evaluate_integrated_real_world_loading_enrichment(
      integrated_results = results,
      loading_scope = loading_scope,
      loading_top_n = loading_top_n,
      feature_sets = feature_sets,
      enrichment_id_map = enrichment_id_map,
      alpha = alpha,
      enrichment_backend = enrichment_backend,
      msigdb_db_species = msigdb_db_species,
      msigdb_species = msigdb_species,
      msigdb_collection = msigdb_collection,
      msigdb_subcollection = msigdb_subcollection,
      enrichment_id_type = enrichment_id_type,
      enrichment_score = enrichment_score,
      term2gene = enrichment_term2gene,
      universe = enrichment_universe,
      min_feature_overlap_warn = min_feature_overlap_warn,
      min_feature_overlap_fraction = min_feature_overlap_fraction,
      enrichment_extra_args = enrichment_extra_args
    )
  }

  timing <- extract_integrated_timing(results)
  inference <- list()

  if (run_inference_effective) {
    if ("isvdx" %in% inference_methods) {
      isvdx_fit <- if (!is.null(results$isvdx) && !is.null(results$isvdx$fit)) results$isvdx$fit else NULL
      time_fit <- 0
      if (is.null(isvdx_fit)) {
        tm_fit <- system.time({
          isvdx_fit <- isvdx(
            Y = train_inputs$Y,
            K = K,
            formula = xformula,
            data = train_inputs$xdata,
            zformula = zformula,
            zdata = train_inputs$zdata,
            family = family,
            graphData = graphData,
            graphDistFn = graphDistFn,
            numThreads = numThreads,
            modelType = svdx_model_type,
            ...
          )
        })
        time_fit <- unname(tm_fit["elapsed"])
      }

      tm_inf <- system.time({
        inf <- inferV(
          fit = isvdx_fit,
          Ytest = test_inputs$Y,
          newdata = test_inputs$xdata,
          newzdata = test_inputs$zdata,
          zformula = zformula,
          family = family,
          alpha = alpha,
          z_adjustment = z_adjustment,
          return_design = return_design,
          verbose = verbose
        )
      })

      inference_enrichment <- NULL
      if (do_inference_enrichment) {
        dataset_names <- names(test_inputs$Y) %||% paste0("dataset_", seq_along(test_inputs$Y))
        inference_enrichment <- dplyr::bind_rows(lapply(dataset_names, function(ds) {
          if (is.null(enrichment_id_map) || is.null(enrichment_id_map[[ds]])) return(NULL)
          ep_args <- c(
            list(
              inf = inf,
              dataset = ds,
              feature_sets = if (!is.null(feature_sets)) feature_sets[[ds]] else NULL,
              alpha = alpha,
              scope = inference_scope,
              universe = enrichment_universe,
              enrichment_backend = enrichment_backend,
              msigdb_db_species = msigdb_db_species,
              msigdb_species = msigdb_species,
              msigdb_collection = msigdb_collection,
              msigdb_subcollection = msigdb_subcollection,
              enrichment_id_type = enrichment_id_type,
              enrichment_score = enrichment_score,
              term2gene = enrichment_term2gene,
              enrichment_id_map = if (!is.null(enrichment_id_map)) enrichment_id_map[[ds]] else NULL,
              enrichment_extra_args = enrichment_extra_args,
              min_feature_overlap_warn = min_feature_overlap_warn,
              min_feature_overlap_fraction = min_feature_overlap_fraction
            )
          )
          enr <- do.call(enrichPathways, ep_args)
          if (is.null(enr) || is.null(enr$results)) return(NULL)
          enr$results |>
            dplyr::mutate(dataset = ds, .before = 1)
        }))
      }

      inference$isvdx <- list(
        fit = isvdx_fit,
        inf = inf,
        inf_summary = summary(inf, alpha = alpha),
        inference_enrichment = inference_enrichment,
        timing = tibble::tibble(
          method = "isvdx",
          time_fit = time_fit,
          time_inference = unname(tm_inf["elapsed"]),
          time_total = time_fit + unname(tm_inf["elapsed"])
        )
      )
    }

    if ("isvdx_adjusted" %in% inference_methods) {
      if (.benchmark_family_has_count(family)) {
        stop("isvdx_adjusted inference is currently defined only for Gaussian integrated real-world runs.", call. = FALSE)
      }
      adjusted_zformula <- if (isTRUE(isvdx_adjusted_use_z)) zformula else rep(list(~1), length(Y_method_train))
      isvdx_adjusted_fit <- if (!is.null(results$isvdx_adjusted) && !is.null(results$isvdx_adjusted$fit)) results$isvdx_adjusted$fit else NULL
      time_fit <- 0
      if (is.null(isvdx_adjusted_fit)) {
        tm_fit <- system.time({
          isvdx_adjusted_fit <- isvdx(
            Y = Y_method_train, K = K, formula = xformula, data = train_inputs$xdata,
            zformula = adjusted_zformula, zdata = train_inputs$zdata, family = family,
            graphData = graphData, graphDistFn = graphDistFn, numThreads = numThreads,
            modelType = svdx_model_type, ...
          )
        })
        time_fit <- unname(tm_fit["elapsed"])
      }
      tm_inf <- system.time({
        inf <- inferV(
          fit = isvdx_adjusted_fit, Ytest = Y_method_test, newdata = test_inputs$xdata,
          newzdata = test_inputs$zdata, zformula = adjusted_zformula, family = family,
          alpha = alpha, z_adjustment = z_adjustment, return_design = return_design,
          verbose = verbose
        )
      })
      inference_enrichment <- NULL
      if (do_inference_enrichment) {
        dataset_names <- names(test_inputs$Y) %||% paste0("dataset_", seq_along(test_inputs$Y))
        inference_enrichment <- dplyr::bind_rows(lapply(dataset_names, function(ds) {
          if (is.null(enrichment_id_map) || is.null(enrichment_id_map[[ds]])) return(NULL)
          ep_args <- c(list(
            inf = inf, dataset = ds,
            feature_sets = if (!is.null(feature_sets)) feature_sets[[ds]] else NULL,
            alpha = alpha, scope = inference_scope, universe = enrichment_universe,
            enrichment_backend = enrichment_backend, msigdb_db_species = msigdb_db_species,
            msigdb_species = msigdb_species, msigdb_collection = msigdb_collection,
            msigdb_subcollection = msigdb_subcollection, enrichment_id_type = enrichment_id_type,
            enrichment_score = enrichment_score, term2gene = enrichment_term2gene,
            enrichment_id_map = if (!is.null(enrichment_id_map)) enrichment_id_map[[ds]] else NULL,
            enrichment_extra_args = enrichment_extra_args,
            min_feature_overlap_warn = min_feature_overlap_warn,
            min_feature_overlap_fraction = min_feature_overlap_fraction
          ))
          enr <- do.call(enrichPathways, ep_args)
          if (is.null(enr) || is.null(enr$results)) return(NULL)
          enr$results |> dplyr::mutate(dataset = ds, .before = 1)
        }))
      }
      inference$isvdx_adjusted <- list(
        fit = isvdx_adjusted_fit, inf = inf, inf_summary = summary(inf, alpha = alpha),
        inference_enrichment = inference_enrichment,
        timing = tibble::tibble(
          method = "isvdx_adjusted", time_fit = time_fit,
          time_inference = unname(tm_inf["elapsed"]),
          time_total = time_fit + unname(tm_inf["elapsed"])
        )
      )
    }

    if ("separate_svdx" %in% inference_methods) {
      dataset_names <- names(train_inputs$Y) %||% paste0("dataset_", seq_along(train_inputs$Y))
      sep_fits <- list()
      sep_infs <- list()
      time_fit <- numeric(length(dataset_names)); names(time_fit) <- dataset_names
      time_inf <- numeric(length(dataset_names)); names(time_inf) <- dataset_names

      trained_sep <- if (!is.null(results$separate_svdx) && !is.null(results$separate_svdx$fit)) results$separate_svdx$fit else NULL

      for (ds in dataset_names) {
        fam_d <- if (is.null(family)) NULL else family[[ds]]
        if (!is.null(trained_sep) && !is.null(trained_sep[[ds]]) && !is.null(trained_sep[[ds]]$fit)) {
          sep_fits[[ds]] <- trained_sep[[ds]]$fit
          time_fit[ds] <- 0
        } else {
          tm_fit <- system.time({
            sep_fits[[ds]] <- svdx(
              Y = train_inputs$Y[[ds]],
              K = K,
              formula = xformula,
              data = train_inputs$xdata,
              zformula = zformula[[ds]],
              zdata = train_inputs$zdata[[ds]],
              family = fam_d,
              graphData = graphData,
              graphDistFn = graphDistFn,
              numThreads = numThreads,
              modelType = svdx_model_type,
              ...
            )
          })
          time_fit[ds] <- unname(tm_fit["elapsed"])
        }

        tm_inf <- system.time({
          sep_infs[[ds]] <- inferV(
            fit = sep_fits[[ds]],
            Ytest = test_inputs$Y[[ds]],
            newdata = test_inputs$xdata,
            newzdata = test_inputs$zdata[[ds]],
            zformula = zformula[[ds]],
            family = fam_d,
            alpha = alpha,
            z_adjustment = z_adjustment,
            return_design = return_design,
            verbose = verbose
          )
        })
        time_inf[ds] <- unname(tm_inf["elapsed"])
      }

      inference_enrichment <- NULL
      if (do_inference_enrichment) {
        inference_enrichment <- dplyr::bind_rows(lapply(dataset_names, function(ds) {
          if (is.null(enrichment_id_map) || is.null(enrichment_id_map[[ds]])) return(NULL)
          ep_args <- c(
            list(
              inf = sep_infs[[ds]],
              feature_sets = if (!is.null(feature_sets)) feature_sets[[ds]] else NULL,
              alpha = alpha,
              scope = inference_scope,
              universe = enrichment_universe,
              enrichment_backend = enrichment_backend,
              msigdb_db_species = msigdb_db_species,
              msigdb_species = msigdb_species,
              msigdb_collection = msigdb_collection,
              msigdb_subcollection = msigdb_subcollection,
              enrichment_id_type = enrichment_id_type,
              enrichment_score = enrichment_score,
              term2gene = enrichment_term2gene,
              enrichment_id_map = if (!is.null(enrichment_id_map)) enrichment_id_map[[ds]] else NULL,
              enrichment_extra_args = enrichment_extra_args,
              min_feature_overlap_warn = min_feature_overlap_warn,
              min_feature_overlap_fraction = min_feature_overlap_fraction
            )
          )
          enr <- do.call(enrichPathways, ep_args)
          if (is.null(enr) || is.null(enr$results)) return(NULL)
          enr$results |>
            dplyr::mutate(dataset = ds, .before = 1)
        }))
      }

      inference$separate_svdx <- list(
        fit = sep_fits,
        inf = sep_infs,
        inf_summary = lapply(sep_infs, summary, alpha = alpha),
        inference_enrichment = inference_enrichment,
        timing = tibble::tibble(
          method = "separate_svdx",
          time_fit = sum(time_fit, na.rm = TRUE),
          time_inference = sum(time_inf, na.rm = TRUE),
          time_total = sum(time_fit + time_inf, na.rm = TRUE)
        ),
        timing_by_dataset = tibble::tibble(
          dataset = dataset_names,
          time_fit = as.numeric(time_fit),
          time_inference = as.numeric(time_inf),
          time_total = as.numeric(time_fit + time_inf)
        )
      )
    }
  }

  prediction_metrics <- list()
  if (has_test && !is.null(results$isvdx) && !is.null(results$isvdx$fit)) {
    prediction_metrics$isvdx <- tryCatch(
      .benchmark_isvdx_prediction_metrics(
        fit = results$isvdx$fit,
        Ytrain = train_inputs$Y,
        Ytest = test_inputs$Y,
        xdata_test = test_inputs$xdata,
        zdata_test = test_inputs$zdata,
        zformula = zformula,
        method = "isvdx",
        numThreads = numThreads
      ),
      error = function(e) {
        warning("Could not compute ISVDX test prediction metrics: ", e$message, call. = FALSE)
        NULL
      }
    )
  }
  if (has_test && !is.null(results$separate_svdx) && !is.null(results$separate_svdx$fit)) {
    prediction_metrics$separate_svdx <- tryCatch(
      .benchmark_separate_svdx_prediction_metrics(
        fits = results$separate_svdx$fit,
        Ytrain = train_inputs$Y,
        Ytest = test_inputs$Y,
        xdata_test = test_inputs$xdata,
        zdata_test = test_inputs$zdata,
        zformula = zformula,
        method = "separate_svdx",
        numThreads = numThreads
      ),
      error = function(e) {
        warning("Could not compute separate SVDX test prediction metrics: ", e$message, call. = FALSE)
        NULL
      }
    )
  }
  prediction_metrics <- dplyr::bind_rows(prediction_metrics)
  if (nrow(prediction_metrics) == 0L) prediction_metrics <- NULL

  out <- list(
    Y = train_inputs$Y,
    xdata = train_inputs$xdata,
    zdata = train_inputs$zdata,
    fdata = fdata,
    feature_sets = feature_sets,
    enrichment_id_map = enrichment_id_map,
    zformula = zformula,
    results = results,
    feature_scores = feature_scores,
    loading_enrichment = loading_enrichment,
    inference = inference,
    alpha = alpha,
    split = split,
    timing = timing,
    standardization = std,
    block_scale = rep(1, length(Y)),
    block_scale_method = "none",
    xformula = xformula,
    family = family,
    modelType = modelType,
    K = K,
    methods = methods,
    inference_methods = inference_methods,
    run_factorization = run_factorization,
    run_inference = run_inference_effective,
    run_inference_requested = run_inference,
    has_test = has_test,
    prediction_metrics = prediction_metrics,
    target_sets = target_sets
  )

  .validate_integrated_real_world_output(out)
  class(out) <- "realWorldIntegratedBenchmark"

  if (!is.null(target_sets)) {
    out$target_auc <- do.call(
      evaluate_integrated_target_set_auc,
      c(list(rw_int = out, target_sets = target_sets), target_set_auc_args)
    )
  }

  out
}


plot_real_world_integrated_U <- function(
    rw,
    method,
    component = 1,
    x,
    color = NULL,
    dataset = NULL,
    connect = FALSE,
    line_group = NULL,
    line_alpha = 0.55,
    line_width = 0.4,
    point_alpha = 0.75
) {
  .validate_real_world_integrated_benchmark(rw)
  res <- .rwbench_get_integrated_method_result(rw, method)

  if (!x %in% colnames(rw$xdata)) {
    stop("x column not found in rw$xdata: ", x, call. = FALSE)
  }
  if (!is.null(color)) {
    missing_color <- setdiff(as.character(color), colnames(rw$xdata))
    if (length(missing_color) > 0L) {
      stop("color column(s) not found in rw$xdata: ",
           paste(missing_color, collapse = ", "), call. = FALSE)
    }
  }
  if (!is.null(line_group)) {
    missing_line_group <- setdiff(as.character(line_group), colnames(rw$xdata))
    if (length(missing_line_group) > 0L) {
      stop("line_group column(s) not found in rw$xdata: ",
           paste(missing_line_group, collapse = ", "), call. = FALSE)
    }
  }

  U <- res$U
  if (is.null(U)) stop("Selected method does not contain a U matrix.", call. = FALSE)

  if (is.list(U)) {
    if (is.null(dataset)) {
      stop("This method has dataset-specific U. Please provide dataset.", call. = FALSE)
    }
    U <- U[[dataset]]
    if (is.null(U)) {
      stop("Dataset not found in method-specific U.", call. = FALSE)
    }
  }

  U <- as.data.frame(U)
  colnames(U) <- colnames(U) %||% paste0("component_", seq_len(ncol(U)))

  k <- if (is.numeric(component)) component else match(component, colnames(U))
  if (is.na(k) || k < 1 || k > ncol(U)) stop("Invalid component.", call. = FALSE)

  df <- cbind(rw$xdata, component_value = U[[k]])

  color_label <- if (is.null(color)) NULL else paste(as.character(color), collapse = " x ")
  color_info <- .benchmark_plot_interaction_var(df, color, name = ".benchmark_plot_color")
  df <- color_info$data
  color_aes <- color_info$var

  if (!is.null(line_group) && identical(as.character(line_group), as.character(color))) {
    line_group_aes <- color_aes
  } else {
    line_info <- .benchmark_plot_interaction_var(df, line_group, name = ".benchmark_plot_line_group")
    df <- line_info$data
    line_group_aes <- line_info$var
  }

  aes_args <- ggplot2::aes(x = .data[[x]], y = component_value)
  if (!is.null(color_aes)) aes_args$colour <- rlang::expr(.data[[!!color_aes]])

  p <- ggplot2::ggplot(df, aes_args)

  if (isTRUE(connect)) {
    line_var <- line_group_aes %||% color_aes
    df_line <- df
    if (!is.null(line_var)) {
      df_line <- df_line[order(df_line[[line_var]], df_line[[x]]), , drop = FALSE]
    } else {
      df_line <- df_line[order(df_line[[x]]), , drop = FALSE]
    }
    line_aes <- aes_args
    if (!is.null(line_var)) {
      line_aes$group <- rlang::expr(.data[[!!line_var]])
    } else {
      line_aes$group <- rlang::expr(1)
    }
    p <- p + ggplot2::geom_line(
      data = df_line,
      mapping = line_aes,
      alpha = line_alpha,
      linewidth = line_width
    )
  }

  p +
    ggplot2::geom_point(alpha = point_alpha) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = x,
      y = paste(method, "component", k),
      color = color_label
    )
}


plot_real_world_integrated_U_trajectory <- function(
    rw,
    method,
    component = 1,
    x,
    group = NULL,
    color = group,
    dataset = NULL,
    line_alpha = 0.55,
    line_width = 0.4,
    point_alpha = 0.75
) {
  plot_real_world_integrated_U(
    rw = rw,
    method = method,
    component = component,
    x = x,
    color = color,
    dataset = dataset,
    connect = !is.null(group),
    line_group = group,
    line_alpha = line_alpha,
    line_width = line_width,
    point_alpha = point_alpha
  )
}

plot_real_world_integrated_U_embedding <- function(
    rw,
    method,
    components = c(1, 2),
    color = NULL,
    dataset = NULL
) {
  .validate_real_world_integrated_benchmark(rw)
  res <- .rwbench_get_integrated_method_result(rw, method)

  U <- res$U
  if (is.null(U)) stop("Selected method does not contain a U matrix.", call. = FALSE)

  if (is.list(U)) {
    if (is.null(dataset)) {
      stop("This method has dataset-specific U. Please provide dataset.", call. = FALSE)
    }
    U <- U[[dataset]]
    if (is.null(U)) {
      stop("Dataset not found in method-specific U.", call. = FALSE)
    }
  }

  U <- as.data.frame(U)
  colnames(U) <- colnames(U) %||% paste0("component_", seq_len(ncol(U)))

  if (length(components) < 2L) {
    stop("components must contain at least two entries.", call. = FALSE)
  }
  k1 <- if (is.numeric(components[1])) components[1] else match(components[1], colnames(U))
  k2 <- if (is.numeric(components[2])) components[2] else match(components[2], colnames(U))
  if (any(is.na(c(k1, k2))) || any(c(k1, k2) < 1) || any(c(k1, k2) > ncol(U))) {
    stop("Invalid components.", call. = FALSE)
  }

  df <- cbind(
    rw$xdata,
    comp1 = U[[k1]],
    comp2 = U[[k2]]
  )

  aes_args <- ggplot2::aes(x = comp1, y = comp2)
  if (!is.null(color)) aes_args$colour <- rlang::expr(.data[[!!color]])

  ggplot2::ggplot(df, aes_args) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = paste(method, "component", k1),
      y = paste(method, "component", k2),
      color = color
    )
}


plot_real_world_integrated_V_heatmap <- function(
    rw,
    method,
    dataset,
    component = NULL,
    n_features = 50,
    features = NULL,
    order_by = NULL,
    annotate = NULL,
    scale_rows = TRUE,
    row_labels = NULL,
    fdata = rw$fdata %||% NULL
) {
  .validate_real_world_integrated_benchmark(rw)
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required.", call. = FALSE)
  }

  res <- .rwbench_get_integrated_method_result(rw, method)

  d <- if (is.numeric(dataset)) dataset else match(dataset, names(rw$Y))
  if (is.na(d)) stop("Unknown dataset.", call. = FALSE)

  if (is.null(res$V_list)) stop("Selected method does not contain V_list.", call. = FALSE)
  if (is.null(res$feature_score_list)) stop("Selected method does not contain feature_score_list.", call. = FALSE)
  V <- res$V_list[[d]]
  scores <- res$feature_score_list[[d]]
  if (is.null(V)) stop("Selected dataset does not contain a V matrix.", call. = FALSE)
  if (is.null(scores)) stop("Selected dataset does not contain feature scores.", call. = FALSE)
  feature_names <- rownames(V) %||% rownames(rw$Y[[d]]) %||% paste0("feature_", seq_len(nrow(V)))
  dataset_name <- names(rw$Y)[d]

  names(scores) <- feature_names

  if (!is.null(component)) {
    k <- if (is.numeric(component)) component else match(component, colnames(V))
    if (is.na(k) || k < 1 || k > ncol(V)) stop("Invalid component.", call. = FALSE)
    scores <- abs(V[, k])
    names(scores) <- feature_names
  }

  if (is.null(features)) {
    features <- names(sort(scores, decreasing = TRUE))[
      seq_len(min(n_features, length(scores)))
    ]
  }
  resolved_features <- .rwbench_resolve_features_and_labels(
    features = features,
    feature_names = feature_names,
    fdata = fdata,
    row_labels = row_labels,
    dataset = dataset_name
  )
  top_features <- resolved_features$feature
  feature_labels <- resolved_features$label

  mat <- rw$Y[[d]][top_features, , drop = FALSE]
  rownames(mat) <- feature_labels

  if (scale_rows) {
    mat <- t(scale(t(mat)))
    mat[!is.finite(mat)] <- 0
  }

  meta <- rw$xdata

  if (!is.null(order_by)) {
    missing_order <- setdiff(order_by, colnames(meta))
    if (length(missing_order) > 0L) {
      stop("order_by columns not found in xdata: ", paste(missing_order, collapse = ", "), call. = FALSE)
    }
    ord <- do.call(order, meta[order_by])
    mat <- mat[, ord, drop = FALSE]
    meta <- meta[ord, , drop = FALSE]
  }

  ann <- NULL
  if (!is.null(annotate)) {
    missing_annotate <- setdiff(annotate, colnames(meta))
    if (length(missing_annotate) > 0L) {
      stop("annotate columns not found in xdata: ", paste(missing_annotate, collapse = ", "), call. = FALSE)
    }
    ann <- meta[, annotate, drop = FALSE]
    rownames(ann) <- colnames(mat)
  }

  pheatmap::pheatmap(
    mat,
    annotation_col = ann,
    show_rownames = TRUE,
    show_colnames = FALSE,
    main = paste(method, names(rw$Y)[d], "top loading features")
  )
}

plot_real_world_isvdx_U_terms <- function(
    rw,
    component = 1,
    terms = NULL,
    x = "t",
    group = NULL,
    fill = "typical",
    n_grid = 200,
    ...
) {
  res <- rw$results$isvdx
  if (is.null(res)) {
    stop("rw does not contain an isvdx fit.", call. = FALSE)
  }

  fit <- res$fit
  if (is.null(fit)) {
    stop("isvdx result does not contain a fit object.", call. = FALSE)
  }

  plotUTerms_fn <- .benchmark_get_svdx_internal("plotUTerms")

  plotUTerms_fn(
    fit = fit,
    component = component,
    terms = terms,
    x = x,
    group = group,
    fill = fill,
    n_grid = n_grid,
    ...
  )
}
