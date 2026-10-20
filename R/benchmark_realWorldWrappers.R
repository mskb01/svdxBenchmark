# ============================================================
# Real-world across-method benchmarking wrappers
# ============================================================
# No ground truth assumed.
#
# Produces:
#   A) V/loadings-based summaries:
#      - feature scores
#      - enrichment from loadings
#      - top-feature heatmaps
#
#   B) U/component summaries:
#      - fitted component plots against metadata
#
#   C) svdx-specific:
#      - inferV()
#      - inference enrichment
# ============================================================

#######USAGE######


#
# rw <- run_real_world_method_benchmark(
#   Y = Y,
#   data = meta,
#   xformula = ~ group * t + severity * t,
#   zformula = ~ batch,
#   K = 6,
#   feature_sets = feature_sets,
#   methods = c("svdx", "svd", "pca", "nmf", "mixomics_pls", "feature_lm")
# )
#
# rw$loading_enrichment
#
# plot_real_world_U(rw, method = "PCA", component = 1, x = "t", color = "group")
#
# plot_real_world_V_heatmap(
#   rw,
#   method = "svdx",
#   n_features = 50,
#   order_by = c("group", "t"),
#   annotate = c("group", "batch")
# )
#
# svdx_inf <- run_real_world_svdx_inference(
#   rw,
#   feature_sets = feature_sets,
#   alpha = 0.1
# )
#
# plot_union_volcano(svdx_inf$inf)
# plot_enrichment_dotplot(svdx_inf$enrichment)
# plot_real_world_svdx_U(
#   rw,
#   component = 1,
#   x = "t",
#   group = "group",
#   fill = "typical",
#   n_grid = 200
# )
# plot_real_world_svdx_U( #2D
#   rw,
#   component = 1,
#   x = "t",
#   y = "severity",
#   fixed = list(group = "1"),
#   fill = "typical",
#   n_grid = 100
# )
# plot_real_world_svdx_U_terms( #term-wise interpol
#   rw,
#   component = 1,
#   terms = c("group:t", "severity:t"),
#   x = "t",
#   group = "group",
#   fill = "typical"
# )
# So the real-world wrapper should have two U plot routes:
#
#   plot_real_world_U()
#
# for generic methods using observed fitted scores only.


# ------------------------------------------------------------
# Feature label helpers for real-world plotting
# ------------------------------------------------------------
# These are intentionally defined here as fallbacks, even though
# benchmark_utils.R also defines them, so real-world plotting remains
# robust when wrapper files are sourced after/without the newest utility file.
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

if (!exists(".benchmark_family_has_count", mode = "function")) {
  .benchmark_family_has_count <- function(family) {
    if (is.null(family)) return(FALSE)
    any(unlist(family, use.names = FALSE) %in% c("poisson", "quasipoisson"), na.rm = TRUE)
  }
}


.safe_fisher_enrichment <- function(discovered, in_set) {
  fn <- .benchmark_get_svdx_internal(".safe_fisher_enrichment")
  fn(discovered, in_set)
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
#
# plot_real_world_svdx_U()
# plot_real_world_svdx_U_terms()
#
# plot_real_world_U_embedding(
#   rw,
#   method = "PCA",
#   components = c(1, 2),
#   color = "group"
# )
# plot_real_world_svdx_U_embedding( #svdx observed scores
#   rw,
#   components = c(1, 2),
#   color = "group"
# )
# plot_real_world_svdx_U_embedding( #interpolated
#   rw,
#   components = c(1, 2),
#   vary = "t",
#   group = "group",
#   color = "group",
#   fill = "typical",
#   n_grid = 200
# )
# plot_real_world_svdx_U_embedding( #term specific
#   rw,
#   components = c(1, 2),
#   include = "severity:t",
#   vary = "t",
#   color = "group",
#   group = "group",
#   fill = "typical"
# )



# ------------------------------------------------------------
# Internal validation helpers for real-world wrappers
# ------------------------------------------------------------

.validate_real_world_graph_inputs <- function(graphData, graphDistFn, n_samples) {
  if (!is.null(graphData)) {
    if (!is.data.frame(graphData)) {
      stop("graphData must be a data.frame when supplied.", call. = FALSE)
    }
    if (nrow(graphData) < 1L) {
      stop("graphData must contain at least one row.", call. = FALSE)
    }
    if (ncol(graphData) != n_samples) {
      stop("ncol(graphData) must equal ncol(Y).", call. = FALSE)
    }
  }
  if (!is.null(graphDistFn) && !is.function(graphDistFn)) {
    stop("graphDistFn must be a function when supplied.", call. = FALSE)
  }
  invisible(TRUE)
}

.validate_feature_sets_basic <- function(feature_sets, feature_names = NULL) {
  if (is.null(feature_sets)) return(invisible(TRUE))
  if (length(feature_sets) == 0L) {
    stop("feature_sets must contain at least one feature set.", call. = FALSE)
  }
  if (!is.null(feature_names)) {
    fs_features <- unique(unlist(feature_sets, use.names = FALSE))
    if (length(fs_features) > 0L && length(intersect(fs_features, feature_names)) == 0L) {
      warning(
        "No overlap detected between supplied feature sets and benchmark feature names.",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

.validate_real_world_method_benchmark_output <- function(x) {
  required <- c("benchmark_results", "raw_fits", "feature_scores", "loading_enrichment", "methods")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop(
      "Internal error: real-world benchmark output is missing required fields: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.validate_real_world_method_benchmark <- function(rw_bench) {
  if (!inherits(rw_bench, "realWorldMethodBenchmark")) {
    stop("rw_bench must inherit from 'realWorldMethodBenchmark'.", call. = FALSE)
  }
  invisible(TRUE)
}

.get_real_world_method_result <- function(rw_bench, method) {
  .validate_real_world_method_benchmark(rw_bench)
  if (is.null(rw_bench$benchmark_results) || !is.list(rw_bench$benchmark_results)) {
    stop("rw_bench$benchmark_results must be a named/list object.", call. = FALSE)
  }

  # Prefer exact list-name lookup. Real-world benchmark_results are
  # method-keyed lists, and some method result objects may carry display
  # labels or legacy method fields that do not exactly match the requested key.
  if (!is.null(names(rw_bench$benchmark_results)) && method %in% names(rw_bench$benchmark_results)) {
    res <- rw_bench$benchmark_results[[method]]
    if (is.null(res)) {
      stop("Method not found in rw_bench$benchmark_results: ", method, call. = FALSE)
    }
    return(res)
  }

  idx <- which(vapply(
    rw_bench$benchmark_results,
    function(z) identical(z$method, method),
    logical(1)
  ))
  if (length(idx) == 0L) {
    available <- names(rw_bench$benchmark_results)
    available <- available[!is.na(available) & nzchar(available)]
    msg <- paste0("Method not found in rw_bench$benchmark_results: ", method)
    if (length(available)) {
      msg <- paste0(msg, ". Available methods: ", paste(available, collapse = ", "))
    }
    stop(msg, call. = FALSE)
  }
  if (length(idx) > 1L) {
    stop("Multiple benchmark results found for method: ", method, call. = FALSE)
  }
  rw_bench$benchmark_results[[idx[1]]]
}


# ------------------------------------------------------------
# Canonical real-world plotting helpers
# ------------------------------------------------------------
# These helpers intentionally use rwbench-specific names and are defined
# unconditionally.  This prevents stale helper definitions from older sourced
# benchmark files from silently changing plotting behavior.

.rwbench_get_method_result <- function(rw_bench, method) {
  .validate_real_world_method_benchmark(rw_bench)
  if (is.null(rw_bench$benchmark_results) || !is.list(rw_bench$benchmark_results)) {
    stop("rw_bench$benchmark_results must be a named/list object.", call. = FALSE)
  }

  if (!is.null(names(rw_bench$benchmark_results)) && method %in% names(rw_bench$benchmark_results)) {
    res <- rw_bench$benchmark_results[[method]]
    if (is.null(res)) {
      stop("Method not found in rw_bench$benchmark_results: ", method, call. = FALSE)
    }
    return(res)
  }

  idx <- which(vapply(
    rw_bench$benchmark_results,
    function(z) identical(z$method, method),
    logical(1)
  ))

  if (length(idx) == 0L) {
    available <- names(rw_bench$benchmark_results)
    available <- available[!is.na(available) & nzchar(available)]
    msg <- paste0("Method not found in rw_bench$benchmark_results: ", method)
    if (length(available)) {
      msg <- paste0(msg, ". Available methods: ", paste(available, collapse = ", "))
    }
    stop(msg, call. = FALSE)
  }
  if (length(idx) > 1L) {
    stop("Multiple benchmark results found for method: ", method, call. = FALSE)
  }
  rw_bench$benchmark_results[[idx[1]]]
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

  if (is.null(labs)) {
    return(data.frame(feature = feature_names, label = feature_names, stringsAsFactors = FALSE))
  }

  if (is.vector(labs) && !is.list(labs)) {
    if (is.null(names(labs))) {
      if (length(labs) != length(feature_names)) {
        stop("Unnamed row_labels/fdata labels must have the same length as feature names.", call. = FALSE)
      }
      return(data.frame(feature = feature_names, label = as.character(labs), stringsAsFactors = FALSE))
    }
    return(data.frame(feature = as.character(names(labs)), label = as.character(labs), stringsAsFactors = FALSE))
  }

  if (is.data.frame(labs)) {
    if (ncol(labs) < 2L) {
      stop("row_labels as a data frame must have at least two columns: feature and display label.", call. = FALSE)
    }
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

  label_map <- .rwbench_fdata_row_labels(
    feature_names = feature_names,
    fdata = fdata,
    row_labels = row_labels,
    dataset = dataset
  )

  resolved <- requested
  unresolved <- !(requested %in% feature_names)
  if (any(unresolved)) {
    dup_labels <- unique(label_map$label[duplicated(label_map$label)])
    requested_labels <- intersect(requested[unresolved], label_map$label)
    ambiguous <- intersect(requested_labels, dup_labels)
    if (length(ambiguous)) {
      stop(
        "Feature label(s) map to multiple features: ", paste(ambiguous, collapse = ", "),
        ". Please specify original feature IDs or unique labels.",
        call. = FALSE
      )
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
  hits <- match(resolved, label_map$feature)
  use <- !is.na(hits)
  labels[use] <- label_map$label[hits[use]]

  list(feature = resolved, label = labels)
}

# ------------------------------------------------------------
# Feature scores from fitted method results
# ------------------------------------------------------------

extract_method_feature_scores <- function(bench_results) {
  dplyr::bind_rows(lapply(bench_results, function(res) {
    feature_names <- rownames(res$V) %||% paste0("feature_", seq_len(nrow(res$V)))

    tibble::tibble(
      method = res$method,
      feature = feature_names,
      feature_score = res$feature_score
    )
  }))
}


# ------------------------------------------------------------
# Enrichment from loadings / feature scores
# ------------------------------------------------------------

enrich_method_loadings <- function(
    bench_results,
    feature_sets,
    top_n = 100,
    alpha = NULL
) {
  dplyr::bind_rows(lapply(bench_results, function(res) {
    if (is.null(res$V) || is.null(res$feature_score) || length(res$feature_score) == 0L) {
      return(NULL)
    }
    feature_names <- rownames(res$V) %||% paste0("feature_", seq_len(nrow(res$V)))
    scores <- res$feature_score
    names(scores) <- feature_names

    ranked <- sort(scores, decreasing = TRUE)
    selected <- names(ranked)[seq_len(min(top_n, length(ranked)))]

    universe <- feature_names

    out <- lapply(names(feature_sets), function(set_name) {
      set_features <- intersect(feature_sets[[set_name]], universe)

      discovered <- universe %in% selected
      in_set <- universe %in% set_features

      sf <- .safe_fisher_enrichment(discovered, in_set)

      tibble::tibble(
        method = res$method,
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
}



.make_real_world_train_test_split <- function(
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
        stop("Character train/test columns require colnames(Y).", call. = FALSE)
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

.enrich_feature_lm_from_pvalues <- function(p_mat, feature_sets, alpha = 0.1) {
  if (is.null(feature_sets) || is.null(p_mat)) return(NULL)
  feature_names <- rownames(p_mat) %||% paste0("feature_", seq_len(nrow(p_mat)))
  p_min <- apply(p_mat, 1, min, na.rm = TRUE)
  q_min <- .safe_bh(p_min)
  selected <- feature_names[which(q_min <= alpha)]
  universe <- feature_names

  dplyr::bind_rows(lapply(names(feature_sets), function(set_name) {
    set_features <- intersect(feature_sets[[set_name]], universe)
    discovered <- universe %in% selected
    in_set <- universe %in% set_features
    sf <- .safe_fisher_enrichment(discovered, in_set)
    tibble::tibble(
      method = "feature_lm",
      pathway = set_name,
      set_size = length(set_features),
      n_selected = length(selected),
      overlap = sum(discovered & in_set),
      odds_ratio = sf$odds_ratio,
      p_value = sf$p_value
    )
  })) |>
    dplyr::mutate(q_value = .safe_bh(p_value)) |>
    dplyr::arrange(q_value, p_value)
}


.benchmark_validate_prepared_single <- function(prepared) {
  if (!is.list(prepared) || is.null(prepared$Y) || is.null(prepared$xdata) || is.null(prepared$xformula)) {
    stop("Prepared single object must contain Y, xdata, and xformula.", call. = FALSE)
  }
  if (isTRUE(prepared$integrated)) {
    stop("Prepared integrated objects should be passed to run_real_world_integrated_benchmark().", call. = FALSE)
  }
  if (!is.matrix(prepared$Y) && !is.data.frame(prepared$Y)) {
    stop("prepared$Y must be a feature-by-sample matrix-like object.", call. = FALSE)
  }
  if (ncol(as.matrix(prepared$Y)) != nrow(as.data.frame(prepared$xdata))) {
    stop("ncol(prepared$Y) must equal nrow(prepared$xdata).", call. = FALSE)
  }
  if (!is.null(prepared$zdata) && ncol(as.matrix(prepared$Y)) != nrow(as.data.frame(prepared$zdata))) {
    stop("ncol(prepared$Y) must equal nrow(prepared$zdata) when zdata is supplied.", call. = FALSE)
  }
  invisible(TRUE)
}

# ------------------------------------------------------------
# Run all real-world methods
# ------------------------------------------------------------

run_real_world_single_benchmark <- function(
    Y,
    data = NULL,
    xformula = NULL,
    zformula = NULL,
    zdata = data,
    K = 10,
    methods = c("svdx", "svdx_adjusted", "svd", "pca", "nmf", "mixomics_pls", "spca", "sparse_nmf", "mefisto_duplicate_view"),
    inference_methods = c("svdx", "svdx_adjusted", "feature_lm"),
    run_factorization = TRUE,
    run_inference = TRUE,
    train_fraction = 0.70,
    train_cols = NULL,
    test_cols = NULL,
    seed = 1,
    alpha = 0.1,
    family = NULL,
    fit_family = NULL,
    inference_family = "gaussian",
    modelType = NULL,
    graphData = NULL,
    graphDistFn = NULL,
    rowOffsets = NULL,
    colOffsets = NULL,
    standardize = TRUE,
    adjust_for_z = TRUE,
    svdx_adjusted_use_z = FALSE,
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
    mefisto_covariates = NULL,
    z_adjustment = c("fitted", "reestimate"),
    numThreads = 4,
    svdx_args = list(),
    ...
) {
  inference_scope <- match.arg(inference_scope)
  loading_scope <- .benchmark_normalize_loading_scope(loading_scope)
  z_adjustment <- match.arg(z_adjustment)
  do_inference_enrichment <- isTRUE(inference_enrichment)

  dots <- list(...)
  if (!is.list(svdx_args)) {
    stop("svdx_args must be a list.", call. = FALSE)
  }
  if (length(svdx_args) > 0L &&
      (is.null(names(svdx_args)) || any(names(svdx_args) == ""))) {
    stop("All elements of svdx_args must be named.", call. = FALSE)
  }
  if (length(dots) > 0L &&
      (is.null(names(dots)) || any(names(dots) == ""))) {
    stop("All arguments supplied through '...' must be named.", call. = FALSE)
  }
  duplicated_svdx_args <- intersect(names(svdx_args), names(dots))
  if (length(duplicated_svdx_args) > 0L) {
    stop(
      "Arguments cannot be supplied through both svdx_args and '...': ",
      paste(duplicated_svdx_args, collapse = ", "),
      call. = FALSE
    )
  }

  protected_svdx_args <- c(
    "Y", "K", "xformula", "formula", "zformula", "data", "zdata",
    "numThreads", "modelType", "family", "rowOffsets", "colOffsets",
    "graphData", "graphDistFn"
  )
  conflicting_svdx_args <- intersect(
    c(names(svdx_args), names(dots)),
    protected_svdx_args
  )
  if (length(conflicting_svdx_args) > 0L) {
    stop(
      "svdx_args/'...' cannot override real-world wrapper-managed arguments: ",
      paste(unique(conflicting_svdx_args), collapse = ", "),
      call. = FALSE
    )
  }

  forwarded_svdx_args <- c(dots, svdx_args)
  if (!is.list(target_set_auc_args)) {
    stop("target_set_auc_args must be a list.", call. = FALSE)
  }

  # Support prepared single-dataset workflow objects, mirroring
  # runSVDXWorkflow()/runSVDXReport() semantics.  A prepared object is
  # expected to be a list containing at least Y/xdata/xformula and, when
  # available, zdata/zformula/family/fdata/modelType.
  if (is.list(Y) && !is.matrix(Y) && !is.data.frame(Y) && !is.null(Y$Y) && !is.null(Y$xdata)) {
    prepared <- Y
    .benchmark_validate_prepared_single(prepared)
    Y <- prepared$Y
    if (is.null(data)) data <- prepared$xdata
    if (is.null(xformula)) xformula <- prepared$xformula
    if (is.null(zdata)) zdata <- prepared$zdata
    if (is.null(zformula)) zformula <- prepared$zformula
    if (is.null(family)) family <- prepared$family
    if (is.null(fdata)) fdata <- prepared$fdata
    if (is.null(enrichment_id_map)) enrichment_id_map <- prepared$enrichment_id_map
    if (is.null(target_sets)) target_sets <- prepared$target_sets
    if (is.null(feature_sets)) feature_sets <- prepared$feature_sets %||% prepared$target_sets
    if (is.null(modelType)) modelType <- prepared$modelType
    if (is.null(rowOffsets)) rowOffsets <- prepared$rowOffsets
    if (is.null(colOffsets)) colOffsets <- prepared$colOffsets
  }

  do_loading_enrichment <- if (is.null(loading_enrichment)) {
    !is.null(feature_sets) || !is.null(enrichment_term2gene) || !is.null(msigdb_collection)
  } else {
    isTRUE(loading_enrichment)
  }

  if (is.null(data)) stop("data must be supplied or available as prepared$xdata.", call. = FALSE)
  if (is.null(xformula)) stop("xformula must be supplied or available as prepared$xformula.", call. = FALSE)

  Y <- as.matrix(Y)
  data <- as.data.frame(data)
  if (!is.null(zdata)) zdata <- as.data.frame(zdata)

  if (ncol(Y) != nrow(data)) {
    stop("ncol(Y) must equal nrow(data).", call. = FALSE)
  }
  if (!is.null(zformula) && !is.null(zdata) && ncol(Y) != nrow(zdata)) {
    stop("ncol(Y) must equal nrow(zdata) when zformula is supplied.", call. = FALSE)
  }
  if (is.null(colnames(Y))) {
    warning("Y has no column names; downstream sample labels may be unavailable.", call. = FALSE)
  }
  if (is.null(rownames(Y))) {
    warning("Y has no row names; feature-level enrichment may be unavailable or less informative.", call. = FALSE)
  }

  if (.benchmark_family_has_count(family)) {
    if (any(!is.finite(Y)) || any(Y < 0) || any(Y != floor(Y))) {
      stop("Count-family real-world benchmarking requires a finite non-negative integer matrix Y.", call. = FALSE)
    }
  }

  if (is.null(rowOffsets)) rowOffsets <- rep(0, nrow(Y))
  if (is.null(colOffsets)) colOffsets <- rep(0, ncol(Y))
  if (!is.numeric(rowOffsets) || length(rowOffsets) != nrow(Y) || any(!is.finite(rowOffsets))) {
    stop("rowOffsets must be NULL or a finite numeric vector of length nrow(Y).", call. = FALSE)
  }
  if (!is.numeric(colOffsets) || length(colOffsets) != ncol(Y) || any(!is.finite(colOffsets))) {
    stop("colOffsets must be NULL or a finite numeric vector of length ncol(Y).", call. = FALSE)
  }

  .validate_real_world_graph_inputs(graphData, graphDistFn, ncol(Y))
  .validate_feature_sets_basic(feature_sets, rownames(Y))
  .benchmark_validate_fdata_matches_features(fdata, rownames(Y), context = "Y")
  if (is.null(enrichment_id_map)) {
    enrichment_id_map <- .benchmark_fdata_to_enrichment_id_map(fdata, required = FALSE)
  }

  enrichment_backend <- match.arg(enrichment_backend)
  enrichment_id_type <- match.arg(enrichment_id_type)

  if ("feature_lm" %in% methods) {
    warning(
      "feature_lm is an inference baseline and is excluded from real-world factorization benchmarking. ",
      "Use inference_methods = 'feature_lm' instead.",
      call. = FALSE
    )
    methods <- setdiff(methods, "feature_lm")
  }

  valid_methods <- c("svdx", "svdx_adjusted", "glm_pca", "svd", "pca", "nmf", "mixomics_pls", "spca", "sparse_nmf", "mefisto_duplicate_view")
  invalid_methods <- setdiff(methods, valid_methods)
  if (length(invalid_methods) > 0L) {
    stop(
      "Unknown real-world methods: ", paste(invalid_methods, collapse = ", "),
      ". Valid options are: ", paste(valid_methods, collapse = ", "),
      call. = FALSE
    )
  }

  valid_inference_methods <- c("svdx", "svdx_adjusted", "feature_lm")
  invalid_inference_methods <- setdiff(inference_methods, valid_inference_methods)
  if (length(invalid_inference_methods) > 0L) {
    stop(
      "Unknown inference_methods: ", paste(invalid_inference_methods, collapse = ", "),
      ". Valid options are: ", paste(valid_inference_methods, collapse = ", "),
      call. = FALSE
    )
  }

  svdx_model_type <- .benchmark_complete_model_type(
    modelType = modelType,
    formula = xformula,
    data = data
  )

  split <- .make_real_world_train_test_split(
    n_samples = ncol(Y),
    train_fraction = train_fraction,
    train_cols = train_cols,
    test_cols = test_cols,
    seed = seed,
    sample_names = colnames(Y)
  )
  has_test <- length(split$test) > 0L
  run_inference_effective <- isTRUE(run_inference) && has_test
  if (isTRUE(run_inference) && !has_test && isTRUE(verbose)) {
    message("Skipping real-world inference because the test split is empty (e.g., train_fraction = 1).")
  }

  count_family <- .benchmark_family_has_count(family)
  if (isTRUE(standardize) && count_family) {
    message("standardize=TRUE ignored for poisson/quasipoisson family; using raw count data for count-native methods.")
    standardize <- FALSE
  }

  if (standardize) {
    std <- standardize_features(Y)
    Y_fit <- std$Y
  } else {
    std <- NULL
    Y_fit <- Y
  }

  # Count-native methods must receive raw counts. Gaussian comparison
  # methods receive a transformed matrix, matching the simulation chain.
  Y_count <- Y
  Y_train <- Y_fit[, split$train, drop = FALSE]
  Y_test <- Y_fit[, split$test, drop = FALSE]
  Y_count_train <- Y_count[, split$train, drop = FALSE]
  Y_count_test <- Y_count[, split$test, drop = FALSE]
  data_train <- data[split$train, , drop = FALSE]
  data_test <- data[split$test, , drop = FALSE]
  zdata_train <- if (!is.null(zdata)) zdata[split$train, , drop = FALSE] else NULL
  zdata_test <- if (!is.null(zdata)) zdata[split$test, , drop = FALSE] else NULL
  rowOffsets_train <- rowOffsets
  rowOffsets_test <- rowOffsets
  colOffsets_train <- colOffsets[split$train]
  colOffsets_test <- colOffsets[split$test]

  Y_baseline_input <- if (count_family) log1p(Y_count_train) else Y_train
  Y_method_train <- if (adjust_for_z && !is.null(zformula)) {
    adjust_Y_by_formula(Y_baseline_input, zdata_train %||% data_train, zformula)
  } else {
    Y_baseline_input
  }
  Y_baseline_nuisance <- Y_baseline_input - Y_method_train

  Y_method_test <- NULL
  if (has_test) {
    Y_baseline_test_input <- if (count_family) log1p(Y_count_test) else Y_test
    Y_method_test <- if (adjust_for_z && !is.null(zformula)) {
      adjust_Y_by_formula(Y_baseline_test_input, zdata_test %||% data_test, zformula)
    } else {
      Y_baseline_test_input
    }
  }

  raw_fits <- list()

  if (run_factorization) {
    if ("svdx" %in% methods) {
      raw_fits$svdx <- do.call(
        fit_svdx_wrapper,
        c(
          list(
            Y = if (count_family) Y_count_train else Y_train,
            K = K,
            xformula = xformula,
            zformula = zformula,
            data = data_train,
            zdata = zdata_train,
            numThreads = numThreads,
            family = fit_family,
            modelType = svdx_model_type,
            graphData = graphData,
            graphDistFn = graphDistFn,
            rowOffsets = rowOffsets_train,
            colOffsets = colOffsets_train
          ),
          forwarded_svdx_args
        )
      )
    }

    if ("svdx_adjusted" %in% methods) {
      if (count_family) {
        stop("svdx_adjusted is currently defined only for Gaussian real-world runs.", call. = FALSE)
      }
      adjusted_zformula <- if (isTRUE(svdx_adjusted_use_z)) zformula else ~1
      raw_fits$svdx_adjusted <- do.call(
        fit_svdx_wrapper,
        c(
          list(
            Y = Y_method_train,
            K = K,
            xformula = xformula,
            zformula = adjusted_zformula,
            data = data_train,
            zdata = zdata_train,
            numThreads = numThreads,
            family = fit_family,
            modelType = svdx_model_type,
            graphData = graphData,
            graphDistFn = graphDistFn,
            rowOffsets = rowOffsets_train,
            colOffsets = colOffsets_train
          ),
          forwarded_svdx_args
        )
      )
      if (!is.null(raw_fits$svdx_adjusted)) raw_fits$svdx_adjusted$method <- "svdx_adjusted"
    }

    if ("glm_pca" %in% methods) {
      if (!count_family) {
        stop("glm_pca is available only for poisson/quasipoisson real-world runs.", call. = FALSE)
      }
      glm_x <- NULL
      if (!is.null(zformula)) {
        glm_x <- stats::model.matrix(zformula, data = zdata_train %||% data_train)
        if ("(Intercept)" %in% colnames(glm_x)) {
          glm_x <- glm_x[, colnames(glm_x) != "(Intercept)", drop = FALSE]
        }
        if (ncol(glm_x) == 0L) glm_x <- NULL
      }
      raw_fits$glm_pca <- fit_glmpca_wrapper(
        Y = Y_count_train,
        K = K,
        X = glm_x,
        colOffsets = colOffsets_train,
        seed = seed
      )
    }

    if ("svd" %in% methods) raw_fits$svd <- fit_svd(Y_method_train, K = K)
    if ("pca" %in% methods) raw_fits$pca <- fit_pca(Y_method_train, K = K)
    if ("nmf" %in% methods) raw_fits$nmf <- fit_nmf(Y_method_train, K = K, seed = seed)

    if ("mixomics_pls" %in% methods) {
      raw_fits$mixomics_pls <- fit_mixomics_pls_formula(
        Y = Y_method_train,
        meta = data_train,
        K = K,
        xformula = xformula
      )
    }

    if ("spca" %in% methods) {
      if (is.null(spca_keepX_grid)) {
        spca_keepX_grid <- round(c(0.01, 0.10, 0.25, 0.50, 0.75) * nrow(Y_method_train))
      }
      raw_fits <- c(
        raw_fits,
        fit_mixomics_spca_grid(Y = Y_method_train, K = K, keepX_grid = spca_keepX_grid)
      )
    }

    if ("sparse_nmf" %in% methods) {
      raw_fits <- c(
        raw_fits,
        fit_sparse_nmf_grid(
          Y = Y_method_train,
          K = K,
          lambda_grid = sparse_nmf_lambda_grid,
          seed = seed
        )
      )
    }

    if ("mefisto_duplicate_view" %in% methods) {
      raw_fits$mefisto_duplicate_view <- fit_mefisto_duplicate_view_wrapper(
        Y = Y_method_train,
        data = data_train,
        K = K,
        covariates = mefisto_covariates,
        xformula = xformula,
        seed = seed
      )
    }
  }

  bench_results <- raw_fits |>
    lapply(as_benchmark_from_fitter) |>
    Filter(Negate(is.null), x = _)

  feature_scores <- extract_method_feature_scores(bench_results)

  loading_enrichment <- NULL
  if (isTRUE(do_loading_enrichment)) {
    loading_enrichment <- evaluate_real_world_loading_enrichment(
      bench_results = bench_results,
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

  inference <- list()

  if (run_inference_effective) {
    if ("svdx" %in% inference_methods) {
      svdx_fit <- if (!is.null(raw_fits$svdx) && !is.null(raw_fits$svdx$fit)) raw_fits$svdx$fit else NULL
      time_fit <- 0

      if (is.null(svdx_fit)) {
        tm_fit <- system.time({
          svdx_fit <- do.call(
            svdx,
            c(
              list(
                Y = if (count_family) Y_count_train else Y_train,
                formula = xformula,
                data = data_train,
                zformula = zformula,
                zdata = zdata_train,
                K = K,
                numThreads = numThreads,
                family = fit_family,
                modelType = svdx_model_type,
                graphData = graphData,
                graphDistFn = graphDistFn,
                rowOffsets = rowOffsets_train,
                colOffsets = colOffsets_train
              ),
              forwarded_svdx_args
            )
          )
        })
        time_fit <- unname(tm_fit["elapsed"])
      }

      tm_inf <- system.time({
        inf <- inferV(
          fit = svdx_fit,
          Ytest = Y_test,
          newdata = data_test,
          newzdata = zdata_test %||% data_test,
          zformula = zformula,
          family = inference_family,
          rowOffsets = rowOffsets_test,
          colOffsets = colOffsets_test,
          alpha = alpha,
          z_adjustment = z_adjustment
        )
      })

      inference_enrichment <- NULL
      if (do_inference_enrichment) {
        ep_args <- c(
          list(
            inf = inf,
            feature_sets = feature_sets,
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
            enrichment_id_map = enrichment_id_map,
            enrichment_extra_args = enrichment_extra_args,
            min_feature_overlap_warn = min_feature_overlap_warn,
            min_feature_overlap_fraction = min_feature_overlap_fraction
          )
        )
        inference_enrichment <- do.call(enrichPathways, ep_args)
      }

      inference$svdx <- list(
        fit = svdx_fit,
        inf = inf,
        inf_summary = summary(inf, alpha = alpha),
        inference_enrichment = inference_enrichment,
        timing = tibble::tibble(
          method = "svdx",
          time_fit = time_fit,
          time_inference = unname(tm_inf["elapsed"]),
          time_total = time_fit + unname(tm_inf["elapsed"])
        )
      )
    }

    if ("svdx_adjusted" %in% inference_methods) {
      if (count_family) {
        stop("svdx_adjusted inference is currently defined only for Gaussian real-world runs.", call. = FALSE)
      }
      adjusted_zformula <- if (isTRUE(svdx_adjusted_use_z)) zformula else ~1
      svdx_adjusted_fit <- if (!is.null(raw_fits$svdx_adjusted) && !is.null(raw_fits$svdx_adjusted$fit)) raw_fits$svdx_adjusted$fit else NULL
      time_fit <- 0
      if (is.null(svdx_adjusted_fit)) {
        tm_fit <- system.time({
          svdx_adjusted_fit <- do.call(
            svdx,
            c(
              list(
                Y = Y_method_train,
                formula = xformula,
                data = data_train,
                zformula = adjusted_zformula,
                zdata = zdata_train,
                K = K,
                numThreads = numThreads,
                family = fit_family,
                modelType = svdx_model_type,
                graphData = graphData,
                graphDistFn = graphDistFn,
                rowOffsets = rowOffsets_train,
                colOffsets = colOffsets_train
              ),
              forwarded_svdx_args
            )
          )
        })
        time_fit <- unname(tm_fit["elapsed"])
      }
      tm_inf <- system.time({
        inf <- inferV(
          fit = svdx_adjusted_fit,
          Ytest = Y_method_test,
          newdata = data_test,
          newzdata = zdata_test %||% data_test,
          zformula = adjusted_zformula,
          family = inference_family,
          rowOffsets = rowOffsets_test,
          colOffsets = colOffsets_test,
          alpha = alpha,
          z_adjustment = z_adjustment
        )
      })
      inference_enrichment <- NULL
      if (do_inference_enrichment) {
        ep_args <- c(list(
          inf = inf, feature_sets = feature_sets, alpha = alpha, scope = inference_scope,
          universe = enrichment_universe, enrichment_backend = enrichment_backend,
          msigdb_db_species = msigdb_db_species, msigdb_species = msigdb_species,
          msigdb_collection = msigdb_collection, msigdb_subcollection = msigdb_subcollection,
          enrichment_id_type = enrichment_id_type, enrichment_score = enrichment_score,
          term2gene = enrichment_term2gene, enrichment_id_map = enrichment_id_map,
          enrichment_extra_args = enrichment_extra_args,
          min_feature_overlap_warn = min_feature_overlap_warn,
          min_feature_overlap_fraction = min_feature_overlap_fraction
        ))
        inference_enrichment <- do.call(enrichPathways, ep_args)
      }
      inference$svdx_adjusted <- list(
        fit = svdx_adjusted_fit,
        inf = inf,
        inf_summary = summary(inf, alpha = alpha),
        inference_enrichment = inference_enrichment,
        timing = tibble::tibble(
          method = "svdx_adjusted",
          time_fit = time_fit,
          time_inference = unname(tm_inf["elapsed"]),
          time_total = time_fit + unname(tm_inf["elapsed"])
        )
      )
    }

    if ("feature_lm" %in% inference_methods) {
      lm_fit <- fit_feature_lm_formula(
        Y = Y_method_train,
        meta = data_train,
        xformula = xformula
      )
      inference$feature_lm <- list(
        fit = lm_fit,
        coef_mat = lm_fit$coef_mat %||% NULL,
        p_mat = lm_fit$p_mat %||% NULL,
        inference_enrichment = if (do_inference_enrichment) {
          .enrich_feature_lm_from_pvalues(lm_fit$p_mat, feature_sets, alpha = alpha)
        } else NULL,
        timing = tibble::tibble(
          method = "feature_lm",
          time_fit = lm_fit$time %||% NA_real_,
          time_inference = NA_real_,
          time_total = lm_fit$time %||% NA_real_
        )
      )
    }
  }

  prediction_metrics <- NULL
  svdx_pred_fit <- NULL
  if (has_test && !is.null(raw_fits$svdx) && !is.null(raw_fits$svdx$fit)) {
    svdx_pred_fit <- raw_fits$svdx$fit
  } else if (has_test && !is.null(inference$svdx) && !is.null(inference$svdx$fit)) {
    svdx_pred_fit <- inference$svdx$fit
  }
  if (has_test && !is.null(svdx_pred_fit)) {
    prediction_metrics <- tryCatch(
      .benchmark_svdx_single_prediction_metrics(
        fit = svdx_pred_fit,
        Ytrain = if (count_family) Y_count_train else Y_train,
        Ytest = if (count_family) Y_count_test else Y_test,
        xdata_test = data_test,
        zdata_test = zdata_test %||% data_test,
        zformula = zformula,
        family = family,
        rowOffsets = rowOffsets_test,
        colOffsets = colOffsets_test,
        method = "svdx",
        dataset = NA_character_,
        numThreads = numThreads
      ),
      error = function(e) {
        warning("Could not compute SVDX test prediction metrics: ", e$message, call. = FALSE)
        NULL
      }
    )
  }

  out <- list(
    Y = if (count_family) Y_count_train else Y_train,
    Y_method = Y_method_train,
    Y_baseline_nuisance = Y_baseline_nuisance,
    data = data_train,
    zdata = zdata_train,
    fdata = fdata,
    feature_sets = feature_sets,
    enrichment_id_map = enrichment_id_map,
    raw_fits = raw_fits,
    benchmark_results = bench_results,
    feature_scores = feature_scores,
    loading_enrichment = loading_enrichment,
    inference = inference,
    split = split,
    standardization = std,
    xformula = xformula,
    zformula = zformula,
    family = family,
    fit_family = fit_family,
    inference_family = inference_family,
    alpha = alpha,
    rowOffsets = rowOffsets_train,
    colOffsets = colOffsets_train,
    test_rowOffsets = rowOffsets_test,
    test_colOffsets = colOffsets_test,
    svdx_args = svdx_args,
    modelType = svdx_model_type,
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

  .validate_real_world_method_benchmark_output(out)
  class(out) <- "realWorldMethodBenchmark"

  if (!is.null(target_sets)) {
    out$target_auc <- do.call(
      evaluate_target_set_auc,
      c(list(rw_bench = out, target_sets = target_sets), target_set_auc_args)
    )
  }

  out
}

run_real_world_method_benchmark <- function(...) {
  run_real_world_single_benchmark(...)
}


# ------------------------------------------------------------
# SVDX-specific inference layer for real-world benchmark
# ------------------------------------------------------------

run_real_world_svdx_inference <- function(
    rw_bench,
    Ytest = NULL,
    newdata = NULL,
    newzdata = newdata,
    zformula = rw_bench$zformula,
    family = NULL,
    alpha = 0.1,
    feature_sets = NULL,
    fdata = rw_bench$fdata %||% NULL,
    enrichment_id_map = NULL,
    inference_enrichment = !is.null(feature_sets) || !is.null(fdata) || !is.null(enrichment_id_map),
    inference_scope = c("all", "component", "union", "global"),
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
    z_adjustment = c("fitted", "reestimate"),
    ...
) {
  inference_scope <- match.arg(inference_scope)
  z_adjustment <- match.arg(z_adjustment)

  .validate_real_world_method_benchmark(rw_bench)
  fit_family <- family
  infer_family <- fit_family %||% "gaussian"
  infer_family <- match.arg(infer_family, choices = c("gaussian", "poisson", "quasipoisson"))
  .validate_feature_sets_basic(feature_sets, rownames(rw_bench$Y))
  .benchmark_validate_fdata_matches_features(fdata, rownames(rw_bench$Y), context = "rw_bench$Y")
  if (is.null(enrichment_id_map)) {
    enrichment_id_map <- .benchmark_fdata_to_enrichment_id_map(fdata, required = FALSE)
  }
  enrichment_backend <- match.arg(enrichment_backend)
  enrichment_id_type <- match.arg(enrichment_id_type)

  svdx_raw <- rw_bench$raw_fits$svdx
  if (is.null(svdx_raw)) {
    stop("rw_bench does not contain an svdx fit.", call. = FALSE)
  }

  fit <- svdx_raw$fit
  if (is.null(fit)) {
    stop("svdx raw fit does not contain a fit object.", call. = FALSE)
  }

  if (is.null(Ytest)) {
    Ytest <- rw_bench$Y
  }

  if (is.null(newdata)) {
    newdata <- rw_bench$data
  }

  if (is.null(newzdata)) {
    newzdata <- newdata
  }

  inf <- inferV(
    fit = fit,
    Ytest = Ytest,
    newdata = newdata,
    newzdata = newzdata,
    zformula = zformula,
    family = infer_family,
    alpha = alpha,
    z_adjustment = z_adjustment,
    ...
  )

  enrichment <- NULL
  if (isTRUE(inference_enrichment)) {
    ep_args <- c(
      list(
        inf = inf,
        feature_sets = feature_sets,
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
        enrichment_id_map = enrichment_id_map,
        enrichment_extra_args = enrichment_extra_args,
        min_feature_overlap_warn = min_feature_overlap_warn,
        min_feature_overlap_fraction = min_feature_overlap_fraction
      )
    )
    enrichment <- do.call(enrichPathways, ep_args)
  }

  out <- list(
    fit = fit,
    inf = inf,
    inf_summary = summary(inf, alpha = alpha),
    enrichment = enrichment,
    alpha = alpha,
    family = family
  )

  class(out) <- "realWorldSVDXInference"
  out
}


# ------------------------------------------------------------
# U plotting for real-world benchmark methods
# ------------------------------------------------------------


.benchmark_plot_interaction_var <- function(data, vars, name = ".benchmark_plot_group") {
  if (is.null(vars)) return(list(data = data, var = NULL))
  vars <- as.character(vars)
  missing <- setdiff(vars, colnames(data))
  if (length(missing) > 0L) {
    stop(
      "Plotting variable(s) not found in data: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(vars) == 1L) return(list(data = data, var = vars))

  out_name <- name
  while (out_name %in% colnames(data)) out_name <- paste0(out_name, "_")
  data[[out_name]] <- do.call(
    interaction,
    c(unname(data[vars]), list(drop = TRUE, sep = " - "))
  )
  list(data = data, var = out_name)
}

plot_real_world_U <- function(
    rw_bench,
    method,
    component = 1,
    x,
    color = NULL,
    facet = NULL,
    connect = FALSE,
    line_group = NULL,
    line_alpha = 0.55,
    line_width = 0.4,
    point_alpha = 0.7
) {
  res <- .rwbench_get_method_result(rw_bench, method)
  if (is.null(res$U)) {
    stop("Selected method does not contain a U matrix.", call. = FALSE)
  }
  if (!x %in% colnames(rw_bench$data)) {
    stop("x column not found in rw_bench$data: ", x, call. = FALSE)
  }
  if (!is.null(color)) {
    missing_color <- setdiff(as.character(color), colnames(rw_bench$data))
    if (length(missing_color) > 0L) {
      stop("color column(s) not found in rw_bench$data: ",
           paste(missing_color, collapse = ", "), call. = FALSE)
    }
  }
  if (!is.null(facet) && !facet %in% colnames(rw_bench$data)) {
    stop("facet column not found in rw_bench$data: ", facet, call. = FALSE)
  }
  if (!is.null(line_group)) {
    missing_line_group <- setdiff(as.character(line_group), colnames(rw_bench$data))
    if (length(missing_line_group) > 0L) {
      stop("line_group column(s) not found in rw_bench$data: ",
           paste(missing_line_group, collapse = ", "), call. = FALSE)
    }
  }

  U <- as.data.frame(res$U)
  colnames(U) <- colnames(U) %||% paste0("component_", seq_len(ncol(U)))

  k <- if (is.numeric(component)) component else match(component, colnames(U))
  if (is.na(k) || k < 1 || k > ncol(U)) {
    stop("Invalid component.", call. = FALSE)
  }

  df <- cbind(rw_bench$data, component_value = U[[k]])

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
    if (!is.null(line_var)) {
      keep_groups <- names(which(vapply(split(df_line[[x]], df_line[[line_var]]), function(xx) length(unique(xx)) > 1L, logical(1))))
      df_line <- df_line[df_line[[line_var]] %in% keep_groups, , drop = FALSE]
    }
    if (nrow(df_line) > 1L && (is.null(line_var) || length(unique(df_line[[line_var]])) > 0L)) {
      p <- p + ggplot2::geom_line(
        data = df_line,
        mapping = line_aes,
        alpha = line_alpha,
        linewidth = line_width
      )
    }
  }

  p <- p +
    ggplot2::geom_point(alpha = point_alpha) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = x,
      y = paste(method, "component", k),
      color = color_label,
      title = paste(method, "sample scores")
    )

  if (!is.null(facet)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet)))
  }

  p
}


plot_real_world_U_trajectory <- function(
    rw_bench,
    method,
    component = 1,
    x,
    group = NULL,
    color = group,
    facet = NULL,
    line_alpha = 0.55,
    line_width = 0.4,
    point_alpha = 0.7
) {
  plot_real_world_U(
    rw_bench = rw_bench,
    method = method,
    component = component,
    x = x,
    color = color,
    facet = facet,
    connect = !is.null(group),
    line_group = group,
    line_alpha = line_alpha,
    line_width = line_width,
    point_alpha = point_alpha
  )
}

plot_real_world_U_embedding <- function(
    rw_bench,
    method,
    components = c(1, 2),
    color = NULL,
    facet = NULL,
    point_alpha = 0.7
) {
  res <- .rwbench_get_method_result(rw_bench, method)
  if (is.null(res$U)) {
    stop("Selected method does not contain a U matrix.", call. = FALSE)
  }

  U <- as.data.frame(res$U)
  colnames(U) <- colnames(U) %||% paste0("component_", seq_len(ncol(U)))

  k1 <- if (is.numeric(components[1])) components[1] else match(components[1], colnames(U))
  k2 <- if (is.numeric(components[2])) components[2] else match(components[2], colnames(U))

  if (any(is.na(c(k1, k2))) || any(c(k1, k2) < 1) || any(c(k1, k2) > ncol(U))) {
    stop("Invalid components.", call. = FALSE)
  }

  c1 <- colnames(U)[k1]
  c2 <- colnames(U)[k2]

  df <- cbind(rw_bench$data, U)

  aes_args <- ggplot2::aes(x = .data[[c1]], y = .data[[c2]])
  if (!is.null(color)) {
    aes_args$colour <- rlang::expr(.data[[!!color]])
  }

  p <- ggplot2::ggplot(df, aes_args) +
    ggplot2::geom_point(alpha = point_alpha) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = paste(method, c1),
      y = paste(method, c2),
      color = color,
      title = paste(method, "component embedding")
    )

  if (!is.null(facet)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet)))
  }

  p
}

# ------------------------------------------------------------
# V heatmap for any benchmark method
# ------------------------------------------------------------

plot_real_world_V_heatmap <- function(
    rw_bench,
    method,
    Y = rw_bench$Y,
    data = rw_bench$data,
    n_features = 50,
    features = NULL,
    cols = NULL,
    order_by = NULL,
    annotate = NULL,
    scale_rows = TRUE,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_feature_names = TRUE,
    show_sample_names = FALSE,
    row_labels = NULL,
    fdata = rw_bench$fdata %||% NULL
) {
  .require_namespace("ComplexHeatmap")

  res <- .rwbench_get_method_result(rw_bench, method)
  if (is.null(res$V)) {
    stop("Selected method does not contain a V matrix.", call. = FALSE)
  }
  if (is.null(res$feature_score)) {
    stop("Selected method does not contain feature scores.", call. = FALSE)
  }

  feature_names <- rownames(Y) %||% paste0("feature_", seq_len(nrow(Y)))

  if (is.null(features)) {
    idx <- order(res$feature_score, decreasing = TRUE, na.last = NA)
    features <- feature_names[head(idx, n_features)]
  }

  resolved_features <- .rwbench_resolve_features_and_labels(
    features = features,
    feature_names = feature_names,
    fdata = fdata,
    row_labels = row_labels
  )
  features <- resolved_features$feature
  feature_labels <- resolved_features$label

  fidx <- match(features, feature_names)
  if (any(is.na(fidx))) {
    stop("Some requested features were not found in rownames(Y).", call. = FALSE)
  }

  prepare_heatmap_matrix <- .benchmark_get_svdx_internal(".prepare_heatmap_matrix")
  prep <- prepare_heatmap_matrix(
    Y = Y,
    data = data,
    features = fidx,
    cols = cols,
    order_by = order_by,
    scale_rows = scale_rows
  )

  mat <- prep$Y
  rownames(mat) <- feature_labels

  if (!is.null(order_by)) {
    missing_order <- setdiff(order_by, colnames(data))
    if (length(missing_order) > 0L) {
      stop("order_by columns not found in data: ", paste(missing_order, collapse = ", "), call. = FALSE)
    }
  }

  ha <- NULL
  if (!is.null(annotate)) {
    missing_annotate <- setdiff(annotate, colnames(prep$data))
    if (length(missing_annotate) > 0L) {
      stop("annotate columns not found in data: ", paste(missing_annotate, collapse = ", "), call. = FALSE)
    }
    ha <- ComplexHeatmap::HeatmapAnnotation(
      df = prep$data[, annotate, drop = FALSE]
    )
  }

  ComplexHeatmap::Heatmap(
    mat,
    name = ifelse(scale_rows, "z", "Y"),
    top_annotation = ha,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    show_row_names = show_feature_names,
    show_column_names = show_sample_names,
    column_title = paste(method, "top feature heatmap")
  )
}


# ------------------------------------------
# SVDX specifics
# ------------------------------------------
plot_real_world_svdx_U <- function(
    rw_bench,
    component = 1,
    include = NULL,
    x,
    y = NULL,
    group = NULL,
    facet = NULL,
    fixed = list(),
    fill = c("error", "typical"),
    n_grid = 100,
    observed = FALSE,
    point_alpha = 0.6
) {
  svdx_raw <- rw_bench$raw_fits$svdx
  if (is.null(svdx_raw)) {
    stop("rw_bench does not contain an svdx fit.", call. = FALSE)
  }
  if (is.null(svdx_raw$fit)) {
    stop("svdx raw fit does not contain a fit object.", call. = FALSE)
  }

  plotU_fn <- .benchmark_get_svdx_internal("plotU")
  plotU_fn(
    fit = svdx_raw$fit,
    data = rw_bench$data,
    component = component,
    include = include,
    x = x,
    y = y,
    group = group,
    facet = facet,
    fixed = fixed,
    fill = fill,
    n_grid = n_grid,
    observed = observed,
    point_alpha = point_alpha
  )
}

plot_real_world_svdx_U_terms <- function(
    rw_bench,
    component = 1,
    terms = NULL,
    x,
    y = NULL,
    group = NULL,
    facet = NULL,
    fixed = list(),
    fill = c("error", "typical"),
    include_intercept = FALSE,
    n_grid = 100
) {
  svdx_raw <- rw_bench$raw_fits$svdx
  if (is.null(svdx_raw)) {
    stop("rw_bench does not contain an svdx fit.", call. = FALSE)
  }
  if (is.null(svdx_raw$fit)) {
    stop("svdx raw fit does not contain a fit object.", call. = FALSE)
  }

  plotUTerms_fn <- .benchmark_get_svdx_internal("plotUTerms")
  plotUTerms_fn(
    fit = svdx_raw$fit,
    data = rw_bench$data,
    component = component,
    terms = terms,
    x = x,
    y = y,
    group = group,
    facet = facet,
    fixed = fixed,
    fill = fill,
    include_intercept = include_intercept,
    n_grid = n_grid
  )
}


plot_real_world_svdx_U_embedding <- function(
    rw_bench,
    components = c(1, 2),
    include = NULL,
    vary = NULL,
    color = NULL,
    group = NULL,
    facet = NULL,
    fixed = list(),
    fill = c("error", "typical"),
    n_grid = 100,
    observed = is.null(vary),
    point_alpha = 0.7
) {
  svdx_raw <- rw_bench$raw_fits$svdx
  if (is.null(svdx_raw)) {
    stop("rw_bench does not contain an svdx fit.", call. = FALSE)
  }
  if (is.null(svdx_raw$fit)) {
    stop("svdx raw fit does not contain a fit object.", call. = FALSE)
  }

  plotUEmbedding_fn <- .benchmark_get_svdx_internal("plotUEmbedding")
  plotUEmbedding_fn(
    fit = svdx_raw$fit,
    data = rw_bench$data,
    components = components,
    include = include,
    vary = vary,
    color = color,
    group = group,
    facet = facet,
    fixed = fixed,
    fill = fill,
    n_grid = n_grid,
    observed = observed,
    point_alpha = point_alpha
  )
}
