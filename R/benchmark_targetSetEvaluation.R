# ============================================================
# Target-set recovery evaluation for real-world benchmarks
# ============================================================
# These helpers evaluate whether method loadings prioritize user-specified
# feature sets, for example circadian genes, interferon genes, bird guilds,
# or any other feature collection. Unlike ordinary enrichment, this returns
# AUCs based on loading-derived feature scores and supports component-rank
# weighting.
# ============================================================

`%||%` <- function(x, y) if (!is.null(x)) x else y

# ------------------------------------------------------------
# MSigDB target-set discovery and loading
# ------------------------------------------------------------

.benchmark_msigdb_call <- function(
    db_species,
    species,
    collection = NULL,
    subcollection = NULL
) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Package 'msigdbr' is required for MSigDB target-set helpers.", call. = FALSE)
  }

  args <- list(
    db_species = as.character(db_species),
    species = as.character(species)
  )
  if (!is.null(collection)) args$collection <- as.character(collection)
  if (!is.null(subcollection)) args$subcollection <- as.character(subcollection)

  do.call(msigdbr::msigdbr, args)
}

.benchmark_msig_col <- function(x, candidates) {
  hit <- intersect(candidates, names(x))
  if (length(hit) == 0L) return(NULL)
  hit[[1L]]
}

#' Search MSigDB gene sets by partial name/description
#'
#' @param query Character string such as "circadian".
#' @param db_species MSigDB database species, e.g. "MM" or "HS".
#' @param species Organism species, e.g. "Mus musculus" or "Homo sapiens".
#' @param collection Optional MSigDB collection, e.g. "M5".
#' @param subcollection Optional MSigDB subcollection, e.g. "BP".
#' @param ignore_case Logical; use case-insensitive matching.
#' @return Data frame of matching gene sets.
search_msigdb_sets <- function(
    query,
    db_species,
    species,
    collection = NULL,
    subcollection = NULL,
    ignore_case = TRUE
) {
  if (missing(query) || length(query) != 1L || !nzchar(as.character(query))) {
    stop("query must be a non-empty character string.", call. = FALSE)
  }

  msig <- .benchmark_msigdb_call(
    db_species = db_species,
    species = species,
    collection = collection,
    subcollection = subcollection
  )

  name_col <- .benchmark_msig_col(msig, c("gs_name", "gs_id"))
  desc_col <- .benchmark_msig_col(msig, c("gs_description", "gs_exact_source"))
  coll_col <- .benchmark_msig_col(msig, c("gs_collection", "gs_cat", "collection"))
  sub_col  <- .benchmark_msig_col(msig, c("gs_subcollection", "gs_subcat", "subcollection"))

  if (is.null(name_col)) {
    stop("MSigDB result does not contain a gene-set name column.", call. = FALSE)
  }

  haystack <- as.character(msig[[name_col]])
  if (!is.null(desc_col)) {
    haystack <- paste(haystack, as.character(msig[[desc_col]] %||% ""))
  }

  hit <- grepl(as.character(query), haystack, ignore.case = ignore_case)
  out <- msig[hit, , drop = FALSE]
  if (nrow(out) == 0L) {
    return(data.frame(
      gs_name = character(),
      gs_collection = character(),
      gs_subcollection = character(),
      gs_description = character(),
      n_genes = integer(),
      stringsAsFactors = FALSE
    ))
  }

  gene_col <- .benchmark_msig_col(out, c("gene_symbol", "ncbi_gene", "ensembl_gene"))
  split_key <- as.character(out[[name_col]])
  n_genes <- if (!is.null(gene_col)) {
    stats::aggregate(
      as.character(out[[gene_col]]),
      by = list(gs_name = split_key),
      FUN = function(z) length(unique(z[!is.na(z) & nzchar(z)]))
    )
  } else {
    stats::aggregate(
      seq_len(nrow(out)),
      by = list(gs_name = split_key),
      FUN = length
    )
  }
  names(n_genes)[2L] <- "n_genes"

  meta <- data.frame(
    gs_name = as.character(out[[name_col]]),
    gs_collection = if (!is.null(coll_col)) as.character(out[[coll_col]]) else NA_character_,
    gs_subcollection = if (!is.null(sub_col)) as.character(out[[sub_col]]) else NA_character_,
    gs_description = if (!is.null(desc_col)) as.character(out[[desc_col]]) else NA_character_,
    stringsAsFactors = FALSE
  )
  meta <- unique(meta)
  meta <- merge(meta, n_genes, by = "gs_name", all.x = TRUE, sort = FALSE)
  meta[order(meta$gs_collection, meta$gs_subcollection, meta$gs_name), , drop = FALSE]
}

#' Load selected MSigDB target sets as a named list
get_msigdb_target_sets <- function(
    set_names,
    db_species,
    species,
    collection = NULL,
    subcollection = NULL,
    id_type = c("symbol", "entrez", "ensembl"),
    ignore_case = FALSE
) {
  id_type <- match.arg(id_type)
  if (missing(set_names) || length(set_names) == 0L) {
    stop("set_names must contain at least one MSigDB gene-set name.", call. = FALSE)
  }

  msig <- .benchmark_msigdb_call(
    db_species = db_species,
    species = species,
    collection = collection,
    subcollection = subcollection
  )

  name_col <- .benchmark_msig_col(msig, c("gs_name", "gs_id"))
  if (is.null(name_col)) stop("MSigDB result does not contain a gene-set name column.", call. = FALSE)

  id_col <- switch(
    id_type,
    symbol = .benchmark_msig_col(msig, c("gene_symbol")),
    entrez = .benchmark_msig_col(msig, c("ncbi_gene", "entrez_gene")),
    ensembl = .benchmark_msig_col(msig, c("ensembl_gene", "ensembl"))
  )
  if (is.null(id_col)) {
    stop("MSigDB result does not contain the requested identifier column for id_type = '", id_type, "'.", call. = FALSE)
  }

  names_in <- as.character(msig[[name_col]])
  wanted <- as.character(set_names)
  if (isTRUE(ignore_case)) {
    keep <- tolower(names_in) %in% tolower(wanted)
  } else {
    keep <- names_in %in% wanted
  }

  msig <- msig[keep, , drop = FALSE]
  if (nrow(msig) == 0L) {
    stop("None of the requested set_names were found in MSigDB result.", call. = FALSE)
  }

  sets <- split(as.character(msig[[id_col]]), as.character(msig[[name_col]]))
  sets <- lapply(sets, function(z) unique(z[!is.na(z) & nzchar(z)]))
  sets[names(sets) %in% unique(as.character(msig[[name_col]]))]
}

# ------------------------------------------------------------
# Loading score engine
# ------------------------------------------------------------

.benchmark_target_rank_weights <- function(K, a = 1, scheme = c("power", "exponential")) {
  scheme <- match.arg(scheme)
  K <- as.integer(K)
  if (!is.finite(K) || K < 1L) stop("K must be a positive integer.", call. = FALSE)
  a <- as.numeric(a)
  if (length(a) != 1L || !is.finite(a) || a < 0) stop("a must be a finite non-negative scalar.", call. = FALSE)

  k <- seq_len(K)
  w <- switch(
    scheme,
    power = 1 / (k ^ a),
    exponential = exp(-a * (k - 1))
  )
  w / sum(w)
}

.benchmark_target_loading_score <- function(V, a = 1, scheme = c("power", "exponential")) {
  if (is.null(V)) stop("V must not be NULL.", call. = FALSE)
  V <- as.matrix(V)
  if (nrow(V) < 1L || ncol(V) < 1L) stop("V must have at least one row and one column.", call. = FALSE)
  scheme <- match.arg(scheme)

  w <- .benchmark_target_rank_weights(ncol(V), a = a, scheme = scheme)
  score <- as.vector(abs(V) %*% w)
  names(score) <- rownames(V) %||% paste0("feature_", seq_len(nrow(V)))
  score
}

.benchmark_target_safe_auc <- function(score, truth) {
  if (exists(".bench_safe_auc_score", mode = "function", inherits = TRUE)) {
    return(.bench_safe_auc_score(score, truth))
  }

  truth <- as.logical(truth)
  score <- as.numeric(score)
  ok <- is.finite(score) & !is.na(truth)
  if (sum(ok) == 0 || length(unique(truth[ok])) < 2L) return(NA_real_)

  if (requireNamespace("pROC", quietly = TRUE)) {
    roc_obj <- pROC::roc(response = truth[ok], predictor = score[ok], quiet = TRUE, direction = "<")
    return(as.numeric(pROC::auc(roc_obj)))
  }

  r <- rank(score[ok], ties.method = "average")
  y <- truth[ok]
  n1 <- sum(y)
  n0 <- sum(!y)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

.benchmark_target_normalize_target_sets <- function(target_sets) {
  if (is.null(target_sets) || length(target_sets) == 0L) {
    stop("target_sets must be a non-empty named list of feature sets.", call. = FALSE)
  }
  if (!is.list(target_sets) || is.data.frame(target_sets)) {
    stop("target_sets must be a named list.", call. = FALSE)
  }
  if (is.null(names(target_sets)) || any(!nzchar(names(target_sets)))) {
    names(target_sets) <- paste0("target_set_", seq_along(target_sets))
  }
  lapply(target_sets, function(z) unique(as.character(unlist(z, use.names = FALSE))))
}

.benchmark_target_expand_sets <- function(target_sets, combine_sets = c("none", "union", "both"), union_name = "union") {
  combine_sets <- match.arg(combine_sets)
  target_sets <- .benchmark_target_normalize_target_sets(target_sets)

  if (identical(combine_sets, "none")) {
    set_mode <- rep("individual", length(target_sets))
    names(set_mode) <- names(target_sets)
    return(list(sets = target_sets, set_mode = set_mode))
  }

  union_set <- unique(as.character(unlist(target_sets, use.names = FALSE)))
  union_set <- union_set[!is.na(union_set) & nzchar(union_set)]
  union_list <- list(union_set)
  names(union_list) <- union_name

  if (identical(combine_sets, "union")) {
    set_mode <- "union"
    names(set_mode) <- names(union_list)
    return(list(sets = union_list, set_mode = set_mode))
  }

  sets <- c(target_sets, union_list)
  set_mode <- c(rep("individual", length(target_sets)), "union")
  names(set_mode) <- names(sets)
  list(sets = sets, set_mode = set_mode)
}


.benchmark_target_enrichment_id_map <- function(fdata = NULL, dataset = NULL, enrichment_id_map = NULL) {
  if (!is.null(enrichment_id_map)) return(enrichment_id_map)
  if (is.null(fdata)) return(NULL)
  fd <- if (exists(".benchmark_resolve_fdata_for_dataset", mode = "function", inherits = TRUE)) {
    .benchmark_resolve_fdata_for_dataset(fdata, dataset = dataset, required = FALSE)
  } else {
    fdata
  }
  if (is.null(fd)) return(NULL)

  if (exists(".benchmark_fdata_to_enrichment_id_map", mode = "function", inherits = TRUE)) {
    return(.benchmark_fdata_to_enrichment_id_map(fd, required = FALSE))
  }

  NULL
}

.benchmark_target_standardize_id_map <- function(id_map, feature_ids) {
  if (is.null(id_map)) return(NULL)
  id_map <- as.data.frame(id_map, stringsAsFactors = FALSE)
  if (nrow(id_map) == 0L) return(NULL)

  feature_col <- intersect(c("feature", "feature_id", "id", "rowname"), names(id_map))
  id_col <- intersect(c("enrichment_id", "gene", "gene_id", "symbol", "entrez", "ensembl"), names(id_map))

  if (length(feature_col) == 0L && !is.null(rownames(id_map))) {
    id_map$feature <- rownames(id_map)
    feature_col <- "feature"
  }
  if (length(feature_col) == 0L && ncol(id_map) >= 2L) feature_col <- names(id_map)[1L]
  if (length(id_col) == 0L && ncol(id_map) >= 2L) id_col <- setdiff(names(id_map), feature_col)[1L]

  if (length(feature_col) == 0L || length(id_col) == 0L) return(NULL)

  out <- data.frame(
    feature = as.character(id_map[[feature_col[1L]]]),
    enrichment_id = as.character(id_map[[id_col[1L]]]),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$feature) & nzchar(out$feature) & !is.na(out$enrichment_id) & nzchar(out$enrichment_id), , drop = FALSE]
  out <- out[out$feature %in% feature_ids, , drop = FALSE]
  out
}

.benchmark_target_membership <- function(
    feature_ids,
    target_set,
    fdata = NULL,
    dataset = NULL,
    enrichment_id_map = NULL,
    target_id_type = c("auto", "feature", "enrichment_id")
) {
  target_id_type <- match.arg(target_id_type)
  feature_ids <- as.character(feature_ids)
  target_set <- unique(as.character(target_set))
  target_set <- target_set[!is.na(target_set) & nzchar(target_set)]

  feature_truth <- feature_ids %in% target_set

  map <- .benchmark_target_standardize_id_map(
    .benchmark_target_enrichment_id_map(fdata = fdata, dataset = dataset, enrichment_id_map = enrichment_id_map),
    feature_ids = feature_ids
  )

  enrichment_truth <- rep(FALSE, length(feature_ids))
  if (!is.null(map) && nrow(map) > 0L) {
    hit_features <- unique(map$feature[map$enrichment_id %in% target_set])
    enrichment_truth <- feature_ids %in% hit_features
  }

  if (identical(target_id_type, "feature")) return(feature_truth)
  if (identical(target_id_type, "enrichment_id")) return(enrichment_truth)

  # Auto: prefer the ID space with greater target overlap.
  if (sum(enrichment_truth, na.rm = TRUE) > sum(feature_truth, na.rm = TRUE)) enrichment_truth else feature_truth
}

.benchmark_target_get_single_result <- function(rw_bench, method) {
  if (!is.null(rw_bench$benchmark_results[[method]])) return(rw_bench$benchmark_results[[method]])

  idx <- which(vapply(rw_bench$benchmark_results, function(z) identical(z$method, method), logical(1)))
  if (length(idx) == 0L) stop("Method not found in rw_bench$benchmark_results: ", method, call. = FALSE)
  rw_bench$benchmark_results[[idx[1L]]]
}

.benchmark_target_get_integrated_result <- function(rw_int, method) {
  if (!is.null(rw_int$results[[method]])) return(rw_int$results[[method]])

  idx <- which(vapply(rw_int$results, function(z) identical(z$method, method), logical(1)))
  if (length(idx) == 0L) stop("Method not found in rw_int$results: ", method, call. = FALSE)
  rw_int$results[[idx[1L]]]
}

# ------------------------------------------------------------
# Real-world single-data target-set AUC
# ------------------------------------------------------------

#' Evaluate target-set recovery from loading scores in single real-world benchmarks
#'
#' @param rw_bench A `realWorldMethodBenchmark` object.
#' @param target_sets Named list of target feature/enrichment-ID sets.
#' @param methods Methods to evaluate; defaults to all benchmark_results.
#' @param a_grid Non-negative component-rank decay parameters.
#' @param scheme Weighting scheme: "power" uses 1/k^a, "exponential" uses exp(-a*(k-1)).
#' @param fdata Optional svdx_fdata; defaults to rw_bench$fdata.
#' @param combine_sets "none" evaluates each set separately, "union" evaluates the union only, and "both" evaluates individual sets plus their union.
#' @param union_name Name used for the union target set.
#' @param target_id_type "auto", "feature", or "enrichment_id".
#' @return Data frame of target-set AUC results.
evaluate_target_set_auc <- function(
    rw_bench,
    target_sets,
    methods = NULL,
    a_grid = c(0, 0.5, 1, 2, 4),
    scheme = c("power", "exponential"),
    fdata = NULL,
    combine_sets = c("none", "union", "both"),
    union_name = "union",
    target_id_type = c("auto", "feature", "enrichment_id")
) {
  if (is.null(rw_bench$benchmark_results)) {
    stop("rw_bench must contain benchmark_results.", call. = FALSE)
  }
  scheme <- match.arg(scheme)
  combine_sets <- match.arg(combine_sets)
  target_id_type <- match.arg(target_id_type)
  expanded_sets <- .benchmark_target_expand_sets(
    target_sets,
    combine_sets = combine_sets,
    union_name = union_name
  )
  target_sets <- expanded_sets$sets
  set_mode <- expanded_sets$set_mode
  fdata <- fdata %||% rw_bench$fdata
  enrichment_id_map <- rw_bench$enrichment_id_map %||% NULL
  methods <- methods %||% names(rw_bench$benchmark_results)

  out <- list()
  ii <- 0L
  for (method in methods) {
    res <- .benchmark_target_get_single_result(rw_bench, method)
    V <- res$V
    if (is.null(V)) next
    V <- as.matrix(V)
    feature_ids <- rownames(V) %||% paste0("feature_", seq_len(nrow(V)))

    for (a in a_grid) {
      score <- .benchmark_target_loading_score(V, a = a, scheme = scheme)
      for (set_name in names(target_sets)) {
        truth <- .benchmark_target_membership(
          feature_ids = feature_ids,
          target_set = target_sets[[set_name]],
          fdata = fdata,
          enrichment_id_map = enrichment_id_map,
          target_id_type = target_id_type
        )
        n_matched <- sum(truth, na.rm = TRUE)
        ii <- ii + 1L
        out[[ii]] <- data.frame(
          method = method,
          target_set = set_name,
          a = as.numeric(a),
          scheme = scheme,
          auc = .benchmark_target_safe_auc(score, truth),
          n_features = length(feature_ids),
          n_target_features = n_matched,
          target_set_size = n_matched,
          target_set_size_raw = length(target_sets[[set_name]]),
          target_set_size_matched = n_matched,
          set_mode = unname(set_mode[[set_name]] %||% "individual"),
          K_returned = ncol(V),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(out) == 0L) return(data.frame())
  do.call(rbind, out)
}

# Backward/readability alias for real-world single workflows.
evaluate_real_world_target_set_auc <- evaluate_target_set_auc

# ------------------------------------------------------------
# Real-world integrated target-set AUC
# ------------------------------------------------------------

evaluate_integrated_target_set_auc <- function(
    rw_int,
    target_sets,
    methods = NULL,
    datasets = NULL,
    a_grid = c(0, 0.5, 1, 2, 4),
    scheme = c("power", "exponential"),
    fdata = NULL,
    combine_sets = c("none", "union", "both"),
    union_name = "union",
    target_id_type = c("auto", "feature", "enrichment_id")
) {
  if (is.null(rw_int$results)) stop("rw_int must contain results.", call. = FALSE)
  scheme <- match.arg(scheme)
  combine_sets <- match.arg(combine_sets)
  target_id_type <- match.arg(target_id_type)
  fdata <- fdata %||% rw_int$fdata
  enrichment_id_map <- rw_int$enrichment_id_map %||% NULL
  methods <- methods %||% names(rw_int$results)
  datasets <- datasets %||% names(rw_int$Y)

  out <- list()
  ii <- 0L
  for (method in methods) {
    res <- .benchmark_target_get_integrated_result(rw_int, method)
    V_list <- res$V_list %||% res$V
    if (is.null(V_list)) next
    if (!is.list(V_list) || is.matrix(V_list) || is.data.frame(V_list)) {
      V_list <- list(dataset = V_list)
    }

    for (ds in intersect(datasets, names(V_list))) {
      V <- as.matrix(V_list[[ds]])
      feature_ids <- rownames(V) %||% paste0("feature_", seq_len(nrow(V)))

      ts_ds <- target_sets
      if (is.list(target_sets) && !is.null(target_sets[[ds]]) && is.list(target_sets[[ds]])) {
        ts_ds <- target_sets[[ds]]
      }
      expanded_sets <- .benchmark_target_expand_sets(
        ts_ds,
        combine_sets = combine_sets,
        union_name = union_name
      )
      ts_ds <- expanded_sets$sets
      set_mode <- expanded_sets$set_mode

      for (a in a_grid) {
        score <- .benchmark_target_loading_score(V, a = a, scheme = scheme)
        for (set_name in names(ts_ds)) {
          truth <- .benchmark_target_membership(
            feature_ids = feature_ids,
            target_set = ts_ds[[set_name]],
            fdata = fdata,
            dataset = ds,
            enrichment_id_map = if (is.list(enrichment_id_map) && !is.data.frame(enrichment_id_map)) enrichment_id_map[[ds]] else enrichment_id_map,
            target_id_type = target_id_type
          )
          n_matched <- sum(truth, na.rm = TRUE)
          ii <- ii + 1L
          out[[ii]] <- data.frame(
            dataset = ds,
            method = method,
            target_set = set_name,
            a = as.numeric(a),
            scheme = scheme,
            auc = .benchmark_target_safe_auc(score, truth),
            n_features = length(feature_ids),
            n_target_features = n_matched,
            target_set_size = n_matched,
            target_set_size_raw = length(ts_ds[[set_name]]),
            target_set_size_matched = n_matched,
            set_mode = unname(set_mode[[set_name]] %||% "individual"),
            K_returned = ncol(V),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (length(out) == 0L) return(data.frame())
  do.call(rbind, out)
}

evaluate_real_world_integrated_target_set_auc <- evaluate_integrated_target_set_auc

# ------------------------------------------------------------
# Plotting
# ------------------------------------------------------------

plot_target_set_auc_curve <- function(
    target_auc,
    facet_by = c("target_set", "dataset_target_set"),
    y_limits = c(0.5, 1)
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.", call. = FALSE)
  facet_by <- match.arg(facet_by)
  df <- as.data.frame(target_auc)
  if (!all(c("method", "target_set", "a", "auc") %in% names(df))) {
    stop("target_auc must contain method, target_set, a, and auc columns.", call. = FALSE)
  }
  if (identical(facet_by, "dataset_target_set") && "dataset" %in% names(df)) {
    df$facet_label <- paste(df$dataset, df$target_set, sep = ": ")
  } else {
    df$facet_label <- df$target_set
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$a, y = .data$auc, color = .data$method, group = .data$method)) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.8, na.rm = TRUE) +
    ggplot2::facet_wrap(~ facet_label) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Component rank penalty a",
      y = "Relevant subset prioritization score",
      color = "Method",
      title = "Relevant subset prioritization by component-rank penalty"
    )
  if (!is.null(y_limits)) p <- p + ggplot2::coord_cartesian(ylim = y_limits)
  p
}

plot_target_set_auc <- function(
    target_auc,
    a = 1,
    facet_by = c("target_set", "dataset_target_set"),
    y_limits = c(0.5, 1)
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.", call. = FALSE)
  facet_by <- match.arg(facet_by)
  df <- as.data.frame(target_auc)
  if (!all(c("method", "target_set", "a", "auc") %in% names(df))) {
    stop("target_auc must contain method, target_set, a, and auc columns.", call. = FALSE)
  }
  df <- df[abs(df$a - a) == min(abs(df$a - a), na.rm = TRUE), , drop = FALSE]
  if (identical(facet_by, "dataset_target_set") && "dataset" %in% names(df)) {
    df$facet_label <- paste(df$dataset, df$target_set, sep = ": ")
  } else {
    df$facet_label <- df$target_set
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$method, y = .data$auc, fill = .data$method)) +
    ggplot2::geom_col(na.rm = TRUE) +
    ggplot2::facet_wrap(~ facet_label) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = NULL,
      y = "Relevant subset prioritization score",
      fill = "Method",
      title = paste0("Relevant subset prioritization at a = ", a)
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  if (!is.null(y_limits)) p <- p + ggplot2::coord_cartesian(ylim = y_limits)
  p
}
