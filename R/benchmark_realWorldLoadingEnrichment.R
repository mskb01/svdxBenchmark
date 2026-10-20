# benchmark_realWorldLoadingEnrichment.R
#
# Real-world loading-based enrichment helpers.
#
# These helpers provide method-comparable loading enrichment for real-world
# benchmark outputs. They reuse the SVDX enrichment backend by constructing a
# minimal svdxInference-compatible proxy from benchmark loading scores.

.benchmark_normalize_loading_scope <- function(loading_scope = c("all", "component", "union", "global")) {
  loading_scope <- match.arg(loading_scope, several.ok = TRUE)
  loading_scope[loading_scope == "global"] <- "union"
  if ("all" %in% loading_scope) {
    loading_scope <- c("component", "union")
  }
  unique(loading_scope)
}

.benchmark_loading_component_names <- function(V) {
  V <- as.matrix(V)
  cn <- colnames(V)
  if (!is.null(cn) && length(cn) == ncol(V) && all(nzchar(cn))) {
    return(as.character(cn))
  }
  paste0("component_", seq_len(ncol(V)))
}

.benchmark_top_loading_mask <- function(score, top_n) {
  score <- as.numeric(score)
  top_n <- min(as.integer(top_n), length(score))
  out <- rep(FALSE, length(score))
  if (!length(score) || is.na(top_n) || top_n <= 0L) return(out)
  ord <- order(score, decreasing = TRUE, na.last = NA)
  if (!length(ord)) return(out)
  out[ord[seq_len(min(top_n, length(ord)))]] <- TRUE
  out
}

.benchmark_make_loading_enrichment_proxy <- function(V, loading_top_n = 100L, alpha = 0.5) {
  V <- as.matrix(V)
  if (!nrow(V) || !ncol(V)) return(NULL)

  feature_names <- rownames(V) %||% paste0("feature_", seq_len(nrow(V)))
  component_names <- .benchmark_loading_component_names(V)
  absV <- abs(V)

  feature_results <- dplyr::bind_rows(lapply(seq_len(ncol(V)), function(k) {
    score <- absV[, k]
    selected <- .benchmark_top_loading_mask(score, loading_top_n)
    tibble::tibble(
      feature = feature_names,
      component = component_names[[k]],
      coef = as.numeric(V[, k]),
      se = NA_real_,
      stat = as.numeric(V[, k]),
      pval = ifelse(selected, 0, 1),
      qval = ifelse(selected, 0, 1)
    )
  }))

  union_score <- apply(absV, 1L, max, na.rm = TRUE)
  union_score[!is.finite(union_score)] <- NA_real_
  best_idx <- apply(absV, 1L, function(x) {
    if (all(is.na(x))) return(NA_integer_)
    which.max(x)
  })
  best_component <- rep(NA_character_, length(best_idx))
  ok <- !is.na(best_idx)
  best_component[ok] <- component_names[best_idx[ok]]
  best_stat <- rep(NA_real_, length(best_idx))
  best_stat[ok] <- V[cbind(which(ok), best_idx[ok])]
  selected_union <- .benchmark_top_loading_mask(union_score, loading_top_n)

  union_results <- tibble::tibble(
    feature = feature_names,
    best_component = best_component,
    best_coef = best_stat,
    best_stat = best_stat,
    max_abs_stat = as.numeric(union_score),
    min_p = ifelse(selected_union, 0, 1),
    union_p = ifelse(selected_union, 0, 1),
    union_q = ifelse(selected_union, 0, 1),
    global_stat = as.numeric(union_score),
    global_p = ifelse(selected_union, 0, 1),
    global_q = ifelse(selected_union, 0, 1),
    significant = selected_union,
    global_significant = selected_union
  )

  out <- list(
    feature_results = feature_results,
    union_results = union_results,
    alpha = alpha,
    family = NULL,
    call = match.call()
  )
  class(out) <- c("svdxInference", "list")
  out
}

.benchmark_enrichment_id_map_for_dataset <- function(enrichment_id_map, dataset = NULL) {
  if (is.null(enrichment_id_map)) return(NULL)
  if (!is.null(dataset) && is.list(enrichment_id_map) && !is.data.frame(enrichment_id_map)) {
    return(enrichment_id_map[[dataset]] %||% NULL)
  }
  enrichment_id_map
}

.benchmark_feature_sets_for_dataset <- function(feature_sets, dataset = NULL) {
  if (is.null(feature_sets)) return(NULL)
  if (!is.null(dataset) && is.list(feature_sets) && !is.null(feature_sets[[dataset]]) && is.list(feature_sets[[dataset]])) {
    return(feature_sets[[dataset]])
  }
  feature_sets
}

.benchmark_enrichment_universe_for_dataset <- function(universe, dataset = NULL) {
  if (is.null(universe)) return(NULL)
  if (!is.null(dataset) && is.list(universe) && !is.data.frame(universe)) {
    return(universe[[dataset]] %||% NULL)
  }
  universe
}

.benchmark_call_loading_enrichment <- function(
    proxy_inf,
    loading_scope = c("all", "component", "union", "global"),
    method,
    dataset = NULL,
    feature_sets = NULL,
    enrichment_id_map = NULL,
    alpha = 0.5,
    enrichment_backend = c("auto", "internal", "clusterProfiler", "fgsea"),
    msigdb_db_species = NULL,
    msigdb_species = NULL,
    msigdb_collection = NULL,
    msigdb_subcollection = NULL,
    enrichment_id_type = c("symbol", "entrez", "ensembl"),
    enrichment_score = "auto",
    term2gene = NULL,
    universe = NULL,
    min_feature_overlap_warn = 5L,
    min_feature_overlap_fraction = 0.01,
    enrichment_extra_args = list(),
    ...
) {
  if (is.null(proxy_inf)) return(NULL)
  loading_scope <- .benchmark_normalize_loading_scope(loading_scope)
  enrichment_backend <- match.arg(enrichment_backend)
  enrichment_id_type <- match.arg(enrichment_id_type)
  score_arg <- if (identical(enrichment_score, "auto")) "max_abs_stat" else enrichment_score

  calls <- list()

  if ("union" %in% loading_scope) {
    calls[["union"]] <- tryCatch(
      enrichPathways(
        proxy_inf,
        feature_sets = feature_sets,
        alpha = alpha,
        scope = "union",
        component = NULL,
        universe = universe,
        enrichment_backend = enrichment_backend,
        msigdb_db_species = msigdb_db_species,
        msigdb_species = msigdb_species,
        msigdb_collection = msigdb_collection,
        msigdb_subcollection = msigdb_subcollection,
        enrichment_id_type = enrichment_id_type,
        enrichment_score = score_arg,
        term2gene = term2gene,
        enrichment_id_map = enrichment_id_map,
        enrichment_extra_args = enrichment_extra_args,
        min_feature_overlap_warn = min_feature_overlap_warn,
        min_feature_overlap_fraction = min_feature_overlap_fraction,
        ...
      ),
      error = function(e) {
        message("Skipping loading enrichment for method '", method, "' scope 'union': ", conditionMessage(e))
        NULL
      }
    )
  }

  if ("component" %in% loading_scope) {
    comps <- unique(as.character(proxy_inf$feature_results$component))
    comps <- comps[!is.na(comps)]
    for (cc in comps) {
      calls[[paste0("component_", cc)]] <- tryCatch(
        enrichPathways(
          proxy_inf,
          feature_sets = feature_sets,
          alpha = alpha,
          scope = "component",
          component = cc,
          universe = universe,
          enrichment_backend = enrichment_backend,
          msigdb_db_species = msigdb_db_species,
          msigdb_species = msigdb_species,
          msigdb_collection = msigdb_collection,
          msigdb_subcollection = msigdb_subcollection,
          enrichment_id_type = enrichment_id_type,
          enrichment_score = score_arg,
          term2gene = term2gene,
          enrichment_id_map = enrichment_id_map,
          enrichment_extra_args = enrichment_extra_args,
          min_feature_overlap_warn = min_feature_overlap_warn,
          min_feature_overlap_fraction = min_feature_overlap_fraction,
          ...
        ),
        error = function(e) {
          message("Skipping loading enrichment for method '", method, "' component '", cc, "': ", conditionMessage(e))
          NULL
        }
      )
    }
  }

  calls <- calls[!vapply(calls, is.null, logical(1))]
  if (!length(calls)) return(NULL)

  results <- dplyr::bind_rows(lapply(calls, function(x) x$results %||% NULL))
  if (!nrow(results)) return(NULL)
  results$method <- method
  if (!is.null(dataset)) results$dataset <- dataset

  raw <- lapply(calls, function(x) x$raw %||% NULL)
  out <- list(
    backend_requested = unique(unlist(lapply(calls, function(x) x$backend_requested %||% NA_character_))),
    backend_resolved = unique(unlist(lapply(calls, function(x) x$backend_resolved %||% NA_character_))),
    backend = unique(unlist(lapply(calls, function(x) x$backend %||% x$backend_resolved))),
    params = list(
      loading_scope = loading_scope,
      method = method,
      dataset = dataset,
      enrichment_id_map = enrichment_id_map,
      loading_enrichment = TRUE
    ),
    results = results,
    raw = raw
  )
  class(out) <- c("svdxEnrichment", "list")
  out
}

evaluate_real_world_loading_enrichment <- function(
    bench_results,
    loading_scope = c("all", "component", "union", "global"),
    loading_top_n = 100,
    feature_sets = NULL,
    enrichment_id_map = NULL,
    alpha = 0.5,
    enrichment_backend = c("auto", "internal", "clusterProfiler", "fgsea"),
    msigdb_db_species = NULL,
    msigdb_species = NULL,
    msigdb_collection = NULL,
    msigdb_subcollection = NULL,
    enrichment_id_type = c("symbol", "entrez", "ensembl"),
    enrichment_score = "auto",
    term2gene = NULL,
    universe = NULL,
    min_feature_overlap_warn = 5L,
    min_feature_overlap_fraction = 0.01,
    enrichment_extra_args = list(),
    ...
) {
  loading_scope <- .benchmark_normalize_loading_scope(loading_scope)
  enrichment_backend <- match.arg(enrichment_backend)
  enrichment_id_type <- match.arg(enrichment_id_type)

  objects <- lapply(bench_results, function(res) {
    if (is.null(res$V)) return(NULL)
    proxy <- .benchmark_make_loading_enrichment_proxy(
      V = res$V,
      loading_top_n = loading_top_n,
      alpha = alpha
    )
    .benchmark_call_loading_enrichment(
      proxy_inf = proxy,
      loading_scope = loading_scope,
      method = res$method %||% NA_character_,
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
      term2gene = term2gene,
      universe = universe,
      min_feature_overlap_warn = min_feature_overlap_warn,
      min_feature_overlap_fraction = min_feature_overlap_fraction,
      enrichment_extra_args = enrichment_extra_args,
      ...
    )
  })
  names(objects) <- vapply(bench_results, function(x) x$method %||% NA_character_, character(1))
  objects <- objects[!vapply(objects, is.null, logical(1))]
  if (!length(objects)) return(NULL)

  out <- list(
    backend_requested = unique(unlist(lapply(objects, function(x) x$backend_requested %||% NA_character_))),
    backend_resolved = unique(unlist(lapply(objects, function(x) x$backend_resolved %||% NA_character_))),
    backend = unique(unlist(lapply(objects, function(x) x$backend %||% x$backend_resolved))),
    params = list(loading_scope = loading_scope, loading_top_n = loading_top_n),
    results = dplyr::bind_rows(lapply(objects, function(x) x$results %||% NULL)),
    raw = lapply(objects, function(x) x$raw %||% NULL)
  )
  names(out$raw) <- names(objects)
  class(out) <- c("svdxEnrichment", "list")
  out
}

evaluate_integrated_real_world_loading_enrichment <- function(
    integrated_results,
    loading_scope = c("all", "component", "union", "global"),
    loading_top_n = 100,
    feature_sets = NULL,
    enrichment_id_map = NULL,
    alpha = 0.5,
    enrichment_backend = c("auto", "internal", "clusterProfiler", "fgsea"),
    msigdb_db_species = NULL,
    msigdb_species = NULL,
    msigdb_collection = NULL,
    msigdb_subcollection = NULL,
    enrichment_id_type = c("symbol", "entrez", "ensembl"),
    enrichment_score = "auto",
    term2gene = NULL,
    universe = NULL,
    min_feature_overlap_warn = 5L,
    min_feature_overlap_fraction = 0.01,
    enrichment_extra_args = list(),
    ...
) {
  loading_scope <- .benchmark_normalize_loading_scope(loading_scope)
  enrichment_backend <- match.arg(enrichment_backend)
  enrichment_id_type <- match.arg(enrichment_id_type)

  objects <- list()
  for (res in integrated_results) {
    if (is.null(res$V_list) || !length(res$V_list)) next
    for (d in seq_along(res$V_list)) {
      dataset <- names(res$V_list)[d] %||% paste0("dataset_", d)
      V <- res$V_list[[d]]
      if (is.null(V)) next
      fs_ds <- .benchmark_feature_sets_for_dataset(feature_sets, dataset = dataset)
      map_ds <- .benchmark_enrichment_id_map_for_dataset(enrichment_id_map, dataset = dataset)
      univ_ds <- .benchmark_enrichment_universe_for_dataset(universe, dataset = dataset)
      proxy <- .benchmark_make_loading_enrichment_proxy(
        V = V,
        loading_top_n = loading_top_n,
        alpha = alpha
      )
      obj <- .benchmark_call_loading_enrichment(
        proxy_inf = proxy,
        loading_scope = loading_scope,
        method = res$method %||% NA_character_,
        dataset = dataset,
        feature_sets = fs_ds,
        enrichment_id_map = map_ds,
        alpha = alpha,
        enrichment_backend = enrichment_backend,
        msigdb_db_species = msigdb_db_species,
        msigdb_species = msigdb_species,
        msigdb_collection = msigdb_collection,
        msigdb_subcollection = msigdb_subcollection,
        enrichment_id_type = enrichment_id_type,
        enrichment_score = enrichment_score,
        term2gene = term2gene,
        universe = univ_ds,
        min_feature_overlap_warn = min_feature_overlap_warn,
        min_feature_overlap_fraction = min_feature_overlap_fraction,
        enrichment_extra_args = enrichment_extra_args,
        ...
      )
      if (!is.null(obj)) objects[[paste(res$method %||% "method", dataset, sep = "__")]] <- obj
    }
  }

  if (!length(objects)) return(NULL)
  out <- list(
    backend_requested = unique(unlist(lapply(objects, function(x) x$backend_requested %||% NA_character_))),
    backend_resolved = unique(unlist(lapply(objects, function(x) x$backend_resolved %||% NA_character_))),
    backend = unique(unlist(lapply(objects, function(x) x$backend %||% x$backend_resolved))),
    params = list(loading_scope = loading_scope, loading_top_n = loading_top_n),
    results = dplyr::bind_rows(lapply(objects, function(x) x$results %||% NULL)),
    raw = lapply(objects, function(x) x$raw %||% NULL)
  )
  class(out) <- c("svdxEnrichment", "list")
  out
}
