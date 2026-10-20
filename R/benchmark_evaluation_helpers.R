# ============================================================
# Benchmark evaluation helpers
# ============================================================
# Shared utilities for single-dataset and integrated benchmarks.
# These functions avoid assumptions about component/module matching.
# Feature recovery is always evaluated against the union of true
# simulated signal features.
# ============================================================

`%||%` <- function(x, y) if (!is.null(x)) x else y

.bench_feature_names <- function(Y) {
  rownames(Y) %||% paste0("feature_", seq_len(nrow(Y)))
}

.bench_resolve_feature_names <- function(Y, idx_or_names) {
  feature_names <- .bench_feature_names(Y)
  if (is.null(idx_or_names)) return(character())
  if (is.character(idx_or_names)) return(intersect(idx_or_names, feature_names))
  idx <- as.integer(idx_or_names)
  idx <- idx[is.finite(idx) & idx >= 1 & idx <= length(feature_names)]
  feature_names[idx]
}

.bench_feature_set_logical <- function(Y, features) {
  feature_names <- .bench_feature_names(Y)
  feature_set_names <- .bench_resolve_feature_names(Y, features)
  out <- feature_names %in% feature_set_names
  names(out) <- feature_names

  out
}

.bench_signal_logical <- function(Y, signal_features) {
  .bench_feature_set_logical(
    Y = Y,
    features = signal_features
  )
}

.bench_batch_logical <- function(Y, batch_features) {
  if (length(batch_features) == 0) {
    return(rep(FALSE, nrow(Y)))
  }
  if (is.list(batch_features)) {
    batch_features <- unique(unlist(batch_features, use.names = FALSE))
  }
  .bench_feature_set_logical(
    Y = Y,
    features = batch_features
  )
}

.bench_safe_auc_score <- function(score, truth) {
  truth <- as.logical(truth)
  score <- as.numeric(score)
  ok <- is.finite(score) & !is.na(truth)
  if (sum(ok) == 0 || length(unique(truth[ok])) < 2) return(NA_real_)

  if (requireNamespace("pROC", quietly = TRUE)) {
    roc_obj <- pROC::roc(
      response = truth[ok],
      predictor = score[ok],
      quiet = TRUE,
      direction = "<"
    )
    return(as.numeric(pROC::auc(roc_obj)))
  }

  r <- rank(score[ok], ties.method = "average")
  y <- truth[ok]
  n1 <- sum(y)
  n0 <- sum(!y)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

.bench_safe_bh <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p) & !is.na(p)
  out[ok] <- stats::p.adjust(p[ok], method = "BH")
  out
}

.bench_safe_cor_vec <- function(x, y, method = "pearson") {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  stats::cor(x[ok], y[ok], method = method)
}

.bench_safe_frob <- function(A, B) {
  if (is.null(A) || is.null(B) || !identical(dim(A), dim(B))) return(NA_real_)
  d <- A - B
  if (!all(is.finite(d))) return(NA_real_)
  sqrt(sum(d^2))
}

.bench_safe_rel_frob <- function(A, B, eps = 1e-8) {
  num <- sqrt(sum((A - B)^2, na.rm = TRUE))
  den <- sqrt(sum(B^2, na.rm = TRUE))
  if (!is.finite(den) || den < eps) return(NA_real_)
  num / den
}

.bench_safe_fisher_enrichment <- function(discovered, in_set) {
  discovered <- as.logical(discovered)
  in_set <- as.logical(in_set)
  ok <- !is.na(discovered) & !is.na(in_set)
  discovered <- discovered[ok]
  in_set <- in_set[ok]

  tab <- matrix(
    c(
      sum(discovered & in_set),
      sum(discovered & !in_set),
      sum(!discovered & in_set),
      sum(!discovered & !in_set)
    ),
    nrow = 2,
    byrow = TRUE
  )

  if (any(!is.finite(tab)) || sum(tab) == 0) {
    return(list(odds_ratio = NA_real_, p_value = NA_real_))
  }

  ft <- tryCatch(stats::fisher.test(tab, alternative = "greater"), error = function(e) NULL)
  if (is.null(ft)) {
    return(list(odds_ratio = NA_real_, p_value = NA_real_))
  }

  list(
    odds_ratio = unname(ft$estimate %||% NA_real_),
    p_value = unname(ft$p.value %||% NA_real_)
  )
}

.bench_add_enrichment_value_transforms <- function(df) {
  if (!is.data.frame(df)) return(df)
  if ("p_value" %in% names(df) && !"neg_log10_p_value" %in% names(df)) {
    df$neg_log10_p_value <- -log10(pmax(as.numeric(df$p_value), .Machine$double.xmin))
  }
  if ("q_value" %in% names(df) && !"neg_log10_q_value" %in% names(df)) {
    df$neg_log10_q_value <- -log10(pmax(as.numeric(df$q_value), .Machine$double.xmin))
  }
  df
}

.bench_module_feature_sets_single <- function(Y, truth) {
  feature_names <- .bench_feature_names(Y)
  module_features <- truth$module_features
  if (is.null(module_features)) return(list())

  out <- lapply(seq_along(module_features), function(m) {
    mf <- module_features[[m]]
    if (is.character(mf)) intersect(mf, feature_names) else feature_names[as.integer(mf)]
  })

  module_type <- truth$module_type %||% rep("module", length(out))
  names(out) <- paste0("module_", seq_along(out), "_", module_type)
  out
}

.bench_module_metadata_single <- function(Y, truth) {
  feature_names <- .bench_feature_names(Y)
  module_features <- truth$module_features
  module_type <- truth$module_type %||% rep(NA_character_, length(module_features))

  dplyr::bind_rows(lapply(seq_along(module_features), function(m) {
    mf <- module_features[[m]]
    module_names <- if (is.character(mf)) mf else feature_names[as.integer(mf)]
    module_names <- intersect(module_names, feature_names)
    tibble::tibble(
      module = paste0("module_", m),
      module_type = module_type[m] %||% NA_character_,
      feature = module_names,
      module_size = length(module_names)
    )
  }))
}

.bench_module_feature_sets_integrated <- function(Yd, truth, dataset) {
  feature_names <- .bench_feature_names(Yd)
  module_features <- truth$module_features[[dataset]]
  if (is.null(module_features)) return(list())

  out <- lapply(seq_along(module_features), function(m) {
    mf <- module_features[[m]]
    if (is.character(mf)) intersect(mf, feature_names) else feature_names[as.integer(mf)]
  })

  module_type <- truth$module_type %||% rep("module", length(out))
  names(out) <- paste0("module_", seq_along(out), "_", module_type)
  out
}

.bench_module_metadata_integrated <- function(Yd, truth, dataset) {
  feature_names <- .bench_feature_names(Yd)
  module_features <- truth$module_features[[dataset]]
  module_type <- truth$module_type %||% rep(NA_character_, length(module_features))

  dplyr::bind_rows(lapply(seq_along(module_features), function(m) {
    mf <- module_features[[m]]
    module_names <- if (is.character(mf)) mf else feature_names[as.integer(mf)]
    module_names <- intersect(module_names, feature_names)
    tibble::tibble(
      module = paste0("module_", m),
      module_type = module_type[m] %||% NA_character_,
      feature = module_names,
      module_size = length(module_names)
    )
  }))
}

.bench_feature_score_tbl <- function(method, Y, score, dataset = NULL) {
  feature_names <- .bench_feature_names(Y)
  score <- as.numeric(score)
  if (length(score) != length(feature_names)) {
    stop("feature score length must equal number of rows/features in Y.", call. = FALSE)
  }

  out <- tibble::tibble(
    method = method,
    feature = feature_names,
    feature_score = score
  )

  if (!is.null(dataset)) out <- dplyr::mutate(out, dataset = dataset, .after = method)
  out
}

# .bench_evaluate_feature_score <- function(method, Y, score, signal_features, dataset = NULL) {
#   signal <- .bench_signal_logical(Y, signal_features)
#   out <- tibble::tibble(
#     method = method,
#     loading_feature_auc_signal_vs_null = .bench_safe_auc_score(score, signal)
#   )
#   if (!is.null(dataset)) out <- dplyr::mutate(out, dataset = dataset, .after = method)
#   out
# }

.bench_evaluate_feature_score <- function(
  method,
  Y,
  score,
  signal_features,
  batch_features = NULL,
  batch_only_features = NULL,
  dataset = NULL
) {

  signal <- .bench_signal_logical(Y, signal_features)
  out <- tibble::tibble(
    method = method,
    loading_feature_auc_signal_vs_non_signal =
      .bench_safe_auc_score(score, signal)
  )
  if (!is.null(batch_features)) {
    batch <- .bench_batch_logical(Y, batch_features)
    out <- dplyr::mutate(
      out,
      loading_feature_auc_batch_vs_non_batch =
        .bench_safe_auc_score(score, batch)
    )
  }
  if (!is.null(batch_only_features)) {
    batch_only <- .bench_batch_logical(Y, batch_only_features)
    keep <- signal | batch_only
    out <- dplyr::mutate(
      out,
      loading_feature_auc_signal_vs_batch_only =
        if (sum(signal[keep]) < 1 || sum(batch_only[keep]) < 1) {
          NA_real_
        } else {
          .bench_safe_auc_score(score[keep], signal[keep])
        }
    )
  }

  if (!is.null(dataset)) {
    out <- dplyr::mutate(
      out,
      dataset = dataset,
      .after = method
    )
  }

  out
}

.bench_evaluate_module_auc <- function(method, Y, V, module_metadata, dataset = NULL) {
  feature_names <- .bench_feature_names(Y)
  V <- as.matrix(V)

  if (nrow(V) != length(feature_names)) {
    stop("V must have one row per feature in Y.", call. = FALSE)
  }

  rownames(V) <- feature_names
  K_use <- ncol(V)

  out <- dplyr::bind_rows(lapply(split(module_metadata, module_metadata$module), function(mm) {
    module <- unique(mm$module)
    module_type <- unique(mm$module_type)
    module_features <- intersect(mm$feature, feature_names)
    module_size <- length(module_features)

    auc_by_component <- vapply(seq_len(K_use), function(k) {
      truth <- feature_names %in% module_features
      score_k <- abs(V[, k])
      .bench_safe_auc_score(score_k, truth)
    }, numeric(1))

    best_component <- if (all(is.na(auc_by_component))) {
      NA_integer_
    } else {
      which.max(replace(auc_by_component, is.na(auc_by_component), -Inf))
    }

    tibble::tibble(
      method = method,
      module = module,
      module_type = module_type,
      module_size = module_size,
      best_component = best_component,
      auc_module_vs_background = if (is.na(best_component)) {
        NA_real_
      } else {
        auc_by_component[[best_component]]
      }
    )
  }))

  if (!is.null(dataset)) {
    out <- dplyr::mutate(out, dataset = dataset, .after = method)
  }

  out
}

.bench_evaluate_module_recall_curves <- function(
    method,
    Y,
    V,
    module_metadata,
    top_frac = seq(0.01, 0.50, by = 0.01),
    dataset = NULL
) {
  feature_names <- .bench_feature_names(Y)
  V <- as.matrix(V)

  if (nrow(V) != length(feature_names)) {
    stop("V must have one row per feature in Y.", call. = FALSE)
  }

  rownames(V) <- feature_names
  K_use <- ncol(V)
  p <- length(feature_names)

  out <- dplyr::bind_rows(
    lapply(split(module_metadata, module_metadata$module), function(mm) {
      module <- unique(mm$module)
      module_type <- unique(mm$module_type)
      module_features <- intersect(mm$feature, feature_names)
      module_size <- length(module_features)

      if (module_size < 1 || K_use < 1) {
        return(dplyr::bind_rows(lapply(top_frac, function(frac) {
          tibble::tibble(
            method = method,
            module = module,
            module_type = module_type,
            best_component = NA_integer_,
            top_frac = frac,
            n_selected = max(1L, ceiling(frac * p)),
            module_size = module_size,
            n_recovered = NA_integer_,
            recall = NA_real_
          )
        })))
      }

      ranked_by_component <- lapply(seq_len(K_use), function(k) {
        names(sort(abs(V[, k]), decreasing = TRUE, na.last = NA))
      })

      dplyr::bind_rows(lapply(top_frac, function(frac) {
        n_selected <- max(1L, ceiling(frac * p))

        recovered_by_component <- vapply(seq_len(K_use), function(k) {
          ranked_k <- ranked_by_component[[k]]
          selected <- ranked_k[seq_len(min(n_selected, length(ranked_k)))]
          sum(module_features %in% selected)
        }, numeric(1))

        recall_by_component <- recovered_by_component / module_size

        best_component <- if (all(is.na(recall_by_component))) {
          NA_integer_
        } else {
          which.max(replace(recall_by_component, is.na(recall_by_component), -Inf))
        }

        tibble::tibble(
          method = method,
          module = module,
          module_type = module_type,
          best_component = best_component,
          top_frac = frac,
          n_selected = n_selected,
          module_size = module_size,
          n_recovered = if (is.na(best_component)) {
            NA_integer_
          } else {
            as.integer(recovered_by_component[[best_component]])
          },
          recall = if (is.na(best_component)) {
            NA_real_
          } else {
            recall_by_component[[best_component]]
          }
        )
      }))
    })
  )

  if (!is.null(dataset)) {
    out <- dplyr::mutate(out, dataset = dataset, .after = method)
  }

  out
}

.bench_summarise_module_recall_curves <- function(recall_curves, recall_at = c(0.05, 0.10, 0.20)) {
  if (is.null(recall_curves) || nrow(recall_curves) == 0) return(tibble::tibble())

  group_cols <- intersect(c("method", "dataset", "module", "module_type"), names(recall_curves))

  module_summary <- recall_curves |>
    dplyr::arrange(dplyr::across(dplyr::all_of(group_cols)), top_frac) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      aurc = {
        x <- top_frac
        y <- recall
        ok <- is.finite(x) & is.finite(y)
        x <- x[ok]
        y <- y[ok]
        if (length(x) < 2) NA_real_ else {
          area <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
          area / (max(x) - min(x))
        }
      },
      final_recall = recall[which.max(top_frac)],
      .groups = "drop"
    )

  recall_at_tbl <- dplyr::bind_rows(lapply(recall_at, function(frac0) {
    recall_curves |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
      dplyr::slice_min(abs(top_frac - frac0), n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        recall_frac = paste0("recall_at_", gsub("\\.", "", as.character(frac0 * 100)), "pct"),
        recall_value = recall
      ) |>
      dplyr::select(dplyr::all_of(group_cols), recall_frac, recall_value)
  })) |>
    tidyr::pivot_wider(names_from = recall_frac, values_from = recall_value)

  module_summary |>
    dplyr::left_join(recall_at_tbl, by = group_cols)
}

.bench_evaluate_loading_enrichment <- function(
    method,
    Y,
    score,
    module_metadata,
    top_n = 100,
    dataset = NULL
) {
  feature_names <- .bench_feature_names(Y)
  score <- as.numeric(score)
  names(score) <- feature_names
  ranked <- names(sort(score, decreasing = TRUE, na.last = NA))
  selected <- ranked[seq_len(min(top_n, length(ranked)))]
  discovered <- feature_names %in% selected

  out <- module_metadata |>
    dplyr::group_by(module, module_type) |>
    dplyr::summarise(
      set_size = dplyr::n(),
      n_selected = length(selected),
      overlap = sum(discovered & (feature_names %in% feature)),
      odds_ratio = {
        sf <- .bench_safe_fisher_enrichment(discovered, feature_names %in% feature)
        sf$odds_ratio
      },
      p_value = {
        sf <- .bench_safe_fisher_enrichment(discovered, feature_names %in% feature)
        sf$p_value
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(
      q_value = .bench_safe_bh(p_value),
      method = method,
      .before = 1
    )

  out <- .bench_add_enrichment_value_transforms(out)
  if (!is.null(dataset)) out <- dplyr::mutate(out, dataset = dataset, .after = method)
  out
}

.bench_extract_union_inference <- function(
    inf,
    alpha = inf$alpha %||% 0.1,
    inference_scope = NULL
) {
  if (!is.null(inference_scope)) {
    inference_scope <- match.arg(inference_scope, c("global", "union"))
  }

  # Prefer a package-provided union_results table. If absent, construct
  # feature-level union p-values from component-wise p-values using
  # Bonferroni correction across tested components per feature:
  #   union_p = min(1, K_feature * min_k p_jk)
  # followed by BH over features.
  if (!is.null(inf$union_results)) {
    ur <- as.data.frame(inf$union_results)
    if (!"feature" %in% names(ur)) {
      stop("inf$union_results must contain a feature column.", call. = FALSE)
    }
    if (!"union_p" %in% names(ur)) {
      if ("min_p" %in% names(ur)) ur$union_p <- ur$min_p
      else if ("p_value" %in% names(ur)) ur$union_p <- ur$p_value
    }
    if (!"union_q" %in% names(ur) && "union_p" %in% names(ur)) {
      ur$union_q <- .bench_safe_bh(ur$union_p)
    }
    if ("global_p" %in% names(ur) && !"global_q" %in% names(ur)) {
      ur$global_q <- .bench_safe_bh(ur$global_p)
    }

    if (is.null(inference_scope)) {
      # Legacy behavior for callers that have not selected a scope explicitly.
      ur$score <- .bench_inference_score_from_union_tbl(ur)
      if ("global_q" %in% names(ur)) {
        ur$discovered <- ur$global_q < alpha
      } else {
        ur$discovered <- ur$union_q < alpha
      }
    } else {
      p_col <- paste0(inference_scope, "_p")
      q_col <- paste0(inference_scope, "_q")
      missing_cols <- setdiff(c(p_col, q_col), names(ur))
      if (length(missing_cols)) {
        stop(
          "Requested inference scope '", inference_scope,
          "' is unavailable: missing ", paste(missing_cols, collapse = ", "), ".",
          call. = FALSE
        )
      }

      p_scope <- suppressWarnings(as.numeric(ur[[p_col]]))
      q_scope <- suppressWarnings(as.numeric(ur[[q_col]]))
      if (!any(is.finite(p_scope) & !is.na(p_scope)) ||
          !any(is.finite(q_scope) & !is.na(q_scope))) {
        stop(
          "Requested inference scope '", inference_scope,
          "' is unavailable: p/q values are all missing or non-finite.",
          call. = FALSE
        )
      }

      ur$score <- -log10(pmax(p_scope, .Machine$double.xmin))
      ur$discovered <- !is.na(q_scope) & is.finite(q_scope) & q_scope <= alpha
      ur$scope <- inference_scope
      ur$component <- inference_scope
    }
    return(tibble::as_tibble(ur))
  }

  if (identical(inference_scope, "global")) {
    stop(
      "Requested inference scope 'global' is unavailable because inf$union_results is absent.",
      call. = FALSE
    )
  }

  fr <- inf$feature_results
  if (is.null(fr)) stop("inf must contain union_results or feature_results.", call. = FALSE)
  fr <- as.data.frame(fr)
  .bench_check_inference_feature_results(fr)

  p_col <- .bench_inference_p_col(fr)
  stat_col <- .bench_inference_stat_col(fr)

  ur <- fr |>
    dplyr::filter(!is.na(.data[[p_col]])) |>
    dplyr::group_by(feature) |>
    dplyr::summarise(
      K_tested = dplyr::n_distinct(component),
      min_p = min(.data[[p_col]], na.rm = TRUE),
      union_p = pmin(1, K_tested * min_p),
      max_abs_stat = if (!is.null(stat_col)) max(abs(.data[[stat_col]]), na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      union_q = .bench_safe_bh(union_p),
      score = ifelse(is.finite(max_abs_stat), max_abs_stat, -log10(pmax(union_p, .Machine$double.xmin))),
      discovered = union_q <= alpha
    )

  if (identical(inference_scope, "union")) {
    ur <- ur |>
      dplyr::mutate(
        score = -log10(pmax(union_p, .Machine$double.xmin)),
        discovered = !is.na(union_q) & is.finite(union_q) & union_q <= alpha,
        scope = "union",
        component = "union"
      )
  }

  ur
}

.bench_check_inference_feature_results <- function(fr) {
  if (!"feature" %in% names(fr)) stop("feature_results must contain feature.", call. = FALSE)
  if (!"component" %in% names(fr)) stop("feature_results must contain component.", call. = FALSE)
  invisible(TRUE)
}

.bench_inference_p_col <- function(fr) {
  candidates <- c("p_value", "p", "pval", "pval_raw")
  hit <- intersect(candidates, names(fr))
  if (length(hit) == 0) stop("Could not find p-value column in feature_results.", call. = FALSE)
  hit[1]
}

.bench_inference_stat_col <- function(fr) {
  candidates <- c("stat", "statistic", "z", "t", "wald_stat", "coef")
  hit <- intersect(candidates, names(fr))
  if (length(hit) == 0) return(NULL)
  hit[1]
}


.bench_inference_score_from_union_tbl <- function(union_tbl) {
  if (is.null(union_tbl) || !is.data.frame(union_tbl)) {
    return(numeric(0))
  }
  n <- nrow(union_tbl)
  score <- rep(NA_real_, n)

  if ("global_p" %in% names(union_tbl)) {
    gp <- suppressWarnings(as.numeric(union_tbl$global_p))
    ok <- is.finite(gp) & !is.na(gp)
    score[ok] <- -log10(pmax(gp[ok], .Machine$double.xmin))
  }

  if ("score" %in% names(union_tbl)) {
    sc <- suppressWarnings(as.numeric(union_tbl$score))
    ok <- (!is.finite(score) | is.na(score)) & is.finite(sc) & !is.na(sc)
    score[ok] <- sc[ok]
  }

  if ("union_p" %in% names(union_tbl)) {
    up <- suppressWarnings(as.numeric(union_tbl$union_p))
    ok <- (!is.finite(score) | is.na(score)) & is.finite(up) & !is.na(up)
    score[ok] <- -log10(pmax(up[ok], .Machine$double.xmin))
  }

  if ("max_abs_stat" %in% names(union_tbl)) {
    st <- suppressWarnings(abs(as.numeric(union_tbl$max_abs_stat)))
    ok <- (!is.finite(score) | is.na(score)) & is.finite(st) & !is.na(st)
    score[ok] <- st[ok]
  }

  score
}

.bench_inference_score_matrix <- function(Y, union_tbl) {
  feature_names <- .bench_feature_names(Y)
  u <- union_tbl |>
    dplyr::filter(.data$feature %in% feature_names)

  score <- rep(NA_real_, length(feature_names))
  names(score) <- feature_names
  if ("score" %in% names(u)) {
    score[u$feature] <- suppressWarnings(as.numeric(u$score))
  } else {
    score[u$feature] <- .bench_inference_score_from_union_tbl(u)
  }

  V <- matrix(score, ncol = 1L)
  rownames(V) <- feature_names
  colnames(V) <- "global"
  V
}

.bench_evaluate_inference_modules <- function(
    method,
    Y,
    union_tbl,
    module_metadata,
    top_frac = seq(0.01, 0.50, by = 0.01),
    dataset = NULL
) {
  V <- .bench_inference_score_matrix(Y, union_tbl)
  module_auc <- .bench_evaluate_module_auc(
    method = method,
    Y = Y,
    V = V,
    module_metadata = module_metadata,
    dataset = dataset
  )
  module_recall_curves <- .bench_evaluate_module_recall_curves(
    method = method,
    Y = Y,
    V = V,
    module_metadata = module_metadata,
    top_frac = top_frac,
    dataset = dataset
  )
  module_recall_summary <- .bench_summarise_module_recall_curves(module_recall_curves)

  list(
    module_auc = module_auc,
    module_recall_curves = module_recall_curves,
    module_recall_summary = module_recall_summary
  )
}

.bench_evaluate_union_inference <- function(
    method,
    Y,
    union_tbl,
    signal_features,
    alpha = 0.1,
    dataset = NULL
) {
  feature_names <- .bench_feature_names(Y)
  signal <- .bench_signal_logical(Y, signal_features)

  u <- union_tbl |>
    dplyr::filter(.data$feature %in% feature_names)

  score <- rep(NA_real_, length(feature_names)); names(score) <- feature_names
  discovered <- rep(FALSE, length(feature_names)); names(discovered) <- feature_names

  if ("score" %in% names(u)) {
    score[u$feature] <- suppressWarnings(as.numeric(u$score))
  } else {
    score[u$feature] <- .bench_inference_score_from_union_tbl(u)
  }
  if ("discovered" %in% names(u)) {
    discovered[u$feature] <- !is.na(u$discovered) & as.logical(u$discovered)
  } else if ("global_q" %in% names(u)) {
    discovered[u$feature] <- !is.na(u$global_q) & u$global_q <= alpha
  } else {
    discovered[u$feature] <- !is.na(u$union_q) & u$union_q <= alpha
  }

  n_disc <- sum(discovered, na.rm = TRUE)
  true_disc <- sum(discovered & signal, na.rm = TRUE)
  false_disc <- sum(discovered & !signal, na.rm = TRUE)

  out <- tibble::tibble(
    method = method,
    inference_auc_signal_vs_null = .bench_safe_auc_score(score, signal),
    n_discoveries = n_disc,
    true_discoveries = true_disc,
    false_discoveries = false_disc,
    power = true_disc / sum(signal, na.rm = TRUE),
    type1_error = false_disc / sum(!signal, na.rm = TRUE),
    empirical_fdr = ifelse(n_disc > 0, false_disc / n_disc, NA_real_),
    #min_union_p = min(u$union_p, na.rm = TRUE),
    #min_union_q = min(u$union_q, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE)
  )

  if (!is.null(dataset)) out <- dplyr::mutate(out, dataset = dataset, .after = method)
  out
}

.bench_union_feature_table <- function(method, Y, union_tbl, dataset = NULL) {
  feature_names <- .bench_feature_names(Y)
  out <- union_tbl |>
    dplyr::filter(.data$feature %in% feature_names) |>
    dplyr::mutate(method = method, .before = 1)
  if (!is.null(dataset)) out <- dplyr::mutate(out, dataset = dataset, .after = method)
  tibble::as_tibble(out)
}

.bench_evaluate_inference_enrichment <- function(
    method,
    Y,
    union_tbl,
    module_metadata,
    alpha = 0.1,
    dataset = NULL
) {
  feature_names <- .bench_feature_names(Y)
  u <- union_tbl |>
    dplyr::filter(.data$feature %in% feature_names)

  discovered <- rep(FALSE, length(feature_names)); names(discovered) <- feature_names
  if ("discovered" %in% names(u)) {
    discovered[u$feature] <- !is.na(u$discovered) & as.logical(u$discovered)
  } else {
    discovered[u$feature] <- !is.na(u$union_q) & u$union_q <= alpha
  }

  out <- module_metadata |>
    dplyr::group_by(module, module_type) |>
    dplyr::summarise(
      set_size = dplyr::n(),
      #n_significant = sum(discovered, na.rm = TRUE),
      n_discoveries = sum(discovered, na.rm = TRUE),
      overlap = sum(discovered & (feature_names %in% feature)),
      odds_ratio = {
        sf <- .bench_safe_fisher_enrichment(discovered, feature_names %in% feature)
        sf$odds_ratio
      },
      p_value = {
        sf <- .bench_safe_fisher_enrichment(discovered, feature_names %in% feature)
        sf$p_value
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(
      q_value = .bench_safe_bh(p_value),
      method = method,
      scope = if ("scope" %in% names(u)) as.character(u$scope[[1]]) else "union",
      component = if ("component" %in% names(u)) as.character(u$component[[1]]) else "union",
      .before = 1
    )

  out <- .bench_add_enrichment_value_transforms(out)
  if (!is.null(dataset)) out <- dplyr::mutate(out, dataset = dataset, .after = method)
  out
}



.bench_evaluate_loading_module_type_enrichment <- function(
    method,
    Y,
    score,
    module_metadata,
    top_n = 100,
    dataset = NULL
) {
  feature_names <- .bench_feature_names(Y)
  score <- as.numeric(score)
  names(score) <- feature_names
  ranked <- names(sort(score, decreasing = TRUE, na.last = NA))
  selected <- ranked[seq_len(min(top_n, length(ranked)))]
  discovered <- feature_names %in% selected

  out <- module_metadata |>
    dplyr::group_by(module_type) |>
    dplyr::summarise(
      set_size = dplyr::n_distinct(feature),
      n_selected = length(selected),
      overlap = sum(discovered & (feature_names %in% unique(feature))),
      odds_ratio = {
        sf <- .bench_safe_fisher_enrichment(discovered, feature_names %in% unique(feature))
        sf$odds_ratio
      },
      p_value = {
        sf <- .bench_safe_fisher_enrichment(discovered, feature_names %in% unique(feature))
        sf$p_value
      },
      .groups = "drop"
    ) |>
    dplyr::mutate(
      q_value = .bench_safe_bh(p_value),
      model_type = .data$module_type,
      method = method,
      .before = 1
    ) |>
    dplyr::relocate(model_type, .after = module_type)

  out <- .bench_add_enrichment_value_transforms(out)
  if (!is.null(dataset)) out <- dplyr::mutate(out, dataset = dataset, .after = method)
  out
}

.bench_evaluate_inference_module_type_enrichment <- function(
    method,
    Y,
    union_tbl,
    module_metadata,
    alpha = 0.1,
    dataset = NULL
) {
  feature_names <- .bench_feature_names(Y)
  u <- union_tbl |>
    dplyr::filter(.data$feature %in% feature_names)

  if (!"scope" %in% names(u)) {
    u$scope <- if ("global_q" %in% names(u)) "global" else "union"
  }
  if (!"component" %in% names(u)) {
    u$component <- u$scope
  }

  groups <- unique(u[c("scope", "component")])
  if (!nrow(groups)) {
    groups <- data.frame(scope = "union", component = "union", stringsAsFactors = FALSE)
  }

  out_list <- lapply(seq_len(nrow(groups)), function(ii) {
    sc <- groups$scope[[ii]]
    comp <- groups$component[[ii]]
    ug <- u[u$scope == sc & u$component == comp, , drop = FALSE]

    discovered <- rep(FALSE, length(feature_names))
    names(discovered) <- feature_names
    if ("discovered" %in% names(ug)) {
      discovered[ug$feature] <- !is.na(ug$discovered) & as.logical(ug$discovered)
    } else if ("global_q" %in% names(ug)) {
      discovered[ug$feature] <- !is.na(ug$global_q) & ug$global_q <= alpha
    } else {
      discovered[ug$feature] <- !is.na(ug$union_q) & ug$union_q <= alpha
    }

    module_metadata |>
      dplyr::group_by(module_type) |>
      dplyr::summarise(
        set_size = dplyr::n_distinct(feature),
        n_discoveries = sum(discovered, na.rm = TRUE),
        overlap = sum(discovered & (feature_names %in% unique(feature)), na.rm = TRUE),
        odds_ratio = {
          sf <- .bench_safe_fisher_enrichment(discovered, feature_names %in% unique(feature))
          sf$odds_ratio
        },
        p_value = {
          sf <- .bench_safe_fisher_enrichment(discovered, feature_names %in% unique(feature))
          sf$p_value
        },
        .groups = "drop"
      ) |>
      dplyr::mutate(
        q_value = .bench_safe_bh(p_value),
        model_type = .data$module_type,
        method = method,
        scope = sc,
        component = comp,
        .before = 1
      ) |>
      dplyr::relocate(model_type, .after = module_type)
  })

  out <- dplyr::bind_rows(out_list)
  out <- .bench_add_enrichment_value_transforms(out)
  if (!is.null(dataset)) out <- dplyr::mutate(out, dataset = dataset, .after = method)
  out
}

# ------------------------------------------------------------
# Module-type summaries
# ------------------------------------------------------------

.bench_module_type_id_cols <- function(df) {
  intersect(
    c(
      "method", "dataset", "module_type", "seed",
      "family", "count_regime_mean",
      "ng", "confounding", "confounding_name",
      "batch_to_biology_ratio", "batch_to_biology_ratio_name",
      "time_imbalance"
    ),
    names(df)
  )
}

# Extract module-level AUC table from either a benchmark object or a data frame.
.bench_get_module_auc_table <- function(x) {
  if (is.data.frame(x)) return(x)
  if (is.list(x) && is.data.frame(x$factorization_module_auc)) {
    return(x$factorization_module_auc)
  }
  NULL
}

# Extract module-level recall summary table from either a benchmark object or a data frame.
.bench_get_module_recall_summary_table <- function(x) {
  if (is.data.frame(x)) return(x)
  if (is.list(x) && is.data.frame(x$factorization_module_recall_summary)) {
    return(x$factorization_module_recall_summary)
  }
  if (is.list(x) && is.data.frame(x$factorization_module_recall)) {
    return(x$factorization_module_recall)
  }
  NULL
}

summarize_module_type_auc <- function(module_auc) {
  module_auc <- .bench_get_module_auc_table(module_auc)
  if (!is.data.frame(module_auc) || nrow(module_auc) == 0L) {
    return(data.frame())
  }
  missing <- setdiff(c("method", "module", "module_type", "auc_module_vs_background"), names(module_auc))
  if (length(missing) > 0L) {
    stop("module_auc is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  group_cols <- .bench_module_type_id_cols(module_auc)
  module_auc |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n_modules = dplyr::n_distinct(.data$module),
      mean_auc = mean(.data$auc_module_vs_background, na.rm = TRUE),
      sd_auc = stats::sd(.data$auc_module_vs_background, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(model_type = .data$module_type) |>
    dplyr::relocate(model_type, .after = module_type)
}

summarize_module_type_recall <- function(module_recall_summary, metric = "aurc") {
  module_recall_summary <- .bench_get_module_recall_summary_table(module_recall_summary)
  if (!is.data.frame(module_recall_summary) || nrow(module_recall_summary) == 0L) {
    return(data.frame())
  }
  missing <- setdiff(c("method", "module", "module_type", metric), names(module_recall_summary))
  if (length(missing) > 0L) {
    stop("module_recall_summary is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  group_cols <- .bench_module_type_id_cols(module_recall_summary)
  module_recall_summary |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n_modules = dplyr::n_distinct(.data$module),
      mean_recall = mean(.data[[metric]], na.rm = TRUE),
      sd_recall = stats::sd(.data[[metric]], na.rm = TRUE),
      recall_metric = metric,
      .groups = "drop"
    ) |>
    dplyr::mutate(model_type = .data$module_type) |>
    dplyr::relocate(model_type, .after = module_type)
}


# ------------------------------------------------------------
# Additive benchmark metrics: reconstruction by truth feature set
# ------------------------------------------------------------

.bench_truth_feature_sets <- function(truth, p, dataset = NULL, index = NULL) {
  get_truth <- function(name) {
    x <- truth[[name]]
    if (is.null(dataset)) return(x)
    if (!is.null(names(x)) && dataset %in% names(x)) return(x[[dataset]])
    if (!is.null(index) && length(x) >= index) return(x[[index]])
    NULL
  }

  signal <- get_truth("signal_features") %||% integer(0)
  null <- get_truth("null_features") %||% setdiff(seq_len(p), signal)
  batch_only <- get_truth("batch_only_features") %||% integer(0)
  clean_null <- get_truth("clean_null_features") %||% integer(0)
  modules <- get_truth("module_features") %||% list()

  base <- list(
    allFeatures = seq_len(p),
    signalFeatures = signal,
    nonSignalFeatures = null,
    batchOnlyFeatures = batch_only,
    cleanNullFeatures = clean_null
  )
  list(base = base, modules = modules)
}

.bench_reconstruction_error_rows <- function(method, Yhat, truth_matrix, truth,
                                             k_scope, dataset = NULL, index = NULL,
                                             reconstruction_source = "training_factorization") {
  Yhat <- as.matrix(Yhat)
  truth_matrix <- as.matrix(truth_matrix)
  if (!all(dim(Yhat) == dim(truth_matrix))) {
    stop("Yhat and truth_matrix must have identical dimensions.", call. = FALSE)
  }

  sets <- .bench_truth_feature_sets(truth, nrow(Yhat), dataset = dataset, index = index)
  one_row <- function(feature_set, idx, module = NA_character_) {
    idx <- as.integer(idx)
    idx <- idx[is.finite(idx) & idx >= 1L & idx <= nrow(Yhat)]
    err <- if (length(idx) > 0L) {
      sum((Yhat[idx, , drop = FALSE] - truth_matrix[idx, , drop = FALSE])^2) /
        (length(idx) * ncol(truth_matrix))
    } else NA_real_
    tibble::tibble(
      method = method,
      dataset = dataset %||% NA_character_,
      k_scope = k_scope,
      feature_set = feature_set,
      module = module,
      reconstruction_source = reconstruction_source,
      reconstruction_error = err,
      n_features_evaluated = length(idx)
    )
  }

  out <- dplyr::bind_rows(lapply(names(sets$base), function(nm) one_row(nm, sets$base[[nm]])))
  if (length(sets$modules) > 0L) {
    module_names <- names(sets$modules)
    if (is.null(module_names) || any(!nzchar(module_names))) module_names <- paste0("module_", seq_along(sets$modules))
    out <- dplyr::bind_rows(out, lapply(seq_along(sets$modules), function(i) {
      one_row("signalModule", sets$modules[[i]], module_names[[i]])
    }))
  }
  out
}

.bench_loading_auc_alias_rows <- function(metrics, common_k_name) {
  if (is.null(metrics) || nrow(metrics) == 0L) return(tibble::tibble())
  common_col <- paste0("loading_feature_auc_signal_vs_non_signal_at_", common_k_name)
  batch_common_col <- paste0("loading_feature_auc_batch_vs_non_batch_at_", common_k_name)
  allk <- metrics |>
    dplyr::transmute(
      method = .data$method,
      dataset = if ("dataset" %in% names(metrics)) .data$dataset else NA_character_,
      k_scope = "allK",
      allSignalFeaturesPooled_auc = .data$loading_feature_auc_signal_vs_non_signal,
      allBatchFeaturesPooled_auc = if ("loading_feature_auc_batch_vs_non_batch" %in% names(metrics)) .data$loading_feature_auc_batch_vs_non_batch else NA_real_
    )
  common <- metrics |>
    dplyr::transmute(
      method = .data$method,
      dataset = if ("dataset" %in% names(metrics)) .data$dataset else NA_character_,
      k_scope = common_k_name,
      allSignalFeaturesPooled_auc = if (common_col %in% names(metrics)) .data[[common_col]] else NA_real_,
      allBatchFeaturesPooled_auc = if (batch_common_col %in% names(metrics)) .data[[batch_common_col]] else NA_real_
    )
  dplyr::bind_rows(allk, common)
}
