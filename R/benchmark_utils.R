#benchmark utils - 

`%||%` <- function(x, y) if (!is.null(x)) x else y

.require_namespace <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.", call. = FALSE)
  }
}

.safe_auc_score <- function(score, truth) {
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
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

.safe_bh <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p) & !is.na(p)
  out[ok] <- stats::p.adjust(p[ok], method = "BH")
  out
}


.safe_fisher_enrichment <- function(discovered, in_set) {
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
  list(
    odds_ratio = unname(ft$estimate %||% NA_real_),
    p_value = unname(ft$p.value %||% NA_real_)
  )
}

.benchmark_find_svdx_function <- function(name, required = FALSE) {
  pkg_opt <- getOption("svdxBenchmark.svdx_package", NULL)
  candidates <- unique(stats::na.omit(c(pkg_opt, "mySvdx", "svdx", "mSvdx", "SVDX")))
  for (pkg in candidates) {
    if (requireNamespace(pkg, quietly = TRUE) &&
        exists(name, envir = asNamespace(pkg), inherits = FALSE)) {
      return(getFromNamespace(name, pkg))
    }
  }

  # Last-resort fallback for interactive/source() workflows.  Package
  # namespace lookup is intentionally preferred so stale functions in the
  # user/global environment do not shadow the installed SVDX implementation.
  if (exists(name, mode = "function", inherits = TRUE)) {
    return(get(name, mode = "function", inherits = TRUE))
  }

  if (isTRUE(required)) {
    stop("Could not find SVDX package function `", name, "`. ",
         "Set options(svdxBenchmark.svdx_package = '<packageName>') ",
         "or load/install the SVDX package.", call. = FALSE)
  }
  NULL
}

.safe_cor_vec <- function(x, y, method = "pearson") {
  x <- as.numeric(x); y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  stats::cor(x[ok], y[ok], method = method)
}

.safe_rel_frob <- function(A, B, eps = 1e-8) {
  num <- sqrt(sum((A - B)^2, na.rm = TRUE))
  den <- sqrt(sum(B^2, na.rm = TRUE))
  if (!is.finite(den) || den < eps) return(NA_real_)
  num / den
}

.benchmark_complete_model_type <- function(
    modelType = NULL,
    formula = NULL,
    data = NULL,
    default_numeric = "cubic"
) {
  if (is.null(formula) || is.null(data)) {
    return(modelType)
  }

  data <- as.data.frame(data)
  vars <- all.vars(formula)
  vars <- intersect(vars, colnames(data))

  out <- modelType %||% list()

  for (v in vars) {
    if (!is.null(out[[v]])) next

    x <- data[[v]]
    if (is.numeric(x)) {
      rng <- range(x, na.rm = TRUE)
      if (length(rng) == 2L && all(is.finite(rng))) {
        out[[v]] <- list(default_numeric, rng)
      }
    }
  }

  if (length(out) == 0L) NULL else out
}


.benchmark_adjust_Y_list_by_formula <- function(
    Y,
    zdata = NULL,
    zformula = NULL
) {
  if (is.null(zformula)) return(Y)
  if (is.null(zdata)) return(Y)

  if (!is.list(Y)) {
    stop("Y must be a list of matrices.", call. = FALSE)
  }
  if (!is.list(zdata)) {
    stop("zdata must be a list when adjusting integrated baseline matrices.", call. = FALSE)
  }
  if (!is.list(zformula)) {
    stop("zformula must be a list when adjusting integrated baseline matrices.", call. = FALSE)
  }

  dataset_names <- names(Y) %||% paste0("dataset_", seq_along(Y))
  names(Y) <- dataset_names
  if (is.null(names(zdata))) names(zdata) <- dataset_names
  if (is.null(names(zformula))) names(zformula) <- dataset_names

  out <- lapply(dataset_names, function(ds) {
    Yd <- Y[[ds]]
    zd <- zdata[[ds]]
    zf <- zformula[[ds]]

    if (is.null(zf) || is.null(zd)) {
      return(Yd)
    }

    adjust_Y_by_formula(Yd, zd, zf)
  })
  stats::setNames(out, dataset_names)
}



.benchmark_get_svdx_internal <- function(name) {
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

  stop(
    "Could not find SVDX package helper `", name, "`. ",
    "Set options(svdxBenchmark.svdx_package = '<packageName>') ",
    "or load/install the SVDX package.",
    call. = FALSE
  )
}

.benchmark_fdata_to_enrichment_id_map <- function(fdata, required = TRUE) {
  if (is.null(fdata)) return(NULL)
  fn <- .benchmark_get_svdx_internal(".svdx_fdata_to_enrichment_id_map")
  fn(fdata, required = required)
}

.benchmark_validate_fdata_matches_features <- function(fdata, features, context = "Y") {
  if (is.null(fdata)) return(invisible(TRUE))
  fn <- .benchmark_get_svdx_internal(".svdx_validate_fdata_matches_features")
  fn(fdata, features, context = context)
  invisible(TRUE)
}


.benchmark_resolve_fdata_for_dataset <- function(fdata, dataset = NULL, required = FALSE) {
  if (is.null(fdata)) return(NULL)
  fn <- .benchmark_get_svdx_internal(".svdx_resolve_fdata_for_dataset")
  fn(fdata, dataset = dataset, required = required)
}

.benchmark_fdata_to_row_labels <- function(fdata, dataset = NULL) {
  if (is.null(fdata)) return(NULL)
  fdata_use <- .benchmark_resolve_fdata_for_dataset(fdata, dataset = dataset, required = FALSE)
  if (is.null(fdata_use)) return(NULL)
  fn <- .benchmark_get_svdx_internal(".svdx_fdata_to_row_labels")
  fn(fdata_use)
}

.benchmark_normalize_row_labels <- function(row_labels = NULL, fdata = NULL, dataset = NULL) {
  if (!is.null(row_labels)) return(row_labels)
  .benchmark_fdata_to_row_labels(fdata, dataset = dataset)
}

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

.benchmark_resolve_feature_ids_and_labels <- function(features, feature_names, row_labels = NULL) {
  feature_names <- as.character(feature_names)
  if (is.null(features)) {
    features <- feature_names
  }
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
