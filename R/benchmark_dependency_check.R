# Dependency helpers for reproducible benchmark runs.

benchmark_dependency_sets <- function(
    include_optional = FALSE,
    include_reports = TRUE
) {
  svdx_pkg <- if (requireNamespace("mySvdx", quietly = TRUE)) {
    "mySvdx"
  } else {
    "svdx"
  }
  core <- c(svdx_pkg, "dplyr", "tibble")
  evaluation <- c("limma", "pROC")
  reports <- if (isTRUE(include_reports)) c("ggplot2") else character(0)
  optional <- if (isTRUE(include_optional)) c("mixOmics", "MOFA2") else character(0)

  list(
    core = core,
    evaluation = evaluation,
    reports = reports,
    optional_competitors = optional
  )
}

check_benchmark_dependencies <- function(
    packages = NULL,
    requested_methods = NULL,
    include_optional = FALSE,
    include_reports = TRUE,
    attach = FALSE,
    stop_on_missing = FALSE
) {
  sets <- benchmark_dependency_sets(
    include_optional = include_optional,
    include_reports = include_reports
  )

  if (is.null(packages)) {
    packages <- unname(unlist(sets, use.names = FALSE))
  }
  if (!is.null(requested_methods) && "glm_pca" %in% requested_methods) {
    packages <- c(packages, "glmpca")
  }
  nmf_methods <- c("nmf", "sparse_nmf", "concat_NMF", "concat_sNMF")
  needs_sklearn <- !is.null(requested_methods) && any(requested_methods %in% nmf_methods)
  if (isTRUE(needs_sklearn)) {
    packages <- c(packages, "reticulate")
  }
  packages <- unique(packages)

  set_name <- rep(NA_character_, length(packages))
  for (nm in names(sets)) {
    set_name[packages %in% sets[[nm]]] <- nm
  }
  set_name[packages == "glmpca"] <- "optional_competitors"
  set_name[packages == "reticulate"] <- "optional_competitors"

  installed <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  version <- rep(NA_character_, length(packages))
  lib_path <- rep(NA_character_, length(packages))

  for (i in seq_along(packages)) {
    if (!installed[[i]]) next
    version[[i]] <- as.character(utils::packageVersion(packages[[i]]))
    lib_path[[i]] <- tryCatch(
      dirname(system.file(package = packages[[i]])),
      error = function(e) NA_character_
    )
  }

  if (isTRUE(attach)) {
    for (pkg in packages[installed]) {
      suppressPackageStartupMessages(
        library(pkg, character.only = TRUE)
      )
    }
  }

  attached <- paste0("package:", packages) %in% search()
  out <- data.frame(
    package = packages,
    set = set_name,
    installed = installed,
    attached = attached,
    version = version,
    lib_path = lib_path,
    stringsAsFactors = FALSE
  )

  missing <- out$package[!out$installed]
  if (length(missing) > 0L && isTRUE(stop_on_missing)) {
    stop(
      "Missing benchmark package dependencies: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  if (isTRUE(needs_sklearn) &&
      requireNamespace("reticulate", quietly = TRUE) &&
      !reticulate::py_module_available("sklearn")) {
    msg <- "Python module 'sklearn' is required for NMF benchmark methods."
    if (isTRUE(stop_on_missing)) {
      stop(msg, call. = FALSE)
    } else {
      warning(msg, call. = FALSE)
    }
  }

  out
}

load_benchmark_dependencies <- function(
    packages = NULL,
    requested_methods = NULL,
    include_optional = FALSE,
    include_reports = TRUE,
    stop_on_missing = TRUE
) {
  check_benchmark_dependencies(
    packages = packages,
    requested_methods = requested_methods,
    include_optional = include_optional,
    include_reports = include_reports,
    attach = TRUE,
    stop_on_missing = stop_on_missing
  )
}
