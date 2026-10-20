#!/usr/bin/env Rscript

# Tiny end-to-end smoke benchmark for package development.
# Run from the benchmark package root with:
#   Rscript inst/scripts/run_smoke_benchmark.R

options(svdxBenchmark.svdx_package = Sys.getenv("SVDX_PACKAGE", "mySvdx"))

.message <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), "|", paste0(..., collapse = ""), "\n")
}

.script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  marker <- "--file="
  hit <- grep(paste0("^", marker), args, value = TRUE)
  if (length(hit) == 0L) return(NULL)
  normalizePath(sub(paste0("^", marker), "", hit[[1]]), mustWork = FALSE)
}

.benchmark_root <- function() {
  env_root <- Sys.getenv("SVDX_BENCHMARK_PATH", unset = NA_character_)
  if (!is.na(env_root) && nzchar(env_root)) {
    return(normalizePath(env_root, mustWork = FALSE))
  }

  sp <- .script_path()
  if (!is.null(sp)) {
    return(normalizePath(file.path(dirname(sp), "..", ".."), mustWork = FALSE))
  }

  normalizePath(getwd(), mustWork = FALSE)
}

.source_benchmark_package <- function(root) {
  r_dir <- file.path(root, "R")
  files <- file.path(r_dir, c(
    "benchmark_dependency_check.R",
    "benchmark_utils.R",
    "benchmark_simulation_svdx_interactions.R",
    "benchmark_simulation_isvdx_interactions.R",
    "benchmark_integrated_preprocess.R",
    "benchmark_evaluation_helpers.R",
    "benchmark_shared_evaluation.R",
    "benchmark_targetSetEvaluation.R",
    "benchmark_single_fitters.R",
    "benchmark_integrated_mefisto.R",
    "benchmark_integrated_fitters.R",
    "benchmark_single_evaluation.R",
    "benchmark_integrated_evaluation.R",
    "benchmark_single_grid.R",
    "benchmark_integrated_grid.R"
  ))

  missing <- files[!file.exists(files)]
  if (length(missing) > 0L) {
    stop(
      "Cannot source benchmark package; missing files: ",
      paste(basename(missing), collapse = ", "),
      call. = FALSE
    )
  }

  invisible(lapply(files, sys.source, envir = .GlobalEnv))
}

.load_benchmark_package <- function(root) {
  use_installed <- identical(tolower(Sys.getenv("SVDX_BENCHMARK_USE_INSTALLED")), "true")
  if (isTRUE(use_installed) && requireNamespace("mySvdxBenchmark", quietly = TRUE)) {
    suppressPackageStartupMessages(library(mySvdxBenchmark))
    return("installed")
  }

  .source_benchmark_package(root)
  "sourced"
}

root <- .benchmark_root()
.message("Benchmark root: ", root)
.message("R home: ", R.home())
.message("R version: ", R.version.string)
.message("R executable: ", commandArgs(FALSE)[[1]])
.message("Library paths: ", paste(.libPaths(), collapse = " | "))

dependency_helper <- file.path(root, "R", "benchmark_dependency_check.R")
if (file.exists(dependency_helper)) {
  sys.source(dependency_helper, envir = .GlobalEnv)
}

suppressPackageStartupMessages({
  library(mySvdx)
})

load_mode <- .load_benchmark_package(root)
.message("Benchmark package load mode: ", load_mode)

deps <- load_benchmark_dependencies(
  packages = c("mySvdx", "dplyr", "tibble"),
  include_optional = FALSE,
  include_reports = FALSE,
  stop_on_missing = TRUE
)
print(deps, row.names = FALSE)

.message("Running tiny single-dataset benchmark")
single <- run_single_benchmark_once(
  seed = 1,
  p = 60,
  G = 2,
  ng = 8,
  B = 2,
  K = 2,
  n_modules = 2,
  zformula = ~ batch,
  factorization_methods = c("svdx", "svd", "pca"),
  inference_methods = character(0),
  run_factorization = TRUE,
  run_inference = FALSE,
  top_frac = c(0.10, 0.20),
  enrichment_top_n = 10,
  numThreads = 2,
  keep_fits = FALSE
)

.message("Single metrics rows: ", nrow(single$factorization$metrics))
.message("Single module AUC rows: ", nrow(single$factorization$module_auc))
print(utils::head(single$factorization$metrics), row.names = FALSE)

.message("Running tiny integrated benchmark")
integrated <- run_integrated_benchmark_once(
  seed = 1,
  D = 2,
  p = c(50, 40),
  G = 2,
  ng = 8,
  B = c(2, 2),
  K = 2,
  n_modules = 2,
  pi_signal = c(0.35, 0.30),
  sigma = c(1, 1),
  module_sd = c(1.5, 1.5),
  batch_feature_fraction = c(0.10, 0.10),
  factorization_methods = c("isvdx", "separate_svdx", "concat_PCA"),
  inference_methods = character(0),
  run_factorization = TRUE,
  run_inference = FALSE,
  top_frac = c(0.10, 0.20),
  enrichment_top_n = 10,
  numThreads = 2,
  keep_fits = FALSE
)

.message("Integrated metrics rows: ", nrow(integrated$factorization$metrics))
.message("Integrated module AUC rows: ", nrow(integrated$factorization$module_auc))
print(utils::head(integrated$factorization$metrics), row.names = FALSE)

.message("Smoke benchmark completed")
