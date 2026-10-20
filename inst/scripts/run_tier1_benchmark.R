#!/usr/bin/env Rscript

# Tier 1 benchmark runner: stable core methods only, small grids, result outputs only.
# Run from the benchmark package root with:
#   Rscript inst/scripts/run_tier1_benchmark.R

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
    "benchmark_integrated_grid.R",
    "benchmark_plot_helpers.R",
    "benchmark_single_plots.R",
    "benchmark_integrated_plots.R",
    "benchmark_reports.R"
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

.write_table <- function(x, file) {
  if (is.null(x)) return(invisible(NULL))
  utils::write.csv(x, file = file, row.names = FALSE)
  invisible(file)
}

.env_flag <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) return(default)
  tolower(value) %in% c("1", "true", "yes", "y")
}

.save_result_tables <- function(result, outdir, prefix, write_csv = FALSE) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  result_file <- file.path(outdir, paste0(prefix, "_results_only.rds"))
  saveRDS(result, result_file)

  if (isTRUE(write_csv)) {
    table_names <- names(result)
    table_names <- table_names[vapply(result, is.data.frame, logical(1))]
    for (nm in table_names) {
      .write_table(result[[nm]], file.path(outdir, paste0(prefix, "_", nm, ".csv")))
    }
  }

  result_file
}

root <- .benchmark_root()
outdir <- Sys.getenv("SVDX_BENCHMARK_OUTDIR", unset = file.path(root, "benchmark_results", "tier1"))
write_csv <- .env_flag("SVDX_BENCHMARK_WRITE_CSV", default = FALSE)

.message("Benchmark root: ", root)
.message("Output directory: ", outdir)
.message("Write CSV tables: ", write_csv)
.message("R home: ", R.home())
.message("R version: ", R.version.string)
.message("R executable: ", commandArgs(FALSE)[[1]])
.message("Library paths: ", paste(.libPaths(), collapse = " | "))

suppressPackageStartupMessages({
  library(mySvdx)
})

.source_benchmark_package(root)

deps <- load_benchmark_dependencies(
  packages = c("mySvdx", "dplyr", "tibble"),
  include_optional = FALSE,
  include_reports = FALSE,
  stop_on_missing = TRUE
)
print(deps, row.names = FALSE)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

.message("Running Tier 1 single-dataset grid")
single_grid <- run_single_benchmark_grid(
  seeds = 1:3,
  ng_grid = c(5, 10, 15, 20),
  confounding_grid = c(0, 0.25, 0.5, 0.75, 1),
  batch_to_biology_ratio_grid = c(1, 10),
  time_imbalance_grid = c(0, 0.5, 1),
  keep_runs = FALSE,
  p = 5000,
  G = 4,
  B = 2,
  K = 2,
  n_modules = 4,
  pi_signal = 0.10,
  zformula = ~ batch,
  factorization_methods = c("svdx", "pca", "mixomics_pls"),
  inference_methods = character(0),
  run_factorization = TRUE,
  run_inference = FALSE,
  top_frac = c(0.05, 0.10, 0.20),
  enrichment_top_n = 20,
  numThreads = 2,
  keep_fits = FALSE
)

single_file <- .save_result_tables(
  result = single_grid,
  outdir = outdir,
  prefix = "single_tier1",
  write_csv = write_csv
)
.message("Saved single result object: ", single_file)
.message("Single metric rows: ", nrow(single_grid$factorization_metrics))

.message("Running Tier 1 integrated grid")
integrated_grid <- run_integrated_benchmark_grid(
  seeds = 1:3,
  ng_grid = c(5, 10, 15, 20),
  confounding_grid = list(
    all_clean = c(0, 0),
    all_025 = c(0.25, 0.25),
    all_050 = c(0.5, 0.5),
    all_075 = c(0.75, 0.75),
    all_confounded = c(1, 1)
  ),
  batch_to_biology_ratio_grid = list(
    balanced = c(1, 1),
    all_strong_batch = c(10, 10)
  ),
  time_imbalance_grid = c(0, 0.5, 1),
  keep_runs = FALSE,
  D = 2,
  p = c(5000, 5000),
  G = 4,
  B = c(2, 2),
  K = 2,
  n_modules = 4,
  pi_signal = c(0.10, 0.10),
  sigma = c(1, 1),
  module_sd = c(1.5, 1.5),
  batch_feature_fraction = c(0.10, 0.10),
  factorization_methods = c("isvdx", "separate_svdx", "concat_PCA"),
  inference_methods = character(0),
  run_factorization = TRUE,
  run_inference = FALSE,
  top_frac = c(0.05, 0.10, 0.20),
  enrichment_top_n = 20,
  numThreads = 2,
  keep_fits = FALSE
)

integrated_file <- .save_result_tables(
  result = integrated_grid,
  outdir = outdir,
  prefix = "integrated_tier1",
  write_csv = write_csv
)
.message("Saved integrated result object: ", integrated_file)
.message("Integrated metric rows: ", nrow(integrated_grid$factorization_metrics))

.message("Writing Tier 1 stratified benchmark plots")
report <- runSimulationBenchmarkReports(
  single_grid = single_grid,
  integrated_grid = integrated_grid,
  outdir = file.path(outdir, "plots"),
  figure_format = "pdf",
  overwrite = TRUE,
  verbose = TRUE
)
.message("Saved Tier 1 plot files: ", nrow(report$index))

.message("Tier 1 benchmark completed")
