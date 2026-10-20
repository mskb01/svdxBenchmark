# ============================================================
# Benchmark report assemblers
# ============================================================
# These functions assemble plots/tables from already-created benchmark
# objects. They do not run benchmark analyses.

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}


.benchmark_report_find_svdx_function <- function(name, required = FALSE) {
  if (exists(".benchmark_find_svdx_function", mode = "function", inherits = TRUE)) {
    return(.benchmark_find_svdx_function(name, required = required))
  }
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
  if (isTRUE(required)) stop("Could not find SVDX package function `", name, "`.", call. = FALSE)
  NULL
}

.benchmark_report_message <- function(..., verbose = TRUE) {
  if (!isTRUE(verbose)) return(invisible(NULL))
  cat(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = ""), "\n")
  invisible(NULL)
}

.benchmark_report_prepare_dir <- function(outdir, overwrite = TRUE) {
  if (!is.character(outdir) || length(outdir) != 1L || is.na(outdir) || !nzchar(outdir)) {
    stop("outdir must be a single non-empty character string.", call. = FALSE)
  }
  outdir <- path.expand(outdir)
  if (dir.exists(outdir) && isTRUE(overwrite)) unlink(outdir, recursive = TRUE, force = TRUE)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(outdir)) stop("Failed to create outdir: ", outdir, call. = FALSE)
  normalizePath(outdir, mustWork = TRUE)
}


.benchmark_report_named_list <- function(x, default_prefix = "dataset") {
  if (is.null(x)) return(NULL)
  if (!is.list(x) || inherits(x, "realWorldMethodBenchmark") || inherits(x, "realWorldIntegratedBenchmark")) {
    x <- list(x)
  }
  if (is.null(names(x)) || any(!nzchar(names(x)))) {
    names(x) <- paste0(default_prefix, "_", seq_along(x))
  }
  x
}


.benchmark_report_write_prediction_metrics <- function(outdir, prediction_metrics, filename = "prediction_metrics.tsv") {
  if (is.null(prediction_metrics) || !is.data.frame(prediction_metrics) || nrow(prediction_metrics) == 0L) {
    return(NULL)
  }
  file <- .benchmark_report_data_file(outdir, filename)
  utils::write.table(prediction_metrics, file = file, sep = "\t", row.names = FALSE, quote = FALSE)
  file
}

.benchmark_report_plot_prediction_metrics <- function(prediction_metrics, y_limits = c(0, 100)) {
  if (is.null(prediction_metrics) || !is.data.frame(prediction_metrics) || nrow(prediction_metrics) == 0L) {
    return(NULL)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  df <- as.data.frame(prediction_metrics)
  if (!"explained_variation_percent" %in% names(df)) return(NULL)
  df$explained_variation_percent <- suppressWarnings(as.numeric(df$explained_variation_percent))
  df$test_mse <- suppressWarnings(as.numeric(df$test_mse))
  df <- df[is.finite(df$explained_variation_percent), , drop = FALSE]
  if (!nrow(df)) return(NULL)
  df$dataset <- if ("dataset" %in% names(df)) as.character(df$dataset) else NA_character_
  df$method <- if ("method" %in% names(df)) as.character(df$method) else "svdx"
  df$label <- ifelse(is.finite(df$test_mse), paste0("MSE=", signif(df$test_mse, 3)), "")
  df$x_label <- ifelse(!is.na(df$dataset) & nzchar(df$dataset), paste(df$method, df$dataset, sep = "\n"), df$method)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x_label, y = explained_variation_percent)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.3, size = 3) +
    ggplot2::labs(
      x = NULL,
      y = "% of unseen test-data variation explained by SVDX predictions"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  if (!is.null(y_limits) && length(y_limits) == 2L) {
    ymax <- max(y_limits[2], max(df$explained_variation_percent, na.rm = TRUE) * 1.12, na.rm = TRUE)
    ymin <- min(y_limits[1], min(df$explained_variation_percent, na.rm = TRUE), na.rm = TRUE)
    p <- p + ggplot2::coord_cartesian(ylim = c(ymin, ymax))
  }
  p
}

.benchmark_report_compute_target_auc <- function(rw_obj, target_sets = NULL, integrated = FALSE, target_set_auc_args = list()) {
  if (!is.null(rw_obj$target_auc)) return(rw_obj$target_auc)
  target_sets <- target_sets %||% rw_obj$target_sets
  if (is.null(target_sets)) return(NULL)
  if (!is.list(target_set_auc_args)) {
    stop("target_set_auc_args must be a list.", call. = FALSE)
  }
  if (isTRUE(integrated)) {
    do.call(
      evaluate_integrated_target_set_auc,
      c(list(rw_int = rw_obj, target_sets = target_sets), target_set_auc_args)
    )
  } else {
    do.call(
      evaluate_target_set_auc,
      c(list(rw_bench = rw_obj, target_sets = target_sets), target_set_auc_args)
    )
  }
}

.benchmark_report_is_list_of <- function(x, class_name) {
  is.list(x) && !inherits(x, class_name) && length(x) > 0L &&
    all(vapply(x, inherits, logical(1), what = class_name))
}

.benchmark_report_safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nzchar(x), x, "dataset")
}

.benchmark_report_validate_figure_format <- function(figure_format = "pdf") {
  match.arg(figure_format, c("pdf", "png"))
}

.benchmark_report_scope_from_stem <- function(stem_l, default = "component") {
  if (grepl("(^|_)global($|_)", stem_l)) return("global")
  if (grepl("(^|_)union($|_)", stem_l)) return("union")
  if (grepl("(^|_)component($|_)", stem_l) || grepl("component_[^_]+", stem_l)) return("component")
  default
}

.benchmark_report_subdir_for_stem <- function(stem) {
  stem0 <- sub("\\.(png|pdf)$", "", as.character(stem), ignore.case = TRUE)
  stem0 <- basename(stem0)
  stem_l <- tolower(stem0)

  if (grepl("^u_", stem_l) || grepl("_(u|u_embedding)_", stem_l) ||
      grepl("_u_", stem_l) || grepl("u_embedding", stem_l)) {
    if (grepl("embedding", stem_l)) return(file.path("U", "embeddings"))
    return(file.path("U", "component_plots"))
  }

  if (grepl("v_loading_heatmap", stem_l)) {
    sc <- .benchmark_report_scope_from_stem(stem_l, default = "component")
    return(file.path("V", "loadings", sc, "heatmap"))
  }
  if (grepl("v_loadings", stem_l)) {
    return(file.path("V", "loadings", "component", "loading_bars"))
  }
  if (grepl("feature_profiles_loading", stem_l)) {
    sc <- .benchmark_report_scope_from_stem(stem_l, default = "component")
    return(file.path("V", "loadings", sc, "feature_profiles"))
  }
  if (grepl("loading_enrichment", stem_l)) {
    sc <- .benchmark_report_scope_from_stem(stem_l, default = "union")
    if (identical(sc, "global")) sc <- "union"
    return(file.path("V", "loadings", sc, "enrichment"))
  }

  if (grepl("v_inference_heatmap", stem_l)) {
    sc <- .benchmark_report_scope_from_stem(stem_l, default = "component")
    return(file.path("V", "inference", sc, "heatmap"))
  }
  if (grepl("v_inference_dotplot", stem_l)) {
    sc <- .benchmark_report_scope_from_stem(stem_l, default = "component")
    return(file.path("V", "inference", sc, "dotplot"))
  }
  if (grepl("feature_profiles_inference", stem_l)) {
    sc <- .benchmark_report_scope_from_stem(stem_l, default = "component")
    return(file.path("V", "inference", sc, "feature_profiles"))
  }
  if (grepl("inference_volcano", stem_l)) {
    sc <- .benchmark_report_scope_from_stem(stem_l, default = "component")
    return(file.path("V", "inference", sc, "volcano"))
  }
  if (grepl("inference_enrichment", stem_l)) {
    sc <- .benchmark_report_scope_from_stem(stem_l, default = "union")
    return(file.path("V", "inference", sc, "enrichment"))
  }
  if (grepl("v_loading_inference_scatter", stem_l)) {
    return(file.path("V", "inference", "loading_inference_scatter"))
  }

  if (grepl("relevant_subset_prioritization", stem_l) || grepl("target_auc", stem_l)) {
    return(file.path("V", "inference", "target_auc"))
  }

  if (grepl("prediction|explained_variation|test_mse", stem_l)) return("prediction")
  if (grepl("runtime", stem_l)) return(file.path("simulation", "runtime"))
  if (grepl("enrichment", stem_l)) return(file.path("simulation", "enrichment"))
  if (grepl("metrics|auc|recall", stem_l)) return(file.path("simulation", "metrics"))

  "misc"
}

.benchmark_report_data_file <- function(outdir, filename) {
  data_dir <- file.path(outdir, "data")
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
  file.path(data_dir, filename)
}


.benchmark_report_plot_file <- function(outdir, stem, figure_format = "pdf", subdir = NULL) {
  figure_format <- .benchmark_report_validate_figure_format(figure_format)
  stem <- sub("\\.(png|pdf)$", "", stem, ignore.case = TRUE)
  subdir <- subdir %||% .benchmark_report_subdir_for_stem(stem)
  plot_dir <- if (is.null(subdir) || !nzchar(subdir)) outdir else file.path(outdir, subdir)
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  file.path(plot_dir, paste0(stem, ".", figure_format))
}

.benchmark_report_save_plot <- function(p, file, width = 10, height = 7,
                                        figure_format = NULL, verbose = TRUE) {
  if (is.null(p)) return(NULL)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .benchmark_report_message("Skipping plot; ggplot2 is not available: ", basename(file), verbose = verbose)
    return(NULL)
  }
  if (!is.null(figure_format)) {
    figure_format <- .benchmark_report_validate_figure_format(figure_format)
    file <- sub("\\.(png|pdf)$", paste0(".", figure_format), file, ignore.case = TRUE)
  }
  ext <- tolower(tools::file_ext(file))

  # ComplexHeatmap objects are not ggplot objects and must be drawn with
  # ComplexHeatmap::draw(), matching the package report helper behavior.
  if (inherits(p, "Heatmap") || inherits(p, "HeatmapList") || inherits(p, "AdditiveUnit")) {
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      .benchmark_report_message("Skipping heatmap plot; ComplexHeatmap is not available: ", basename(file), verbose = verbose)
      return(NULL)
    }
    ok <- tryCatch({
      if (identical(ext, "png")) {
        grDevices::png(file, width = width, height = height, units = "in", res = 300)
      } else {
        grDevices::cairo_pdf(file, width = width, height = height)
      }
      on.exit(grDevices::dev.off(), add = TRUE)
      ComplexHeatmap::draw(p)
      TRUE
    }, error = function(e) {
      .benchmark_report_message("Skipping plot: ", basename(file), " | ", e$message, verbose = verbose)
      FALSE
    })
    return(if (isTRUE(ok)) file else NULL)
  }

  # pheatmap() returns a grid/gtable-backed object, not a ggplot. Save it
  # explicitly so generic integrated real-world V heatmaps can be written.
  if (inherits(p, "pheatmap")) {
    if (!requireNamespace("grid", quietly = TRUE)) {
      .benchmark_report_message("Skipping heatmap plot; grid is not available: ", basename(file), verbose = verbose)
      return(NULL)
    }
    ok <- tryCatch({
      if (identical(ext, "png")) {
        grDevices::png(file, width = width, height = height, units = "in", res = 300)
      } else {
        grDevices::cairo_pdf(file, width = width, height = height)
      }
      on.exit(grDevices::dev.off(), add = TRUE)
      grid::grid.newpage()
      grid::grid.draw(p$gtable)
      TRUE
    }, error = function(e) {
      .benchmark_report_message("Skipping plot: ", basename(file), " | ", e$message, verbose = verbose)
      FALSE
    })
    return(if (isTRUE(ok)) file else NULL)
  }

  device <- switch(
    ext,
    pdf = grDevices::cairo_pdf,
    png = "png",
    NULL
  )
  ok <- tryCatch({
    ggplot2::ggsave(
      filename = file,
      plot = p,
      width = width,
      height = height,
      units = "in",
      dpi = 300,
      device = device
    )
    TRUE
  }, error = function(e) {
    .benchmark_report_message("Skipping plot: ", basename(file), " | ", e$message, verbose = verbose)
    FALSE
  })
  if (isTRUE(ok)) file else NULL
}

.benchmark_report_safe_plot <- function(expr, file, width = 10, height = 7,
                                        figure_format = NULL, verbose = TRUE) {
  p <- tryCatch(expr, error = function(e) e)
  if (inherits(p, "error")) {
    .benchmark_report_message("Skipping plot: ", basename(file), " | ", p$message, verbose = verbose)
    return(NULL)
  }
  .benchmark_report_save_plot(
    p,
    file,
    width = width,
    height = height,
    figure_format = figure_format,
    verbose = verbose
  )
}

.benchmark_report_safe_plot_collection <- function(expr, outdir, stem, width = 10, height = 7,
                                                   figure_format = "pdf", verbose = TRUE) {
  obj <- tryCatch(expr, error = function(e) e)
  if (inherits(obj, "error")) {
    file <- .benchmark_report_plot_file(outdir, stem, figure_format)
    .benchmark_report_message("Skipping plot: ", basename(file), " | ", obj$message, verbose = verbose)
    return(character(0))
  }
  if (is.null(obj)) return(character(0))

  if (inherits(obj, "ggplot") || inherits(obj, "Heatmap") || inherits(obj, "HeatmapList") || inherits(obj, "AdditiveUnit")) {
    file <- .benchmark_report_save_plot(
      obj,
      .benchmark_report_plot_file(outdir, stem, figure_format),
      width = width,
      height = height,
      figure_format = figure_format,
      verbose = verbose
    )
    return(file %||% character(0))
  }

  if (!is.list(obj)) {
    file <- .benchmark_report_plot_file(outdir, stem, figure_format)
    .benchmark_report_message("Skipping plot: ", basename(file), " | returned object is not a plot or list of plots.", verbose = verbose)
    return(character(0))
  }

  nms <- names(obj)
  if (is.null(nms)) nms <- rep("", length(obj))
  out <- character(0)
  for (ii in seq_along(obj)) {
    term_plot <- obj[[ii]]
    if (is.null(term_plot)) next
    term_name <- nms[[ii]]
    if (is.na(term_name) || !nzchar(term_name)) term_name <- paste0("term_", ii)
    file <- .benchmark_report_save_plot(
      term_plot,
      .benchmark_report_plot_file(outdir, paste0(stem, "_", .benchmark_report_safe_name(term_name)), figure_format),
      width = width,
      height = height,
      figure_format = figure_format,
      verbose = verbose
    )
    if (!is.null(file)) out <- c(out, file)
  }
  out
}

.benchmark_report_methods <- function(x, integrated = FALSE) {
  nms <- if (isTRUE(integrated)) names(x$results) else names(x$benchmark_results)
  nms <- nms[!is.na(nms) & nzchar(nms)]
  nms
}

.benchmark_report_numeric_or_factor <- function(data, var) {
  if (!var %in% colnames(data)) return("missing")
  x <- data[[var]]
  if (is.numeric(x) || inherits(x, c("Date", "POSIXct", "POSIXlt", "POSIXt"))) return("numeric")
  if (is.factor(x) || is.character(x) || is.logical(x)) return("factor")
  "other"
}

.benchmark_report_formula_vars <- function(xformula) {
  if (is.null(xformula)) return(character(0))
  if (inherits(xformula, "formula")) return(all.vars(xformula))
  if (is.list(xformula)) {
    vars <- unlist(lapply(xformula, .benchmark_report_formula_vars), use.names = FALSE)
    return(unique(vars[!is.na(vars) & nzchar(vars)]))
  }
  character(0)
}


.benchmark_report_svdx_model_vars <- function(fit, fallback_formula = NULL) {
  vars <- tryCatch(
    all.vars(fit$fitAtOptCV$trainModel$mod$formula),
    error = function(e) character(0)
  )
  if (length(vars) == 0L) vars <- .benchmark_report_formula_vars(fallback_formula)
  unique(vars[!is.na(vars) & nzchar(vars)])
}

.benchmark_report_svdx_process_vars <- function(vars, model_vars, what = "plotting", verbose = TRUE) {
  if (is.null(vars)) return(NULL)
  vars <- unique(as.character(vars))
  keep <- intersect(vars, model_vars)
  dropped <- setdiff(vars, model_vars)
  if (length(dropped) > 0L) {
    .benchmark_report_message(
      "SVDX process ", what, ": ignoring variable(s) not in the fitted process model: ",
      paste(dropped, collapse = ", "),
      verbose = verbose
    )
  }
  if (length(keep) == 0L) NULL else keep
}

.benchmark_report_plotU_call_args <- function(fit, data, component, x, group = NULL,
                                              facet = NULL, fill = "typical",
                                              xformula = NULL, plot_args = list()) {
  args <- list(fit = fit, data = data, component = component, x = x,
               group = group, facet = facet, fill = fill)
  if (!"observed" %in% names(plot_args)) {
    model_vars <- .benchmark_report_formula_vars(xformula)
    args$observed <- !x %in% model_vars
  }
  c(args, plot_args)
}


.benchmark_report_svdx_u_plot_specs <- function(data, xaxes = NULL, u_plot_specs = NULL) {
  data <- as.data.frame(data)
  if (!is.null(u_plot_specs)) {
    specs <- lapply(u_plot_specs, function(spec) {
      var <- .benchmark_report_spec_value(spec, "var", .benchmark_report_spec_value(spec, "x", NULL))
      if (is.null(var)) return(NULL)
      if (!var %in% colnames(data)) return(NULL)
      spec$var <- var
      spec
    })
    specs <- Filter(Negate(is.null), specs)
    if (length(specs) > 0L) return(specs)
  }
  xaxes <- unique(stats::na.omit(xaxes))
  lapply(xaxes[xaxes %in% colnames(data)], function(xa) list(var = xa))
}

.benchmark_report_add_connected_u_layer <- function(p, connect = FALSE, group = NULL,
                                                    line_group = NULL) {
  if (!isTRUE(connect) || is.null(p) || !inherits(p, "ggplot")) return(p)
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("rlang", quietly = TRUE)) {
    return(p)
  }
  line_group <- line_group %||% group
  if (!is.null(line_group)) {
    group_info <- .benchmark_plot_interaction_var(
      p$data, line_group, name = ".benchmark_plot_line_group"
    )
    p + ggplot2::geom_line(
      data = group_info$data,
      mapping = ggplot2::aes(group = .data[[group_info$var]]),
      linewidth = 0.6
    )
  } else {
    p + ggplot2::geom_line(ggplot2::aes(group = 1), linewidth = 0.6)
  }
}

.benchmark_report_default_u_plot_specs <- function(data, xformula = NULL, max_vars = 6L) {
  vars <- character(0)
  if (!is.null(xformula)) vars <- all.vars(xformula)
  vars <- intersect(vars, colnames(data))
  if (length(vars) == 0L) vars <- colnames(data)
  vars <- vars[seq_len(min(length(vars), max_vars))]
  lapply(vars, function(v) list(var = v, type = .benchmark_report_numeric_or_factor(data, v)))
}

.benchmark_report_has_spec_field <- function(spec, name) {
  is.list(spec) && name %in% names(spec)
}

.benchmark_report_spec_value <- function(spec, name, default = NULL) {
  if (.benchmark_report_has_spec_field(spec, name)) spec[[name]] else default
}

.benchmark_report_prepare_u_plot_specs <- function(data, xformula = NULL, u_plot_specs = NULL,
                                                  xaxis = NULL, max_vars = 6L) {
  data <- as.data.frame(data)
  if (is.null(u_plot_specs)) {
    if (!is.null(xaxis) && length(xaxis) == 1L && nzchar(xaxis)) {
      u_plot_specs <- list(list(var = xaxis, type = .benchmark_report_numeric_or_factor(data, xaxis)))
    } else {
      u_plot_specs <- .benchmark_report_default_u_plot_specs(data, xformula = xformula, max_vars = max_vars)
    }
  }
  if (!is.list(u_plot_specs) || (length(u_plot_specs) > 0L && !is.list(u_plot_specs[[1L]]))) {
    stop("u_plot_specs must be a list of plot-specification lists.", call. = FALSE)
  }
  lapply(seq_along(u_plot_specs), function(i) {
    spec <- u_plot_specs[[i]]
    var <- .benchmark_report_spec_value(spec, "var", .benchmark_report_spec_value(spec, "x", NULL))
    if (is.null(var) && !is.null(xaxis) && length(xaxis) == 1L && nzchar(xaxis)) {
      spec$var <- xaxis
      var <- xaxis
    }
    if (!is.null(var) && !.benchmark_report_has_spec_field(spec, "type")) {
      spec$type <- .benchmark_report_numeric_or_factor(data, var)
    }
    spec
  })
}

.benchmark_report_u_plot_one <- function(rw, method, component, spec, integrated = FALSE,
                                         dataset = NULL, color = NULL,
                                         group = NULL, line_group = NULL,
                                         verbose = TRUE) {
  var <- .benchmark_report_spec_value(spec, "var", .benchmark_report_spec_value(spec, "x", NULL))
  if (is.null(var)) return(NULL)
  data_use <- if (isTRUE(integrated)) rw$xdata else rw$data
  typ <- .benchmark_report_spec_value(spec, "type", .benchmark_report_numeric_or_factor(data_use, var))
  color <- .benchmark_report_spec_value(spec, "color", color)
  group <- .benchmark_report_spec_value(spec, "group", group)
  line_group_default <- if (.benchmark_report_has_spec_field(spec, "group")) group else line_group
  line_group <- .benchmark_report_spec_value(spec, "line_group", line_group_default)
  connect <- if (.benchmark_report_has_spec_field(spec, "connect")) {
    isTRUE(spec$connect)
  } else {
    identical(typ, "numeric") && !is.null(line_group)
  }

  if (isTRUE(integrated)) {
    plot_real_world_integrated_U(
      rw, method = method, component = component, x = var,
      color = color, dataset = dataset,
      connect = connect, line_group = line_group
    )
  } else {
    plot_real_world_U(
      rw, method = method, component = component, x = var,
      color = color,
      connect = connect, line_group = line_group
    )
  }
}

.benchmark_report_save_u_panels <- function(rw, outdir, methods, components, u_plot_specs,
                                            integrated = FALSE, dataset = NULL,
                                            color = NULL, group = NULL, line_group = NULL, plots_per_method = 5L,
                                            verbose = TRUE,
                                            figure_format = "pdf") {
  figure_format <- .benchmark_report_validate_figure_format(figure_format)
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(character(0))
  files <- character(0)
  for (method in methods) {
    method_components <- .benchmark_report_filter_components_for_method(
      rw, method = method, components = components, integrated = integrated, dataset = dataset
    )
    if (!length(method_components)) {
      .benchmark_report_message("Skipping U panels for ", method, ": no requested components are available.", verbose = verbose)
      next
    }
    for (spec_i in seq_along(u_plot_specs)) {
      spec <- u_plot_specs[[spec_i]]
      var <- .benchmark_report_spec_value(spec, "var", .benchmark_report_spec_value(spec, "x", paste0("spec", spec_i)))
      for (component in method_components) {
        file <- .benchmark_report_plot_file(outdir, sprintf("U_%s_%s_component_%s", method, var, component), figure_format)
        files <- c(files, .benchmark_report_safe_plot(
          .benchmark_report_u_plot_one(rw, method = method, component = component, spec = spec,
                                       integrated = integrated, dataset = dataset,
                                       color = color, group = group, line_group = line_group, verbose = verbose),
          file = file, width = 8, height = 5, verbose = verbose
        ))
      }
    }
  }
  files[!is.na(files)]
}

.benchmark_report_method_v_result <- function(rw, method, integrated = FALSE, dataset = NULL) {
  res <- tryCatch(
    if (isTRUE(integrated)) .rwbench_get_integrated_method_result(rw, method) else .rwbench_get_method_result(rw, method),
    error = function(e) NULL
  )
  if (is.null(res)) return(NULL)

  if (isTRUE(integrated)) {
    datasets <- names(rw$Y) %||% names(res$V_list) %||% names(res$feature_score_list)
    d <- if (is.numeric(dataset)) dataset else match(dataset, datasets)
    if (is.na(d) || d < 1L) return(NULL)
    V <- res$V_list[[d]] %||% NULL
    scores <- res$feature_score_list[[d]] %||% NULL
    Y <- rw$Y[[d]] %||% NULL
    dataset_name <- datasets[[d]] %||% as.character(dataset)
  } else {
    V <- res$V %||% NULL
    scores <- res$feature_score %||% NULL
    Y <- rw$Y %||% NULL
    dataset_name <- NULL
  }

  if (is.null(V)) return(NULL)
  V <- tryCatch(as.matrix(V), error = function(e) NULL)
  if (is.null(V) || length(dim(V)) != 2L || nrow(V) < 1L || ncol(V) < 1L) return(NULL)

  feature_names <- rownames(V) %||% rownames(Y) %||% paste0("feature_", seq_len(nrow(V)))
  if (length(feature_names) != nrow(V)) feature_names <- paste0("feature_", seq_len(nrow(V)))
  rownames(V) <- feature_names

  if (!is.null(scores)) {
    scores <- as.numeric(scores)
    if (length(scores) == nrow(V)) names(scores) <- feature_names else scores <- NULL
  }
  if (is.null(scores)) {
    scores <- apply(abs(V), 1L, max, na.rm = TRUE)
    names(scores) <- feature_names
  }

  list(V = V, feature_score = scores, Y = Y, dataset = dataset_name, feature_names = feature_names)
}


.benchmark_report_write_tsv <- function(df, file) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(NULL)
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(df, file = file, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  file
}

.benchmark_report_loading_table <- function(
    rw,
    method,
    component = NULL,
    integrated = FALSE,
    dataset = NULL
) {
  vr <- .benchmark_report_method_v_result(rw, method, integrated = integrated, dataset = dataset)
  if (is.null(vr) || is.null(vr$V)) return(NULL)
  V <- as.matrix(vr$V)
  if (!nrow(V) || !ncol(V)) return(NULL)

  feature_id <- rownames(V) %||% vr$feature_names %||% paste0("feature_", seq_len(nrow(V)))
  feature_name <- .benchmark_report_generic_feature_labels(feature_id, fdata = rw$fdata, dataset = dataset)
  component_names <- colnames(V)
  if (is.null(component_names) || length(component_names) != ncol(V) || any(!nzchar(component_names))) {
    component_names <- paste0("component_", seq_len(ncol(V)))
  }

  if (is.null(component)) {
    ranking_statistic <- vr$feature_score
    if (is.null(ranking_statistic) || length(ranking_statistic) != nrow(V)) {
      ranking_statistic <- apply(abs(V), 1L, max, na.rm = TRUE)
    }
    ranking_statistic <- as.numeric(ranking_statistic)
    best_idx <- apply(abs(V), 1L, function(x) if (all(is.na(x))) NA_integer_ else which.max(x))
    best_component <- rep(NA_character_, nrow(V))
    best_loading <- rep(NA_real_, nrow(V))
    ok <- !is.na(best_idx)
    best_component[ok] <- component_names[best_idx[ok]]
    best_loading[ok] <- V[cbind(which(ok), best_idx[ok])]
    out <- data.frame(
      method = method,
      feature_id = feature_id,
      feature_name = feature_name,
      ranking_statistic = ranking_statistic,
      best_component = best_component,
      best_loading = best_loading,
      stringsAsFactors = FALSE
    )
    if (!is.null(dataset)) out$dataset <- dataset
    Vdf <- as.data.frame(V, stringsAsFactors = FALSE)
    names(Vdf) <- paste0("loading_", .benchmark_report_safe_name(component_names))
    out <- cbind(out, Vdf)
  } else {
    component <- suppressWarnings(as.integer(component))
    if (is.na(component) || component < 1L || component > ncol(V)) return(NULL)
    loading <- as.numeric(V[, component])
    out <- data.frame(
      method = method,
      feature_id = feature_id,
      feature_name = feature_name,
      component = component_names[[component]],
      loading = loading,
      ranking_statistic = abs(loading),
      stringsAsFactors = FALSE
    )
    if (!is.null(dataset)) out$dataset <- dataset
  }

  ord <- order(out$ranking_statistic, decreasing = TRUE, na.last = TRUE)
  out <- out[ord, , drop = FALSE]
  out$rank <- seq_len(nrow(out))
  out
}

.benchmark_report_write_loading_tables <- function(
    rw,
    outdir,
    methods,
    components = NULL,
    integrated = FALSE,
    datasets = NULL,
    verbose = TRUE
) {
  files <- character(0)
  datasets <- if (isTRUE(integrated)) datasets %||% names(rw$Y) else NULL

  for (method in methods) {
    ds_values <- if (isTRUE(integrated)) datasets else list(NULL)
    if (!isTRUE(integrated)) ds_values <- list(NULL)
    for (ds in ds_values) {
      method_components <- .benchmark_report_filter_loading_components_for_method(
        rw, method = method, components = components, integrated = integrated, dataset = ds
      )
      vr <- .benchmark_report_method_v_result(rw, method, integrated = integrated, dataset = ds)
      if (is.null(vr)) next

      method_tag <- .benchmark_report_safe_name(method)
      ds_tag <- if (!is.null(ds)) paste0("_", .benchmark_report_safe_name(ds)) else ""

      union_tbl <- .benchmark_report_loading_table(rw, method, component = NULL, integrated = integrated, dataset = ds)
      if (!is.null(union_tbl)) {
        f <- file.path(outdir, "V", "loadings", "union",
                       paste0(method_tag, ds_tag, "_feature_loadings_union.tsv"))
        files <- c(files, .benchmark_report_write_tsv(union_tbl, f))
      }

      for (cc in method_components) {
        comp_tbl <- .benchmark_report_loading_table(rw, method, component = cc, integrated = integrated, dataset = ds)
        if (is.null(comp_tbl)) next
        comp_tag <- .benchmark_report_safe_name(unique(comp_tbl$component)[1])
        f <- file.path(outdir, "V", "loadings", "component",
                       paste0(method_tag, ds_tag, "_feature_loadings_", comp_tag, ".tsv"))
        files <- c(files, .benchmark_report_write_tsv(comp_tbl, f))
      }
    }
  }
  files[!is.na(files)]
}

.benchmark_report_inference_q_column <- function(df, scope) {
  candidates <- switch(scope,
    global = c("global_q", "global_padj", "q_value", "qval", "padj", "p.adjust", "FDR"),
    union = c("union_q", "q_value", "qval", "padj", "p.adjust", "FDR"),
    component = c("q_value", "qval", "q", "padj", "p.adjust", "FDR"),
    character(0)
  )
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) hit[[1]] else NULL
}

.benchmark_report_inference_p_column <- function(df, scope) {
  candidates <- switch(scope,
    global = c("global_p", "p_value", "pval", "p", "pval_raw"),
    union = c("union_p", "p_value", "pval", "p", "pval_raw"),
    component = c("p_value", "pval", "p", "pval_raw"),
    character(0)
  )
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) hit[[1]] else NULL
}

.benchmark_report_prepare_inference_table <- function(df, scope, rw, method, dataset = NULL) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(NULL)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (!"feature" %in% names(df)) return(NULL)

  df$feature_name <- .benchmark_report_generic_feature_labels(df$feature, fdata = rw$fdata, dataset = dataset)
  if (!"method" %in% names(df)) df$method <- method
  if (!is.null(dataset) && !"dataset" %in% names(df)) df$dataset <- dataset

  q_col <- .benchmark_report_inference_q_column(df, scope)
  if (is.null(q_col)) {
    p_col <- .benchmark_report_inference_p_column(df, scope)
    if (!is.null(p_col)) {
      df$export_q_value <- stats::p.adjust(as.numeric(df[[p_col]]), method = "BH")
      q_col <- "export_q_value"
    }
  }
  if (!is.null(q_col)) {
    alpha <- rw$alpha %||% 0.1
    q_value <- as.numeric(df[[q_col]])
    df$significant_fdr <- !is.na(q_value) & q_value <= alpha
    df$fdr_threshold <- alpha
    df <- df[order(q_value, na.last = TRUE), , drop = FALSE]
  }
  df$rank <- seq_len(nrow(df))
  df
}

.benchmark_report_inference_objects <- function(rw, method) {
  obj <- rw$inference[[method]] %||% NULL
  if (is.null(obj) || is.null(obj$inf)) return(list())
  inf <- obj$inf
  is_inference_object <- function(x) {
    is.list(x) && (!is.null(x$feature_results) || !is.null(x$union_results))
  }
  if (is_inference_object(inf)) return(list(`__all__` = inf))
  if (is.list(inf) && length(inf)) {
    keep <- vapply(inf, is_inference_object, logical(1))
    if (any(keep)) return(inf[keep])
  }
  list()
}

.benchmark_report_write_inference_tables <- function(
    rw,
    outdir,
    methods,
    inference_scope = c("all", "component", "union", "global"),
    verbose = TRUE
) {
  scopes <- .benchmark_report_expand_inference_scope(match.arg(inference_scope))
  files <- character(0)

  for (method in methods) {
    inf_objects <- .benchmark_report_inference_objects(rw, method)
    if (!length(inf_objects)) next
    method_tag <- .benchmark_report_safe_name(method)

    for (obj_name in names(inf_objects)) {
      inf <- inf_objects[[obj_name]]
      dataset_hint <- if (!identical(obj_name, "__all__")) obj_name else NULL

      if ("union" %in% scopes || "global" %in% scopes) {
        ur <- tryCatch(.bench_extract_union_inference(inf), error = function(e) NULL)
        if (!is.null(ur) && nrow(ur)) {
          ds_values <- if ("dataset" %in% names(ur)) unique(as.character(ur$dataset)) else dataset_hint
          if (is.null(ds_values) || !length(ds_values)) ds_values <- NA_character_
          for (ds in ds_values) {
            ur_ds <- if (!is.na(ds) && "dataset" %in% names(ur)) ur[as.character(ur$dataset) == ds, , drop = FALSE] else ur
            ds_use <- if (is.na(ds)) dataset_hint else ds
            ds_tag <- if (!is.null(ds_use)) paste0("_", .benchmark_report_safe_name(ds_use)) else ""

            if ("union" %in% scopes) {
              tbl <- .benchmark_report_prepare_inference_table(ur_ds, "union", rw, method, ds_use)
              if (!is.null(tbl)) {
                f <- file.path(outdir, "V", "inference", "union",
                               paste0(method_tag, ds_tag, "_feature_inference_union.tsv"))
                files <- c(files, .benchmark_report_write_tsv(tbl, f))
              }
            }
            if ("global" %in% scopes && any(c("global_q", "global_p", "global_stat") %in% names(ur_ds))) {
              tbl <- .benchmark_report_prepare_inference_table(ur_ds, "global", rw, method, ds_use)
              if (!is.null(tbl)) {
                f <- file.path(outdir, "V", "inference", "global",
                               paste0(method_tag, ds_tag, "_feature_inference_global.tsv"))
                files <- c(files, .benchmark_report_write_tsv(tbl, f))
              }
            }
          }
        }
      }

      if ("component" %in% scopes) {
        fr <- inf$feature_results %||% NULL
        if (!is.null(fr) && is.data.frame(fr) && all(c("feature", "component") %in% names(fr))) {
          fr <- as.data.frame(fr, stringsAsFactors = FALSE)
          ds_values <- if ("dataset" %in% names(fr)) unique(as.character(fr$dataset)) else dataset_hint
          if (is.null(ds_values) || !length(ds_values)) ds_values <- NA_character_
          for (ds in ds_values) {
            fr_ds <- if (!is.na(ds) && "dataset" %in% names(fr)) fr[as.character(fr$dataset) == ds, , drop = FALSE] else fr
            ds_use <- if (is.na(ds)) dataset_hint else ds
            ds_tag <- if (!is.null(ds_use)) paste0("_", .benchmark_report_safe_name(ds_use)) else ""
            comps <- unique(fr_ds$component)
            comps <- comps[!is.na(comps)]
            for (comp in comps) {
              comp_rows <- fr_ds[as.character(fr_ds$component) == as.character(comp), , drop = FALSE]
              tbl <- .benchmark_report_prepare_inference_table(comp_rows, "component", rw, method, ds_use)
              if (is.null(tbl)) next
              comp_tag <- .benchmark_report_safe_name(comp)
              f <- file.path(outdir, "V", "inference", "component",
                             paste0(method_tag, ds_tag, "_feature_inference_", comp_tag, ".tsv"))
              files <- c(files, .benchmark_report_write_tsv(tbl, f))
            }
          }
        }
      }
    }
  }
  files[!is.na(files)]
}

.benchmark_report_filter_loading_components_for_method <- function(
    rw,
    method,
    components,
    integrated = FALSE,
    dataset = NULL
) {
  vr <- .benchmark_report_method_v_result(rw, method, integrated = integrated, dataset = dataset)
  if (is.null(vr) || is.null(vr$V)) return(integer(0))
  K <- ncol(vr$V)
  if (!is.finite(K) || K < 1L) return(integer(0))
  if (is.null(components)) return(seq_len(K))
  components <- suppressWarnings(as.integer(components))
  components <- components[!is.na(components) & components >= 1L]
  if (!length(components)) return(integer(0))
  components[components <= K]
}

.benchmark_report_generic_loading_features <- function(
    rw,
    method,
    component = NULL,
    n_features = 50L,
    integrated = FALSE,
    dataset = NULL
) {
  vr <- .benchmark_report_method_v_result(rw, method, integrated = integrated, dataset = dataset)
  if (is.null(vr)) return(character(0))
  V <- vr$V
  if (!is.null(component)) {
    component <- suppressWarnings(as.integer(component))
    if (is.na(component) || component < 1L || component > ncol(V)) return(character(0))
    score <- abs(V[, component])
    names(score) <- rownames(V)
  } else {
    score <- vr$feature_score
  }
  score <- score[is.finite(score)]
  if (!length(score)) return(character(0))
  names(sort(score, decreasing = TRUE))[seq_len(min(as.integer(n_features), length(score)))]
}

.benchmark_report_generic_feature_labels <- function(features, fdata = NULL, dataset = NULL) {
  features <- as.character(features)
  labels <- features
  if (is.null(fdata) || !length(features)) return(labels)

  # Respect the semantic feature/label columns recorded by as_svdx_fdata().
  # For integrated analyses, first resolve the dataset-specific fdata object.
  resolve_fdata <- .benchmark_report_find_svdx_function(".svdx_resolve_fdata_for_dataset", required = FALSE)
  to_row_labels <- .benchmark_report_find_svdx_function(".svdx_fdata_to_row_labels", required = FALSE)
  fd_use <- fdata
  if (!is.null(resolve_fdata)) {
    fd_use <- tryCatch(
      resolve_fdata(fdata, dataset = dataset, required = FALSE),
      error = function(e) fdata
    )
  } else if (is.list(fdata) && !is.data.frame(fdata) && !inherits(fdata, "svdx_fdata")) {
    if (!is.null(dataset) && dataset %in% names(fdata)) fd_use <- fdata[[dataset]]
  }

  if (!is.null(to_row_labels)) {
    row_labels <- tryCatch(to_row_labels(fd_use), error = function(e) NULL)
    if (!is.null(row_labels) && all(c("feature", "label") %in% colnames(row_labels))) {
      idx <- match(features, as.character(row_labels$feature))
      repl <- as.character(row_labels$label[idx])
      ok <- !is.na(idx) & !is.na(repl) & nzchar(repl)
      labels[ok] <- repl[ok]
      return(labels)
    }
  }

  # Backward-compatible fallback for plain, non-svdx fdata frames.
  fd <- tryCatch(as.data.frame(fd_use), error = function(e) NULL)
  if (is.null(fd) || !nrow(fd)) return(labels)
  key <- NULL
  for (nm in c("feature", "feature_id", "id", "gene", "gene_id", "name")) {
    if (nm %in% names(fd)) { key <- nm; break }
  }
  lab_col <- NULL
  for (nm in c("label", "feature_label", "gene_name", "symbol", "term_name", "name")) {
    if (nm %in% names(fd) && !identical(nm, key)) { lab_col <- nm; break }
  }
  if (is.null(key) && !is.null(rownames(fd))) {
    idx <- match(features, rownames(fd))
  } else if (!is.null(key)) {
    idx <- match(features, as.character(fd[[key]]))
  } else {
    return(labels)
  }
  if (!is.null(lab_col)) {
    repl <- as.character(fd[[lab_col]][idx])
    ok <- !is.na(idx) & !is.na(repl) & nzchar(repl)
    labels[ok] <- repl[ok]
  }
  labels
}

.benchmark_report_plot_generic_loading_bars <- function(
    rw,
    method,
    component,
    n_features = 50L,
    integrated = FALSE,
    dataset = NULL,
    fdata = NULL
) {
  vr <- .benchmark_report_method_v_result(rw, method, integrated = integrated, dataset = dataset)
  if (is.null(vr)) stop("Selected method does not contain a usable V/loading matrix.", call. = FALSE)
  V <- vr$V
  component <- suppressWarnings(as.integer(component))
  if (is.na(component) || component < 1L || component > ncol(V)) stop("Invalid component.", call. = FALSE)
  score <- V[, component]
  ord <- order(abs(score), decreasing = TRUE, na.last = NA)
  ord <- head(ord, min(as.integer(n_features), length(ord)))
  df <- data.frame(
    feature = rownames(V)[ord],
    loading = as.numeric(score[ord]),
    abs_loading = abs(as.numeric(score[ord])),
    stringsAsFactors = FALSE
  )
  df$label <- .benchmark_report_generic_feature_labels(df$feature, fdata = fdata, dataset = vr$dataset)
  label_lookup <- stats::setNames(as.character(df$label), df$feature)
  df$feature <- factor(df$feature, levels = rev(df$feature))
  title <- paste(method, "component", component, "top loadings")
  if (!is.null(vr$dataset)) title <- paste(method, vr$dataset, "component", component, "top loadings")
  ggplot2::ggplot(df, ggplot2::aes(x = feature, y = loading)) +
    ggplot2::geom_col() +
    ggplot2::scale_x_discrete(labels = label_lookup) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = NULL, y = "Loading", title = title)
}

.benchmark_report_plot_generic_feature_profiles <- function(
    Y,
    data,
    features,
    x,
    method,
    scope,
    component = NULL,
    dataset = NULL,
    fdata = NULL,
    top_n = 5L,
    summarize = c("mean", "median", "none"),
    show_se = TRUE,
    se_width = 0.15,
    scale = "zscore",
    color = "feature",
    facet = "feature",
    line = TRUE,
    points = TRUE
) {
  summarize <- match.arg(summarize)
  if (is.null(x) || !x %in% colnames(data)) stop("x column not found in metadata.", call. = FALSE)
  x_values <- data[[x]]
  x_is_discrete <- is.factor(x_values) || is.character(x_values) || is.logical(x_values)
  Y <- as.matrix(Y)
  features <- as.character(features)
  features <- features[features %in% rownames(Y)]
  features <- head(features, min(as.integer(top_n), length(features)))
  if (!length(features)) stop("No requested loading-selected features were found in Y.", call. = FALSE)
  mat <- Y[features, , drop = FALSE]
  if (!is.null(scale) && identical(scale, "zscore")) {
    mat <- t(scale(t(mat)))
    mat[!is.finite(mat)] <- 0
  }
  labels <- .benchmark_report_generic_feature_labels(features, fdata = fdata, dataset = dataset)
  rownames(mat) <- labels
  df <- data.frame(
    sample = rep(seq_len(ncol(mat)), each = nrow(mat)),
    feature = rep(rownames(mat), times = ncol(mat)),
    value = as.vector(mat),
    stringsAsFactors = FALSE
  )
  meta <- as.data.frame(data)
  df$x <- meta[[x]][df$sample]
  for (nm in colnames(meta)) {
    if (!nm %in% names(df)) df[[nm]] <- meta[[nm]][df$sample]
  }

  color_var <- if (!is.null(color) && length(color) == 1L && color %in% names(df)) color else "feature"
  facet_var <- if (!is.null(facet) && length(facet) == 1L && facet %in% names(df)) facet else NULL

  title <- paste(method, "loading-selected feature profiles", scope)
  if (!is.null(component)) title <- paste(title, "component", component)
  if (!is.null(dataset)) title <- paste(title, dataset)

  group_vars <- unique(c("feature", color_var, facet_var))
  group_vars <- group_vars[!is.na(group_vars) & nzchar(group_vars) & group_vars %in% names(df)]
  df$.profile_group <- interaction(df[group_vars], drop = TRUE, lex.order = TRUE)

  if (identical(summarize, "none")) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = value, colour = .data[[color_var]], group = .profile_group))
    if (isTRUE(line)) p <- p + ggplot2::geom_line(alpha = 0.45)
    if (isTRUE(points)) p <- p + ggplot2::geom_point(alpha = 0.55)
    p <- p + ggplot2::theme_bw() +
      ggplot2::labs(x = x, y = if (identical(scale, "zscore")) "Scaled expression" else "Expression", colour = color_var, title = title)
    if (!is.null(facet_var)) p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)), scales = "free_y")
    if (isTRUE(x_is_discrete)) {
      p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
    }
    return(p)
  }

  center_fun <- if (identical(summarize, "median")) stats::median else mean
  by_cols <- unique(c("feature", "x", color_var, facet_var))
  by_cols <- by_cols[by_cols %in% names(df)]
  centers <- stats::aggregate(df$value, by = df[by_cols], FUN = center_fun, na.rm = TRUE)
  names(centers)[ncol(centers)] <- "center"
  ses <- stats::aggregate(df$value, by = df[by_cols], FUN = function(z) stats::sd(z, na.rm = TRUE) / sqrt(sum(is.finite(z))))
  names(ses)[ncol(ses)] <- "se"
  pd <- merge(centers, ses, by = by_cols, all.x = TRUE)
  pd$.profile_group <- interaction(pd[group_vars], drop = TRUE, lex.order = TRUE)

  p <- ggplot2::ggplot(pd, ggplot2::aes(x = x, y = center, colour = .data[[color_var]], group = .profile_group)) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = x, y = paste0(summarize, if (identical(scale, "zscore")) " scaled expression" else " expression"), colour = color_var, title = title)
  if (isTRUE(line)) p <- p + ggplot2::geom_line()
  if (isTRUE(points)) p <- p + ggplot2::geom_point()
  if (isTRUE(show_se)) {
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = center - se, ymax = center + se), width = se_width, alpha = 0.7)
  }
  if (!is.null(facet_var)) p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)), scales = "free_y")
  if (isTRUE(x_is_discrete)) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  p
}


.benchmark_report_resolve_annotation_colors <- function(df, reference_df = NULL) {
  resolver <- .benchmark_report_find_svdx_function(".resolve_annotation_colors", required = FALSE)
  if (!is.null(resolver)) {
    return(resolver(df, reference_df = reference_df))
  }
  if (is.null(df) || ncol(df) == 0L) return(NULL)
  if (is.null(reference_df)) reference_df <- df
  out <- list()
  for (nm in colnames(df)) {
    x <- df[[nm]]
    ref <- if (nm %in% colnames(reference_df)) reference_df[[nm]] else x
    if (is.numeric(x) || is.integer(x)) {
      rng <- range(ref, na.rm = TRUE)
      if (length(rng) != 2L || any(!is.finite(rng))) next
      if (identical(rng[1], rng[2])) {
        delta <- if (rng[1] == 0) 0.5 else abs(rng[1]) * 0.05
        rng <- rng + c(-delta, delta)
      }
      if (requireNamespace("circlize", quietly = TRUE)) {
        out[[nm]] <- circlize::colorRamp2(rng, grDevices::hcl.colors(2L, palette = "BluYl"))
      }
      next
    }
    lev <- if (is.factor(ref)) levels(ref) else sort(unique(as.character(ref[!is.na(ref)])))
    if (length(lev) == 0L) next
    out[[nm]] <- stats::setNames(grDevices::hcl.colors(length(lev), palette = "Dark 3"), lev)
  }
  if (length(out) == 0L) NULL else out
}

.benchmark_report_plot_generic_loading_heatmap <- function(
    rw,
    method,
    component = NULL,
    n_features = 50L,
    features = NULL,
    integrated = FALSE,
    dataset = NULL,
    order_by = NULL,
    annotate = NULL,
    scale_rows = TRUE,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_feature_names = TRUE,
    show_sample_names = FALSE,
    fdata = NULL
) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required for generic loading heatmaps.", call. = FALSE)
  }
  prepare_heatmap_matrix <- .benchmark_report_find_svdx_function(".prepare_heatmap_matrix", required = TRUE)
  vr <- .benchmark_report_method_v_result(rw, method = method, integrated = integrated, dataset = dataset)
  if (is.null(vr)) stop("Selected method does not contain a usable V/loading matrix.", call. = FALSE)

  Y <- vr$Y
  if (is.null(Y)) stop("Selected method does not contain an associated Y matrix.", call. = FALSE)
  Y <- as.matrix(Y)
  feature_names <- rownames(Y) %||% vr$feature_names %||% paste0("feature_", seq_len(nrow(Y)))
  if (length(feature_names) != nrow(Y)) feature_names <- paste0("feature_", seq_len(nrow(Y)))
  rownames(Y) <- feature_names

  if (is.null(features)) {
    features <- .benchmark_report_generic_loading_features(
      rw,
      method = method,
      component = component,
      n_features = n_features,
      integrated = integrated,
      dataset = dataset
    )
  }
  features <- as.character(features)
  features <- features[features %in% feature_names]
  features <- head(features, min(as.integer(n_features), length(features)))
  if (!length(features)) stop("No requested loading-selected features were found in Y.", call. = FALSE)
  fidx <- match(features, feature_names)

  data <- if (isTRUE(integrated)) as.data.frame(rw$xdata) else as.data.frame(rw$data)
  if (nrow(data) != ncol(Y)) {
    stop("Generic loading heatmap received inconsistent data: nrow(metadata) != ncol(Y).", call. = FALSE)
  }
  if (is.null(annotate) && !is.null(order_by)) annotate <- order_by

  prep <- prepare_heatmap_matrix(
    Y = Y,
    data = data,
    features = fidx,
    cols = NULL,
    order_by = order_by,
    scale_rows = scale_rows
  )
  mat <- prep$Y
  rownames(mat) <- .benchmark_report_generic_feature_labels(features, fdata = fdata, dataset = vr$dataset)

  ha <- NULL
  if (!is.null(annotate)) {
    missing_annotate <- setdiff(annotate, colnames(prep$data))
    if (length(missing_annotate) > 0L) {
      stop("annotate variable(s) not found in data: ", paste(missing_annotate, collapse = ", "), call. = FALSE)
    }
    ann_df <- prep$data[, annotate, drop = FALSE]
    ann_ref_df <- data[, annotate, drop = FALSE]
    ha <- ComplexHeatmap::HeatmapAnnotation(
      df = ann_df,
      col = .benchmark_report_resolve_annotation_colors(ann_df, reference_df = ann_ref_df)
    )
  }

  title <- paste(method, "top loading features")
  if (!is.null(vr$dataset)) title <- paste(method, vr$dataset, "top loading features")
  if (!is.null(component)) title <- paste(title, "component", component)

  ComplexHeatmap::Heatmap(
    mat,
    name = ifelse(scale_rows, "z", "Y"),
    top_annotation = ha,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    show_row_names = show_feature_names,
    show_column_names = show_sample_names,
    column_title = title
  )
}

.benchmark_report_save_generic_loading_analyses <- function(
    rw,
    outdir,
    methods,
    components,
    integrated = FALSE,
    datasets = NULL,
    xaxis = NULL,
    features_top_n = 50L,
    feature_profiles_top_n = 5L,
    loading_heatmaps = TRUE,
    loading_feature_profiles = TRUE,
    feature_profile_summarize = c("mean", "median", "none"),
    feature_profile_show_se = TRUE,
    plot_args = list(),
    figure_format = "pdf",
    verbose = TRUE
) {
  feature_profile_summarize <- match.arg(feature_profile_summarize)
  figure_format <- .benchmark_report_validate_figure_format(figure_format)
  plot_args <- plot_args %||% list()
  heatmap_order_by <- plot_args$heatmap_order_by %||% NULL
  heatmap_annotate <- plot_args$heatmap_annotate %||% heatmap_order_by
  profile_color <- plot_args$feature_profile_color %||% "feature"
  profile_facet <- plot_args$feature_profile_facet %||% "feature"
  profile_se_width <- plot_args$feature_profile_se_width %||% 0.15
  profile_scale <- plot_args$feature_profile_scale %||% "zscore"
  files <- character(0)

  if (isTRUE(integrated)) {
    datasets <- datasets %||% names(rw$Y)
    meta <- as.data.frame(rw$xdata)
    for (method in methods) {
      for (ds in datasets) {
        method_components <- .benchmark_report_filter_loading_components_for_method(
          rw, method = method, components = components, integrated = TRUE, dataset = ds
        )
        if (!length(method_components)) {
          .benchmark_report_message("Skipping integrated benchmark loading analyses for ", method, " / ", ds, ": no requested components are available.", verbose = verbose)
          next
        }
        for (component in method_components) {
          stem <- sprintf("benchmark_v_loadings_%s_%s_component_%s", .benchmark_report_safe_name(method), .benchmark_report_safe_name(ds), component)
          files <- c(files, .benchmark_report_safe_plot(
            .benchmark_report_plot_generic_loading_bars(rw, method = method, component = component, integrated = TRUE, dataset = ds, n_features = features_top_n, fdata = rw$fdata),
            .benchmark_report_plot_file(outdir, stem, figure_format),
            width = 8, height = 6, verbose = verbose
          ))
        }

        if (isTRUE(loading_heatmaps)) {
          stem <- sprintf("benchmark_v_loading_heatmap_%s_%s_union", .benchmark_report_safe_name(method), .benchmark_report_safe_name(ds))
          files <- c(files, .benchmark_report_safe_plot(
            .benchmark_report_plot_generic_loading_heatmap(rw, method = method, dataset = ds, component = NULL, n_features = features_top_n, integrated = TRUE, order_by = heatmap_order_by, annotate = heatmap_annotate, fdata = rw$fdata),
            .benchmark_report_plot_file(outdir, stem, figure_format),
            width = 10, height = 9, verbose = verbose
          ))
          for (component in method_components) {
            stem <- sprintf("benchmark_v_loading_heatmap_%s_%s_component_%s", .benchmark_report_safe_name(method), .benchmark_report_safe_name(ds), component)
            files <- c(files, .benchmark_report_safe_plot(
              .benchmark_report_plot_generic_loading_heatmap(rw, method = method, dataset = ds, component = component, n_features = features_top_n, integrated = TRUE, order_by = heatmap_order_by, annotate = heatmap_annotate, fdata = rw$fdata),
              .benchmark_report_plot_file(outdir, stem, figure_format),
              width = 10, height = 9, verbose = verbose
            ))
          }
        }

        if (isTRUE(loading_feature_profiles) && !is.null(xaxis) && xaxis %in% colnames(meta)) {
          Y_ds <- rw$Y[[ds]]
          union_feats <- .benchmark_report_generic_loading_features(rw, method = method, component = NULL, n_features = feature_profiles_top_n, integrated = TRUE, dataset = ds)
          if (length(union_feats)) {
            stem <- sprintf("benchmark_feature_profiles_loading_%s_%s_union", .benchmark_report_safe_name(method), .benchmark_report_safe_name(ds))
            files <- c(files, .benchmark_report_safe_plot(
              .benchmark_report_plot_generic_feature_profiles(Y_ds, meta, union_feats, x = xaxis, method = method, scope = "union", dataset = ds, fdata = rw$fdata, top_n = feature_profiles_top_n, summarize = feature_profile_summarize, show_se = feature_profile_show_se, se_width = profile_se_width, scale = profile_scale, color = profile_color, facet = profile_facet, line = TRUE, points = TRUE),
              .benchmark_report_plot_file(outdir, stem, figure_format),
              width = 9, height = 6, verbose = verbose
            ))
          }
          for (component in method_components) {
            feats <- .benchmark_report_generic_loading_features(rw, method = method, component = component, n_features = feature_profiles_top_n, integrated = TRUE, dataset = ds)
            if (!length(feats)) next
            stem <- sprintf("benchmark_feature_profiles_loading_%s_%s_component_%s", .benchmark_report_safe_name(method), .benchmark_report_safe_name(ds), component)
            files <- c(files, .benchmark_report_safe_plot(
              .benchmark_report_plot_generic_feature_profiles(Y_ds, meta, feats, x = xaxis, method = method, scope = "component", component = component, dataset = ds, fdata = rw$fdata, top_n = feature_profiles_top_n, summarize = feature_profile_summarize, show_se = feature_profile_show_se, se_width = profile_se_width, scale = profile_scale, color = profile_color, facet = profile_facet, line = TRUE, points = TRUE),
              .benchmark_report_plot_file(outdir, stem, figure_format),
              width = 9, height = 6, verbose = verbose
            ))
          }
        }
      }
    }
  } else {
    meta <- as.data.frame(rw$data)
    for (method in methods) {
      method_components <- .benchmark_report_filter_loading_components_for_method(
        rw, method = method, components = components, integrated = FALSE
      )
      if (!length(method_components)) {
        .benchmark_report_message("Skipping benchmark loading analyses for ", method, ": no requested components are available.", verbose = verbose)
        next
      }
      for (component in method_components) {
        stem <- sprintf("benchmark_v_loadings_%s_component_%s", .benchmark_report_safe_name(method), component)
        files <- c(files, .benchmark_report_safe_plot(
          .benchmark_report_plot_generic_loading_bars(rw, method = method, component = component, integrated = FALSE, n_features = features_top_n, fdata = rw$fdata),
          .benchmark_report_plot_file(outdir, stem, figure_format),
          width = 8, height = 6, verbose = verbose
        ))
      }

      if (isTRUE(loading_heatmaps)) {
        stem <- sprintf("benchmark_v_loading_heatmap_%s_union", .benchmark_report_safe_name(method))
        files <- c(files, .benchmark_report_safe_plot(
          .benchmark_report_plot_generic_loading_heatmap(rw, method = method, component = NULL, n_features = features_top_n, integrated = FALSE, order_by = heatmap_order_by, annotate = heatmap_annotate, fdata = rw$fdata),
          .benchmark_report_plot_file(outdir, stem, figure_format),
          width = 10, height = 9, verbose = verbose
        ))
        for (component in method_components) {
          feats <- .benchmark_report_generic_loading_features(rw, method = method, component = component, n_features = features_top_n, integrated = FALSE)
          if (!length(feats)) next
          stem <- sprintf("benchmark_v_loading_heatmap_%s_component_%s", .benchmark_report_safe_name(method), component)
          files <- c(files, .benchmark_report_safe_plot(
            .benchmark_report_plot_generic_loading_heatmap(rw, method = method, features = feats, component = component, n_features = features_top_n, integrated = FALSE, order_by = heatmap_order_by, annotate = heatmap_annotate, fdata = rw$fdata),
            .benchmark_report_plot_file(outdir, stem, figure_format),
            width = 10, height = 9, verbose = verbose
          ))
        }
      }

      if (isTRUE(loading_feature_profiles) && !is.null(xaxis) && xaxis %in% colnames(meta)) {
        union_feats <- .benchmark_report_generic_loading_features(rw, method = method, component = NULL, n_features = feature_profiles_top_n, integrated = FALSE)
        if (length(union_feats)) {
          stem <- sprintf("benchmark_feature_profiles_loading_%s_union", .benchmark_report_safe_name(method))
          files <- c(files, .benchmark_report_safe_plot(
            .benchmark_report_plot_generic_feature_profiles(rw$Y, meta, union_feats, x = xaxis, method = method, scope = "union", fdata = rw$fdata, top_n = feature_profiles_top_n, summarize = feature_profile_summarize, show_se = feature_profile_show_se, se_width = profile_se_width, scale = profile_scale, color = profile_color, facet = profile_facet, line = TRUE, points = TRUE),
            .benchmark_report_plot_file(outdir, stem, figure_format),
            width = 9, height = 6, verbose = verbose
          ))
        }
        for (component in method_components) {
          feats <- .benchmark_report_generic_loading_features(rw, method = method, component = component, n_features = feature_profiles_top_n, integrated = FALSE)
          if (!length(feats)) next
          stem <- sprintf("benchmark_feature_profiles_loading_%s_component_%s", .benchmark_report_safe_name(method), component)
          files <- c(files, .benchmark_report_safe_plot(
            .benchmark_report_plot_generic_feature_profiles(rw$Y, meta, feats, x = xaxis, method = method, scope = "component", component = component, fdata = rw$fdata, top_n = feature_profiles_top_n, summarize = feature_profile_summarize, show_se = feature_profile_show_se, se_width = profile_se_width, scale = profile_scale, color = profile_color, facet = profile_facet, line = TRUE, points = TRUE),
            .benchmark_report_plot_file(outdir, stem, figure_format),
            width = 9, height = 6, verbose = verbose
          ))
        }
      }
    }
  }
  files[!is.na(files)]
}


.benchmark_report_get_svdx_package_function <- function(name, required = TRUE) {
  .benchmark_report_find_svdx_function(name, required = required)
}

.benchmark_report_first_xaxis <- function(xformula = NULL, data = NULL, u_plot_specs = NULL) {
  if (!is.null(u_plot_specs) && length(u_plot_specs) > 0L) {
    for (sp in u_plot_specs) {
      v <- sp$var %||% sp$x
      if (!is.null(v) && length(v) == 1L && nzchar(v)) return(v)
    }
  }
  vars <- if (!is.null(xformula)) all.vars(xformula) else character(0)
  if (!is.null(data)) vars <- intersect(vars, colnames(as.data.frame(data)))
  if (length(vars) > 0L) vars[[1L]] else NULL
}

.benchmark_report_nonnull_args <- function(args) {
  args[!vapply(args, is.null, logical(1))]
}

.benchmark_report_subset_args <- function(args, allowed) {
  if (is.null(args) || !length(args)) return(list())
  if (is.null(names(args))) return(list())
  args[names(args) %in% allowed]
}

.benchmark_report_svdx_plot_args <- function(plot_args, target = c("plotU", "plotUTerms", "plotUEmbedding")) {
  target <- match.arg(target)
  allowed <- switch(
    target,
    plotU = c("include", "y", "fixed", "n_grid", "observed", "point_alpha"),
    plotUTerms = c("terms", "y", "fixed", "include_intercept", "n_grid"),
    plotUEmbedding = c("include", "vary", "fixed", "n_grid", "observed", "point_alpha")
  )
  if (is.null(plot_args) || !length(plot_args)) return(list())
  if (is.null(names(plot_args))) return(list())
  if ("yaxis" %in% names(plot_args) && !"y" %in% names(plot_args)) {
    plot_args$y <- plot_args$yaxis
  }
  .benchmark_report_subset_args(plot_args, allowed)
}

.benchmark_report_row_offset_for_Y <- function(Y, row_offset = NULL) {
  if (is.null(Y)) return(NULL)
  if (!is.null(row_offset)) {
    if (length(row_offset) != nrow(Y)) {
      stop("row_offset length must equal nrow(Y) for package report delegation.", call. = FALSE)
    }
    return(as.numeric(row_offset))
  }
  rep(0, nrow(Y))
}

.benchmark_report_col_offset_for_Y <- function(Y, col_offset = NULL) {
  if (is.null(Y)) return(NULL)
  if (!is.null(col_offset)) {
    if (length(col_offset) != ncol(Y)) {
      stop("col_offset length must equal ncol(Y) for package report delegation.", call. = FALSE)
    }
    return(as.numeric(col_offset))
  }
  rep(0, ncol(Y))
}

.benchmark_report_row_offset_for_Y_list <- function(Y, row_offset = NULL) {
  if (is.null(Y)) return(NULL)
  if (!is.list(Y)) stop("Integrated package report delegation requires Y to be a list.", call. = FALSE)
  if (!is.null(row_offset)) {
    if (!is.list(row_offset)) stop("row_offset must be a list for integrated reports.", call. = FALSE)
    out <- row_offset
  } else {
    out <- lapply(Y, function(y) rep(0, nrow(y)))
  }
  if (is.null(names(out)) && !is.null(names(Y))) names(out) <- names(Y)
  for (nm in names(Y) %||% seq_along(Y)) {
    yy <- Y[[nm]]
    rr <- out[[nm]]
    if (length(rr) != nrow(yy)) {
      stop("Each integrated row_offset length must equal nrow(Y[[dataset]]).", call. = FALSE)
    }
    out[[nm]] <- as.numeric(rr)
  }
  out
}

.benchmark_report_col_offset_for_Y_list <- function(Y, col_offset = NULL) {
  if (is.null(Y)) return(NULL)
  if (!is.list(Y)) stop("Integrated package report delegation requires Y to be a list.", call. = FALSE)
  if (!is.null(col_offset)) {
    if (!is.list(col_offset)) stop("col_offset must be a list for integrated reports.", call. = FALSE)
    out <- col_offset
  } else {
    out <- lapply(Y, function(y) rep(0, ncol(y)))
  }
  if (is.null(names(out)) && !is.null(names(Y))) names(out) <- names(Y)
  for (nm in names(Y) %||% seq_along(Y)) {
    yy <- Y[[nm]]
    cc <- out[[nm]]
    if (length(cc) != ncol(yy)) {
      stop("Each integrated col_offset length must equal ncol(Y[[dataset]]).", call. = FALSE)
    }
    out[[nm]] <- as.numeric(cc)
  }
  out
}

.benchmark_report_collect_files <- function(outdir) {
  if (!dir.exists(outdir)) return(character(0))
  list.files(outdir, recursive = TRUE, full.names = TRUE)
}


.benchmark_report_extract_svdx_fit <- function(rw) {
  raw <- rw$raw_fits$svdx %||% NULL
  if (!is.null(raw$fit)) return(raw$fit)
  inf <- rw$inference$svdx %||% NULL
  if (!is.null(inf$fit)) return(inf$fit)
  NULL
}

.benchmark_report_extract_svdx_inf <- function(rw) {
  inf <- rw$inference$svdx %||% NULL
  inf$inf %||% NULL
}

.benchmark_report_extract_svdx_inf_summary <- function(rw) {
  inf <- rw$inference$svdx %||% NULL
  inf$inf_summary %||% NULL
}

.benchmark_report_as_enrichment_object <- function(enrichment, method = NULL) {
  if (is.null(enrichment)) return(NULL)
  if (inherits(enrichment, "svdxEnrichment") && is.list(enrichment) && !is.null(enrichment$results)) {
    if (!is.null(method) && is.data.frame(enrichment$results) && "method" %in% names(enrichment$results)) {
      enrichment$results <- enrichment$results[enrichment$results$method == method, , drop = FALSE]
    }
    return(enrichment)
  }
  if (is.list(enrichment) && !is.data.frame(enrichment) && !is.null(enrichment$results)) {
    if (!is.null(method) && is.data.frame(enrichment$results) && "method" %in% names(enrichment$results)) {
      enrichment$results <- enrichment$results[enrichment$results$method == method, , drop = FALSE]
    }
    return(enrichment)
  }
  if (is.data.frame(enrichment)) {
    enrichment <- as.data.frame(enrichment)
    if (!is.null(method) && "method" %in% names(enrichment)) {
      enrichment <- enrichment[enrichment$method == method, , drop = FALSE]
    }
    return(list(results = enrichment))
  }
  NULL
}

.benchmark_report_extract_svdx_loading_enrichment <- function(rw) {
  load_enr <- rw$loading_enrichment %||% NULL
  # `loading_enrichment` may be a bare tibble/data.frame for single-dataset
  # loading enrichment. Do not probe `$svdx` on data frames/tibbles: tibble
  # emits an "Unknown or uninitialised column" warning, which report writing
  # converts to a skipped plot. Only dataset/method-index list objects should
  # be unwrapped by `svdx`.
  if (!is.null(load_enr) && !is.data.frame(load_enr) && is.list(load_enr)) {
    load_enr_svdx <- load_enr[["svdx"]]
    if (!is.null(load_enr_svdx)) load_enr <- load_enr_svdx
  }
  .benchmark_report_as_enrichment_object(load_enr, method = "svdx")
}

.benchmark_report_extract_svdx_inference_enrichment <- function(rw) {
  inf_enr <- rw$inference$svdx$inference_enrichment %||% NULL
  .benchmark_report_as_enrichment_object(inf_enr)
}

.benchmark_report_extract_svdx_enrichment <- function(rw) {
  # Backward-compatible workflow object field: prefer inference enrichment,
  # then fall back to loading enrichment if inference enrichment is absent.
  .benchmark_report_extract_svdx_inference_enrichment(rw) %||%
    .benchmark_report_extract_svdx_loading_enrichment(rw)
}

.benchmark_report_make_svdx_workflow <- function(rw, verbose = TRUE) {
  if (!inherits(rw, "realWorldMethodBenchmark")) return(NULL)
  fit <- .benchmark_report_extract_svdx_fit(rw)
  if (is.null(fit)) {
    .benchmark_report_message("Skipping SVDX-specific analyses: no stored SVDX fit found.", verbose = verbose)
    return(NULL)
  }
  if (is.null(rw$Y) || is.null(rw$data)) {
    .benchmark_report_message("Skipping SVDX-specific analyses: benchmark object lacks Y/data.", verbose = verbose)
    return(NULL)
  }
  if (ncol(as.matrix(rw$Y)) != nrow(as.data.frame(rw$data))) {
    .benchmark_report_message("Skipping SVDX-specific analyses: ncol(Y) does not equal nrow(data).", verbose = verbose)
    return(NULL)
  }
  zdata <- rw$zdata %||% rw$data
  if (!is.null(zdata) && nrow(as.data.frame(zdata)) != ncol(as.matrix(rw$Y))) {
    .benchmark_report_message("Skipping SVDX-specific analyses: nrow(zdata) does not equal ncol(Y).", verbose = verbose)
    return(NULL)
  }

  summary_svdx_fit_fn <- .benchmark_report_find_svdx_function("summary_svdx_fit", required = TRUE)
  fit_summary <- summary_svdx_fit_fn(
    fit = fit,
    Y = as.matrix(rw$Y)
  )
  wf <- list(
    split = rw$split,
    Y = as.matrix(rw$Y),
    data = as.data.frame(rw$data),
    zdata = if (!is.null(zdata)) as.data.frame(zdata) else NULL,
    Ytrain = as.matrix(rw$Y),
    Ytest = NULL,
    data_train = as.data.frame(rw$data),
    data_test = NULL,
    zdata_train = if (!is.null(zdata)) as.data.frame(zdata) else NULL,
    zdata_test = NULL,
    fit = fit,
    fit_summary = fit_summary,
    inf = .benchmark_report_extract_svdx_inf(rw),
    inf_summary = .benchmark_report_extract_svdx_inf_summary(rw),
    enrichment = .benchmark_report_extract_svdx_enrichment(rw),
    fdata = rw$fdata,
    xformula = rw$xformula,
    zformula = rw$zformula,
    modelType = rw$modelType,
    family = rw$family,
    fit_family = rw$family,
    inference_family = rw$family,
    fit_description = "benchmark-stored svdx fit",
    alpha = rw$alpha %||% 0.1,
    call = NULL
  )
  class(wf) <- "svdxWorkflow"
  wf
}

.benchmark_report_extract_isvdx_fit <- function(rw) {
  res <- rw$results$isvdx %||% NULL
  if (!is.null(res$fit)) return(res$fit)
  inf <- rw$inference$isvdx %||% NULL
  if (!is.null(inf$fit)) return(inf$fit)
  NULL
}

.benchmark_report_extract_isvdx_inf <- function(rw) {
  inf <- rw$inference$isvdx %||% NULL
  inf$inf %||% NULL
}

.benchmark_report_extract_isvdx_inf_summary <- function(rw) {
  inf <- rw$inference$isvdx %||% NULL
  inf$inf_summary %||% NULL
}

.benchmark_report_as_integrated_enrichment_object <- function(enrichment, method = NULL) {
  if (is.null(enrichment)) return(NULL)
  if (inherits(enrichment, "svdxEnrichment") && is.list(enrichment) && !is.null(enrichment$results)) {
    if (!is.null(method)) {
      if (is.data.frame(enrichment$results) && "method" %in% names(enrichment$results)) {
        enrichment$results <- enrichment$results[enrichment$results$method == method, , drop = FALSE]
      } else if (is.list(enrichment$results)) {
        enrichment$results <- lapply(enrichment$results, function(x) {
          if (is.data.frame(x) && "method" %in% names(x)) x[x$method == method, , drop = FALSE] else x
        })
      }
    }
    return(enrichment)
  }
  if (is.list(enrichment) && !is.data.frame(enrichment) && !is.null(enrichment$results)) {
    if (!is.null(method)) {
      if (is.data.frame(enrichment$results) && "method" %in% names(enrichment$results)) {
        enrichment$results <- enrichment$results[enrichment$results$method == method, , drop = FALSE]
      } else if (is.list(enrichment$results)) {
        enrichment$results <- lapply(enrichment$results, function(x) {
          if (is.data.frame(x) && "method" %in% names(x)) x[x$method == method, , drop = FALSE] else x
        })
      }
    }
    return(enrichment)
  }
  if (is.data.frame(enrichment)) {
    enrichment <- as.data.frame(enrichment)
    if (!is.null(method) && "method" %in% names(enrichment)) {
      enrichment <- enrichment[enrichment$method == method, , drop = FALSE]
    }
    if ("dataset" %in% names(enrichment)) {
      return(list(results = split(enrichment, enrichment$dataset)))
    }
    return(list(results = enrichment))
  }
  NULL
}

.benchmark_report_extract_isvdx_loading_enrichment <- function(rw) {
  load_enr <- rw$loading_enrichment %||% NULL
  .benchmark_report_as_integrated_enrichment_object(load_enr, method = "isvdx")
}

.benchmark_report_extract_isvdx_inference_enrichment <- function(rw) {
  inf_enr <- rw$inference$isvdx$inference_enrichment %||% NULL
  .benchmark_report_as_integrated_enrichment_object(inf_enr)
}

.benchmark_report_extract_isvdx_enrichment <- function(rw) {
  # Backward-compatible workflow object field: prefer inference enrichment,
  # then fall back to loading enrichment if inference enrichment is absent.
  .benchmark_report_extract_isvdx_inference_enrichment(rw) %||%
    .benchmark_report_extract_isvdx_loading_enrichment(rw)
}

.benchmark_report_make_isvdx_workflow <- function(rw, verbose = TRUE) {
  if (!inherits(rw, "realWorldIntegratedBenchmark")) return(NULL)
  fit <- .benchmark_report_extract_isvdx_fit(rw)
  if (is.null(fit)) {
    .benchmark_report_message("Skipping ISVDX-specific analyses: no stored ISVDX fit found.", verbose = verbose)
    return(NULL)
  }
  if (is.null(rw$Y) || !is.list(rw$Y) || is.null(rw$xdata)) {
    .benchmark_report_message("Skipping ISVDX-specific analyses: benchmark object lacks Y list/xdata.", verbose = verbose)
    return(NULL)
  }
  for (ds in names(rw$Y) %||% seq_along(rw$Y)) {
    if (ncol(as.matrix(rw$Y[[ds]])) != nrow(as.data.frame(rw$xdata))) {
      .benchmark_report_message("Skipping ISVDX-specific analyses: ncol(Y[[", ds, "]]) does not equal nrow(xdata).", verbose = verbose)
      return(NULL)
    }
  }
  fit_summary <- summary(
    fit,
    Y = rw$Y
  )
  wf <- list(
    split = rw$split,
    Y = rw$Y,
    data = as.data.frame(rw$xdata),
    zdata = rw$zdata,
    Ytrain = rw$Y,
    Ytest = NULL,
    data_train = as.data.frame(rw$xdata),
    data_test = NULL,
    zdata_train = rw$zdata,
    zdata_test = NULL,
    fit = fit,
    fit_summary = fit_summary,
    inf = .benchmark_report_extract_isvdx_inf(rw),
    inf_summary = .benchmark_report_extract_isvdx_inf_summary(rw),
    enrichment = .benchmark_report_extract_isvdx_enrichment(rw),
    fdata = rw$fdata,
    xformula = rw$xformula,
    zformula = rw$zformula,
    modelType = rw$modelType,
    family = rw$family,
    fit_family = rw$family,
    inference_family = rw$family,
    fit_description = "benchmark-stored isvdx fit",
    alpha = rw$alpha %||% 0.1,
    call = NULL
  )
  class(wf) <- "isvdxWorkflow"
  wf
}

.benchmark_report_save_package_workflow_plot <- function(wf, type, file, ..., width = 10, height = 7,
                                                         figure_format = NULL, verbose = TRUE) {
  .benchmark_report_safe_plot(
    plot(wf, type = type, ...),
    file = file,
    width = width,
    height = height,
    figure_format = figure_format,
    verbose = verbose
  )
}

.benchmark_report_component_count_from_fit <- function(fit, default = 1L) {
  U <- tryCatch(fit$fitAtOptCV$umat, error = function(e) NULL)
  if (!is.null(U)) return(ncol(as.matrix(U)))
  V <- tryCatch(fit$fitAtOptCV$vmat, error = function(e) NULL)
  if (!is.null(V)) return(ncol(as.matrix(V)))
  default
}


.benchmark_report_svdx_fn <- function(name, required = TRUE) {
  .benchmark_report_find_svdx_function(name, required = required)
}

.benchmark_report_expand_inference_scope <- function(inference_scope) {
  inference_scope <- match.arg(inference_scope, c("all", "component", "union", "global"))
  if (identical(inference_scope, "all")) c("component", "union", "global") else inference_scope
}


.benchmark_report_component_count_from_U <- function(U, dataset = NULL, default = 0L) {
  if (is.null(U)) return(as.integer(default))
  if (is.list(U) && !is.data.frame(U)) {
    if (!is.null(dataset) && !is.null(U[[dataset]])) {
      U <- U[[dataset]]
    } else if (length(U) > 0L) {
      U <- U[[1L]]
    } else {
      return(as.integer(default))
    }
  }
  U <- tryCatch(as.matrix(U), error = function(e) NULL)
  if (is.null(U) || length(dim(U)) != 2L) return(as.integer(default))
  as.integer(ncol(U))
}

.benchmark_report_method_component_count <- function(rw, method, integrated = FALSE, dataset = NULL, default = 0L) {
  res <- tryCatch(
    if (isTRUE(integrated)) .rwbench_get_integrated_method_result(rw, method) else .rwbench_get_method_result(rw, method),
    error = function(e) NULL
  )
  if (is.null(res)) return(as.integer(default))
  .benchmark_report_component_count_from_U(res$U, dataset = dataset, default = default)
}

.benchmark_report_filter_components_for_method <- function(rw, method, components, integrated = FALSE, dataset = NULL) {
  K <- .benchmark_report_method_component_count(rw, method, integrated = integrated, dataset = dataset, default = 0L)
  if (!is.finite(K) || K < 1L) return(integer(0))
  if (is.null(components)) return(seq_len(K))
  components <- suppressWarnings(as.integer(components))
  components <- components[!is.na(components) & components >= 1L]
  if (!length(components)) return(integer(0))
  components[components <= K]
}

.benchmark_report_inference_components <- function(inf, components = NULL) {
  fr <- inf$feature_results %||% NULL
  if (is.null(fr) || !is.data.frame(fr) || !"component" %in% names(fr)) return(integer(0))
  comps <- unique(fr$component)
  comps <- comps[!is.na(comps)]
  if (!length(comps)) return(integer(0))
  n <- length(comps)
  if (is.null(components)) return(seq_len(n))
  idx <- suppressWarnings(as.integer(components))
  idx <- idx[!is.na(idx) & idx >= 1L & idx <= n]
  unique(idx)
}

.benchmark_report_loading_features_from_fit <- function(fit, component, n_features = 5L, dataset = NULL) {
  selector <- .benchmark_report_svdx_fn(".select_features_from_loadings", required = FALSE)
  if (is.null(selector)) return(NULL)
  mat <- tryCatch({
    get_vmat <- .benchmark_report_svdx_fn(".get_vmat", required = FALSE)
    if (is.null(get_vmat)) NULL else get_vmat(fit, dataset = dataset)
  }, error = function(e) NULL)
  idx <- tryCatch(
    selector(fit, component = component, n_features = n_features, dataset = dataset),
    error = function(e) integer(0)
  )
  if (is.null(mat) || !length(idx)) return(NULL)
  rownames(mat)[idx]
}

.benchmark_report_loading_union_features_from_fit <- function(fit, components, n_features = 50L, dataset = NULL) {
  union_selector <- .benchmark_report_svdx_fn(".report_loading_union_features", required = FALSE)
  if (is.null(union_selector)) return(NULL)

  out <- tryCatch(
    union_selector(fit, components = components, n_features = n_features, dataset = dataset),
    error = function(e) NULL
  )
  if (is.null(out) || !length(out)) return(NULL)
  as.character(out)
}

.benchmark_report_normalize_enrichment_for_dotplot <- function(enrichment) {
  if (is.null(enrichment)) return(NULL)
  if (inherits(enrichment, "svdxEnrichment")) return(enrichment)
  if (is.list(enrichment) && !is.null(enrichment$results)) {
    enrichment$results <- .benchmark_report_normalize_enrichment_for_dotplot(enrichment$results)
    return(enrichment)
  }
  if (!is.data.frame(enrichment)) return(enrichment)
  df <- as.data.frame(enrichment)
  if (!nrow(df)) return(df)
  need <- c("term_name", "n_overlap", "score", "pvalue", "padj")
  if (all(need %in% names(df))) return(df)

  normalizer <- .benchmark_report_svdx_fn(".normalize_enrichment_result", required = FALSE)
  if (!is.null(normalizer) && all(c("pathway", "p_value", "q_value", "odds_ratio") %in% names(df))) {
    out <- tryCatch(normalizer(df, backend = "internal"), error = function(e) NULL)
    if (!is.null(out)) return(out)
  }

  pick <- function(candidates) {
    hit <- candidates[candidates %in% names(df)][1]
    if (is.na(hit)) NULL else hit
  }
  term_col <- pick(c("term_name", "pathway", "term", "set", "category", "target_set"))
  overlap_col <- pick(c("n_overlap", "overlap", "overlap_size", "Count", "count"))
  p_col <- pick(c("pvalue", "p_value", "pval", "p", "P.Value"))
  q_col <- pick(c("padj", "q_value", "qval", "q", "p.adjust", "FDR"))
  score_col <- pick(c("score", "odds_ratio", "NES", "statistic", "enrichment_score"))
  if (!is.null(term_col) && !"term_name" %in% names(df)) df$term_name <- as.character(df[[term_col]])
  if (!is.null(overlap_col) && !"n_overlap" %in% names(df)) df$n_overlap <- suppressWarnings(as.integer(df[[overlap_col]]))
  if (!is.null(p_col) && !"pvalue" %in% names(df)) df$pvalue <- suppressWarnings(as.numeric(df[[p_col]]))
  if (!"padj" %in% names(df)) {
    if (!is.null(q_col)) {
      df$padj <- suppressWarnings(as.numeric(df[[q_col]]))
    } else if ("pvalue" %in% names(df)) {
      df$padj <- stats::p.adjust(df$pvalue, method = "BH")
    }
  }
  if (!"score" %in% names(df)) {
    if (!is.null(score_col)) {
      df$score <- suppressWarnings(as.numeric(df[[score_col]]))
    } else if ("padj" %in% names(df)) {
      df$score <- -log10(pmax(df$padj, .Machine$double.xmin))
    }
  }
  df
}

.benchmark_report_enrichment_dotplot_object <- function(enrichment, context = "enrichment", verbose = TRUE) {
  df <- .benchmark_report_normalize_enrichment_for_dotplot(enrichment)
  if (inherits(df, "svdxEnrichment") || (is.list(df) && !is.null(df$results))) {
    df <- df$results
  }
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) {
    .benchmark_report_message("Skipping ", context, ": no enrichment rows available after normalization.", verbose = verbose)
    return(NULL)
  }
  df <- as.data.frame(df)
  required <- c("term_name", "n_overlap", "score", "pvalue", "padj")
  missing <- setdiff(required, names(df))
  if (length(missing)) {
    .benchmark_report_message(
      "Skipping ", context, ": missing enrichment plot columns: ",
      paste(missing, collapse = ", "), ".",
      verbose = verbose
    )
    return(NULL)
  }

  df$term_name <- as.character(df$term_name)
  df$n_overlap <- suppressWarnings(as.numeric(df$n_overlap))
  df$score <- suppressWarnings(as.numeric(df$score))
  df$pvalue <- suppressWarnings(as.numeric(df$pvalue))
  df$padj <- suppressWarnings(as.numeric(df$padj))

  keep <- stats::complete.cases(df[, required, drop = FALSE]) &
    nzchar(df$term_name) &
    is.finite(df$n_overlap) &
    is.finite(df$score) &
    is.finite(df$pvalue) &
    is.finite(df$padj)
  df <- df[keep, , drop = FALSE]
  if (!nrow(df)) {
    .benchmark_report_message(
      "Skipping ", context, ": no rows remain after coercing enrichment plot columns.",
      verbose = verbose
    )
    return(NULL)
  }

  structure(list(results = df), class = c("svdxEnrichment", "list"))
}


.benchmark_report_safe_stem_value <- function(x, fallback = "NA") {
  if (length(x) == 0L || is.na(x) || !nzchar(as.character(x))) return(fallback)
  out <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  out <- gsub("_+", "_", out)
  out <- gsub("^_|_$", "", out)
  if (!nzchar(out)) fallback else out
}

.benchmark_report_enrichment_split_values <- function(df) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(data.frame())
  scope <- if ("scope" %in% names(df)) as.character(df$scope) else rep("all", nrow(df))
  scope[is.na(scope) | !nzchar(scope)] <- "all"
  component <- if ("component" %in% names(df)) as.character(df$component) else rep(NA_character_, nrow(df))
  component[is.na(component) | !nzchar(component)] <- NA_character_
  out <- unique(data.frame(scope = scope, component = component, stringsAsFactors = FALSE))
  ord <- order(out$scope, is.na(out$component), out$component, na.last = TRUE)
  out[ord, , drop = FALSE]
}

.benchmark_report_plot_enrichment_by_scope_component <- function(
    enrichment,
    prefix,
    outdir,
    plot_enrichment_dotplot_fn,
    enrichment_top_n = 20L,
    figure_format = "pdf",
    verbose = TRUE
) {
  files <- character(0)
  if (is.null(enrichment) || is.null(plot_enrichment_dotplot_fn)) return(files)

  enr_obj <- .benchmark_report_enrichment_dotplot_object(
    enrichment,
    context = prefix,
    verbose = verbose
  )
  if (is.null(enr_obj)) return(files)
  df <- as.data.frame(enr_obj$results)
  if (!nrow(df)) return(files)

  split_values <- .benchmark_report_enrichment_split_values(df)
  if (!nrow(split_values)) return(files)

  for (ii in seq_len(nrow(split_values))) {
    scope <- split_values$scope[ii]
    component <- split_values$component[ii]

    keep <- rep(TRUE, nrow(df))
    if ("scope" %in% names(df)) {
      df_scope <- as.character(df$scope)
      df_scope[is.na(df_scope) | !nzchar(df_scope)] <- "all"
      keep <- keep & df_scope == scope
    }
    if ("component" %in% names(df)) {
      df_component <- as.character(df$component)
      df_component[is.na(df_component) | !nzchar(df_component)] <- NA_character_
      if (is.na(component)) {
        keep <- keep & is.na(df_component)
      } else {
        keep <- keep & !is.na(df_component) & df_component == component
      }
    }
    plot_rows <- df[keep, , drop = FALSE]
    if (!nrow(plot_rows)) next

    stem <- paste0(prefix, "_", .benchmark_report_safe_stem_value(scope))
    if (!is.na(component)) {
      stem <- paste0(stem, "_", .benchmark_report_safe_stem_value(component))
    }
    context <- if (!is.na(component)) {
      sprintf("%s [%s / %s]", prefix, scope, component)
    } else {
      sprintf("%s [%s]", prefix, scope)
    }

    plot_obj <- structure(list(results = plot_rows), class = c("svdxEnrichment", "list"))
    files <- c(files, .benchmark_report_safe_plot(
      plot_enrichment_dotplot_fn(plot_obj, top_n = enrichment_top_n),
      .benchmark_report_plot_file(outdir, stem, figure_format),
      width = 9,
      height = 6,
      verbose = verbose
    ))
  }

  files[!is.na(files)]
}

.benchmark_report_inference_component_features <- function(inf, component, n_features = 50L) {
  fr <- inf$feature_results %||% NULL
  if (is.null(fr) || !is.data.frame(fr) || !"feature" %in% names(fr) || !"component" %in% names(fr)) {
    return(character(0))
  }
  comp_vals <- unique(fr$component)
  component_use <- component
  if (!component_use %in% comp_vals) {
    comp_chr <- as.character(component)
    comp_vals_chr <- as.character(comp_vals)
    if (comp_chr %in% comp_vals_chr) component_use <- comp_vals[match(comp_chr, comp_vals_chr)]
  }
  fr <- fr[fr$component == component_use, , drop = FALSE]
  if (nrow(fr) == 0L) return(character(0))
  if ("qval" %in% names(fr)) {
    ord <- order(fr$qval, fr$pval %||% fr$qval, na.last = NA)
  } else if ("pval" %in% names(fr)) {
    ord <- order(fr$pval, na.last = NA)
  } else if ("stat" %in% names(fr)) {
    ord <- order(abs(fr$stat), decreasing = TRUE, na.last = NA)
  } else {
    ord <- seq_len(nrow(fr))
  }
  unique(as.character(fr$feature[ord]))[seq_len(min(n_features, length(unique(as.character(fr$feature[ord])))))]
}

.benchmark_report_svdx_inference_enrichment_results <- function(enrichment) {
  if (is.null(enrichment)) return(NULL)
  if (is.data.frame(enrichment)) return(enrichment)
  if (is.list(enrichment) && !is.null(enrichment$results)) return(enrichment$results)
  NULL
}

.benchmark_report_isvdx_dataset_inf <- function(inf, dataset) {
  if (is.null(inf)) return(NULL)
  if (!is.null(inf$dataset_results) && !is.null(inf$dataset_results[[dataset]])) {
    return(inf$dataset_results[[dataset]])
  }
  if (!is.null(inf$results) && !is.null(inf$results[[dataset]])) {
    return(inf$results[[dataset]])
  }
  NULL
}

.benchmark_report_isvdx_enrichment_results <- function(enrichment, dataset) {
  if (is.null(enrichment)) return(NULL)
  if (is.data.frame(enrichment) && "dataset" %in% names(enrichment)) {
    out <- enrichment[enrichment$dataset == dataset, , drop = FALSE]
    return(if (nrow(out) > 0L) out else NULL)
  }
  if (is.list(enrichment) && !is.null(enrichment$results)) {
    if (is.data.frame(enrichment$results) && "dataset" %in% names(enrichment$results)) {
      out <- enrichment$results[enrichment$results$dataset == dataset, , drop = FALSE]
      return(if (nrow(out) > 0L) out else NULL)
    }
    if (is.list(enrichment$results) && !is.null(enrichment$results[[dataset]])) {
      return(enrichment$results[[dataset]])
    }
  }
  NULL
}

.benchmark_report_call_svdx_workflow_plots <- function(
    rw,
    outdir,
    xaxes = NULL,
    components = 1:3,
    color = NULL,
    group = NULL,
    line_group = NULL,
    facet = NULL,
    u_plot_specs = NULL,
    features_top_n = 50L,
    feature_profiles_top_n = 5L,
    enrichment_top_n = 20L,
    volcano_labels_top_n = 25L,
    volcano_label_features = NULL,
    inference_scope = c("all", "component", "union", "global"),
    loading_heatmaps = TRUE,
    loading_feature_profiles = TRUE,
    loading_enrichment = TRUE,
    inference_volcanoes = TRUE,
    inference_heatmaps = TRUE,
    inference_feature_profiles = TRUE,
    inference_dotplots = TRUE,
    inference_enrichment = TRUE,
    feature_profile_summarize = c("mean", "median", "none"),
    feature_profile_show_se = TRUE,
    plot_args = list(),
    figure_format = "pdf",
    verbose = TRUE
) {
  inference_scope <- match.arg(inference_scope)
  feature_profile_summarize <- match.arg(feature_profile_summarize)
  files <- character(0)

  fit <- .benchmark_report_extract_svdx_fit(rw)
  if (is.null(fit)) {
    .benchmark_report_message("Skipping SVDX-specific analyses: no stored SVDX fit found.", verbose = verbose)
    return(files)
  }
  data <- as.data.frame(rw$data)
  Y <- as.matrix(rw$Y)
  zdata <- if (!is.null(rw$zdata)) as.data.frame(rw$zdata) else NULL
  if (ncol(Y) != nrow(data)) {
    stop("SVDX-specific report bridge received inconsistent data: ncol(rw$Y) != nrow(rw$data).", call. = FALSE)
  }
  if (!is.null(zdata) && nrow(zdata) != ncol(Y)) {
    stop("SVDX-specific report bridge received inconsistent zdata: nrow(rw$zdata) != ncol(rw$Y).", call. = FALSE)
  }

  plotU_fn <- .benchmark_report_svdx_fn("plotU")
  plotUTerms_fn <- .benchmark_report_svdx_fn("plotUTerms")
  plotUEmbedding_fn <- .benchmark_report_svdx_fn("plotUEmbedding")
  plotVLoadings_fn <- .benchmark_report_svdx_fn("plotVLoadings")
  plotVHeatmap_fn <- .benchmark_report_svdx_fn("plotVHeatmap")
  plotFeatureProfiles_fn <- .benchmark_report_svdx_fn("plotFeatureProfiles", required = FALSE)
  plot_union_volcano_fn <- .benchmark_report_svdx_fn("plot_union_volcano", required = FALSE)
  plot_component_volcano_fn <- .benchmark_report_svdx_fn("plot_component_volcano", required = FALSE)
  plotVInferenceHeatmap_fn <- .benchmark_report_svdx_fn("plotVInferenceHeatmap", required = FALSE)
  plotVInferenceDotplot_fn <- .benchmark_report_svdx_fn("plotVInferenceDotplot", required = FALSE)
  plotVLoadingInferenceScatter_fn <- .benchmark_report_svdx_fn("plotVLoadingInferenceScatter", required = FALSE)
  plot_enrichment_dotplot_fn <- .benchmark_report_svdx_fn("plot_enrichment_dotplot", required = FALSE)
  report_get_inference_component_name_fn <- .benchmark_report_svdx_fn(".report_get_inference_component_name", required = FALSE)
  report_component_label_fn <- .benchmark_report_svdx_fn(".report_component_label", required = FALSE)

  K_eff <- .benchmark_report_component_count_from_fit(
    fit,
    default = if (is.null(components)) 0L else max(components, na.rm = TRUE)
  )
  if (is.null(components)) {
    components <- seq_len(K_eff)
  } else {
    components <- suppressWarnings(as.integer(components))
    components <- components[!is.na(components) & components >= 1L & components <= K_eff]
    if (length(components) == 0L) components <- seq_len(min(3L, K_eff))
  }
  u_specs <- .benchmark_report_svdx_u_plot_specs(data, xaxes = xaxes, u_plot_specs = u_plot_specs)
  svdx_model_vars <- .benchmark_report_svdx_model_vars(fit, fallback_formula = rw$xformula)
  process_facet <- .benchmark_report_svdx_process_vars(
    facet, svdx_model_vars, what = "faceting", verbose = verbose
  )

  # U plots remain x-axis specific.  U-term plots are term-specific in the
  # current SVDX API, so call plotUTerms() once per component to avoid
  # duplicated output for every requested x-axis.
  for (spec in u_specs) {
    xa <- .benchmark_report_spec_value(spec, "var", .benchmark_report_spec_value(spec, "x", NULL))
    if (is.null(xa) || !xa %in% colnames(data)) next
    spec_group <- .benchmark_report_svdx_process_vars(
      .benchmark_report_spec_value(spec, "group", group),
      svdx_model_vars,
      what = "grouping",
      verbose = verbose
    )
    spec_line_group <- .benchmark_report_svdx_process_vars(
      .benchmark_report_spec_value(spec, "line_group", line_group),
      svdx_model_vars,
      what = "line grouping",
      verbose = verbose
    )
    spec_connect <- if (.benchmark_report_has_spec_field(spec, "connect")) isTRUE(spec$connect) else FALSE
    for (cc in components) {
      p_u <- do.call(plotU_fn, .benchmark_report_plotU_call_args(
        fit = fit, data = data, component = cc, x = xa,
        group = spec_group, facet = process_facet, fill = "typical",
        xformula = rw$xformula, plot_args = .benchmark_report_svdx_plot_args(plot_args, "plotU")
      ))
      p_u <- .benchmark_report_add_connected_u_layer(
        p_u, connect = spec_connect, group = spec_group, line_group = spec_line_group
      )
      files <- c(files, .benchmark_report_safe_plot(
        p_u,
        .benchmark_report_plot_file(outdir, sprintf("svdx_u_%s_component_%s", xa, cc), figure_format),
        width = 8, height = 5, verbose = verbose
      ))
    }
  }

  for (cc in components) {
    files <- c(files, .benchmark_report_safe_plot_collection(
      do.call(plotUTerms_fn, c(list(fit = fit, data = data, component = cc, facet = process_facet, fill = "typical"), .benchmark_report_svdx_plot_args(plot_args, "plotUTerms"))),
      outdir = outdir,
      stem = sprintf("svdx_u_terms_component_%s", cc),
      width = 9, height = 6, figure_format = figure_format, verbose = verbose
    ))
  }

  if (K_eff >= 2L) {
    emb_plot_args <- .benchmark_report_svdx_plot_args(plot_args, "plotUEmbedding")
    emb_observed <- if ("observed" %in% names(emb_plot_args)) {
      isTRUE(emb_plot_args$observed)
    } else {
      is.null(emb_plot_args$vary)
    }
    emb_data <- data
    emb_color <- color
    emb_group <- group
    if (isTRUE(emb_observed)) {
      color_info <- .benchmark_plot_interaction_var(
        emb_data, color, name = ".benchmark_embedding_color"
      )
      emb_data <- color_info$data
      emb_color <- color_info$var
      if (!is.null(group) && identical(as.character(group), as.character(color))) {
        emb_group <- emb_color
      } else {
        group_info <- .benchmark_plot_interaction_var(
          emb_data, group, name = ".benchmark_embedding_group"
        )
        emb_data <- group_info$data
        emb_group <- group_info$var
      }
    } else if (length(color) > 1L || length(group) > 1L) {
      stop(
        "Multivariable color/group aesthetics for predicted SVDX embeddings are not supported by the current plotUEmbedding API. ",
        "Use observed = TRUE or supply scalar embedding aesthetics.",
        call. = FALSE
      )
    }
    emb_args <- c(
      list(
        fit = fit, data = emb_data, components = c(1, 2),
        color = emb_color, group = emb_group, facet = facet, fill = "typical"
      ),
      emb_plot_args
    )
    p_emb <- do.call(plotUEmbedding_fn, emb_args)
    if (!is.null(color) && length(color) > 1L && inherits(p_emb, "ggplot")) {
      p_emb <- p_emb + ggplot2::labs(color = paste(as.character(color), collapse = " x "))
    }
    files <- c(files, .benchmark_report_safe_plot(
      p_emb,
      .benchmark_report_plot_file(outdir, "svdx_u_embedding_components_1_2", figure_format),
      width = 8, height = 6, verbose = verbose
    ))
  }

  if (isTRUE(loading_heatmaps)) {
    union_feats <- .benchmark_report_loading_union_features_from_fit(
      fit,
      components = components,
      n_features = features_top_n
    )
    if (!is.null(union_feats) && length(union_feats)) {
      files <- c(files, .benchmark_report_safe_plot(
        plotVHeatmap_fn(
          object = fit,
          Y = Y,
          data = data,
          zdata = zdata,
          features = union_feats,
          order_by = plot_args$heatmap_order_by %||% NULL,
          annotate = plot_args$heatmap_annotate %||% NULL,
          fdata = rw$fdata
        ),
        .benchmark_report_plot_file(outdir, "svdx_v_loading_heatmap_union", figure_format),
        width = 10, height = 9, verbose = verbose
      ))
    }
  }

  for (cc in components) {
    files <- c(files, .benchmark_report_safe_plot(
      plotVLoadings_fn(fit = fit, component = cc, n = features_top_n, fdata = rw$fdata),
      .benchmark_report_plot_file(outdir, sprintf("svdx_v_loadings_component_%s", cc), figure_format),
      width = 8, height = 6, verbose = verbose
    ))
    if (isTRUE(loading_heatmaps)) {
      files <- c(files, .benchmark_report_safe_plot(
        plotVHeatmap_fn(
          object = fit,
        Y = Y,
        data = data,
        zdata = zdata,
        component = cc,
        n_features = features_top_n,
        order_by = plot_args$heatmap_order_by %||% NULL,
        annotate = plot_args$heatmap_annotate %||% NULL,
          fdata = rw$fdata
        ),
        .benchmark_report_plot_file(outdir, sprintf("svdx_v_loading_heatmap_component_%s", cc), figure_format),
        width = 10, height = 9, verbose = verbose
      ))
    }
  }

  wf <- NULL
  profile_x <- xaxes[1] %||% NULL
  if (!is.null(plotFeatureProfiles_fn) && !is.null(profile_x) && profile_x %in% colnames(data)) {
    wf <- .benchmark_report_make_svdx_workflow(rw, verbose = verbose)
    if (isTRUE(loading_feature_profiles) && !is.null(wf)) {
      feats <- .benchmark_report_loading_union_features_from_fit(
        fit,
        components = components,
        n_features = feature_profiles_top_n
      )
      if (!is.null(feats) && length(feats)) {
        files <- c(files, .benchmark_report_safe_plot(
          plotFeatureProfiles_fn(
            wf,
            x = profile_x,
            features = feats,
            top_n = feature_profiles_top_n,
            color = plot_args$feature_profile_color %||% "feature",
            facet = plot_args$feature_profile_facet %||% "feature",
            summarize = feature_profile_summarize,
            show_se = feature_profile_show_se,
            se_width = plot_args$feature_profile_se_width %||% 0.15,
            scale = plot_args$feature_profile_scale %||% "zscore",
            line = TRUE,
            points = TRUE,
            fdata = rw$fdata
          ),
          .benchmark_report_plot_file(outdir, "svdx_feature_profiles_loading_union", figure_format),
          width = 9, height = 6, verbose = verbose
        ))
      }
      for (cc in components) {
        feats <- .benchmark_report_loading_features_from_fit(fit, component = cc, n_features = feature_profiles_top_n)
        if (is.null(feats) || !length(feats)) next
        files <- c(files, .benchmark_report_safe_plot(
          plotFeatureProfiles_fn(
            wf,
            x = profile_x,
            features = feats,
            top_n = feature_profiles_top_n,
            color = plot_args$feature_profile_color %||% "feature",
            facet = plot_args$feature_profile_facet %||% "feature",
            summarize = feature_profile_summarize,
            show_se = feature_profile_show_se,
            se_width = plot_args$feature_profile_se_width %||% 0.15,
            scale = plot_args$feature_profile_scale %||% "zscore",
            line = TRUE,
            points = TRUE,
            fdata = rw$fdata
          ),
          .benchmark_report_plot_file(outdir, sprintf("svdx_feature_profiles_loading_component_%s", cc), figure_format),
          width = 9, height = 6, verbose = verbose
        ))
      }
    }
  }

  inf <- .benchmark_report_extract_svdx_inf(rw)
  if (!is.null(inf)) {
    scopes <- .benchmark_report_expand_inference_scope(inference_scope)
    if (isTRUE(inference_volcanoes)) {
      if ("union" %in% scopes && !is.null(plot_union_volcano_fn)) {
        files <- c(files, .benchmark_report_safe_plot(
          plot_union_volcano_fn(
            inf,
            scope = "union",
            fdata = rw$fdata,
            label_top_n = volcano_labels_top_n,
            label_features = volcano_label_features
          ),
          .benchmark_report_plot_file(outdir, "svdx_inference_volcano_union", figure_format),
          width = 8, height = 6, verbose = verbose
        ))
      }
      if ("global" %in% scopes && !is.null(plot_union_volcano_fn)) {
        files <- c(files, .benchmark_report_safe_plot(
          plot_union_volcano_fn(
            inf,
            scope = "global",
            fdata = rw$fdata,
            label_top_n = volcano_labels_top_n,
            label_features = volcano_label_features
          ),
          .benchmark_report_plot_file(outdir, "svdx_inference_volcano_global", figure_format),
          width = 8, height = 6, verbose = verbose
        ))
      }
      if ("component" %in% scopes && !is.null(plot_component_volcano_fn) && !is.null(report_get_inference_component_name_fn)) {
        for (component in components) {
          lab <- tryCatch(
            if (!is.null(report_component_label_fn)) report_component_label_fn(component) else paste0("component_", component),
            error = function(e) paste0("component_", component)
          )
          comp_name <- tryCatch(
            report_get_inference_component_name_fn(inf, component),
            error = function(e) {
              .benchmark_report_message("Skipping SVDX component volcano for ", lab, ": ", conditionMessage(e), verbose = verbose)
              NULL
            }
          )
          if (is.null(comp_name)) next
          files <- c(files, .benchmark_report_safe_plot(
            plot_component_volcano_fn(
              inf,
              component = comp_name,
              fdata = rw$fdata,
              label_top_n = volcano_labels_top_n,
              label_features = volcano_label_features
            ),
            .benchmark_report_plot_file(outdir, sprintf("svdx_inference_volcano_%s", lab), figure_format),
            width = 8, height = 6, verbose = verbose
          ))
        }
      }
    }
    if (isTRUE(inference_heatmaps) && !is.null(plotVHeatmap_fn)) {
      if ("union" %in% scopes) {
        files <- c(files, .benchmark_report_safe_plot(
          plotVHeatmap_fn(
            object = inf,
            Y = Y,
            data = data,
            zdata = zdata,
            inference_scope = "union",
            select_by = "q",
            n_features = features_top_n,
            order_by = plot_args$heatmap_order_by %||% NULL,
            annotate = plot_args$heatmap_annotate %||% NULL,
            fdr = rw$alpha %||% 0.1,
            fdata = rw$fdata
          ),
          .benchmark_report_plot_file(outdir, "svdx_v_inference_heatmap_union", figure_format),
          width = 10, height = 8, verbose = verbose
        ))
      }
      if ("global" %in% scopes) {
        files <- c(files, .benchmark_report_safe_plot(
          plotVHeatmap_fn(
            object = inf,
            Y = Y,
            data = data,
            zdata = zdata,
            inference_scope = "global",
            select_by = "q",
            n_features = features_top_n,
            order_by = plot_args$heatmap_order_by %||% NULL,
            annotate = plot_args$heatmap_annotate %||% NULL,
            fdr = rw$alpha %||% 0.1,
            fdata = rw$fdata
          ),
          .benchmark_report_plot_file(outdir, "svdx_v_inference_heatmap_global", figure_format),
          width = 10, height = 8, verbose = verbose
        ))
      }
      if ("component" %in% scopes && !is.null(report_get_inference_component_name_fn)) {
        for (component in components) {
          lab <- tryCatch(
            if (!is.null(report_component_label_fn)) report_component_label_fn(component) else paste0("component_", component),
            error = function(e) paste0("component_", component)
          )
          comp_name <- tryCatch(
            report_get_inference_component_name_fn(inf, component),
            error = function(e) {
              .benchmark_report_message("Skipping SVDX component inference heatmap for ", lab, ": ", conditionMessage(e), verbose = verbose)
              NULL
            }
          )
          if (is.null(comp_name)) next
          files <- c(files, .benchmark_report_safe_plot(
            plotVHeatmap_fn(
              object = inf,
              Y = Y,
              data = data,
              zdata = zdata,
              inference_scope = "component",
              component = comp_name,
              select_by = "q",
              n_features = features_top_n,
              order_by = plot_args$heatmap_order_by %||% NULL,
              annotate = plot_args$heatmap_annotate %||% NULL,
              fdr = rw$alpha %||% 0.1,
              fdata = rw$fdata
            ),
            .benchmark_report_plot_file(outdir, sprintf("svdx_v_inference_heatmap_%s", lab), figure_format),
            width = 10, height = 8, verbose = verbose
          ))
        }
      }
    }
    if (isTRUE(inference_dotplots) && !is.null(plotVInferenceDotplot_fn)) {
      if ("union" %in% scopes) {
        files <- c(files, .benchmark_report_safe_plot(
          plotVInferenceDotplot_fn(inf, n_features = features_top_n, inference_scope = "union", select_by = "q", fdr = rw$alpha %||% 0.1, fdata = rw$fdata),
          .benchmark_report_plot_file(outdir, "svdx_v_inference_dotplot_union", figure_format),
          width = 9, height = 6, verbose = verbose
        ))
      }
      if ("global" %in% scopes) {
        files <- c(files, .benchmark_report_safe_plot(
          plotVInferenceDotplot_fn(inf, n_features = features_top_n, inference_scope = "global", select_by = "q", fdr = rw$alpha %||% 0.1, fdata = rw$fdata),
          .benchmark_report_plot_file(outdir, "svdx_v_inference_dotplot_global", figure_format),
          width = 9, height = 6, verbose = verbose
        ))
      }
      if ("component" %in% scopes && !is.null(report_get_inference_component_name_fn)) {
        for (component in components) {
          lab <- tryCatch(
            if (!is.null(report_component_label_fn)) report_component_label_fn(component) else paste0("component_", component),
            error = function(e) paste0("component_", component)
          )
          comp_name <- tryCatch(
            report_get_inference_component_name_fn(inf, component),
            error = function(e) {
              .benchmark_report_message("Skipping SVDX component inference dotplot for ", lab, ": ", conditionMessage(e), verbose = verbose)
              NULL
            }
          )
          if (is.null(comp_name)) next
          files <- c(files, .benchmark_report_safe_plot(
            plotVInferenceDotplot_fn(inf, inference_scope = "component", component = comp_name, n_features = features_top_n, select_by = "q", fdr = rw$alpha %||% 0.1, fdata = rw$fdata),
            .benchmark_report_plot_file(outdir, sprintf("svdx_v_inference_dotplot_%s", lab), figure_format),
            width = 9, height = 6, verbose = verbose
          ))
        }
      }
    }
    if (isTRUE(inference_feature_profiles) && !is.null(plotFeatureProfiles_fn) && !is.null(profile_x) && profile_x %in% colnames(data)) {
      if ("union" %in% scopes) {
        files <- c(files, .benchmark_report_safe_plot(
          plotFeatureProfiles_fn(inf, Y = Y, data = data, source = "passed", x = profile_x, inference_scope = "union", component = NULL, select_by = "q", top_n = feature_profiles_top_n, color = plot_args$feature_profile_color %||% "feature", facet = plot_args$feature_profile_facet %||% "feature", summarize = feature_profile_summarize, show_se = feature_profile_show_se, se_width = plot_args$feature_profile_se_width %||% 0.15, scale = plot_args$feature_profile_scale %||% "zscore", line = TRUE, points = TRUE, fdata = rw$fdata),
          .benchmark_report_plot_file(outdir, "svdx_feature_profiles_inference_union", figure_format),
          width = 9, height = 6, verbose = verbose
        ))
      }
      if ("global" %in% scopes) {
        files <- c(files, .benchmark_report_safe_plot(
          plotFeatureProfiles_fn(inf, Y = Y, data = data, source = "passed", x = profile_x, inference_scope = "global", component = NULL, select_by = "q", top_n = feature_profiles_top_n, color = plot_args$feature_profile_color %||% "feature", facet = plot_args$feature_profile_facet %||% "feature", summarize = feature_profile_summarize, show_se = feature_profile_show_se, se_width = plot_args$feature_profile_se_width %||% 0.15, scale = plot_args$feature_profile_scale %||% "zscore", line = TRUE, points = TRUE, fdata = rw$fdata),
          .benchmark_report_plot_file(outdir, "svdx_feature_profiles_inference_global", figure_format),
          width = 9, height = 6, verbose = verbose
        ))
      }
      if ("component" %in% scopes && !is.null(report_get_inference_component_name_fn)) {
        for (component in components) {
          lab <- tryCatch(
            if (!is.null(report_component_label_fn)) report_component_label_fn(component) else paste0("component_", component),
            error = function(e) paste0("component_", component)
          )
          comp_name <- tryCatch(
            report_get_inference_component_name_fn(inf, component),
            error = function(e) {
              .benchmark_report_message("Skipping SVDX component inference feature profiles for ", lab, ": ", conditionMessage(e), verbose = verbose)
              NULL
            }
          )
          if (is.null(comp_name)) next
          files <- c(files, .benchmark_report_safe_plot(
            plotFeatureProfiles_fn(inf, Y = Y, data = data, source = "passed", x = profile_x, inference_scope = "component", component = comp_name, select_by = "q", top_n = feature_profiles_top_n, color = plot_args$feature_profile_color %||% "feature", facet = plot_args$feature_profile_facet %||% "feature", summarize = feature_profile_summarize, show_se = feature_profile_show_se, se_width = plot_args$feature_profile_se_width %||% 0.15, scale = plot_args$feature_profile_scale %||% "zscore", line = TRUE, points = TRUE, fdata = rw$fdata),
            .benchmark_report_plot_file(outdir, sprintf("svdx_feature_profiles_inference_%s", lab), figure_format),
            width = 9, height = 6, verbose = verbose
          ))
        }
      }
    }
    if (!is.null(plotVLoadingInferenceScatter_fn)) {
      for (cc in components) {
        files <- c(files, .benchmark_report_safe_plot(
          plotVLoadingInferenceScatter_fn(fit = fit, inf = inf, component = cc),
          .benchmark_report_plot_file(outdir, sprintf("svdx_v_loading_inference_scatter_%s", cc), figure_format),
          width = 8, height = 6, verbose = verbose
        ))
      }
    }
  }

  loading_enr_obj <- .benchmark_report_extract_svdx_loading_enrichment(rw)
  if (isTRUE(loading_enrichment) && !is.null(loading_enr_obj) && !is.null(plot_enrichment_dotplot_fn)) {
    files <- c(files, .benchmark_report_plot_enrichment_by_scope_component(
      loading_enr_obj,
      prefix = "svdx_loading_enrichment",
      outdir = outdir,
      plot_enrichment_dotplot_fn = plot_enrichment_dotplot_fn,
      enrichment_top_n = enrichment_top_n,
      figure_format = figure_format,
      verbose = verbose
    ))
  }

  inference_enr_obj <- .benchmark_report_extract_svdx_inference_enrichment(rw)
  if (isTRUE(inference_enrichment) && !is.null(inference_enr_obj) && !is.null(plot_enrichment_dotplot_fn)) {
    files <- c(files, .benchmark_report_plot_enrichment_by_scope_component(
      inference_enr_obj,
      prefix = "svdx_inference_enrichment",
      outdir = outdir,
      plot_enrichment_dotplot_fn = plot_enrichment_dotplot_fn,
      enrichment_top_n = enrichment_top_n,
      figure_format = figure_format,
      verbose = verbose
    ))
  }

  files[!is.na(files)]
}

.benchmark_report_call_isvdx_workflow_plots <- function(
    rw,
    outdir,
    xaxes = NULL,
    components = 1:3,
    datasets = NULL,
    color = NULL,
    group = NULL,
    line_group = NULL,
    facet = NULL,
    u_plot_specs = NULL,
    features_top_n = 50L,
    feature_profiles_top_n = 5L,
    enrichment_top_n = 20L,
    volcano_labels_top_n = 25L,
    volcano_label_features = NULL,
    inference_scope = c("all", "component", "union", "global"),
    loading_heatmaps = TRUE,
    loading_feature_profiles = TRUE,
    loading_enrichment = TRUE,
    inference_volcanoes = TRUE,
    inference_heatmaps = TRUE,
    inference_feature_profiles = TRUE,
    inference_dotplots = TRUE,
    inference_enrichment = TRUE,
    feature_profile_summarize = c("mean", "median", "none"),
    feature_profile_show_se = TRUE,
    plot_args = list(),
    figure_format = "pdf",
    verbose = TRUE
) {
  inference_scope <- match.arg(inference_scope)
  feature_profile_summarize <- match.arg(feature_profile_summarize)
  files <- character(0)

  fit <- .benchmark_report_extract_isvdx_fit(rw)
  if (is.null(fit)) {
    .benchmark_report_message("Skipping ISVDX-specific analyses: no stored ISVDX fit found.", verbose = verbose)
    return(files)
  }
  if (is.null(rw$Y) || !is.list(rw$Y) || is.null(rw$xdata)) {
    stop("ISVDX-specific report bridge requires rw$Y list and rw$xdata.", call. = FALSE)
  }

  data <- as.data.frame(rw$xdata)
  datasets <- datasets %||% names(rw$Y)
  datasets <- intersect(datasets, names(rw$Y))
  if (length(datasets) == 0L) return(files)
  for (ds in datasets) {
    if (ncol(as.matrix(rw$Y[[ds]])) != nrow(data)) {
      stop("ISVDX-specific report bridge received inconsistent dimensions for dataset: ", ds, call. = FALSE)
    }
  }

  plotU_fn <- .benchmark_report_svdx_fn("plotU")
  plotUTerms_fn <- .benchmark_report_svdx_fn("plotUTerms")
  plotUEmbedding_fn <- .benchmark_report_svdx_fn("plotUEmbedding")
  plotVLoadings_fn <- .benchmark_report_svdx_fn("plotVLoadings")
  plotVHeatmap_fn <- .benchmark_report_svdx_fn("plotVHeatmap")
  plotFeatureProfiles_fn <- .benchmark_report_svdx_fn("plotFeatureProfiles", required = FALSE)
  plot_union_volcano_fn <- .benchmark_report_svdx_fn("plot_union_volcano", required = FALSE)
  plot_component_volcano_fn <- .benchmark_report_svdx_fn("plot_component_volcano", required = FALSE)
  plotVInferenceHeatmap_fn <- .benchmark_report_svdx_fn("plotVInferenceHeatmap", required = FALSE)
  plotVInferenceDotplot_fn <- .benchmark_report_svdx_fn("plotVInferenceDotplot", required = FALSE)
  plotVLoadingInferenceScatter_fn <- .benchmark_report_svdx_fn("plotVLoadingInferenceScatter", required = FALSE)
  plot_enrichment_dotplot_fn <- .benchmark_report_svdx_fn("plot_enrichment_dotplot", required = FALSE)
  report_get_inference_component_name_fn <- .benchmark_report_svdx_fn(".report_get_inference_component_name", required = FALSE)
  report_component_label_fn <- .benchmark_report_svdx_fn(".report_component_label", required = FALSE)

  K_eff <- .benchmark_report_component_count_from_fit(
    fit,
    default = if (is.null(components)) 0L else max(components, na.rm = TRUE)
  )
  if (is.null(components)) {
    components <- seq_len(K_eff)
  } else {
    components <- suppressWarnings(as.integer(components))
    components <- components[!is.na(components) & components >= 1L & components <= K_eff]
    if (length(components) == 0L) components <- seq_len(min(3L, K_eff))
  }
  u_specs <- .benchmark_report_svdx_u_plot_specs(data, xaxes = xaxes, u_plot_specs = u_plot_specs)

  for (spec in u_specs) {
    xa <- .benchmark_report_spec_value(spec, "var", .benchmark_report_spec_value(spec, "x", NULL))
    if (is.null(xa) || !xa %in% colnames(data)) next
    spec_group <- .benchmark_report_spec_value(spec, "group", group)
    spec_line_group <- .benchmark_report_spec_value(spec, "line_group", line_group)
    spec_connect <- if (.benchmark_report_has_spec_field(spec, "connect")) isTRUE(spec$connect) else FALSE
    for (cc in components) {
      p_u <- do.call(plotU_fn, .benchmark_report_plotU_call_args(
        fit = fit, data = data, component = cc, x = xa,
        group = spec_group, facet = facet, fill = "typical",
        xformula = rw$xformula, plot_args = .benchmark_report_svdx_plot_args(plot_args, "plotU")
      ))
      p_u <- .benchmark_report_add_connected_u_layer(
        p_u, connect = spec_connect, group = spec_group, line_group = spec_line_group
      )
      files <- c(files, .benchmark_report_safe_plot(
        p_u,
        .benchmark_report_plot_file(outdir, sprintf("isvdx_u_%s_component_%s", xa, cc), figure_format),
        width = 8, height = 5, verbose = verbose
      ))
    }
  }

  for (cc in components) {
    files <- c(files, .benchmark_report_safe_plot_collection(
      do.call(plotUTerms_fn, c(list(fit = fit, data = data, component = cc, facet = facet, fill = "typical"), .benchmark_report_svdx_plot_args(plot_args, "plotUTerms"))),
      outdir = outdir,
      stem = sprintf("isvdx_u_terms_component_%s", cc),
      width = 9, height = 6, figure_format = figure_format, verbose = verbose
    ))
  }

  if (K_eff >= 2L) {
    files <- c(files, .benchmark_report_safe_plot(
      do.call(plotUEmbedding_fn, c(list(fit = fit, data = data, components = c(1, 2), color = color, group = group, facet = facet, fill = "typical"), .benchmark_report_svdx_plot_args(plot_args, "plotUEmbedding"))),
      .benchmark_report_plot_file(outdir, "isvdx_u_embedding_components_1_2", figure_format),
      width = 8, height = 6, verbose = verbose
    ))
  }

  wf <- NULL
  profile_x <- xaxes[1] %||% NULL
  if (!is.null(plotFeatureProfiles_fn) && !is.null(profile_x) && profile_x %in% colnames(data)) {
    wf <- .benchmark_report_make_isvdx_workflow(rw, verbose = verbose)
  }

  inf <- .benchmark_report_extract_isvdx_inf(rw)
  for (ds in datasets) {
    ds_dir <- file.path(outdir, .benchmark_report_safe_name(ds))
    dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
    Y_ds <- as.matrix(rw$Y[[ds]])
    zdata_ds <- if (!is.null(rw$zdata) && !is.null(rw$zdata[[ds]])) as.data.frame(rw$zdata[[ds]]) else NULL
    if (isTRUE(loading_heatmaps)) {
      union_feats <- .benchmark_report_loading_union_features_from_fit(
        fit,
        components = components,
        n_features = features_top_n,
        dataset = ds
      )
      if (!is.null(union_feats) && length(union_feats)) {
        files <- c(files, .benchmark_report_safe_plot(
          plotVHeatmap_fn(
            object = fit,
            Y = Y_ds,
            data = data,
            zdata = zdata_ds,
            dataset = ds,
            features = union_feats,
            order_by = plot_args$heatmap_order_by %||% NULL,
            annotate = plot_args$heatmap_annotate %||% NULL,
            fdata = rw$fdata
          ),
          .benchmark_report_plot_file(ds_dir, sprintf("isvdx_v_loading_heatmap_%s_union", ds), figure_format),
          width = 10, height = 9, verbose = verbose
        ))
      }
    }
    for (cc in components) {
      files <- c(files, .benchmark_report_safe_plot(
        plotVLoadings_fn(fit = fit, dataset = ds, component = cc, n = features_top_n, fdata = rw$fdata),
        .benchmark_report_plot_file(ds_dir, sprintf("isvdx_v_loadings_%s_component_%s", ds, cc), figure_format),
        width = 8, height = 6, verbose = verbose
      ))
      if (isTRUE(loading_heatmaps)) {
        files <- c(files, .benchmark_report_safe_plot(
          plotVHeatmap_fn(
            object = fit,
            Y = Y_ds,
          data = data,
          zdata = zdata_ds,
          dataset = ds,
          component = cc,
          n_features = features_top_n,
            fdata = rw$fdata
          ),
          .benchmark_report_plot_file(ds_dir, sprintf("isvdx_v_loading_heatmap_%s_component_%s", ds, cc), figure_format),
          width = 10, height = 9, verbose = verbose
        ))
      }
    }

    if (isTRUE(loading_feature_profiles) && !is.null(plotFeatureProfiles_fn) && !is.null(wf) && !is.null(profile_x) && profile_x %in% colnames(data)) {
      feats <- .benchmark_report_loading_union_features_from_fit(
        fit,
        components = components,
        n_features = feature_profiles_top_n,
        dataset = ds
      )
      if (!is.null(feats) && length(feats)) {
        files <- c(files, .benchmark_report_safe_plot(
          plotFeatureProfiles_fn(wf, dataset = ds, x = profile_x, features = feats, top_n = feature_profiles_top_n, color = plot_args$feature_profile_color %||% "feature", facet = plot_args$feature_profile_facet %||% "feature", summarize = feature_profile_summarize, show_se = feature_profile_show_se, se_width = plot_args$feature_profile_se_width %||% 0.15, scale = plot_args$feature_profile_scale %||% "zscore", line = TRUE, points = TRUE, fdata = rw$fdata),
          .benchmark_report_plot_file(ds_dir, sprintf("isvdx_feature_profiles_loading_%s_union", ds), figure_format),
          width = 9, height = 6, verbose = verbose
        ))
      }
      for (cc in components) {
        feats <- .benchmark_report_loading_features_from_fit(fit, component = cc, n_features = feature_profiles_top_n, dataset = ds)
        if (is.null(feats) || !length(feats)) next
        files <- c(files, .benchmark_report_safe_plot(
          plotFeatureProfiles_fn(wf, dataset = ds, x = profile_x, features = feats, top_n = feature_profiles_top_n, color = plot_args$feature_profile_color %||% "feature", facet = plot_args$feature_profile_facet %||% "feature", summarize = feature_profile_summarize, show_se = feature_profile_show_se, se_width = plot_args$feature_profile_se_width %||% 0.15, scale = plot_args$feature_profile_scale %||% "zscore", line = TRUE, points = TRUE, fdata = rw$fdata),
          .benchmark_report_plot_file(ds_dir, sprintf("isvdx_feature_profiles_loading_%s_component_%s", ds, cc), figure_format),
          width = 9, height = 6, verbose = verbose
        ))
      }
    }

    ds_inf <- .benchmark_report_isvdx_dataset_inf(inf, ds)
    if (!is.null(ds_inf)) {
      scopes <- .benchmark_report_expand_inference_scope(inference_scope)
      if (isTRUE(inference_volcanoes)) {
        if ("union" %in% scopes && !is.null(plot_union_volcano_fn)) {
          files <- c(files, .benchmark_report_safe_plot(
            plot_union_volcano_fn(
              ds_inf,
              scope = "union",
              fdata = rw$fdata,
              dataset = ds,
              label_top_n = volcano_labels_top_n,
              label_features = volcano_label_features
            ),
            .benchmark_report_plot_file(ds_dir, sprintf("isvdx_inference_volcano_%s_union", ds), figure_format),
            width = 8, height = 6, verbose = verbose
          ))
        }
        if ("global" %in% scopes && !is.null(plot_union_volcano_fn)) {
          files <- c(files, .benchmark_report_safe_plot(
            plot_union_volcano_fn(
              ds_inf,
              scope = "global",
              fdata = rw$fdata,
              dataset = ds,
              label_top_n = volcano_labels_top_n,
              label_features = volcano_label_features
            ),
            .benchmark_report_plot_file(ds_dir, sprintf("isvdx_inference_volcano_%s_global", ds), figure_format),
            width = 8, height = 6, verbose = verbose
          ))
        }
        if ("component" %in% scopes && !is.null(plot_component_volcano_fn) && !is.null(report_get_inference_component_name_fn)) {
          for (component in components) {
            lab <- tryCatch(
              if (!is.null(report_component_label_fn)) report_component_label_fn(component) else paste0("component_", component),
              error = function(e) paste0("component_", component)
            )
            comp_name <- tryCatch(
              report_get_inference_component_name_fn(ds_inf, component),
              error = function(e) {
                .benchmark_report_message("Skipping ISVDX component volcano for ", ds, " / ", lab, ": ", conditionMessage(e), verbose = verbose)
                NULL
              }
            )
            if (is.null(comp_name)) next
            files <- c(files, .benchmark_report_safe_plot(
              plot_component_volcano_fn(
                ds_inf,
                component = comp_name,
                fdata = rw$fdata,
                dataset = ds,
                label_top_n = volcano_labels_top_n,
                label_features = volcano_label_features
              ),
              .benchmark_report_plot_file(ds_dir, sprintf("isvdx_inference_volcano_%s_%s", ds, lab), figure_format),
              width = 8, height = 6, verbose = verbose
            ))
          }
        }
      }
      if (isTRUE(inference_heatmaps) && !is.null(plotVHeatmap_fn)) {
        if ("union" %in% scopes) {
          files <- c(files, .benchmark_report_safe_plot(
            plotVHeatmap_fn(object = ds_inf, Y = Y_ds, data = data, zdata = zdata_ds, dataset = ds, inference_scope = "union", select_by = "q", n_features = features_top_n, order_by = plot_args$heatmap_order_by %||% NULL, annotate = plot_args$heatmap_annotate %||% NULL, fdr = rw$alpha %||% 0.1, fdata = rw$fdata),
            .benchmark_report_plot_file(ds_dir, sprintf("isvdx_v_inference_heatmap_%s_union", ds), figure_format),
            width = 10, height = 8, verbose = verbose
          ))
        }
        if ("global" %in% scopes) {
          files <- c(files, .benchmark_report_safe_plot(
            plotVHeatmap_fn(object = ds_inf, Y = Y_ds, data = data, zdata = zdata_ds, dataset = ds, inference_scope = "global", select_by = "q", n_features = features_top_n, order_by = plot_args$heatmap_order_by %||% NULL, annotate = plot_args$heatmap_annotate %||% NULL, fdr = rw$alpha %||% 0.1, fdata = rw$fdata),
            .benchmark_report_plot_file(ds_dir, sprintf("isvdx_v_inference_heatmap_%s_global", ds), figure_format),
            width = 10, height = 8, verbose = verbose
          ))
        }
        if ("component" %in% scopes && !is.null(report_get_inference_component_name_fn)) {
          for (component in components) {
            lab <- tryCatch(
              if (!is.null(report_component_label_fn)) report_component_label_fn(component) else paste0("component_", component),
              error = function(e) paste0("component_", component)
            )
            comp_name <- tryCatch(
              report_get_inference_component_name_fn(ds_inf, component),
              error = function(e) {
                .benchmark_report_message("Skipping ISVDX component inference heatmap for ", ds, " / ", lab, ": ", conditionMessage(e), verbose = verbose)
                NULL
              }
            )
            if (is.null(comp_name)) next
            files <- c(files, .benchmark_report_safe_plot(
              plotVHeatmap_fn(object = ds_inf, Y = Y_ds, data = data, zdata = zdata_ds, dataset = ds, inference_scope = "component", component = comp_name, select_by = "q", n_features = features_top_n, order_by = plot_args$heatmap_order_by %||% NULL, annotate = plot_args$heatmap_annotate %||% NULL, fdr = rw$alpha %||% 0.1, fdata = rw$fdata),
              .benchmark_report_plot_file(ds_dir, sprintf("isvdx_v_inference_heatmap_%s_%s", ds, lab), figure_format),
              width = 10, height = 8, verbose = verbose
            ))
          }
        }
      }
      if (isTRUE(inference_dotplots) && !is.null(plotVInferenceDotplot_fn)) {
        if ("union" %in% scopes) {
          files <- c(files, .benchmark_report_safe_plot(
            plotVInferenceDotplot_fn(ds_inf, dataset = ds, n_features = features_top_n, inference_scope = "union", select_by = "q", fdr = rw$alpha %||% 0.1, fdata = rw$fdata),
            .benchmark_report_plot_file(ds_dir, sprintf("isvdx_v_inference_dotplot_%s_union", ds), figure_format),
            width = 9, height = 6, verbose = verbose
          ))
        }
        if ("global" %in% scopes) {
          files <- c(files, .benchmark_report_safe_plot(
            plotVInferenceDotplot_fn(ds_inf, dataset = ds, n_features = features_top_n, inference_scope = "global", select_by = "q", fdr = rw$alpha %||% 0.1, fdata = rw$fdata),
            .benchmark_report_plot_file(ds_dir, sprintf("isvdx_v_inference_dotplot_%s_global", ds), figure_format),
            width = 9, height = 6, verbose = verbose
          ))
        }
        if ("component" %in% scopes && !is.null(report_get_inference_component_name_fn)) {
          for (component in components) {
            lab <- tryCatch(
              if (!is.null(report_component_label_fn)) report_component_label_fn(component) else paste0("component_", component),
              error = function(e) paste0("component_", component)
            )
            comp_name <- tryCatch(
              report_get_inference_component_name_fn(ds_inf, component),
              error = function(e) {
                .benchmark_report_message("Skipping ISVDX component inference dotplot for ", ds, " / ", lab, ": ", conditionMessage(e), verbose = verbose)
                NULL
              }
            )
            if (is.null(comp_name)) next
            files <- c(files, .benchmark_report_safe_plot(
              plotVInferenceDotplot_fn(ds_inf, dataset = ds, inference_scope = "component", component = comp_name, n_features = features_top_n, select_by = "q", fdr = rw$alpha %||% 0.1, fdata = rw$fdata),
              .benchmark_report_plot_file(ds_dir, sprintf("isvdx_v_inference_dotplot_%s_%s", ds, lab), figure_format),
              width = 9, height = 6, verbose = verbose
            ))
          }
        }
      }
      if (isTRUE(inference_feature_profiles) && !is.null(plotFeatureProfiles_fn) && !is.null(profile_x) && profile_x %in% colnames(data)) {
        scopes <- .benchmark_report_expand_inference_scope(inference_scope)
        if ("union" %in% scopes) {
          files <- c(files, .benchmark_report_safe_plot(
            plotFeatureProfiles_fn(ds_inf, Y = Y_ds, data = data, dataset = ds, source = "passed", x = profile_x, inference_scope = "union", component = NULL, select_by = "q", top_n = feature_profiles_top_n, color = plot_args$feature_profile_color %||% "feature", facet = plot_args$feature_profile_facet %||% "feature", summarize = feature_profile_summarize, show_se = feature_profile_show_se, se_width = plot_args$feature_profile_se_width %||% 0.15, scale = plot_args$feature_profile_scale %||% "zscore", line = TRUE, points = TRUE, fdata = rw$fdata),
            .benchmark_report_plot_file(ds_dir, sprintf("isvdx_feature_profiles_inference_%s_union", ds), figure_format),
            width = 9, height = 6, verbose = verbose
          ))
        }
        if ("global" %in% scopes) {
          files <- c(files, .benchmark_report_safe_plot(
            plotFeatureProfiles_fn(ds_inf, Y = Y_ds, data = data, dataset = ds, source = "passed", x = profile_x, inference_scope = "global", component = NULL, select_by = "q", top_n = feature_profiles_top_n, color = plot_args$feature_profile_color %||% "feature", facet = plot_args$feature_profile_facet %||% "feature", summarize = feature_profile_summarize, show_se = feature_profile_show_se, se_width = plot_args$feature_profile_se_width %||% 0.15, scale = plot_args$feature_profile_scale %||% "zscore", line = TRUE, points = TRUE, fdata = rw$fdata),
            .benchmark_report_plot_file(ds_dir, sprintf("isvdx_feature_profiles_inference_%s_global", ds), figure_format),
            width = 9, height = 6, verbose = verbose
          ))
        }
        if ("component" %in% scopes && !is.null(report_get_inference_component_name_fn)) {
          for (component in components) {
            lab <- tryCatch(
              if (!is.null(report_component_label_fn)) report_component_label_fn(component) else paste0("component_", component),
              error = function(e) paste0("component_", component)
            )
            comp_name <- tryCatch(
              report_get_inference_component_name_fn(ds_inf, component),
              error = function(e) {
                .benchmark_report_message("Skipping ISVDX component inference feature profiles for ", ds, " / ", lab, ": ", conditionMessage(e), verbose = verbose)
                NULL
              }
            )
            if (is.null(comp_name)) next
            files <- c(files, .benchmark_report_safe_plot(
              plotFeatureProfiles_fn(ds_inf, Y = Y_ds, data = data, dataset = ds, source = "passed", x = profile_x, inference_scope = "component", component = comp_name, select_by = "q", top_n = feature_profiles_top_n, color = plot_args$feature_profile_color %||% "feature", facet = plot_args$feature_profile_facet %||% "feature", summarize = feature_profile_summarize, show_se = feature_profile_show_se, se_width = plot_args$feature_profile_se_width %||% 0.15, scale = plot_args$feature_profile_scale %||% "zscore", line = TRUE, points = TRUE, fdata = rw$fdata),
              .benchmark_report_plot_file(ds_dir, sprintf("isvdx_feature_profiles_inference_%s_%s", ds, lab), figure_format),
              width = 9, height = 6, verbose = verbose
            ))
          }
        }
      }
      if (!is.null(plotVLoadingInferenceScatter_fn)) {
        for (cc in components) {
          files <- c(files, .benchmark_report_safe_plot(
            plotVLoadingInferenceScatter_fn(fit = fit, inf = ds_inf, dataset = ds, component = cc),
            .benchmark_report_plot_file(ds_dir, sprintf("isvdx_v_loading_inference_scatter_%s_%s", ds, cc), figure_format),
            width = 8, height = 6, verbose = verbose
          ))
        }
      }
    }

    loading_enr_results <- .benchmark_report_isvdx_enrichment_results(.benchmark_report_extract_isvdx_loading_enrichment(rw), ds)
    if (isTRUE(loading_enrichment) && !is.null(loading_enr_results) && !is.null(plot_enrichment_dotplot_fn)) {
      files <- c(files, .benchmark_report_plot_enrichment_by_scope_component(
        list(results = as.data.frame(loading_enr_results)),
        prefix = sprintf("isvdx_loading_enrichment_%s", ds),
        outdir = ds_dir,
        plot_enrichment_dotplot_fn = plot_enrichment_dotplot_fn,
        enrichment_top_n = enrichment_top_n,
        figure_format = figure_format,
        verbose = verbose
      ))
    }

    inference_enr_results <- .benchmark_report_isvdx_enrichment_results(.benchmark_report_extract_isvdx_inference_enrichment(rw), ds)
    if (isTRUE(inference_enrichment) && !is.null(inference_enr_results) && !is.null(plot_enrichment_dotplot_fn)) {
      files <- c(files, .benchmark_report_plot_enrichment_by_scope_component(
        list(results = as.data.frame(inference_enr_results)),
        prefix = sprintf("isvdx_inference_enrichment_%s", ds),
        outdir = ds_dir,
        plot_enrichment_dotplot_fn = plot_enrichment_dotplot_fn,
        enrichment_top_n = enrichment_top_n,
        figure_format = figure_format,
        verbose = verbose
      ))
    }
  }

  files[!is.na(files)]
}

#' Run reports for a single real-world benchmark object.
#'
#' This function assembles method-comparison plots from a completed
#' run_real_world_single_benchmark() result. It does not rerun methods.
#'

.benchmark_report_loading_enrichment_results <- function(rw) {
  load_enr <- rw$loading_enrichment %||% NULL
  if (is.null(load_enr)) return(NULL)
  if (is.data.frame(load_enr)) return(as.data.frame(load_enr))
  if (is.list(load_enr) && !is.null(load_enr$results)) {
    return(as.data.frame(load_enr$results))
  }
  NULL
}

.benchmark_report_plot_real_world_loading_enrichment <- function(
    rw,
    outdir,
    methods = NULL,
    datasets = NULL,
    integrated = FALSE,
    enrichment_top_n = 20L,
    figure_format = "pdf",
    verbose = TRUE
) {
  plot_enrichment_dotplot_fn <- .benchmark_report_svdx_fn("plot_enrichment_dotplot", required = FALSE)
  if (is.null(plot_enrichment_dotplot_fn)) return(character(0))

  enr <- .benchmark_report_loading_enrichment_results(rw)
  if (is.null(enr) || !nrow(enr)) return(character(0))

  if (!"method" %in% names(enr)) {
    .benchmark_report_message("Skipping benchmark loading enrichment plots: no method column found.", verbose = verbose)
    return(character(0))
  }
  if (!"scope" %in% names(enr)) {
    .benchmark_report_message("Skipping benchmark loading enrichment plots: no scope column found.", verbose = verbose)
    return(character(0))
  }

  methods <- methods %||% unique(as.character(enr$method))
  methods <- intersect(as.character(methods), unique(as.character(enr$method)))
  if (!length(methods)) return(character(0))

  if (isTRUE(integrated)) {
    if (!"dataset" %in% names(enr)) {
      .benchmark_report_message("Skipping integrated benchmark loading enrichment plots: no dataset column found.", verbose = verbose)
      return(character(0))
    }
    datasets <- datasets %||% unique(as.character(enr$dataset))
    datasets <- intersect(as.character(datasets), unique(as.character(enr$dataset)))
    if (!length(datasets)) return(character(0))
  } else {
    datasets <- NA_character_
  }

  safe_chr <- function(x, fallback = "NA") {
    x <- as.character(x %||% fallback)
    x[is.na(x) | !nzchar(x)] <- fallback
    .benchmark_report_safe_name(x)
  }

  files <- character(0)
  for (ds in datasets) {
    ds_rows <- if (isTRUE(integrated)) enr[as.character(enr$dataset) == ds, , drop = FALSE] else enr
    if (!nrow(ds_rows)) next
    for (method in methods) {
      method_rows <- ds_rows[as.character(ds_rows$method) == method, , drop = FALSE]
      if (!nrow(method_rows)) next
      scopes <- unique(as.character(method_rows$scope))
      scopes <- scopes[!is.na(scopes) & nzchar(scopes)]
      for (scope in scopes) {
        scope_rows <- method_rows[as.character(method_rows$scope) == scope, , drop = FALSE]
        if (!nrow(scope_rows)) next
        components <- if (identical(scope, "component") && "component" %in% names(scope_rows)) {
          unique(as.character(scope_rows$component))
        } else {
          NA_character_
        }
        components <- components[!is.na(components) & nzchar(components)]
        if (!length(components)) components <- NA_character_
        for (component in components) {
          plot_rows <- scope_rows
          if (!is.na(component) && "component" %in% names(plot_rows)) {
            plot_rows <- plot_rows[as.character(plot_rows$component) == component, , drop = FALSE]
          }
          if (!nrow(plot_rows)) next

          stem <- if (isTRUE(integrated)) {
            if (!is.na(component)) {
              sprintf("benchmark_loading_enrichment_%s_%s_%s_%s", safe_chr(ds), safe_chr(method), safe_chr(scope), safe_chr(component))
            } else {
              sprintf("benchmark_loading_enrichment_%s_%s_%s", safe_chr(ds), safe_chr(method), safe_chr(scope))
            }
          } else {
            if (!is.na(component)) {
              sprintf("benchmark_loading_enrichment_%s_%s_%s", safe_chr(method), safe_chr(scope), safe_chr(component))
            } else {
              sprintf("benchmark_loading_enrichment_%s_%s", safe_chr(method), safe_chr(scope))
            }
          }

          context <- if (isTRUE(integrated)) {
            if (!is.na(component)) {
              sprintf("benchmark loading enrichment [%s / %s / %s / %s]", ds, method, scope, component)
            } else {
              sprintf("benchmark loading enrichment [%s / %s / %s]", ds, method, scope)
            }
          } else {
            if (!is.na(component)) {
              sprintf("benchmark loading enrichment [%s / %s / %s]", method, scope, component)
            } else {
              sprintf("benchmark loading enrichment [%s / %s]", method, scope)
            }
          }
          enr_obj <- .benchmark_report_enrichment_dotplot_object(
            list(results = as.data.frame(plot_rows)),
            context = context,
            verbose = verbose
          )
          if (is.null(enr_obj)) next

          #.benchmark_report_message("Writing ", context, " (", nrow(enr_obj$results), " rows).", verbose = verbose)
          files <- c(files, .benchmark_report_safe_plot(
            plot_enrichment_dotplot_fn(
              enr_obj,
              top_n = enrichment_top_n
            ),
            .benchmark_report_plot_file(outdir, stem, figure_format),
            width = 9,
            height = 6,
            verbose = verbose
          ))
        }
      }
    }
  }
  files[!is.na(files)]
}

#' If `target_sets` is supplied, target-set prioritization results are
#' computed automatically with `evaluate_target_set_auc()`. If the benchmark
#' object already contains `target_auc`, the stored result is used.
runSingleBenchmarkReports <- function(
    rw_single,
    outdir = "benchmark_reports/single_real_world",
    methods = NULL,
    components = NULL,
    u_plot_specs = NULL,
    color = NULL,
    group = NULL,
    line_group = NULL,
    xaxis = NULL,
    facet = NULL,
    target_sets = NULL,
    target_set_auc_args = list(),
    include_svdx_specific = TRUE,
    features_top_n = 50L,
    feature_profiles_top_n = 5L,
    enrichment_top_n = 20L,
    volcano_labels_top_n = 25L,
    volcano_label_features = NULL,
    inference_scope = c("all", "component", "union", "global"),
    loading_heatmaps = TRUE,
    loading_feature_profiles = TRUE,
    loading_enrichment = TRUE,
    inference_volcanoes = TRUE,
    inference_heatmaps = TRUE,
    inference_feature_profiles = TRUE,
    inference_dotplots = TRUE,
    inference_enrichment = TRUE,
    feature_profile_summarize = c("mean", "median", "none"),
    feature_profile_show_se = TRUE,
    svdx_report_args = list(),
    figure_format = "pdf",
    overwrite = TRUE,
    verbose = TRUE
) {
  inference_scope <- match.arg(inference_scope)
  feature_profile_summarize <- match.arg(feature_profile_summarize)
  figure_format <- .benchmark_report_validate_figure_format(figure_format)
  if (!inherits(rw_single, "realWorldMethodBenchmark")) {
    stop("rw_single must inherit from 'realWorldMethodBenchmark'.", call. = FALSE)
  }
  outdir <- .benchmark_report_prepare_dir(outdir, overwrite = overwrite)
  methods <- methods %||% .benchmark_report_methods(rw_single, integrated = FALSE)
  u_plot_specs <- .benchmark_report_prepare_u_plot_specs(
    rw_single$data, xformula = rw_single$xformula, u_plot_specs = u_plot_specs, xaxis = xaxis
  )
  xaxis <- xaxis %||% .benchmark_report_first_xaxis(rw_single$xformula, rw_single$data, u_plot_specs)

  target_auc <- .benchmark_report_compute_target_auc(
    rw_obj = rw_single,
    target_sets = target_sets,
    integrated = FALSE,
    target_set_auc_args = target_set_auc_args
  )

  saveRDS(rw_single, .benchmark_report_data_file(outdir, "rw_single.rds"))
  .benchmark_report_write_prediction_metrics(outdir, rw_single$prediction_metrics)
  files <- character(0)
  if (!is.null(rw_single$prediction_metrics)) {
    files <- c(files, .benchmark_report_safe_plot(
      .benchmark_report_plot_prediction_metrics(rw_single$prediction_metrics),
      .benchmark_report_plot_file(outdir, "test_explained_variation_svdx_predictions", figure_format),
      width = 7, height = 5, verbose = verbose, figure_format = figure_format
    ))
  }

  .benchmark_report_message("Writing single real-world U plots", verbose = verbose)
  files <- c(files, .benchmark_report_save_u_panels(
    rw_single, outdir = outdir, methods = methods, components = components,
    u_plot_specs = u_plot_specs, integrated = FALSE, color = color, group = group, line_group = line_group, verbose = verbose, figure_format = figure_format
  ))

  .benchmark_report_message("Writing single real-world loading/V analyses", verbose = verbose)
  files <- c(files, .benchmark_report_save_generic_loading_analyses(
    rw_single, outdir = outdir, methods = methods, components = components,
    integrated = FALSE, xaxis = xaxis,
    features_top_n = features_top_n, feature_profiles_top_n = feature_profiles_top_n,
    loading_heatmaps = loading_heatmaps, loading_feature_profiles = loading_feature_profiles,
    feature_profile_summarize = feature_profile_summarize,
    feature_profile_show_se = feature_profile_show_se,
    plot_args = svdx_report_args,
    figure_format = figure_format, verbose = verbose
  ))
  files <- c(files, .benchmark_report_write_loading_tables(
    rw_single, outdir = outdir, methods = methods, components = components,
    integrated = FALSE, verbose = verbose
  ))
  files <- c(files, .benchmark_report_write_inference_tables(
    rw_single, outdir = outdir, methods = names(rw_single$inference %||% list()),
    inference_scope = inference_scope, verbose = verbose
  ))

  if (!is.null(target_auc)) {
    files <- c(files, .benchmark_report_safe_plot(
      plot_target_set_auc_curve(target_auc, y_limits = NULL),
      .benchmark_report_plot_file(outdir, "relevant_subset_prioritization_curve", figure_format),
      width = 9, height = 6, verbose = verbose, figure_format = figure_format
    ))
    files <- c(files, .benchmark_report_safe_plot(
      plot_target_set_auc(target_auc, y_limits = NULL),
      .benchmark_report_plot_file(outdir, "relevant_subset_prioritization_summary", figure_format),
      width = 9, height = 6, verbose = verbose
    ))
  }

  if (isTRUE(loading_enrichment)) {
    files <- c(files, .benchmark_report_plot_real_world_loading_enrichment(
      rw_single,
      outdir = outdir,
      methods = NULL,
      integrated = FALSE,
      enrichment_top_n = enrichment_top_n,
      figure_format = figure_format,
      verbose = verbose
    ))
  }

  if (isTRUE(include_svdx_specific) && "svdx" %in% names(rw_single$raw_fits %||% list())) {
    svdx_dir <- file.path(outdir, "SVDX_specific_analyses")
    dir.create(svdx_dir, showWarnings = FALSE, recursive = TRUE)
    files <- c(files, .benchmark_report_call_svdx_workflow_plots(
      rw_single,
      outdir = svdx_dir,
      xaxes = vapply(u_plot_specs, function(sp) .benchmark_report_spec_value(sp, "var", .benchmark_report_spec_value(sp, "x", NA_character_)), character(1)),
      components = components,
      color = color,
      group = group,
      line_group = line_group,
      facet = facet,
      u_plot_specs = u_plot_specs,
      features_top_n = features_top_n,
      feature_profiles_top_n = feature_profiles_top_n,
      enrichment_top_n = enrichment_top_n,
      volcano_labels_top_n = volcano_labels_top_n,
      volcano_label_features = volcano_label_features,
      inference_scope = inference_scope,
      loading_heatmaps = loading_heatmaps,
      loading_feature_profiles = loading_feature_profiles,
      loading_enrichment = loading_enrichment,
      inference_volcanoes = inference_volcanoes,
      inference_heatmaps = inference_heatmaps,
      inference_feature_profiles = inference_feature_profiles,
      inference_dotplots = inference_dotplots,
      inference_enrichment = inference_enrichment,
      feature_profile_summarize = feature_profile_summarize,
      feature_profile_show_se = feature_profile_show_se,
      plot_args = svdx_report_args,
      figure_format = figure_format,
      verbose = verbose
    ))
    files <- c(files, .benchmark_report_write_loading_tables(
      rw_single, outdir = svdx_dir, methods = "svdx", components = components,
      integrated = FALSE, verbose = verbose
    ))
    files <- c(files, .benchmark_report_write_inference_tables(
      rw_single, outdir = svdx_dir,
      methods = intersect("svdx", names(rw_single$inference %||% list())),
      inference_scope = inference_scope, verbose = verbose
    ))
  }

  index <- data.frame(file = files, stringsAsFactors = FALSE)
  utils::write.csv(index, .benchmark_report_data_file(outdir, "report_index.csv"), row.names = FALSE)
  invisible(list(outdir = outdir, files = files, index = index))
}

#' Run reports for an integrated real-world benchmark object.
#'
#' If `target_sets` is supplied, target-set prioritization results are
#' computed automatically with `evaluate_integrated_target_set_auc()`. If the
#' benchmark object already contains `target_auc`, the stored result is used.
runIntegratedBenchmarkReports <- function(
    rw_integrated,
    outdir = "benchmark_reports/integrated_real_world",
    methods = NULL,
    components = NULL,
    u_plot_specs = NULL,
    dataset = NULL,
    color = NULL,
    group = NULL,
    line_group = NULL,
    xaxis = NULL,
    facet = NULL,
    target_sets = NULL,
    target_set_auc_args = list(),
    include_isvdx_specific = TRUE,
    features_top_n = 50L,
    feature_profiles_top_n = 5L,
    enrichment_top_n = 20L,
    volcano_labels_top_n = 25L,
    volcano_label_features = NULL,
    inference_scope = c("all", "component", "union", "global"),
    loading_heatmaps = TRUE,
    loading_feature_profiles = TRUE,
    loading_enrichment = TRUE,
    inference_volcanoes = TRUE,
    inference_heatmaps = TRUE,
    inference_feature_profiles = TRUE,
    inference_dotplots = TRUE,
    inference_enrichment = TRUE,
    feature_profile_summarize = c("mean", "median", "none"),
    feature_profile_show_se = TRUE,
    isvdx_report_args = list(),
    figure_format = "pdf",
    overwrite = TRUE,
    verbose = TRUE
) {
  inference_scope <- match.arg(inference_scope)
  feature_profile_summarize <- match.arg(feature_profile_summarize)
  figure_format <- .benchmark_report_validate_figure_format(figure_format)
  if (!inherits(rw_integrated, "realWorldIntegratedBenchmark")) {
    stop("rw_integrated must inherit from 'realWorldIntegratedBenchmark'.", call. = FALSE)
  }
  outdir <- .benchmark_report_prepare_dir(outdir, overwrite = overwrite)
  methods <- methods %||% .benchmark_report_methods(rw_integrated, integrated = TRUE)
  u_plot_specs <- .benchmark_report_prepare_u_plot_specs(
    rw_integrated$xdata, xformula = rw_integrated$xformula, u_plot_specs = u_plot_specs, xaxis = xaxis
  )
  xaxis <- xaxis %||% .benchmark_report_first_xaxis(rw_integrated$xformula, rw_integrated$xdata, u_plot_specs)

  target_auc <- .benchmark_report_compute_target_auc(
    rw_obj = rw_integrated,
    target_sets = target_sets,
    integrated = TRUE,
    target_set_auc_args = target_set_auc_args
  )

  saveRDS(rw_integrated, .benchmark_report_data_file(outdir, "rw_integrated.rds"))
  .benchmark_report_write_prediction_metrics(outdir, rw_integrated$prediction_metrics)
  files <- character(0)
  if (!is.null(rw_integrated$prediction_metrics)) {
    files <- c(files, .benchmark_report_safe_plot(
      .benchmark_report_plot_prediction_metrics(rw_integrated$prediction_metrics),
      .benchmark_report_plot_file(outdir, "test_explained_variation_svdx_predictions", figure_format),
      width = 8, height = 5.5, verbose = verbose, figure_format = figure_format
    ))
  }

  .benchmark_report_message("Writing integrated real-world U plots", verbose = verbose)
  files <- c(files, .benchmark_report_save_u_panels(
    rw_integrated, outdir = outdir, methods = methods, components = components,
    u_plot_specs = u_plot_specs, integrated = TRUE, dataset = dataset,
    color = color, group = group, line_group = line_group, verbose = verbose, figure_format = figure_format
  ))

  .benchmark_report_message("Writing integrated real-world loading/V analyses", verbose = verbose)
  files <- c(files, .benchmark_report_save_generic_loading_analyses(
    rw_integrated, outdir = outdir, methods = methods, components = components,
    integrated = TRUE, datasets = dataset %||% names(rw_integrated$Y), xaxis = xaxis,
    features_top_n = features_top_n, feature_profiles_top_n = feature_profiles_top_n,
    loading_heatmaps = loading_heatmaps, loading_feature_profiles = loading_feature_profiles,
    feature_profile_summarize = feature_profile_summarize,
    feature_profile_show_se = feature_profile_show_se,
    plot_args = isvdx_report_args,
    figure_format = figure_format, verbose = verbose
  ))
  files <- c(files, .benchmark_report_write_loading_tables(
    rw_integrated, outdir = outdir, methods = methods, components = components,
    integrated = TRUE, datasets = dataset %||% names(rw_integrated$Y), verbose = verbose
  ))
  files <- c(files, .benchmark_report_write_inference_tables(
    rw_integrated, outdir = outdir, methods = names(rw_integrated$inference %||% list()),
    inference_scope = inference_scope, verbose = verbose
  ))

  if (!is.null(target_auc)) {
    files <- c(files, .benchmark_report_safe_plot(
      plot_target_set_auc_curve(target_auc, facet_by = "dataset_target_set", y_limits = NULL),
      .benchmark_report_plot_file(outdir, "relevant_subset_prioritization_curve", figure_format),
      width = 10, height = 7, verbose = verbose
    ))
    files <- c(files, .benchmark_report_safe_plot(
      plot_target_set_auc(target_auc, facet_by = "dataset_target_set", y_limits = NULL),
      .benchmark_report_plot_file(outdir, "relevant_subset_prioritization_summary", figure_format),
      width = 10, height = 7, verbose = verbose
    ))
  }

  if (isTRUE(loading_enrichment)) {
    files <- c(files, .benchmark_report_plot_real_world_loading_enrichment(
      rw_integrated,
      outdir = outdir,
      methods = NULL,
      datasets = dataset,
      integrated = TRUE,
      enrichment_top_n = enrichment_top_n,
      figure_format = figure_format,
      verbose = verbose
    ))
  }

  if (isTRUE(include_isvdx_specific) && any(c("isvdx", "svdx", "separate_svdx") %in% names(rw_integrated$results %||% list()))) {
    isvdx_dir <- file.path(outdir, "ISVDX_specific_analyses")
    dir.create(isvdx_dir, showWarnings = FALSE, recursive = TRUE)
    files <- c(files, .benchmark_report_call_isvdx_workflow_plots(
      rw_integrated,
      outdir = isvdx_dir,
      xaxes = vapply(u_plot_specs, function(sp) .benchmark_report_spec_value(sp, "var", .benchmark_report_spec_value(sp, "x", NA_character_)), character(1)),
      components = components,
      datasets = dataset %||% names(rw_integrated$Y),
      color = color,
      group = group,
      line_group = line_group,
      facet = facet,
      u_plot_specs = u_plot_specs,
      features_top_n = features_top_n,
      feature_profiles_top_n = feature_profiles_top_n,
      enrichment_top_n = enrichment_top_n,
      volcano_labels_top_n = volcano_labels_top_n,
      volcano_label_features = volcano_label_features,
      inference_scope = inference_scope,
      loading_heatmaps = loading_heatmaps,
      loading_feature_profiles = loading_feature_profiles,
      loading_enrichment = loading_enrichment,
      inference_volcanoes = inference_volcanoes,
      inference_heatmaps = inference_heatmaps,
      inference_feature_profiles = inference_feature_profiles,
      inference_dotplots = inference_dotplots,
      inference_enrichment = inference_enrichment,
      feature_profile_summarize = feature_profile_summarize,
      feature_profile_show_se = feature_profile_show_se,
      plot_args = isvdx_report_args,
      figure_format = figure_format,
      verbose = verbose
    ))
    for (ds_export in dataset %||% names(rw_integrated$Y)) {
      ds_export_dir <- file.path(isvdx_dir, .benchmark_report_safe_name(ds_export))
      files <- c(files, .benchmark_report_write_loading_tables(
        rw_integrated, outdir = ds_export_dir,
        methods = intersect("isvdx", names(rw_integrated$results %||% list())),
        components = components, integrated = TRUE,
        datasets = ds_export, verbose = verbose
      ))
    }
    files <- c(files, .benchmark_report_write_inference_tables(
      rw_integrated, outdir = isvdx_dir,
      methods = intersect("isvdx", names(rw_integrated$inference %||% list())),
      inference_scope = inference_scope, verbose = verbose
    ))
  }

  index <- data.frame(file = files, stringsAsFactors = FALSE)
  utils::write.csv(index, .benchmark_report_data_file(outdir, "report_index.csv"), row.names = FALSE)
  invisible(list(outdir = outdir, files = files, index = index))
}

#' Run simulation benchmark reports.
runSimulationBenchmarkReports <- function(
    single_grid = NULL,
    integrated_grid = NULL,
    outdir = "benchmark_reports/simulation",
    figure_format = "pdf",
    overwrite = TRUE,
    include_module_plots = FALSE,
    include_moduleType_plots = TRUE,
    enrichment_metrics = c("odds_ratio", "overlap", "p_value", "q_value", "neg_log10_p_value", "neg_log10_q_value"),
    verbose = TRUE
) {
  figure_format <- .benchmark_report_validate_figure_format(figure_format)
  outdir <- .benchmark_report_prepare_dir(outdir, overwrite = overwrite)
  enrichment_metrics <- intersect(
    as.character(enrichment_metrics),
    c("odds_ratio", "overlap", "p_value", "q_value", "neg_log10_p_value", "neg_log10_q_value")
  )
  files <- character(0)
  if (!is.null(single_grid)) {
    sd <- file.path(outdir, "single")
    dir.create(sd, showWarnings = FALSE, recursive = TRUE)
    saveRDS(single_grid, .benchmark_report_data_file(sd, "single_grid.rds"))
    .benchmark_report_write_prediction_metrics(sd, single_grid$prediction_metrics)
    if (!is.null(single_grid$prediction_metrics)) {
      files <- c(files, .benchmark_report_safe_plot(
        .benchmark_report_plot_prediction_metrics(single_grid$prediction_metrics),
        .benchmark_report_plot_file(sd, "test_explained_variation_svdx_predictions", figure_format),
        width = 7, height = 5, verbose = verbose, figure_format = figure_format
      ))
    }
    fns <- c("plot_single_factorization_metrics", "plot_single_runtime")
    if (isTRUE(include_moduleType_plots)) {
      fns <- c(
        fns,
        "plot_single_module_type_auc",
        "plot_single_module_type_recall",
        "plot_single_inference_module_type_auc",
        "plot_single_inference_module_type_recall"
      )
    }
    if (isTRUE(include_module_plots)) {
      fns <- c(
        fns,
        "plot_single_module_auc",
        "plot_single_module_recall_summary",
        "plot_single_inference_module_auc",
        "plot_single_inference_module_recall_summary"
      )
    }
    if (!is.null(single_grid$inference_metrics) && nrow(single_grid$inference_metrics) > 0L) {
      fns <- c(fns, "plot_single_inference_metrics")
    }
    for (fn_name in fns) if (exists(fn_name, mode = "function")) {
      files <- c(files, .benchmark_report_safe_plot_collection(
        get(fn_name, mode = "function")(single_grid),
        sd,
        fn_name,
        figure_format = figure_format,
        verbose = verbose
      ))
    }
    if (isTRUE(include_module_plots)) {
      if (exists("plot_single_loading_enrichment", mode = "function")) {
        for (metric in enrichment_metrics) {
          files <- c(files, .benchmark_report_safe_plot_collection(
            plot_single_loading_enrichment(single_grid, metric = metric),
            sd,
            paste0("plot_single_loading_enrichment_", metric),
            figure_format = figure_format,
            verbose = verbose
          ))
        }
      }
      if (exists("plot_single_inference_module_enrichment", mode = "function")) {
        for (metric in enrichment_metrics) {
          files <- c(files, .benchmark_report_safe_plot_collection(
            plot_single_inference_module_enrichment(single_grid, metric = metric, exclude_feature_lm = FALSE),
            sd,
            paste0("plot_single_inference_module_enrichment_", metric),
            figure_format = figure_format,
            verbose = verbose
          ))
        }
      }
    }
    if (isTRUE(include_moduleType_plots)) {
      if (exists("plot_single_loading_module_type_enrichment", mode = "function")) {
        for (metric in enrichment_metrics) {
          files <- c(files, .benchmark_report_safe_plot_collection(
            plot_single_loading_module_type_enrichment(single_grid, metric = metric),
            sd,
            paste0("plot_single_loading_module_type_enrichment_", metric),
            figure_format = figure_format,
            verbose = verbose
          ))
        }
      }
      if (exists("plot_single_inference_module_type_enrichment", mode = "function")) {
        for (metric in enrichment_metrics) {
          files <- c(files, .benchmark_report_safe_plot_collection(
            plot_single_inference_module_type_enrichment(single_grid, metric = metric),
            sd,
            paste0("plot_single_inference_module_type_enrichment_", metric),
            figure_format = figure_format,
            verbose = verbose
          ))
        }
      }
    }
  }
  if (!is.null(integrated_grid)) {
    id <- file.path(outdir, "integrated")
    dir.create(id, showWarnings = FALSE, recursive = TRUE)
    saveRDS(integrated_grid, .benchmark_report_data_file(id, "integrated_grid.rds"))
    .benchmark_report_write_prediction_metrics(id, integrated_grid$prediction_metrics)
    if (!is.null(integrated_grid$prediction_metrics)) {
      files <- c(files, .benchmark_report_safe_plot(
        .benchmark_report_plot_prediction_metrics(integrated_grid$prediction_metrics),
        .benchmark_report_plot_file(id, "test_explained_variation_svdx_predictions", figure_format),
        width = 8, height = 5.5, verbose = verbose, figure_format = figure_format
      ))
    }
    fns <- c("plot_integrated_factorization_metrics", "plot_integrated_runtime")
    if (isTRUE(include_moduleType_plots)) {
      fns <- c(
        fns,
        "plot_integrated_module_type_auc",
        "plot_integrated_module_type_recall",
        "plot_integrated_inference_module_type_auc",
        "plot_integrated_inference_module_type_recall"
      )
    }
    if (isTRUE(include_module_plots)) {
      fns <- c(
        fns,
        "plot_integrated_module_auc",
        "plot_integrated_module_recall_summary",
        "plot_integrated_inference_module_auc",
        "plot_integrated_inference_module_recall_summary"
      )
    }
    if (!is.null(integrated_grid$inference_metrics) && nrow(integrated_grid$inference_metrics) > 0L) {
      fns <- c(fns, "plot_integrated_inference_metrics")
    }
    for (fn_name in fns) if (exists(fn_name, mode = "function")) {
      files <- c(files, .benchmark_report_safe_plot_collection(
        get(fn_name, mode = "function")(integrated_grid),
        id,
        fn_name,
        figure_format = figure_format,
        verbose = verbose
      ))
    }
    if (isTRUE(include_module_plots)) {
      if (exists("plot_integrated_loading_enrichment", mode = "function")) {
        for (metric in enrichment_metrics) {
          files <- c(files, .benchmark_report_safe_plot_collection(
            plot_integrated_loading_enrichment(integrated_grid, metric = metric),
            id,
            paste0("plot_integrated_loading_enrichment_", metric),
            figure_format = figure_format,
            verbose = verbose
          ))
        }
      }
      if (exists("plot_integrated_inference_module_enrichment", mode = "function")) {
        for (metric in enrichment_metrics) {
          files <- c(files, .benchmark_report_safe_plot_collection(
            plot_integrated_inference_module_enrichment(integrated_grid, metric = metric),
            id,
            paste0("plot_integrated_inference_module_enrichment_", metric),
            figure_format = figure_format,
            verbose = verbose
          ))
        }
      }
    }
    if (isTRUE(include_moduleType_plots)) {
      if (exists("plot_integrated_loading_module_type_enrichment", mode = "function")) {
        for (metric in enrichment_metrics) {
          files <- c(files, .benchmark_report_safe_plot_collection(
            plot_integrated_loading_module_type_enrichment(integrated_grid, metric = metric),
            id,
            paste0("plot_integrated_loading_module_type_enrichment_", metric),
            figure_format = figure_format,
            verbose = verbose
          ))
        }
      }
      if (exists("plot_integrated_inference_module_type_enrichment", mode = "function")) {
        for (metric in enrichment_metrics) {
          files <- c(files, .benchmark_report_safe_plot_collection(
            plot_integrated_inference_module_type_enrichment(integrated_grid, metric = metric),
            id,
            paste0("plot_integrated_inference_module_type_enrichment_", metric),
            figure_format = figure_format,
            verbose = verbose
          ))
        }
      }
    }
  }
  index <- data.frame(file = files, stringsAsFactors = FALSE)
  utils::write.csv(index, .benchmark_report_data_file(outdir, "report_index.csv"), row.names = FALSE)
  invisible(list(outdir = outdir, files = files, index = index))
}


#' Run reports for one or more real-world benchmark objects.
#'
#' `rw_single` and `rw_integrated` may each be either a single completed
#' benchmark object or a named list of completed benchmark objects.  Each
#' dataset is written to its own subdirectory under `outdir`.
#'
#' `single_target_sets` and `integrated_target_sets` are user-facing
#' convenience inputs. If the benchmark object does not already contain
#' stored `target_auc` results, target-set prioritization is evaluated
#' automatically before plotting.
runRealWorldBenchmarkReports <- function(
    rw_single = NULL,
    rw_integrated = NULL,
    single_target_sets = NULL,
    integrated_target_sets = NULL,
    target_set_auc_args = list(),
    outdir = "benchmark_reports/real_world",
    overwrite = TRUE,
    verbose = TRUE,
    figure_format = "pdf",
    methods = NULL,
    components = 1:5,
    u_plot_specs = NULL,
    color = NULL,
    group = NULL,
    line_group = NULL,
    xaxis = NULL,
    facet = NULL,
    dataset = NULL,
    include_svdx_specific = TRUE,
    include_isvdx_specific = TRUE,
    features_top_n = 50L,
    feature_profiles_top_n = 5L,
    enrichment_top_n = 20L,
    volcano_labels_top_n = 25L,
    volcano_label_features = NULL,
    inference_scope = c("all", "component", "union", "global"),
    loading_heatmaps = TRUE,
    loading_feature_profiles = TRUE,
    loading_enrichment = TRUE,
    inference_volcanoes = TRUE,
    inference_heatmaps = TRUE,
    inference_feature_profiles = TRUE,
    inference_dotplots = TRUE,
    inference_enrichment = TRUE,
    feature_profile_summarize = c("mean", "median", "none"),
    feature_profile_show_se = TRUE,
    svdx_report_args = list(),
    isvdx_report_args = list(),
    single_report_args = list(),
    integrated_report_args = list(),
    ...
) {
  inference_scope <- match.arg(inference_scope)
  feature_profile_summarize <- match.arg(feature_profile_summarize)
  figure_format <- .benchmark_report_validate_figure_format(figure_format)
  if (!is.list(target_set_auc_args)) {
    stop("target_set_auc_args must be a list.", call. = FALSE)
  }
  outdir <- .benchmark_report_prepare_dir(outdir, overwrite = overwrite)
  out <- list(outdir = outdir)

  merge_args <- function(...) {
    parts <- list(...)
    out <- list()
    for (part in parts) {
      if (is.null(part) || !length(part)) next
      nms <- names(part)
      if (is.null(nms)) next
      for (nm in nms[nzchar(nms)]) out[[nm]] <- part[[nm]]
    }
    out
  }

  filter_supported_args <- function(args, fun) {
    if (is.null(args) || !length(args)) return(list())
    args[names(args) %in% names(formals(fun))]
  }

  dots <- list(...)
  common_report_args <- list(
    methods = methods,
    components = components,
    u_plot_specs = u_plot_specs,
    color = color,
    group = group,
    line_group = line_group,
    xaxis = xaxis,
    facet = facet,
    features_top_n = features_top_n,
    feature_profiles_top_n = feature_profiles_top_n,
    enrichment_top_n = enrichment_top_n,
    volcano_labels_top_n = volcano_labels_top_n,
    volcano_label_features = volcano_label_features,
    inference_scope = inference_scope,
    loading_heatmaps = loading_heatmaps,
    loading_feature_profiles = loading_feature_profiles,
    loading_enrichment = loading_enrichment,
    inference_volcanoes = inference_volcanoes,
    inference_heatmaps = inference_heatmaps,
    inference_feature_profiles = inference_feature_profiles,
    inference_dotplots = inference_dotplots,
    inference_enrichment = inference_enrichment,
    feature_profile_summarize = feature_profile_summarize,
    feature_profile_show_se = feature_profile_show_se
  )

  if (!is.null(rw_single)) {
    single_list <- .benchmark_report_named_list(rw_single, default_prefix = "single")
    out$single <- vector("list", length(single_list))
    names(out$single) <- names(single_list)
    single_dots <- filter_supported_args(dots, runSingleBenchmarkReports)
    for (nm in names(single_list)) {
      single_args <- merge_args(
        list(
          rw_single = single_list[[nm]],
          target_sets = single_target_sets,
          target_set_auc_args = target_set_auc_args,
          outdir = file.path(outdir, "single", .benchmark_report_safe_name(nm)),
          overwrite = TRUE,
          verbose = verbose,
          figure_format = figure_format
        ),
        common_report_args,
        list(
          include_svdx_specific = include_svdx_specific,
          svdx_report_args = svdx_report_args
        ),
        single_dots,
        single_report_args
      )
      out$single[[nm]] <- do.call(runSingleBenchmarkReports, single_args)
    }
  }

  if (!is.null(rw_integrated)) {
    integrated_list <- .benchmark_report_named_list(rw_integrated, default_prefix = "integrated")
    out$integrated <- vector("list", length(integrated_list))
    names(out$integrated) <- names(integrated_list)
    integrated_dots <- filter_supported_args(dots, runIntegratedBenchmarkReports)
    for (nm in names(integrated_list)) {
      integrated_args <- merge_args(
        list(
          rw_integrated = integrated_list[[nm]],
          target_sets = integrated_target_sets,
          target_set_auc_args = target_set_auc_args,
          outdir = file.path(outdir, "integrated", .benchmark_report_safe_name(nm)),
          overwrite = TRUE,
          verbose = verbose,
          figure_format = figure_format
        ),
        common_report_args,
        list(
          dataset = dataset,
          include_isvdx_specific = include_isvdx_specific,
          isvdx_report_args = isvdx_report_args
        ),
        integrated_dots,
        integrated_report_args
      )
      out$integrated[[nm]] <- do.call(runIntegratedBenchmarkReports, integrated_args)
    }
  }

  invisible(out)
}

#' Dispatch benchmark reports for any available benchmark objects.
#'
#' For real-world reports, `single_target_sets` and
#' `integrated_target_sets` are forwarded to `runRealWorldBenchmarkReports()`
#' and evaluated automatically using stored benchmark results when available.
runBenchmarkReports <- function(
    single_grid = NULL,
    integrated_grid = NULL,
    rw_single = NULL,
    rw_integrated = NULL,
    single_target_sets = NULL,
    integrated_target_sets = NULL,
    target_set_auc_args = list(),
    outdir = "benchmark_reports",
    overwrite = TRUE,
    verbose = TRUE,
    figure_format = "pdf",
    features_top_n = 50L,
    feature_profiles_top_n = 5L,
    enrichment_top_n = 20L,
    volcano_labels_top_n = 25L,
    volcano_label_features = NULL,
    inference_scope = c("all", "component", "union", "global"),
    loading_heatmaps = TRUE,
    loading_feature_profiles = TRUE,
    loading_enrichment = TRUE,
    inference_volcanoes = TRUE,
    inference_heatmaps = TRUE,
    inference_feature_profiles = TRUE,
    inference_dotplots = TRUE,
    inference_enrichment = TRUE,
    feature_profile_summarize = c("mean", "median", "none"),
    feature_profile_show_se = TRUE,
    include_module_plots = FALSE,
    include_moduleType_plots = TRUE,
    enrichment_metrics = c("odds_ratio", "overlap", "p_value", "q_value", "neg_log10_p_value", "neg_log10_q_value"),
    ...
) {
  inference_scope <- match.arg(inference_scope)
  feature_profile_summarize <- match.arg(feature_profile_summarize)
  figure_format <- .benchmark_report_validate_figure_format(figure_format)
  outdir <- .benchmark_report_prepare_dir(outdir, overwrite = overwrite)
  out <- list(outdir = outdir)
  if (!is.null(single_grid) || !is.null(integrated_grid)) {
    out$simulation <- runSimulationBenchmarkReports(
      single_grid = single_grid, integrated_grid = integrated_grid,
      outdir = file.path(outdir, "simulation"), overwrite = TRUE,
      verbose = verbose, figure_format = figure_format,
      include_module_plots = include_module_plots,
      include_moduleType_plots = include_moduleType_plots,
      enrichment_metrics = enrichment_metrics
    )
  }
  if (!is.null(rw_single) || !is.null(rw_integrated)) {
    out$real_world <- runRealWorldBenchmarkReports(
      rw_single = rw_single,
      rw_integrated = rw_integrated,
      single_target_sets = single_target_sets,
      integrated_target_sets = integrated_target_sets,
      target_set_auc_args = target_set_auc_args,
      outdir = file.path(outdir, "real_world"),
      overwrite = TRUE,
      verbose = verbose,
      figure_format = figure_format,
      features_top_n = features_top_n,
      feature_profiles_top_n = feature_profiles_top_n,
      enrichment_top_n = enrichment_top_n,
      volcano_labels_top_n = volcano_labels_top_n,
      volcano_label_features = volcano_label_features,
      inference_scope = inference_scope,
      loading_heatmaps = loading_heatmaps,
      loading_feature_profiles = loading_feature_profiles,
      loading_enrichment = loading_enrichment,
      inference_volcanoes = inference_volcanoes,
      inference_heatmaps = inference_heatmaps,
      inference_feature_profiles = inference_feature_profiles,
      inference_dotplots = inference_dotplots,
      inference_enrichment = inference_enrichment,
      feature_profile_summarize = feature_profile_summarize,
      feature_profile_show_se = feature_profile_show_se,
      ...
    )
  }
  invisible(out)
}
