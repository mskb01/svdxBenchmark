# ============================================================
# Benchmark plot helpers
# ============================================================
# Central rule:
#   Variables listed in average_over are removed from aggregation identity,
#   split identity, plot keys, titles, and facets. They should not appear in
#   plot names or labels after averaging.
# ============================================================

.bp_check_columns <- function(df, required, name = deparse(substitute(df))) {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(name, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.bp_clean_key <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[nchar(x) == 0] <- "all"
  x
}

.bp_value_label <- function(var, value) {
  if (var %in% c("time_imbalance")) {
    return(paste0("timeimb_", .bp_clean_key(value)))
  }
  if (var %in% c("batch_to_biology_ratio", "batch_to_biology_ratio_name")) {
    return(.bp_clean_key(value))
  }
  if (var %in% c("confounding", "confounding_name")) {
    return(paste0("conf_", .bp_clean_key(value)))
  }
  if (var %in% c("count_regime_mean")) {
    return(paste0("countmean_", .bp_clean_key(value)))
  }
  if (var %in% c("family")) {
    return(paste0("family_", .bp_clean_key(value)))
  }
  .bp_clean_key(value)
}

.bp_key_from_values <- function(values) {
  if (length(values) == 0) return("all")
  paste(vapply(names(values), function(nm) {
    .bp_value_label(nm, values[[nm]])
  }, character(1)), collapse = "__")
}

.bp_key_from_df <- function(df, split_vars) {
  split_vars <- intersect(split_vars, names(df))
  if (length(split_vars) == 0) return("all")
  vals <- lapply(split_vars, function(v) unique(df[[v]])[1])
  names(vals) <- split_vars
  .bp_key_from_values(vals)
}

.bp_title_from_df <- function(df, split_vars, prefix = NULL) {
  split_vars <- intersect(split_vars, names(df))
  pieces <- vapply(split_vars, function(v) {
    paste0(v, " = ", unique(df[[v]])[1])
  }, character(1))
  out <- paste(c(prefix, pieces), collapse = " | ")
  if (identical(out, "")) "Benchmark plot" else out
}

.bp_split <- function(df, split_vars, average_over = character()) {
  split_vars <- setdiff(split_vars, average_over)
  split_vars <- intersect(split_vars, names(df))

  if (length(split_vars) == 0) {
    return(list(all = list(data = df, split_vars = character())))
  }

  key_df <- df[split_vars]
  split_key <- interaction(key_df, drop = TRUE, lex.order = TRUE, sep = "___")
  raw <- split(df, split_key)

  out <- lapply(raw, function(dd) {
    list(data = dd, split_vars = split_vars)
  })

  names(out) <- vapply(out, function(x) .bp_key_from_df(x$data, x$split_vars), character(1))
  out
}

.bp_summarise <- function(
    df,
    value_col,
    id_cols,
    average_over = "seed",
    show_se = TRUE
) {
  .bp_check_columns(df, value_col, "plot data")

  average_over <- intersect(average_over, names(df))
  id_cols <- intersect(id_cols, names(df))
  group_cols <- setdiff(id_cols, average_over)

  if (length(group_cols) == 0) {
    out <- df |>
      dplyr::summarise(
        value = mean(.data[[value_col]], na.rm = TRUE),
        sd = stats::sd(.data[[value_col]], na.rm = TRUE),
        n = sum(is.finite(.data[[value_col]])),
        se = dplyr::if_else(n > 1, sd / sqrt(n), NA_real_),
        .groups = "drop"
      )
  } else {
    out <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
      dplyr::summarise(
        value = mean(.data[[value_col]], na.rm = TRUE),
        sd = stats::sd(.data[[value_col]], na.rm = TRUE),
        n = sum(is.finite(.data[[value_col]])),
        se = dplyr::if_else(n > 1, sd / sqrt(n), NA_real_),
        .groups = "drop"
      )
  }

  if (!show_se) out$se <- NA_real_
  out
}

.bp_metrics_long <- function(df, metric_cols) {
  metric_cols <- intersect(metric_cols, names(df))
  if (length(metric_cols) == 0) {
    stop("None of the requested metric columns are present.", call. = FALSE)
  }
  tidyr::pivot_longer(
    df,
    cols = dplyr::all_of(metric_cols),
    names_to = "metric",
    values_to = "metric_value"
  )
}

.bp_line_plot <- function(
    df,
    facet_var = NULL,
    x_var = "ng",
    title = "Benchmark plot",
    ylab = "Value",
    xlab = x_var,
    log10_x = identical(x_var, "count_regime_mean")
) {
  .bp_check_columns(df, c("method", x_var), "plot data")

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[x_var]],
      y = .data$value,
      color = .data$method,
      group = .data$method
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = ylab,
      color = "Method"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )

  grouping_cols <- "method"
  if (!is.null(facet_var) && facet_var %in% names(df)) {
    grouping_cols <- c(grouping_cols, facet_var)
  }

  n_x <- df |>
    dplyr::filter(is.finite(.data$value)) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_cols))) |>
    dplyr::summarise(
      n_x = dplyr::n_distinct(.data[[x_var]]),
      .groups = "drop"
    ) |>
    dplyr::pull(.data$n_x)

  if (any(n_x >= 2)) {
    p <- p + ggplot2::geom_line(na.rm = TRUE)
  }

  if ("se" %in% names(df) && any(is.finite(df$se))) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$value - .data$se, ymax = .data$value + .data$se),
      width = 0.05,
      alpha = 0.6,
      na.rm = TRUE
    )
  }

  if (isTRUE(log10_x)) {
    positive_x <- df[[x_var]][is.finite(df[[x_var]])]
    if (length(positive_x) > 0L && all(positive_x > 0)) {
      p <- p + ggplot2::scale_x_log10()
    }
  }

  if (!is.null(facet_var) && facet_var %in% names(df)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)))
  }

  p
}

.bp_make_metric_plots <- function(
    df,
    metric_cols,
    id_cols,
    split_vars,
    facet_var,
    x_var = "ng",
    average_over = "seed",
    show_se = TRUE,
    title_prefix = NULL
) {
  long <- .bp_metrics_long(df, metric_cols)

  id_cols <- unique(c(id_cols, "metric"))
  split_vars <- unique(c("metric", split_vars))

  sm <- .bp_summarise(
    long,
    value_col = "metric_value",
    id_cols = id_cols,
    average_over = average_over,
    show_se = show_se
  )

  pieces <- .bp_split(sm, split_vars = split_vars, average_over = average_over)

  plots <- lapply(pieces, function(piece) {
    dd <- piece$data
    .bp_line_plot(
      dd,
      facet_var = facet_var,
      x_var = x_var,
      title = .bp_title_from_df(dd, piece$split_vars, prefix = title_prefix),
      ylab = unique(dd$metric)[1] %||% "Value",
      xlab = x_var
    )
  })

  names(plots) <- names(pieces)
  plots
}

.bp_make_value_plots <- function(
    df,
    value_col,
    id_cols,
    split_vars,
    facet_var,
    x_var = "ng",
    average_over = "seed",
    show_se = TRUE,
    title_prefix = NULL,
    ylab = value_col
) {
  sm <- .bp_summarise(
    df,
    value_col = value_col,
    id_cols = id_cols,
    average_over = average_over,
    show_se = show_se
  )

  pieces <- .bp_split(sm, split_vars = split_vars, average_over = average_over)

  plots <- lapply(pieces, function(piece) {
    dd <- piece$data
    .bp_line_plot(
      dd,
      facet_var = facet_var,
      x_var = x_var,
      title = .bp_title_from_df(dd, piece$split_vars, prefix = title_prefix),
      ylab = ylab,
      xlab = x_var
    )
  })

  names(plots) <- names(pieces)
  plots
}

save_plot_set <- function(plots, output_dir, width = 8, height = 5) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  invisible(lapply(names(plots), function(nm) {
    ggplot2::ggsave(
      filename = file.path(output_dir, paste0(nm, ".pdf")),
      plot = plots[[nm]],
      width = width,
      height = height
    )
  }))
}
