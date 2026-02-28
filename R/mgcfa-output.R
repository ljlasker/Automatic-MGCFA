#' Print Method for MGCFA Results
#'
#' @param x An object returned by \code{mgcfa_auto()}.
#' @param digits Number of digits used for numeric formatting.
#' @param rounding Numeric rounding mode. One of \code{"signif"},
#'   \code{"round"}, or \code{"none"}.
#' @param show_freed_parameters Logical; if \code{TRUE}, print per-step freed
#'   constraints selected by automatic partial discovery.
#' @param verbose Logical; if \code{TRUE}, print full diagnostic detail in
#'   addition to the simplified overview.
#' @param ... Unused.
#'
#' @return The input object invisibly.
#' @export
print.mgcfa_result <- function(
  x,
  digits = 3L,
  rounding = c("signif", "round", "none"),
  show_freed_parameters = FALSE,
  verbose = FALSE,
  ...
) {
  x <- .mgcfa_as_result(x)
  rounding <- match.arg(rounding)
  cat("Mode:", x$mode, "\n")
  cat("Model type:", x$model_type, "\n")
  .mgcfa_print_result_input_summary(x$input_summary %||% NULL, x$input_notes %||% character())
  cat("Steps:", paste(names(x$fits), collapse = ", "), "\n\n")

  overview <- .mgcfa_step_overview(x)
  if (nrow(overview) > 0L) {
    overview <- .mgcfa_format_numeric_df(overview, digits = digits, rounding = rounding)
    print(overview, row.names = FALSE)
  }
  if (!is.null(x$not_applicable_steps) && length(x$not_applicable_steps) > 0L) {
    cat("\nNot-applicable steps:", paste(names(x$not_applicable_steps), collapse = ", "), "\n")
  }
  if (isTRUE(x$stopped_early)) {
    cat("\nStopped early at step:", x$stopped_after_step %||% "unknown", "\n")
    if (!is.null(x$stopped_reason) && nzchar(x$stopped_reason)) {
      cat("Reason:", x$stopped_reason, "\n")
    }
  }

  if (!isTRUE(verbose)) {
    if (isTRUE(show_freed_parameters) && !is.null(x$freed_parameters) && length(x$freed_parameters) > 0L) {
      cat("\nFreed parameters by step\n")
      for (step in names(x$freed_parameters)) {
        vals <- x$freed_parameters[[step]]
        if (length(vals) > 0L) {
          cat(step, ":", paste(vals, collapse = ", "), "\n")
        }
      }
    }
    cat("\nGuidance: use `mgcfa_help_inputs()` for input recipes and `verbose = TRUE` for full diagnostics.\n")
    return(invisible(x))
  }

  cat("\nFull fit table\n")
  fit_tbl <- x$fit_table
  fit_tbl[] <- lapply(fit_tbl, .mgcfa_round_values, digits = digits, rounding = rounding)
  print(fit_tbl)

  if (!is.null(x$lrt_table) && nrow(x$lrt_table) > 0L) {
    cat("\nNested chi-square differences\n")
    lrt_tbl <- x$lrt_table
    if (is.data.frame(lrt_tbl)) {
      lrt_tbl <- .mgcfa_format_numeric_df(lrt_tbl, digits = digits, rounding = rounding)
    }
    print(lrt_tbl, row.names = FALSE)
  }

  if (!is.null(x$practical_change_table) && nrow(x$practical_change_table) > 0L) {
    cat("\nPractical fit changes by step\n")
    pc_tbl <- .mgcfa_format_numeric_df(
      x$practical_change_table,
      digits = digits,
      rounding = rounding
    )
    print(pc_tbl, row.names = FALSE)
  }

  if (!is.null(x$step_failures) && length(x$step_failures) > 0L) {
    failed_steps <- names(Filter(function(z) isTRUE(z$failed), x$step_failures))
    if (length(failed_steps) > 0L) {
      cat("\nFailed invariance steps:", paste(failed_steps, collapse = ", "), "\n")
    }
  }

  failed_outputs <- x$failed_step_outputs
  if (is.null(failed_outputs) || length(failed_outputs) == 0L) {
    failed_outputs <- Filter(function(z) isTRUE(z$failed), x$step_failures %||% list())
  }
  if (length(failed_outputs) > 0L) {
    cat("\nFailed non-partial model outputs by step\n")
    for (step in names(failed_outputs)) {
      sf <- failed_outputs[[step]]
      cat(step, "\n")
      if (!is.null(sf$criterion) && nzchar(as.character(sf$criterion))) {
        cat("  Criterion:", sf$criterion, "\n")
      }
      if (!is.null(sf$threshold) && is.finite(sf$threshold)) {
        cat(
          "  Threshold:",
          .mgcfa_round_values(sf$threshold, digits = digits, rounding = rounding),
          "\n"
        )
      }
      if (!is.null(sf$value) && is.finite(sf$value)) {
        cat(
          "  Value:",
          .mgcfa_round_values(sf$value, digits = digits, rounding = rounding),
          "\n"
        )
      }
      if (!is.null(sf$reason) && nzchar(as.character(sf$reason))) {
        cat("  Reason:", sf$reason, "\n")
      }
      if (!is.null(sf$failed_fit_measures)) {
        fail_vals <- .mgcfa_round_values(
          as.numeric(sf$failed_fit_measures),
          digits = digits,
          rounding = rounding
        )
        names(fail_vals) <- names(sf$failed_fit_measures)
        print(fail_vals)
      }
    }
  }

  single_partial <- x$partial_search
  if (!is.null(single_partial) && isTRUE(single_partial$triggered)) {
    cat("\nAutomatic partial search\n")
    cat("Decision:", single_partial$decision, "\n")
    if (!is.null(single_partial$acceptable_models)) {
      cat("Acceptable candidates found:", nrow(single_partial$acceptable_models), "\n")
    }
    if (!is.null(single_partial$top_models) && nrow(single_partial$top_models) > 0L) {
      top_tbl <- .mgcfa_format_numeric_df(
        single_partial$top_models,
        digits = digits,
        rounding = rounding
      )
      print(top_tbl, row.names = FALSE)
    }
    if (length(single_partial$selected_partial) > 0L) {
      cat("Selected partial constraints:\n")
      cat(paste(single_partial$selected_partial, collapse = ", "), "\n")
    }
    if (!is.null(single_partial$selected_stage_reached) &&
        identical(single_partial$selected_stage_reached, FALSE)) {
      cat("Selected candidate is exploratory (stage not reached: all step constraints freed).\n")
    }
  }

  if (!is.null(x$partial_searches) && length(x$partial_searches) > 0L) {
    cat("\nAutomatic partial searches by step\n")
    for (step in names(x$partial_searches)) {
      s <- x$partial_searches[[step]]
      cat(step, ":", s$decision %||% "unknown", "\n")
      if (!is.null(s$acceptable_models)) {
        cat("  acceptable:", nrow(s$acceptable_models), "\n")
      }
      if (!is.null(s$selected_partial) && length(s$selected_partial) > 0L) {
        cat("  selected constraints:", paste(s$selected_partial, collapse = ", "), "\n")
      }
      if (!is.null(s$selected_stage_reached) && identical(s$selected_stage_reached, FALSE)) {
        cat("  selected: exploratory (stage not reached)\n")
      }
    }
  }

  if (isTRUE(show_freed_parameters) && !is.null(x$freed_parameters) && length(x$freed_parameters) > 0L) {
    cat("\nFreed parameters by step\n")
    for (step in names(x$freed_parameters)) {
      vals <- x$freed_parameters[[step]]
      if (length(vals) > 0L) {
        cat(step, ":", paste(vals, collapse = ", "), "\n")
      }
    }
  }

  invisible(x)
}

#' Tidy MGCFA Fit Results
#'
#' Converts \code{mgcfa_result} objects into a long-format data.frame suitable
#' for reporting and plotting.
#'
#' @param x An object returned by \code{mgcfa_auto()}.
#' @param measures Fit measures to extract.
#' @param include_non_partial Logical; if \code{TRUE}, include failed
#'   non-partial models for failed steps when available.
#' @param baseline_step Step used for delta calculations.
#' @param add_delta Logical; if \code{TRUE}, include changes from
#'   \code{baseline_step}.
#' @param digits Number of digits used for numeric formatting in the returned
#'   table.
#' @param rounding Numeric rounding mode. One of \code{"signif"},
#'   \code{"round"}, or \code{"none"}.
#'
#' @return A data.frame in long format with one row per
#'   step/measure/model-variant.
#' @export
mgcfa_tidy_fit <- function(
  x,
  measures = c("chisq", "df", "cfi", "rmsea", "srmr", "aic", "bic"),
  include_non_partial = TRUE,
  baseline_step = "configural",
  add_delta = TRUE,
  digits = 3L,
  rounding = c("signif", "round", "none")
) {
  x <- .mgcfa_as_result(x)
  rounding <- match.arg(rounding)
  measures <- unique(as.character(measures))
  if (length(measures) == 0L) {
    stop("`measures` must include at least one fit measure name.", call. = FALSE)
  }

  step_order <- names(x$fits)
  if (!(baseline_step %in% step_order)) {
    stop("`baseline_step` must be one of the fitted step names.", call. = FALSE)
  }

  out <- list()
  out_i <- 0L

  for (step_i in seq_along(step_order)) {
    step <- step_order[[step_i]]
    fm <- suppressWarnings(lavaan::fitMeasures(x$fits[[step]], measures))
    freed_terms <- .mgcfa_normalize_terms(x$freed_parameters[[step]] %||% character())
    freed_label <- if (length(freed_terms) > 0L) paste(freed_terms, collapse = ", ") else ""
    for (m in measures) {
      out_i <- out_i + 1L
      out[[out_i]] <- data.frame(
        step = step,
        step_index = step_i,
        measure = m,
        value = if (m %in% names(fm)) as.numeric(fm[[m]]) else NA_real_,
        variant = "selected",
        series = "Selected",
        line_type = "solid",
        model_style = if (length(freed_terms) > 0L) "partial" else "full",
        n_freed = length(freed_terms),
        freed_terms = freed_label,
        stringsAsFactors = FALSE
      )
    }
  }

  if (isTRUE(include_non_partial) && !is.null(x$step_failures) && length(x$step_failures) > 0L) {
    for (step in names(x$step_failures)) {
      sf <- x$step_failures[[step]]
      if (!isTRUE(sf$failed) || is.null(sf$failed_fit_measures)) {
        next
      }
      step_i <- match(step, step_order)
      if (is.na(step_i)) {
        next
      }
      fm <- sf$failed_fit_measures
      for (m in measures) {
        out_i <- out_i + 1L
        out[[out_i]] <- data.frame(
          step = step,
          step_index = step_i,
          measure = m,
          value = if (m %in% names(fm)) as.numeric(fm[[m]]) else NA_real_,
          variant = "non_partial",
          series = "Non-Partial (Failed)",
          line_type = "dashed",
          model_style = "failed_non_partial",
          n_freed = 0L,
          freed_terms = "",
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(out) == 0L) {
    return(data.frame())
  }

  df <- do.call(rbind, out)
  df$step <- factor(df$step, levels = step_order, ordered = TRUE)
  df$variant <- factor(df$variant, levels = c("selected", "non_partial"))
  df$series <- factor(df$series, levels = c("Selected", "Non-Partial (Failed)"))

  if (isTRUE(add_delta)) {
    base_idx <- df$step == baseline_step & df$variant == "selected"
    base_map <- stats::setNames(df$value[base_idx], df$measure[base_idx])
    df$baseline_step <- baseline_step
    df$baseline_value <- as.numeric(base_map[df$measure])
    df$delta <- df$value - df$baseline_value
  }

  ord <- order(df$step_index, df$measure, df$variant)
  rownames(df) <- NULL
  df <- df[ord, , drop = FALSE]

  num_cols <- vapply(df, is.numeric, logical(1L))
  num_cols[match(c("step_index", "n_freed"), names(df), nomatch = 0L)] <- FALSE
  if (any(num_cols)) {
    df[num_cols] <- lapply(df[num_cols], .mgcfa_round_values, digits = digits, rounding = rounding)
  }
  df
}

#' Compare Raw vs Summary MGCFA Results
#'
#' Compares fit indices from a raw-data MGCFA run and a summary-matrix MGCFA run
#' across shared steps.
#'
#' @param raw_result An object returned by \code{mgcfa_auto()} from raw data.
#' @param summary_result An object returned by \code{mgcfa_auto()} from summary
#'   matrices.
#' @param measures Fit measures to compare.
#' @param tolerance Absolute difference tolerance for equivalence flags.
#'
#' @return A data.frame with per-step fit differences and
#'   \code{within_tolerance}. An \code{overall_within_tolerance} attribute is
#'   attached.
#' @export
mgcfa_compare_results <- function(
  raw_result,
  summary_result,
  measures = c("chisq", "df", "cfi", "rmsea", "srmr", "aic", "bic"),
  tolerance = 1e-4
) {
  raw_result <- .mgcfa_as_result(raw_result)
  summary_result <- .mgcfa_as_result(summary_result)
  measures <- unique(as.character(measures))
  if (length(measures) == 0L) {
    stop("`measures` must include at least one fit measure.", call. = FALSE)
  }
  tolerance <- as.numeric(tolerance)[[1L]]
  if (!is.finite(tolerance) || tolerance < 0) {
    stop("`tolerance` must be a non-negative numeric scalar.", call. = FALSE)
  }

  steps <- intersect(names(raw_result$fits), names(summary_result$fits))
  if (length(steps) == 0L) {
    stop("No overlapping steps to compare between the two results.", call. = FALSE)
  }

  out <- list()
  out_i <- 0L
  for (step in steps) {
    fm_raw <- suppressWarnings(lavaan::fitMeasures(raw_result$fits[[step]], measures))
    fm_sum <- suppressWarnings(lavaan::fitMeasures(summary_result$fits[[step]], measures))
    for (m in measures) {
      out_i <- out_i + 1L
      raw_val <- as.numeric(fm_raw[[m]] %||% NA_real_)
      sum_val <- as.numeric(fm_sum[[m]] %||% NA_real_)
      abs_diff <- abs(sum_val - raw_val)
      out[[out_i]] <- data.frame(
        step = step,
        measure = m,
        raw_value = raw_val,
        summary_value = sum_val,
        abs_diff = abs_diff,
        within_tolerance = is.finite(abs_diff) && abs_diff <= tolerance,
        stringsAsFactors = FALSE
      )
    }
  }
  df <- do.call(rbind, out)
  attr(df, "overall_within_tolerance") <- all(df$within_tolerance, na.rm = TRUE)
  df
}

#' Plot MGCFA Fit Evolution
#'
#' Builds a \code{ggplot2} visualization of fit statistics across invariance
#' levels, including failed non-partial models (dashed) and selected models
#' (solid).
#'
#' @param x An object returned by \code{mgcfa_auto()}.
#' @param measures Fit measures to plot.
#' @param include_non_partial Logical; include failed non-partial lines.
#' @param baseline_step Step used for delta plotting.
#' @param plot_type One of \code{"delta"} or \code{"raw"}.
#' @param color_values Named colors for \code{"selected"} and
#'   \code{"non_partial"}. You may also provide \code{"selected_mean"}.
#' @param linetype_values Named linetypes for \code{"selected"} and
#'   \code{"non_partial"}. You may also provide \code{"selected_mean"}.
#' @param point_size Point size.
#' @param line_size Line width.
#' @param facet_scales Facet scale behavior for \code{ggplot2::facet_wrap()}.
#' @param facet_ncol Optional facet column count.
#' @param title Optional custom plot title. If \code{NULL}, a default title is
#'   generated from the selected fit measure.
#'
#' @return A \code{ggplot} object.
#' @export
mgcfa_plot_fit <- function(
  x,
  measures = c("aic"),
  include_non_partial = TRUE,
  baseline_step = "configural",
  plot_type = c("delta", "raw"),
  color_values = c(selected = "#1768AC", non_partial = "#B9473E"),
  linetype_values = c(selected = "solid", non_partial = "dashed"),
  point_size = 2.2,
  line_size = 0.9,
  facet_scales = "free_y",
  facet_ncol = NULL,
  title = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plotting.", call. = FALSE)
  }
  plot_type <- match.arg(plot_type)
  x <- .mgcfa_as_result(x)

  df <- mgcfa_tidy_fit(
    x = x,
    measures = measures,
    include_non_partial = include_non_partial,
    baseline_step = baseline_step,
    add_delta = TRUE,
    rounding = "none"
  )
  if (nrow(df) == 0L) {
    stop("No rows available to plot.", call. = FALSE)
  }

  n_latent <- .mgcfa_n_latent(x)
  step_levels <- levels(df$step)
  step_labels <- .mgcfa_step_label_map(step_levels, n_latent = n_latent, line_break = TRUE)
  baseline_label <- .mgcfa_step_label_map(baseline_step, n_latent = n_latent, line_break = FALSE)[[1L]]

  measure_levels <- unique(as.character(df$measure))
  measure_labels <- stats::setNames(vapply(measure_levels, .mgcfa_measure_label, character(1L)), measure_levels)
  df$measure_label <- measure_labels[as.character(df$measure)]

  df$variant_plot <- as.character(df$variant)
  df$series_plot <- ifelse(df$variant_plot == "selected", "Selected", "Non-Partial (Failed)")

  has_partial_lv_var <- length(.mgcfa_normalize_terms((x$freed_parameters %||% list())[["lv.variances"]] %||% character())) > 0L
  split_mean_color <- isTRUE(n_latent == 1L) &&
    has_partial_lv_var &&
    ("lv.variances" %in% as.character(df$step)) &&
    ("means" %in% as.character(df$step))

  if (split_mean_color) {
    means_idx <- match("means", step_levels)
    is_mean_selected <- df$variant_plot == "selected" & df$step_index >= means_idx
    df$variant_plot[is_mean_selected] <- "selected_mean"
    df$series_plot[df$variant_plot == "selected"] <- "Selected (Variance)"
    df$series_plot[df$variant_plot == "selected_mean"] <- "Selected (Mean)"

    bridge_from_step <- if (!is.na(means_idx) && means_idx > 1L) {
      as.character(step_levels[[means_idx - 1L]])
    } else {
      "lv.variances"
    }
    bridge_rows <- df[df$variant == "selected" & as.character(df$step) == bridge_from_step, , drop = FALSE]
    if (nrow(bridge_rows) > 0L) {
      bridge_rows$variant_plot <- "selected_mean"
      bridge_rows$series_plot <- "Selected (Mean)"
      df <- rbind(df, bridge_rows)
    }
  }

  df$value_plot <- if (identical(plot_type, "delta")) df$delta else df$value
  df <- df[order(df$step_index, df$measure, df$series_plot, df$variant_plot), , drop = FALSE]

  single_measure <- length(unique(as.character(df$measure))) == 1L
  primary_measure <- if (single_measure) unique(as.character(df$measure))[[1L]] else NULL
  primary_label <- if (single_measure) .mgcfa_measure_label(primary_measure) else "Fit Measure"
  y_lab <- if (identical(plot_type, "delta")) {
    if (single_measure) {
      sprintf("Change in %s From %s", primary_label, baseline_label)
    } else {
      sprintf("Change in Value From %s", baseline_label)
    }
  } else {
    if (single_measure) primary_label else "Fit Value"
  }
  title <- title %||% if (single_measure) {
    sprintf("%s Changes Across Measurement Invariance Stages", primary_label)
  } else {
    "Fit Measure Changes Across Measurement Invariance Stages"
  }

  palette_default <- c(
    selected = "#1768AC",
    selected_mean = "#2A9D8F",
    non_partial = "#B9473E"
  )
  linetype_default <- c(
    selected = "solid",
    selected_mean = "solid",
    non_partial = "dashed"
  )

  if (length(color_values) > 0L) {
    palette_default[names(color_values)] <- color_values
  }
  if (length(linetype_values) > 0L) {
    linetype_default[names(linetype_values)] <- linetype_values
  }

  present_variants <- unique(as.character(df$variant_plot))
  present_variants <- intersect(names(palette_default), present_variants)
  legend_labels <- c(
    selected = "Selected",
    selected_mean = "Selected (Mean)",
    non_partial = "Non-Partial (Failed)"
  )
  legend_labels <- legend_labels[present_variants]
  legend_colors <- palette_default[present_variants]
  legend_types <- linetype_default[present_variants]
  series_plot <- measure_label <- variant_plot <- NULL

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = step,
      y = value_plot,
      group = interaction(series_plot, measure_label, drop = TRUE),
      color = variant_plot,
      linetype = variant_plot
    )
  ) +
    ggplot2::geom_line(linewidth = line_size, alpha = 0.95, na.rm = TRUE) +
    ggplot2::geom_point(size = point_size, na.rm = TRUE) +
    ggplot2::scale_x_discrete(labels = step_labels) +
    ggplot2::scale_color_manual(
      values = legend_colors,
      breaks = present_variants,
      labels = legend_labels
    ) +
    ggplot2::scale_linetype_manual(
      values = legend_types,
      breaks = present_variants,
      labels = legend_labels
    ) +
    ggplot2::labs(
      x = "Measurement Invariance Stage",
      y = y_lab,
      color = "Model Line",
      linetype = "Model Line",
      title = title
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_text(angle = 18, hjust = 1)
    )

  if (identical(plot_type, "delta")) {
    p <- p + ggplot2::geom_hline(
      yintercept = 0,
      color = "#7A7A7A",
      linewidth = 0.35
    )
  }
  if (!single_measure) {
    p <- p + ggplot2::facet_wrap(~measure_label, scales = facet_scales, ncol = facet_ncol)
  }

  p
}
