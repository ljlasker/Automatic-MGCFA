#' Run MGCFA with a Simple Interface
#'
#' Convenience wrapper that accepts raw data or summary statistics, prepares the
#' inputs, and runs \code{mgcfa_auto()} with practical defaults.
#'
#' @param model Lavaan syntax for \code{model_type = "custom"}.
#' @param model_type One of \code{"custom"}, \code{"single_factor"},
#'   \code{"efa"}, or \code{"pca"}.
#' @param input Optional object returned by \code{mgcfa_prepare_input()}.
#' @param data Optional raw data.frame.
#' @param group Optional grouping column in \code{data}.
#' @param sample_cov Optional covariance matrix/list.
#' @param sample_cor Optional correlation matrix/list.
#' @param sample_mean Optional means (summary mode).
#' @param sample_nobs Optional sample sizes (summary mode).
#' @param sample_sd Optional SDs (required for correlations unless auto-filled).
#' @param group_labels Optional group labels (summary mode).
#' @param include_steps Invariance steps to run.
#' @param ordered Optional ordered/categorical indicator specification.
#' @param partial_auto_search Behavior for automatic partial search.
#' @param ... Additional arguments passed to \code{mgcfa_auto()}.
#'
#' @return An object of class \code{"mgcfa_result"}.
#' @export
mgcfa_run_simple <- function(
  model,
  model_type = c("custom", "single_factor", "efa", "pca"),
  input = NULL,
  data = NULL,
  group = NULL,
  sample_cov = NULL,
  sample_cor = NULL,
  sample_mean = NULL,
  sample_nobs = NULL,
  sample_sd = NULL,
  group_labels = NULL,
  include_steps = c("configural", "metric", "scalar", "strict", "lv.variances", "lv.covariances", "residual.covariances", "regressions", "means"),
  ordered = NULL,
  partial_auto_search = c("always", "prompt", "never"),
  ...
) {
  model_type <- match.arg(model_type)
  partial_auto_search <- match.arg(partial_auto_search)
  if (is.null(input)) {
    input <- mgcfa_prepare_input(
      data = data,
      group = group,
      sample_cov = sample_cov,
      sample_cor = sample_cor,
      sample_mean = sample_mean,
      sample_nobs = sample_nobs,
      sample_sd = sample_sd,
      group_labels = group_labels
    )
  }
  mgcfa_auto(
    model_type = model_type,
    model = model,
    input = input,
    include_steps = include_steps,
    ordered = ordered,
    partial_auto_search = partial_auto_search,
    ...
  )
}

#' Build a Structured MGCFA Report
#'
#' Creates a compact report object with core tables, diagnostics, and optional
#' plots and exports.
#'
#' @param x An object returned by \code{mgcfa_auto()}.
#' @param measures Fit measures for tidy output.
#' @param plot_measures Measures used in default delta/raw plots.
#' @param include_non_partial Logical; include failed non-partial lines in plots.
#' @param digits Number of digits for printed/exported numeric summaries.
#' @param rounding Numeric rounding mode (\code{"signif"}, \code{"round"}, or
#'   \code{"none"}).
#' @param include_plots Logical; include ggplot objects when \code{ggplot2} is
#'   available.
#' @param output_dir Optional directory path to export CSV tables and an RDS
#'   report object.
#'
#' @return A list of class \code{"mgcfa_report"}.
#' @export
mgcfa_report <- function(
  x,
  measures = c("chisq", "df", "cfi", "rmsea", "srmr", "aic", "bic"),
  plot_measures = c("aic", "bic"),
  include_non_partial = TRUE,
  digits = 3L,
  rounding = c("signif", "round", "none"),
  include_plots = TRUE,
  output_dir = NULL
) {
  x <- .mgcfa_as_result(x)
  rounding <- match.arg(rounding)
  overview <- .mgcfa_format_numeric_df(.mgcfa_step_overview(x), digits = digits, rounding = rounding)
  tidy <- mgcfa_tidy_fit(
    x,
    measures = measures,
    include_non_partial = include_non_partial,
    add_delta = TRUE,
    digits = digits,
    rounding = rounding
  )
  failures <- .mgcfa_failure_summary(x, digits = digits, rounding = rounding)
  decision_trace <- .mgcfa_format_numeric_df(x$decision_trace %||% data.frame(), digits = digits, rounding = rounding)
  practical <- .mgcfa_format_numeric_df(x$practical_change_table %||% data.frame(), digits = digits, rounding = rounding)

  plots <- list()
  if (isTRUE(include_plots) && requireNamespace("ggplot2", quietly = TRUE)) {
    plots$delta <- mgcfa_plot_fit(
      x,
      measures = plot_measures,
      include_non_partial = include_non_partial,
      plot_type = "delta"
    )
    plots$raw <- mgcfa_plot_fit(
      x,
      measures = plot_measures,
      include_non_partial = include_non_partial,
      plot_type = "raw"
    )
  }

  report <- list(
    overview = overview,
    fit_table = .mgcfa_format_numeric_df(as.data.frame(x$fit_table), digits = digits, rounding = rounding),
    lrt_table = .mgcfa_format_numeric_df(x$lrt_table %||% data.frame(), digits = digits, rounding = rounding),
    tidy_fit = tidy,
    decision_trace = decision_trace,
    practical_change_table = practical,
    failures = failures,
    input_summary = x$input_summary %||% NULL,
    input_notes = x$input_notes %||% character(),
    metadata = x$metadata %||% list(),
    plots = plots
  )
  class(report) <- "mgcfa_report"

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    utils::write.csv(report$overview, file.path(output_dir, "overview.csv"), row.names = FALSE)
    utils::write.csv(report$tidy_fit, file.path(output_dir, "tidy_fit.csv"), row.names = FALSE)
    if (nrow(report$decision_trace) > 0L) {
      utils::write.csv(report$decision_trace, file.path(output_dir, "decision_trace.csv"), row.names = FALSE)
    }
    if (nrow(report$practical_change_table) > 0L) {
      utils::write.csv(report$practical_change_table, file.path(output_dir, "practical_change_table.csv"), row.names = FALSE)
    }
    if (nrow(report$failures) > 0L) {
      utils::write.csv(report$failures, file.path(output_dir, "failures.csv"), row.names = FALSE)
    }
    saveRDS(report, file = file.path(output_dir, "mgcfa_report.rds"))
  }
  report
}
