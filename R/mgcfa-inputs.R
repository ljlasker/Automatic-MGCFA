#' Prepare MGCFA Inputs from Raw or Summary Data
#'
#' Builds a standardized input object that can be passed to
#' \code{mgcfa_auto(input = ...)}. This simplifies setup when users have raw
#' data, covariance matrices, or correlation matrices with Ns/means/SDs.
#'
#' @param data Optional raw data.frame.
#' @param group Optional grouping column name in \code{data}.
#' @param sample_cov Optional covariance matrix or list of covariance matrices.
#' @param sample_cor Optional correlation matrix or list of correlation
#'   matrices. Alias for \code{sample_cov} + \code{matrices_are_cor = TRUE}.
#' @param sample_mean Optional vector/list of group means (summary mode).
#' @param sample_nobs Optional sample size scalar/vector/list (summary mode).
#' @param sample_sd Optional SD vector/list when using correlations.
#' @param group_labels Optional explicit group labels (summary mode).
#' @param infer_cor_from_sample_cov Logical; if \code{TRUE}, auto-detects when
#'   \code{sample_cov} appears to contain correlation matrices.
#' @param cor_without_sd_action Action when correlation matrices are supplied
#'   without \code{sample_sd}: \code{"auto_unit"} (fill unit SDs) or
#'   \code{"error"}.
#' @param summary_check_action Validation action for summary matrices:
#'   \code{"warn"}, \code{"error"}, or \code{"none"}.
#'
#' @return An object of class \code{"mgcfa_input"} containing
#'   \code{mgcfa_args}, \code{mode}, and \code{summary}.
#' @export
mgcfa_prepare_input <- function(
  data = NULL,
  group = NULL,
  sample_cov = NULL,
  sample_cor = NULL,
  sample_mean = NULL,
  sample_nobs = NULL,
  sample_sd = NULL,
  group_labels = NULL,
  infer_cor_from_sample_cov = TRUE,
  cor_without_sd_action = c("auto_unit", "error"),
  summary_check_action = c("warn", "error", "none")
) {
  cor_without_sd_action <- match.arg(cor_without_sd_action)
  summary_check_action <- match.arg(summary_check_action)
  using_raw <- !is.null(data)
  using_summary <- !is.null(sample_cov) || !is.null(sample_cor)
  if (identical(using_raw, using_summary)) {
    stop(
      "Provide either raw inputs (`data` + `group`) or summary inputs (`sample_cov`/`sample_cor`).",
      call. = FALSE
    )
  }

  notes <- character()

  ordered_vars <- NULL
  if (using_raw) {
    if (!is.data.frame(data)) {
      stop("`data` must be a data.frame.", call. = FALSE)
    }
    if (!is.character(group) || length(group) != 1L || !nzchar(group) || !(group %in% names(data))) {
      stop("`group` must name a grouping column in `data`.", call. = FALSE)
    }
    g <- as.character(data[[group]])
    g <- g[!is.na(g)]
    grp <- unique(g)
    if (length(grp) < 2L) {
      stop("Raw-data mode requires at least two groups.", call. = FALSE)
    }
    grp_tab <- table(g)
    summary <- list(
      mode = "raw_data",
      n_rows = nrow(data),
      n_columns = ncol(data),
      group = group,
      n_groups = length(grp),
      group_labels = as.character(names(grp_tab)),
      nobs_by_group = as.integer(grp_tab)
    )
    out <- list(
      mode = "raw_data",
      mgcfa_args = list(data = data, group = group),
      summary = summary,
      notes = notes
    )
    class(out) <- "mgcfa_input"
    return(out)
  }

  if (!is.null(sample_cor) && !is.null(sample_cov)) {
    stop("Provide `sample_cov` or `sample_cor`, not both.", call. = FALSE)
  }
  matrices_are_cor <- !is.null(sample_cor)
  sample_cov_use <- sample_cov
  if (!is.null(sample_cor)) {
    sample_cov_use <- sample_cor
  }

  cov_list <- if (is.matrix(sample_cov_use)) list(sample_cov_use) else sample_cov_use
  if (!is.list(cov_list) || length(cov_list) < 1L) {
    stop("`sample_cov`/`sample_cor` must be a matrix or a non-empty list of matrices.", call. = FALSE)
  }
  cov_list <- lapply(cov_list, function(m) {
    if (!is.matrix(m) || nrow(m) != ncol(m)) {
      stop("Each supplied matrix must be square.", call. = FALSE)
    }
    storage.mode(m) <- "double"
    m
  })

  if (!isTRUE(matrices_are_cor) && isTRUE(infer_cor_from_sample_cov)) {
    if (all(vapply(cov_list, .mgcfa_is_correlation_matrix, logical(1L)))) {
      matrices_are_cor <- TRUE
      notes <- c(notes, "Detected correlation matrices from `sample_cov`.")
    }
  }

  if (isTRUE(matrices_are_cor) && is.null(sample_sd)) {
    if (identical(cor_without_sd_action, "error")) {
      stop(
        "Correlation matrices require `sample_sd` (or use `cor_without_sd_action = \"auto_unit\"`).",
        call. = FALSE
      )
    }
    sample_sd <- lapply(cov_list, function(m) {
      nm <- colnames(m) %||% paste0("V", seq_len(ncol(m)))
      stats::setNames(rep(1, ncol(m)), nm)
    })
    notes <- c(notes, "No SDs provided for correlations; unit SDs were inserted automatically.")
  }

  n_groups <- length(cov_list)
  gl <- group_labels
  if (is.null(gl)) {
    gl <- names(cov_list)
    if (is.null(gl) || any(!nzchar(gl))) {
      gl <- paste0("g", seq_len(n_groups))
    }
  }

  if (is.null(sample_nobs)) {
    stop("Summary mode requires `sample_nobs` (scalar, vector, or list).", call. = FALSE)
  }
  .mgcfa_validate_summary_matrices(
    sample_cov = cov_list,
    sample_nobs = sample_nobs,
    sample_mean = sample_mean,
    action = summary_check_action
  )

  summary <- list(
    mode = "summary_matrices",
    n_groups = n_groups,
    n_variables = ncol(cov_list[[1L]]),
    matrix_type = if (isTRUE(matrices_are_cor)) "correlation" else "covariance",
    has_means = !is.null(sample_mean),
    has_sds = !is.null(sample_sd),
    group_labels = as.character(gl)
  )
  out <- list(
    mode = "summary_matrices",
    mgcfa_args = list(
      sample_cov = sample_cov_use,
      sample_mean = sample_mean,
      sample_nobs = sample_nobs,
      group_labels = gl,
      matrices_are_cor = matrices_are_cor,
      sample_sd = sample_sd
    ),
    summary = summary,
    notes = notes
  )
  class(out) <- "mgcfa_input"
  out
}

#' Print Method for Prepared MGCFA Inputs
#'
#' @param x An object returned by \code{mgcfa_prepare_input()}.
#' @param ... Unused.
#'
#' @return The input object invisibly.
#' @export
print.mgcfa_input <- function(x, ...) {
  if (!inherits(x, "mgcfa_input")) {
    stop("`x` must be an object returned by `mgcfa_prepare_input()`.", call. = FALSE)
  }
  cat("Prepared input mode:", x$mode, "\n")
  s <- x$summary %||% list()
  if (identical(x$mode, "raw_data")) {
    cat("Rows:", s$n_rows %||% NA_integer_, " Groups:", s$n_groups %||% NA_integer_, "\n")
    cat("Group variable:", s$group %||% "", "\n")
  } else {
    cat("Groups:", s$n_groups %||% NA_integer_, " Variables:", s$n_variables %||% NA_integer_, "\n")
    cat("Matrix type:", s$matrix_type %||% "", "\n")
    cat("Means provided:", isTRUE(s$has_means), " SDs provided:", isTRUE(s$has_sds), "\n")
  }
  if (!is.null(x$notes) && length(x$notes) > 0L) {
    cat("Notes:\n")
    for (msg in x$notes) {
      cat("- ", msg, "\n", sep = "")
    }
  }
  invisible(x)
}

#' Quick Input Guidance for AutomaticMGCFA
#'
#' Prints concise recipes for common input types.
#'
#' @return Invisibly returns a data.frame of recipes.
#' @export
mgcfa_help_inputs <- function() {
  txt <- data.frame(
    scenario = c(
      "Raw data",
      "Covariance matrices",
      "Correlation matrices + SDs",
      "Correlation matrices without SDs"
    ),
    minimal_requirements = c(
      "data, group",
      "sample_cov, sample_nobs",
      "sample_cor (or sample_cov), sample_sd, sample_nobs",
      "mgcfa_prepare_input(..., sample_cor=..., cor_without_sd_action='auto_unit')"
    ),
    stringsAsFactors = FALSE
  )
  cat("AutomaticMGCFA input recipes\n\n")
  print(txt, row.names = FALSE)
  cat("\nTip: use `mgcfa_prepare_input()` and pass it into `mgcfa_auto(input = prepared)`.\n")
  invisible(txt)
}
