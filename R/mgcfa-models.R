#' Build SEM Syntax for MGCFA
#'
#' Creates lavaan model syntax from a user-supplied model, a single-factor
#' template, or thresholded exploratory loadings from EFA/PCA.
#'
#' @param model_type One of \code{"custom"}, \code{"single_factor"},
#'   \code{"efa"}, or \code{"pca"}.
#' @param model Lavaan syntax used when \code{model_type = "custom"}.
#' @param variables Character vector of indicator variable names.
#' @param data Optional data.frame used for EFA/PCA model construction.
#' @param input_matrix Optional covariance or correlation matrix used for
#'   EFA/PCA model construction.
#' @param matrix_nobs Optional sample size for matrix-based EFA/PCA.
#' @param n_factors Number of factors/components for EFA/PCA model generation.
#' @param loading_threshold Absolute loading cut point used to retain paths.
#' @param rotation Rotation passed to \code{psych::fa()} or
#'   \code{psych::principal()}.
#' @param allow_cross_loadings Logical. If \code{TRUE}, keep all loadings above
#'   threshold. If \code{FALSE}, keep only the strongest loading per variable.
#'
#' @return A character string with lavaan model syntax.
#' @export
mgcfa_build_model <- function(
  model_type = c("custom", "single_factor", "efa", "pca"),
  model = NULL,
  variables = NULL,
  data = NULL,
  input_matrix = NULL,
  matrix_nobs = NULL,
  n_factors = 1L,
  loading_threshold = 0.30,
  rotation = "oblimin",
  allow_cross_loadings = FALSE
) {
  model_type <- match.arg(model_type)

  if (model_type == "custom") {
    if (!is.character(model) || length(model) != 1L || !nzchar(model)) {
      stop("`model` must be a non-empty lavaan syntax string when `model_type = \"custom\"`.", call. = FALSE)
    }
    return(model)
  }

  if (is.null(variables)) {
    if (!is.null(data)) {
      variables <- colnames(data)
    } else if (!is.null(input_matrix)) {
      variables <- colnames(input_matrix)
    }
  }

  if (is.null(variables) || length(variables) < 2L) {
    stop("`variables` must include at least two observed indicators.", call. = FALSE)
  }

  variables <- as.character(variables)

  if (model_type == "single_factor") {
    return(paste0("F1 =~ ", paste(variables, collapse = " + ")))
  }

  if (!requireNamespace("psych", quietly = TRUE)) {
    stop("Package `psych` is required for EFA/PCA model generation.", call. = FALSE)
  }

  n_factors <- as.integer(n_factors)
  if (is.na(n_factors) || n_factors < 1L) {
    stop("`n_factors` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(loading_threshold) || length(loading_threshold) != 1L || loading_threshold < 0) {
    stop("`loading_threshold` must be a non-negative numeric scalar.", call. = FALSE)
  }

  loadings <- .mgcfa_exploratory_loadings(
    model_type = model_type,
    variables = variables,
    data = data,
    input_matrix = input_matrix,
    matrix_nobs = matrix_nobs,
    n_factors = n_factors,
    rotation = rotation
  )

  .mgcfa_syntax_from_loadings(
    loadings = loadings,
    loading_threshold = loading_threshold,
    allow_cross_loadings = isTRUE(allow_cross_loadings)
  )
}

#' Build Summary Inputs from Raw Data
#'
#' Creates group-wise summary inputs compatible with \code{mgcfa_auto()} from
#' raw data. Optionally applies paper-style decimal formatting (truncate or
#' round) to correlations/covariances, SDs, and means.
#'
#' @param data Data frame containing indicators and a grouping variable.
#' @param group Grouping variable name in \code{data}.
#' @param variables Character vector of indicator names.
#' @param summary_profile Summary preset. \code{"raw_equivalent"} (default)
#'   targets reproducibility versus raw-data MGCFA by using covariance matrices
#'   and no decimal formatting. \code{"paper"} uses compact table-style summary
#'   output (correlation matrices + SDs by default, truncated decimals unless
#'   overridden). \code{"custom"} uses the supplied \code{matrix_type},
#'   \code{format}, and digit options directly.
#' @param matrix_type One of \code{"cor"} or \code{"cov"}. If \code{NULL},
#'   defaults are chosen from \code{summary_profile}.
#' @param cor_digits Decimal places for correlation matrices when
#'   \code{matrix_type = "cor"}. Ignored when \code{NULL}.
#' @param sd_digits Decimal places for group SD vectors when
#'   \code{matrix_type = "cor"}. Ignored when \code{NULL}.
#' @param mean_digits Decimal places for group means. Ignored when \code{NULL}.
#' @param cov_digits Decimal places for covariance matrices when
#'   \code{matrix_type = "cov"}. Ignored when \code{NULL}.
#' @param format Decimal formatting method: \code{"none"},
#'   \code{"truncate"}, or \code{"round"}. If \code{NULL}, defaults are chosen
#'   from \code{summary_profile}.
#' @param drop_incomplete Logical; if \code{TRUE}, uses complete cases for
#'   \code{variables} within each group.
#' 
#' @return A list containing \code{sample_cov}, \code{sample_mean},
#'   \code{sample_nobs}, and \code{group_labels}, plus \code{matrices_are_cor}
#'   and \code{sample_sd} (NULL in covariance mode). Also returns
#'   \code{summary_profile} and \code{mgcfa_args}, a ready-to-pass argument list
#'   for \code{mgcfa_auto()} matrix mode.
#' @export
mgcfa_make_summary <- function(
  data,
  group,
  variables,
  summary_profile = c("raw_equivalent", "paper", "custom"),
  matrix_type = NULL,
  cor_digits = NULL,
  sd_digits = NULL,
  mean_digits = NULL,
  cov_digits = NULL,
  format = NULL,
  drop_incomplete = TRUE
) {
  summary_profile <- match.arg(summary_profile)
  if (is.null(matrix_type)) {
    matrix_type <- switch(
      summary_profile,
      raw_equivalent = "cov",
      paper = "cor",
      custom = "cor"
    )
  }
  matrix_type <- match.arg(matrix_type, choices = c("cor", "cov"))

  if (is.null(format)) {
    format <- switch(
      summary_profile,
      raw_equivalent = "none",
      paper = "truncate",
      custom = "none"
    )
  }
  format <- match.arg(format, choices = c("none", "truncate", "round"))

  if (identical(summary_profile, "paper")) {
    if (identical(matrix_type, "cor")) {
      cor_digits <- cor_digits %||% 2L
      sd_digits <- sd_digits %||% 2L
    } else {
      cov_digits <- cov_digits %||% 2L
    }
    mean_digits <- mean_digits %||% 2L
  }

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  if (!is.character(group) || length(group) != 1L || !nzchar(group)) {
    stop("`group` must be a single non-empty column name.", call. = FALSE)
  }
  if (!(group %in% names(data))) {
    stop("`group` must exist in `data`.", call. = FALSE)
  }
  if (is.null(variables) || length(variables) < 2L) {
    stop("`variables` must contain at least two indicator names.", call. = FALSE)
  }
  variables <- as.character(variables)
  if (!all(variables %in% names(data))) {
    stop("All `variables` must exist in `data`.", call. = FALSE)
  }
  if (!all(vapply(data[variables], is.numeric, logical(1L)))) {
    stop("All selected `variables` must be numeric.", call. = FALSE)
  }
  if (!is.logical(drop_incomplete) || length(drop_incomplete) != 1L || is.na(drop_incomplete)) {
    stop("`drop_incomplete` must be TRUE or FALSE.", call. = FALSE)
  }

  g_raw <- data[[group]]
  if (anyNA(g_raw)) {
    stop("`group` contains missing values; please remove or recode them first.", call. = FALSE)
  }
  g <- as.character(g_raw)
  groups <- if (is.factor(g_raw)) {
    lev <- levels(g_raw)
    lev[lev %in% g]
  } else {
    unique(g)
  }
  if (length(groups) < 2L) {
    stop("At least two groups are required to build multi-group summary inputs.", call. = FALSE)
  }

  sample_cov <- vector("list", length(groups))
  sample_sd <- if (identical(matrix_type, "cor")) vector("list", length(groups)) else NULL
  sample_mean <- vector("list", length(groups))
  sample_nobs <- vector("list", length(groups))

  for (i in seq_along(groups)) {
    gi <- groups[[i]]
    dgi <- data[g == gi, variables, drop = FALSE]

    if (isTRUE(drop_incomplete)) {
      keep <- stats::complete.cases(dgi)
      dgi <- dgi[keep, , drop = FALSE]
    }

    n_i <- nrow(dgi)
    if (n_i < 2L) {
      stop(sprintf("Group `%s` has fewer than 2 complete rows for `variables`.", gi), call. = FALSE)
    }

    if (identical(matrix_type, "cor")) {
      Ri <- stats::cor(dgi)
      Ri <- .mgcfa_decimal_matrix(
        Ri,
        digits = cor_digits,
        format = format,
        is_correlation = TRUE
      )

      sdi <- apply(dgi, 2, stats::sd)
      sdi <- .mgcfa_decimal_values(sdi, digits = sd_digits, format = format)
      names(sdi) <- variables

      sample_cov[[i]] <- Ri
      sample_sd[[i]] <- sdi
    } else {
      Si <- stats::cov(dgi)
      Si <- .mgcfa_decimal_matrix(
        Si,
        digits = cov_digits,
        format = format,
        is_correlation = FALSE
      )
      sample_cov[[i]] <- Si
    }

    mi <- colMeans(dgi)
    mi <- .mgcfa_decimal_values(mi, digits = mean_digits, format = format)
    names(mi) <- variables

    sample_mean[[i]] <- mi
    sample_nobs[[i]] <- n_i
  }

  out <- list(
    sample_cov = sample_cov,
    sample_sd = sample_sd,
    sample_mean = sample_mean,
    sample_nobs = sample_nobs,
    group_labels = groups,
    matrices_are_cor = identical(matrix_type, "cor"),
    summary_profile = summary_profile
  )
  out$mgcfa_args <- list(
    sample_cov = out$sample_cov,
    sample_mean = out$sample_mean,
    sample_nobs = out$sample_nobs,
    group_labels = out$group_labels,
    matrices_are_cor = out$matrices_are_cor,
    sample_sd = out$sample_sd
  )
  out
}
