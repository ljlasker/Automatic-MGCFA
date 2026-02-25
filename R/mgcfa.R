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
#' @param matrix_type One of \code{"cor"} or \code{"cov"}.
#' @param cor_digits Decimal places for correlation matrices when
#'   \code{matrix_type = "cor"}. Ignored when \code{NULL}.
#' @param sd_digits Decimal places for group SD vectors when
#'   \code{matrix_type = "cor"}. Ignored when \code{NULL}.
#' @param mean_digits Decimal places for group means. Ignored when \code{NULL}.
#' @param cov_digits Decimal places for covariance matrices when
#'   \code{matrix_type = "cov"}. Ignored when \code{NULL}.
#' @param format Decimal formatting method: \code{"none"},
#'   \code{"truncate"}, or \code{"round"}.
#' @param drop_incomplete Logical; if \code{TRUE}, uses complete cases for
#'   \code{variables} within each group.
#'
#' @return A list containing \code{sample_cov}, \code{sample_mean},
#'   \code{sample_nobs}, and \code{group_labels}, plus \code{matrices_are_cor}
#'   and \code{sample_sd} (NULL in covariance mode).
#' @export
mgcfa_make_summary <- function(
  data,
  group,
  variables,
  matrix_type = c("cor", "cov"),
  cor_digits = NULL,
  sd_digits = NULL,
  mean_digits = NULL,
  cov_digits = NULL,
  format = c("none", "truncate", "round"),
  drop_incomplete = TRUE
) {
  matrix_type <- match.arg(matrix_type)
  format <- match.arg(format)

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

  list(
    sample_cov = sample_cov,
    sample_sd = sample_sd,
    sample_mean = sample_mean,
    sample_nobs = sample_nobs,
    group_labels = groups,
    matrices_are_cor = identical(matrix_type, "cor")
  )
}

#' Run Automatic MGCFA from Raw Data or Matrices
#'
#' Fits a sequential multi-group SEM/CFA invariance workflow with a model
#' supplied directly or generated from single-factor, EFA, or PCA templates.
#'
#' @param model_type One of \code{"custom"}, \code{"single_factor"},
#'   \code{"efa"}, or \code{"pca"}.
#' @param model Lavaan syntax when \code{model_type = "custom"}.
#' @param variables Variables used for model generation in non-custom modes.
#' @param n_factors Number of factors/components for EFA/PCA generation.
#' @param loading_threshold Absolute loading threshold for EFA/PCA/PCA-derived
#'   syntax.
#' @param rotation Rotation for exploratory model generation.
#' @param allow_cross_loadings Logical for EFA/PCA syntax generation.
#' @param data Raw-data input (data.frame) for multi-group SEM.
#' @param group Grouping variable name in \code{data}.
#' @param sample_cov Matrix or list of group covariance/correlation matrices.
#' @param sample_mean Optional vector/list of group means (required for scalar
#'   and stricter steps in matrix mode).
#' @param sample_nobs Scalar, vector, or list of group sample sizes in matrix
#'   mode.
#' @param group_labels Optional explicit group labels.
#' @param matrices_are_cor Logical. If \code{TRUE}, \code{sample_cov} is treated
#'   as correlation matrices and converted using \code{sample_sd}.
#' @param sample_sd Group SD vectors used when \code{matrices_are_cor = TRUE}.
#' @param include_steps Character vector of invariance steps to fit.
#' @param partial Optional named list where each element is a
#'   \code{group.partial} character vector for a step.
#' @param estimator Optional lavaan estimator.
#' @param std.lv_configural \code{std.lv} setting for the configural step.
#' @param std.lv_constrained \code{std.lv} setting for constrained steps.
#' @param orthogonal Passed to lavaan.
#' @param fit_measures Fit statistics extracted from each fitted step.
#' @param partial_failure_criterion Criterion used to decide whether a constrained
#'   invariance step failed relative to the previous step in the fitted
#'   sequence. One of
#'   \code{"chisq_pvalue"}, \code{"delta_cfi"}, \code{"aic_bic_weight"},
#'   \code{"measure_change"}, or \code{"none"}.
#' @param partial_failure_threshold Threshold used with
#'   \code{partial_failure_criterion}. Defaults are \code{0.05} for
#'   \code{"chisq_pvalue"}, \code{0.01} for \code{"delta_cfi"},
#'   \code{0.5} for \code{"aic_bic_weight"}, and \code{0} for
#'   \code{"measure_change"}.
#' @param partial_failure_measure Fit measure used when
#'   \code{partial_failure_criterion = "measure_change"}.
#' @param partial_failure_direction Direction used when
#'   \code{partial_failure_criterion = "measure_change"}.
#' @param partial_auto_search Behavior when a constrained step fails:
#'   \code{"prompt"}, \code{"never"}, or \code{"always"}.
#' @param partial_search_criterion Criterion used during automatic partial
#'   invariance search. One of \code{"chisq_pvalue"}, \code{"delta_cfi"},
#'   \code{"aic_bic_weight"}, or \code{"measure_change"}.
#' @param partial_search_threshold Threshold used for
#'   \code{partial_search_criterion}.
#' @param partial_search_measure Fit measure used when
#'   \code{partial_search_criterion = "measure_change"}. Defaults to
#'   \code{partial_failure_measure}.
#' @param partial_search_direction Direction used when
#'   \code{partial_search_criterion = "measure_change"}. Defaults to
#'   \code{partial_failure_direction}.
#' @param partial_ic_bic_weight Weight on BIC when
#'   \code{partial_failure_criterion} or \code{partial_search_criterion} is
#'   \code{"aic_bic_weight"}. The AIC weight is
#'   \code{1 - partial_ic_bic_weight}.
#' @param partial_search_max_free Maximum number of additional constraints to
#'   free during automatic partial search at the failed step. The algorithm
#'   always retains at least one active equality constraint for that step.
#' @param partial_search_top_n Number of closest candidate models returned from
#'   partial search.
#' @param partial_search_stop_on_accept Logical; stop search at first acceptable
#'   candidate when \code{TRUE}.
#' @param partial_search_rank Ranking for candidate output:
#'   \code{"closest"} or \code{"best"}.
#' @param partial_search_use_best_if_no_pass Logical; if no candidate meets the
#'   threshold, continue invariance testing using the best ranked candidate.
#' @param partial_search_candidate_source Candidate release-term source.
#'   \code{"score"} uses only score-test releasable constraints,
#'   \code{"all"} uses all releasable equality terms detected for the step,
#'   and \code{"auto"} uses \code{"all"} for
#'   \code{partial_search_exhaustive_steps} and \code{"score"} otherwise.
#' @param partial_search_exhaustive_steps Steps that should default to exhaustive
#'   release-term enumeration when \code{partial_search_candidate_source = "auto"}.
#'   Defaults to \code{c("lv.variances", "means")}.
#' @param partial_search_max_models Maximum number of candidate models evaluated
#'   during automatic partial search for a step.
#' @param stop_at_first_unacceptable Logical; if \code{TRUE}, stop fitting
#'   higher invariance stages after the first constrained stage that remains
#'   unacceptable relative to the previous fitted stage.
#'
#' @return An object of class \code{"mgcfa_result"}.
#' @export
mgcfa_auto <- function(
  model_type = c("custom", "single_factor", "efa", "pca"),
  model = NULL,
  variables = NULL,
  n_factors = 1L,
  loading_threshold = 0.30,
  rotation = "oblimin",
  allow_cross_loadings = FALSE,
  data = NULL,
  group = NULL,
  sample_cov = NULL,
  sample_mean = NULL,
  sample_nobs = NULL,
  group_labels = NULL,
  matrices_are_cor = FALSE,
  sample_sd = NULL,
  include_steps = c("configural", "metric", "scalar", "strict", "lv.variances", "means"),
  partial = NULL,
  estimator = NULL,
  std.lv_configural = TRUE,
  std.lv_constrained = FALSE,
  orthogonal = FALSE,
  fit_measures = c("chisq", "df", "npar", "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic", "bic"),
  partial_failure_criterion = c("chisq_pvalue", "delta_cfi", "aic_bic_weight", "measure_change", "none"),
  partial_failure_threshold = NULL,
  partial_failure_measure = "aic",
  partial_failure_direction = c("decrease", "increase"),
  partial_auto_search = c("prompt", "never", "always"),
  partial_search_criterion = NULL,
  partial_search_threshold = NULL,
  partial_search_measure = NULL,
  partial_search_direction = NULL,
  partial_ic_bic_weight = 0.5,
  partial_search_max_free = NULL,
  partial_search_top_n = 5L,
  partial_search_stop_on_accept = TRUE,
  partial_search_rank = c("closest", "best"),
  partial_search_use_best_if_no_pass = TRUE,
  partial_search_candidate_source = c("auto", "score", "all"),
  partial_search_exhaustive_steps = c("lv.variances", "means"),
  partial_search_max_models = 5000L,
  stop_at_first_unacceptable = TRUE
) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package `lavaan` is required for MGCFA.", call. = FALSE)
  }

  model_type <- match.arg(model_type)
  available_steps <- c("configural", "metric", "scalar", "strict", "lv.variances", "means")
  include_steps <- unique(match.arg(include_steps, choices = available_steps, several.ok = TRUE))
  partial_failure_criterion <- match.arg(
    partial_failure_criterion,
    choices = c("chisq_pvalue", "delta_cfi", "aic_bic_weight", "measure_change", "none")
  )
  partial_failure_direction <- match.arg(partial_failure_direction)
  partial_auto_search <- match.arg(partial_auto_search, choices = c("prompt", "never", "always"))
  partial_search_rank <- match.arg(partial_search_rank)
  partial_search_candidate_source <- match.arg(partial_search_candidate_source, choices = c("auto", "score", "all"))

  if (!is.null(partial_search_criterion)) {
    partial_search_criterion <- match.arg(
      partial_search_criterion,
      choices = c("chisq_pvalue", "delta_cfi", "aic_bic_weight", "measure_change")
    )
  }
  if (!is.null(partial_search_direction)) {
    partial_search_direction <- match.arg(partial_search_direction, choices = c("decrease", "increase"))
  }

  partial_failure_measure <- as.character(partial_failure_measure)[[1L]]
  if (!nzchar(partial_failure_measure)) {
    stop("`partial_failure_measure` must be a non-empty fit measure name.", call. = FALSE)
  }
  if (!is.null(partial_search_measure)) {
    partial_search_measure <- as.character(partial_search_measure)[[1L]]
    if (!nzchar(partial_search_measure)) {
      stop("`partial_search_measure` must be NULL or a non-empty fit measure name.", call. = FALSE)
    }
  }
  if (!is.numeric(partial_ic_bic_weight) || length(partial_ic_bic_weight) != 1L || !is.finite(partial_ic_bic_weight) ||
      partial_ic_bic_weight < 0 || partial_ic_bic_weight > 1) {
    stop("`partial_ic_bic_weight` must be a numeric scalar in [0, 1].", call. = FALSE)
  }

  if (partial_failure_criterion != "none" && is.null(partial_failure_threshold)) {
    partial_failure_threshold <- .mgcfa_default_partial_threshold(partial_failure_criterion)
  }
  if (!is.null(partial_failure_threshold)) {
    if (!is.numeric(partial_failure_threshold) || length(partial_failure_threshold) != 1L || !is.finite(partial_failure_threshold)) {
      stop("`partial_failure_threshold` must be a finite numeric scalar.", call. = FALSE)
    }
  }

  if (!is.null(partial_search_threshold) &&
      (!is.numeric(partial_search_threshold) || length(partial_search_threshold) != 1L || !is.finite(partial_search_threshold))) {
    stop("`partial_search_threshold` must be a finite numeric scalar.", call. = FALSE)
  }

  if (!is.null(partial_search_max_free)) {
    partial_search_max_free <- as.integer(partial_search_max_free)
    if (is.na(partial_search_max_free) || partial_search_max_free < 0L) {
      stop("`partial_search_max_free` must be NULL or a non-negative integer.", call. = FALSE)
    }
  }

  partial_search_top_n <- as.integer(partial_search_top_n)
  if (is.na(partial_search_top_n) || partial_search_top_n < 1L) {
    stop("`partial_search_top_n` must be a positive integer.", call. = FALSE)
  }
  if (!is.logical(partial_search_stop_on_accept) || length(partial_search_stop_on_accept) != 1L ||
      is.na(partial_search_stop_on_accept)) {
    stop("`partial_search_stop_on_accept` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(partial_search_use_best_if_no_pass) || length(partial_search_use_best_if_no_pass) != 1L ||
      is.na(partial_search_use_best_if_no_pass)) {
    stop("`partial_search_use_best_if_no_pass` must be TRUE or FALSE.", call. = FALSE)
  }
  partial_search_max_models <- as.integer(partial_search_max_models)
  if (is.na(partial_search_max_models) || partial_search_max_models < 1L) {
    stop("`partial_search_max_models` must be a positive integer.", call. = FALSE)
  }
  partial_search_exhaustive_steps <- unique(match.arg(
    partial_search_exhaustive_steps,
    choices = available_steps,
    several.ok = TRUE
  ))

  if (!is.logical(stop_at_first_unacceptable) || length(stop_at_first_unacceptable) != 1L || is.na(stop_at_first_unacceptable)) {
    stop("`stop_at_first_unacceptable` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(partial) && !is.list(partial)) {
    stop("`partial` must be NULL or a named list of step-specific `group.partial` vectors.", call. = FALSE)
  }

  using_raw <- !is.null(data)
  using_matrices <- !is.null(sample_cov)
  if (identical(using_raw, using_matrices)) {
    stop("Provide either raw `data` (with `group`) or matrix input via `sample_cov`.", call. = FALSE)
  }

  if (using_raw) {
    if (!is.data.frame(data)) {
      stop("`data` must be a data.frame.", call. = FALSE)
    }
    if (!is.character(group) || length(group) != 1L || !nzchar(group) || !(group %in% names(data))) {
      stop("`group` must name a grouping column in `data`.", call. = FALSE)
    }

    if (is.null(variables) && model_type != "custom") {
      is_num <- vapply(data, is.numeric, logical(1L))
      variables <- setdiff(names(data)[is_num], group)
    }

    model_syntax <- mgcfa_build_model(
      model_type = model_type,
      model = model,
      variables = variables,
      data = if (model_type %in% c("efa", "pca")) data[, variables, drop = FALSE] else NULL,
      n_factors = n_factors,
      loading_threshold = loading_threshold,
      rotation = rotation,
      allow_cross_loadings = allow_cross_loadings
    )

    fit_args_base <- list(
      data = data,
      group = group
    )
    if (!is.null(group_labels)) {
      fit_args_base$group.label <- as.character(group_labels)
    }
  } else {
    prepared <- .mgcfa_prepare_matrix_input(
      sample_cov = sample_cov,
      sample_mean = sample_mean,
      sample_nobs = sample_nobs,
      group_labels = group_labels,
      matrices_are_cor = matrices_are_cor,
      sample_sd = sample_sd
    )

    if (is.null(variables) && model_type != "custom") {
      variables <- colnames(prepared$sample_cov[[1L]])
    }

    model_syntax <- mgcfa_build_model(
      model_type = model_type,
      model = model,
      variables = variables,
      input_matrix = if (model_type %in% c("efa", "pca")) prepared$sample_cov[[1L]] else NULL,
      matrix_nobs = if (model_type %in% c("efa", "pca")) prepared$sample_nobs[[1L]] else NULL,
      n_factors = n_factors,
      loading_threshold = loading_threshold,
      rotation = rotation,
      allow_cross_loadings = allow_cross_loadings
    )

    if (is.null(prepared$sample_mean)) {
      no_mean_steps <- c("scalar", "strict", "lv.variances", "means")
      dropped <- intersect(include_steps, no_mean_steps)
      if (length(dropped) > 0L) {
        warning(
          "Dropping steps that require `sample_mean` in matrix mode: ",
          paste(dropped, collapse = ", ")
        )
        include_steps <- setdiff(include_steps, dropped)
      }
    }

    if (length(include_steps) == 0L) {
      stop("No estimable invariance steps remain after input validation.", call. = FALSE)
    }

    fit_args_base <- list(
      sample.cov = prepared$sample_cov,
      sample.nobs = prepared$sample_nobs
    )
    if (!is.null(prepared$sample_mean)) {
      fit_args_base$sample.mean <- prepared$sample_mean
    }
    if (!is.null(prepared$group_labels)) {
      fit_args_base$group.label <- prepared$group_labels
    }
  }

  step_equal <- list(
    configural = NULL,
    metric = c("loadings"),
    scalar = c("loadings", "intercepts"),
    strict = c("loadings", "intercepts", "residuals"),
    "lv.variances" = c("loadings", "intercepts", "residuals", "lv.variances"),
    means = c("loadings", "intercepts", "residuals", "lv.variances", "means")
  )

  fits <- list()
  step_failures <- list()
  partial_searches <- list()
  freed_parameters <- list()
  auto_partial_terms <- character()
  recovered_steps <- character()
  partial_failure <- NULL
  partial_search <- NULL
  recovered_partial <- FALSE
  stopped_early <- FALSE
  stopped_after_step <- NULL
  stopped_reason <- NULL

  for (i in seq_along(include_steps)) {
    step <- include_steps[[i]]
    prev_step <- if (i > 1L) include_steps[[i - 1L]] else NULL
    prev_fit <- if (!is.null(prev_step)) fits[[prev_step]] else NULL
    geq <- step_equal[[step]]
    added_constraints <- .mgcfa_added_constraints(step_equal, step = step, prev_step = prev_step)

    step_partial <- character()
    if (!is.null(partial) && !is.null(partial[[step]])) {
      step_partial <- .mgcfa_normalize_terms(partial[[step]])
    }

    if (!identical(step, "configural") && length(auto_partial_terms) > 0L) {
      step_partial <- unique(c(step_partial, auto_partial_terms))
    }

    args <- c(
      list(
        model = model_syntax,
        std.lv = if (identical(step, "configural")) isTRUE(std.lv_configural) else isTRUE(std.lv_constrained),
        orthogonal = isTRUE(orthogonal)
      ),
      fit_args_base
    )

    if (!is.null(estimator)) {
      args$estimator <- estimator
    }

    if (!is.null(geq)) {
      args$group.equal <- geq
    }

    if (length(step_partial) > 0L) {
      args$group.partial <- step_partial
    }

    fit_or_error <- tryCatch(
      do.call("sem", args, envir = asNamespace("lavaan")),
      error = function(e) e
    )

    if (inherits(fit_or_error, "error")) {
      can_evaluate <- !is.null(prev_fit) && !is.null(geq) && partial_failure_criterion != "none"
      if (can_evaluate) {
        fail_info <- list(
          failed = TRUE,
          reason = fit_or_error$message,
          criterion = partial_failure_criterion,
          threshold = partial_failure_threshold,
          value = NA_real_,
          gap = NA_real_,
          from_step = prev_step,
          to_step = step,
          failed_fit = NULL,
          failed_fit_measures = NULL
        )
        step_failures[[step]] <- fail_info

        do_search <- .mgcfa_should_run_partial_search(
          mode = partial_auto_search,
          reason = sprintf("Step `%s` estimation failed (%s).", step, fit_or_error$message)
        )

        search_criterion <- partial_search_criterion %||% "chisq_pvalue"
        search_threshold <- partial_search_threshold %||% .mgcfa_default_partial_threshold(search_criterion)
        search_measure <- partial_search_measure %||% partial_failure_measure
        search_direction <- partial_search_direction %||% partial_failure_direction
        step_candidate_source <- .mgcfa_resolve_candidate_source(
          step = step,
          source = partial_search_candidate_source,
          exhaustive_steps = partial_search_exhaustive_steps
        )
        step_stop_on_accept <- isTRUE(partial_search_stop_on_accept)
        if (identical(partial_search_candidate_source, "auto") &&
            (step %in% partial_search_exhaustive_steps) &&
            identical(step_candidate_source, "all")) {
          step_stop_on_accept <- FALSE
        }
        if (do_search) {
          search_out <- .mgcfa_search_partial_step(
            step = step,
            group_equal = geq,
            added_constraints = added_constraints,
            model_syntax = model_syntax,
            fit_args_base = fit_args_base,
            estimator = estimator,
            std_lv = isTRUE(std.lv_constrained),
            orthogonal = isTRUE(orthogonal),
            previous_fit = prev_fit,
            base_partial = step_partial,
            criterion = search_criterion,
            threshold = search_threshold,
            eval_measure = search_measure,
            eval_direction = search_direction,
            ic_bic_weight = partial_ic_bic_weight,
            max_free = partial_search_max_free,
            top_n = partial_search_top_n,
            stop_on_accept = step_stop_on_accept,
            rank = partial_search_rank,
            use_best_if_no_pass = isTRUE(partial_search_use_best_if_no_pass),
            candidate_source = step_candidate_source,
            max_models = partial_search_max_models,
            fit_measures = fit_measures
          )
          search_out$triggered <- TRUE
          search_out$decision <- "accepted"
          partial_searches[[step]] <- search_out

          if (!is.null(search_out$selected_fit)) {
            fits[[step]] <- search_out$selected_fit
            if (length(search_out$selected_added_terms) > 0L && isTRUE(search_out$selected_is_acceptable)) {
              auto_partial_terms <- unique(c(auto_partial_terms, search_out$selected_added_terms))
              freed_parameters[[step]] <- search_out$selected_added_terms
            }
            if (isTRUE(search_out$selected_is_acceptable)) {
              recovered_steps <- unique(c(recovered_steps, step))
              recovered_partial <- TRUE
            }
            next
          }
        } else {
          search_out <- list(
            triggered = TRUE,
            decision = if (partial_auto_search == "prompt" && interactive()) "declined" else "not_requested",
            criterion = partial_search_criterion %||% "chisq_pvalue",
            threshold = partial_search_threshold %||%
              .mgcfa_default_partial_threshold(
                partial_search_criterion %||% "chisq_pvalue"
              ),
            criterion_measure = if (identical(search_criterion, "measure_change")) search_measure else NULL,
            criterion_direction = if (identical(search_criterion, "measure_change")) search_direction else NULL,
            criterion_ic_bic_weight = if (identical(search_criterion, "aic_bic_weight")) partial_ic_bic_weight else NULL,
            candidates = data.frame(),
            top_models = data.frame(),
            selected_partial = character(),
            selected_index = NA_integer_,
            selected_fit = NULL,
            selected_added_terms = character()
          )
          partial_searches[[step]] <- search_out
        }

        if (is.null(partial_failure)) {
          partial_failure <- fail_info
        }
        if (is.null(partial_search)) {
          partial_search <- partial_searches[[step]]
        }
      }

      stop(sprintf("MGCFA failed at step `%s`: %s", step, fit_or_error$message), call. = FALSE)
    }

    fits[[step]] <- fit_or_error

    if (!is.null(prev_fit) && !is.null(geq) && partial_failure_criterion != "none") {
      step_eval <- .mgcfa_step_eval(
        previous_fit = prev_fit,
        candidate_fit = fit_or_error,
        criterion = partial_failure_criterion,
        threshold = partial_failure_threshold,
        measure = partial_failure_measure,
        direction = partial_failure_direction,
        ic_bic_weight = partial_ic_bic_weight
      )

      fail_info <- list(
        failed = !isTRUE(step_eval$pass),
        reason = if (isTRUE(step_eval$pass)) NULL else sprintf("Step `%s` failed `%s` criterion.", step, partial_failure_criterion),
        criterion = partial_failure_criterion,
        threshold = step_eval$threshold,
        value = step_eval$value,
        gap = step_eval$gap,
        from_step = prev_step,
        to_step = step,
        failed_fit = if (isTRUE(step_eval$pass)) NULL else fit_or_error,
        failed_fit_measures = if (isTRUE(step_eval$pass)) NULL else lavaan::fitMeasures(fit_or_error, fit_measures),
        details = step_eval$details
      )
      step_failures[[step]] <- fail_info
      if (is.null(partial_failure)) {
        partial_failure <- fail_info
      }

      if (!isTRUE(step_eval$pass)) {
        step_acceptable <- FALSE
        do_search <- .mgcfa_should_run_partial_search(
          mode = partial_auto_search,
          reason = sprintf("Step `%s` failed `%s` criterion.", step, partial_failure_criterion)
        )

        search_criterion <- partial_search_criterion %||% "chisq_pvalue"
        search_threshold <- partial_search_threshold %||% .mgcfa_default_partial_threshold(search_criterion)
        search_measure <- partial_search_measure %||% partial_failure_measure
        search_direction <- partial_search_direction %||% partial_failure_direction
        step_candidate_source <- .mgcfa_resolve_candidate_source(
          step = step,
          source = partial_search_candidate_source,
          exhaustive_steps = partial_search_exhaustive_steps
        )
        step_stop_on_accept <- isTRUE(partial_search_stop_on_accept)
        if (identical(partial_search_candidate_source, "auto") &&
            (step %in% partial_search_exhaustive_steps) &&
            identical(step_candidate_source, "all")) {
          step_stop_on_accept <- FALSE
        }
        if (do_search) {
          search_out <- .mgcfa_search_partial_step(
            step = step,
            group_equal = geq,
            added_constraints = added_constraints,
            model_syntax = model_syntax,
            fit_args_base = fit_args_base,
            estimator = estimator,
            std_lv = isTRUE(std.lv_constrained),
            orthogonal = isTRUE(orthogonal),
            previous_fit = prev_fit,
            base_partial = step_partial,
            criterion = search_criterion,
            threshold = search_threshold,
            eval_measure = search_measure,
            eval_direction = search_direction,
            ic_bic_weight = partial_ic_bic_weight,
            max_free = partial_search_max_free,
            top_n = partial_search_top_n,
            stop_on_accept = step_stop_on_accept,
            rank = partial_search_rank,
            use_best_if_no_pass = isTRUE(partial_search_use_best_if_no_pass),
            candidate_source = step_candidate_source,
            max_models = partial_search_max_models,
            fit_measures = fit_measures
          )
          search_out$triggered <- TRUE
          search_out$decision <- "accepted"
          partial_searches[[step]] <- search_out
          if (is.null(partial_search)) {
            partial_search <- search_out
          }

          if (!is.null(search_out$selected_fit)) {
            fits[[step]] <- search_out$selected_fit
            step_acceptable <- isTRUE(search_out$selected_is_acceptable)
            if (length(search_out$selected_added_terms) > 0L && isTRUE(step_acceptable)) {
              auto_partial_terms <- unique(c(auto_partial_terms, search_out$selected_added_terms))
              freed_parameters[[step]] <- search_out$selected_added_terms
            }
            if (isTRUE(step_acceptable)) {
              recovered_steps <- unique(c(recovered_steps, step))
              recovered_partial <- TRUE
            }
          }
        } else {
          search_out <- list(
            triggered = TRUE,
            decision = if (partial_auto_search == "prompt" && interactive()) "declined" else "not_requested",
            criterion = partial_search_criterion %||% "chisq_pvalue",
            threshold = partial_search_threshold %||% .mgcfa_default_partial_threshold(partial_search_criterion %||% "chisq_pvalue"),
            criterion_measure = if (identical(search_criterion, "measure_change")) search_measure else NULL,
            criterion_direction = if (identical(search_criterion, "measure_change")) search_direction else NULL,
            criterion_ic_bic_weight = if (identical(search_criterion, "aic_bic_weight")) partial_ic_bic_weight else NULL,
            candidates = data.frame(),
            top_models = data.frame(),
            selected_partial = character(),
            selected_index = NA_integer_,
            selected_fit = NULL,
            selected_added_terms = character()
          )
          partial_searches[[step]] <- search_out
          if (is.null(partial_search)) {
            partial_search <- search_out
          }
        }

        if (isTRUE(stop_at_first_unacceptable) && !isTRUE(step_acceptable)) {
          stopped_early <- TRUE
          stopped_after_step <- step
          stopped_reason <- sprintf(
            "Stopped after `%s` because invariance remained unacceptable under `%s` (threshold = %s).",
            step,
            partial_failure_criterion,
            format(step_eval$threshold, trim = TRUE)
          )
          break
        }
      }
    }
  }

  fit_table <- do.call(
    cbind,
    lapply(fits, function(fit) {
      lavaan::fitMeasures(fit, fit_measures)
    })
  )

  lrt_table <- .mgcfa_nested_lrt(fits)
  failed_step_outputs <- Filter(function(z) isTRUE(z$failed), step_failures)

  out <- list(
    call = match.call(),
    mode = if (using_raw) "raw_data" else "summary_matrices",
    model_type = model_type,
    model_syntax = model_syntax,
    fits = fits,
    fit_table = fit_table,
    lrt_table = lrt_table,
    step_failures = step_failures,
    failed_step_outputs = failed_step_outputs,
    partial_searches = partial_searches,
    freed_parameters = freed_parameters,
    recovered_steps = recovered_steps,
    partial_failure = partial_failure,
    partial_search = partial_search,
    recovered_partial = recovered_partial,
    stopped_early = stopped_early,
    stopped_after_step = stopped_after_step,
    stopped_reason = stopped_reason
  )
  class(out) <- "mgcfa_result"
  out
}

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
  cat("Steps:", paste(names(x$fits), collapse = ", "), "\n\n")

  overview <- .mgcfa_step_overview(x)
  if (nrow(overview) > 0L) {
    overview <- .mgcfa_format_numeric_df(overview, digits = digits, rounding = rounding)
    print(overview, row.names = FALSE)
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

    bridge_rows <- df[df$variant == "selected" & as.character(df$step) == "lv.variances", , drop = FALSE]
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

.mgcfa_step_overview <- function(x) {
  if (!inherits(x, "mgcfa_result") || is.null(x$fits) || length(x$fits) == 0L) {
    return(data.frame())
  }
  steps <- names(x$fits)
  out <- vector("list", length(steps))
  for (i in seq_along(steps)) {
    step <- steps[[i]]
    fm <- suppressWarnings(lavaan::fitMeasures(x$fits[[step]], c("cfi", "rmsea", "srmr", "aic", "bic")))
    freed <- .mgcfa_normalize_terms(x$freed_parameters[[step]] %||% character())
    fail_rec <- x$step_failures[[step]]
    status <- if (is.null(fail_rec) || !isTRUE(fail_rec$failed)) {
      "ok"
    } else if (step %in% (x$recovered_steps %||% character())) {
      "recovered_partial"
    } else {
      "failed"
    }
    out[[i]] <- data.frame(
      step = step,
      status = status,
      model_used = if (length(freed) > 0L) "partial" else "full",
      n_freed = length(freed),
      cfi = as.numeric(fm[["cfi"]]),
      rmsea = as.numeric(fm[["rmsea"]]),
      srmr = as.numeric(fm[["srmr"]]),
      aic = as.numeric(fm[["aic"]]),
      bic = as.numeric(fm[["bic"]]),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

.mgcfa_as_result <- function(x) {
  if (inherits(x, "mgcfa_result")) {
    return(x)
  }
  if (inherits(x, "lavaan")) {
    stop(
      "`x` looks like a lavaan fit object. Pass the object returned by `mgcfa_auto()` instead (for example, `out <- mgcfa_auto(...); mgcfa_tidy_fit(out)`).",
      call. = FALSE
    )
  }
  if (is.list(x) && !is.null(x$fits)) {
    if (is.null(x$model_type)) {
      x$model_type <- "unknown"
    }
    if (is.null(x$mode)) {
      x$mode <- "unknown"
    }
    class(x) <- unique(c("mgcfa_result", class(x)))
    return(x)
  }
  stop("`x` must be an object returned by `mgcfa_auto()`.", call. = FALSE)
}

.mgcfa_decimal_values <- function(x, digits = NULL, format = c("none", "truncate", "round")) {
  format <- match.arg(format)
  if (!is.numeric(x) || is.null(digits) || identical(format, "none")) {
    return(x)
  }
  digits <- as.integer(digits)
  if (length(digits) != 1L || is.na(digits) || digits < 0L) {
    stop("`digits` must be NULL or a non-negative integer.", call. = FALSE)
  }

  if (identical(format, "round")) {
    return(round(x, digits = digits))
  }

  scale <- 10^digits
  sign(x) * floor(abs(x) * scale) / scale
}

.mgcfa_decimal_matrix <- function(m, digits = NULL, format = c("none", "truncate", "round"), is_correlation = FALSE) {
  if (!is.matrix(m)) {
    stop("`m` must be a matrix.", call. = FALSE)
  }
  out <- .mgcfa_decimal_values(m, digits = digits, format = format)
  out <- (out + t(out)) / 2
  if (isTRUE(is_correlation)) {
    diag(out) <- 1
    out[out > 1] <- 1
    out[out < -1] <- -1
  }
  dimnames(out) <- dimnames(m)
  out
}

.mgcfa_round_values <- function(x, digits = 3L, rounding = c("signif", "round", "none")) {
  rounding <- match.arg(rounding)
  digits <- as.integer(digits)
  if (is.na(digits) || digits < 1L) {
    stop("`digits` must be a positive integer.", call. = FALSE)
  }
  if (identical(rounding, "none")) {
    return(x)
  }
  if (!is.numeric(x)) {
    return(x)
  }
  if (identical(rounding, "signif")) {
    return(signif(x, digits = digits))
  }
  round(x, digits = digits)
}

.mgcfa_format_numeric_df <- function(df, digits = 3L, rounding = c("signif", "round", "none")) {
  rounding <- match.arg(rounding)
  if (nrow(df) == 0L) {
    return(df)
  }
  is_num <- vapply(df, is.numeric, logical(1L))
  if (any(is_num)) {
    df[is_num] <- lapply(df[is_num], .mgcfa_round_values, digits = digits, rounding = rounding)
  }
  df
}

.mgcfa_n_latent <- function(x) {
  fit <- x$fits[[1L]] %||% NULL
  if (is.null(fit)) {
    return(NA_integer_)
  }
  lv <- tryCatch(
    lavaan::lavNames(fit, "lv"),
    error = function(e) NULL
  )
  if (is.null(lv) || length(lv) == 0L) {
    lv <- tryCatch(
      lavaan::lavInspect(fit, "lv.names"),
      error = function(e) NULL
    )
  }
  if (is.list(lv)) {
    lv <- unique(unlist(lv, use.names = FALSE))
  }
  if (is.null(lv) || length(lv) == 0L) {
    return(NA_integer_)
  }
  as.integer(length(unique(as.character(lv))))
}

.mgcfa_step_label_map <- function(steps, n_latent = NA_integer_, line_break = TRUE) {
  steps <- as.character(steps)
  line_sep <- if (isTRUE(line_break)) "\n" else " "
  latent_variance <- if (isTRUE(!is.na(n_latent) && n_latent == 1L)) {
    paste0("Latent", line_sep, "Variance")
  } else {
    paste0("Latent", line_sep, "Variances")
  }
  latent_mean <- if (isTRUE(!is.na(n_latent) && n_latent == 1L)) {
    paste0("Latent", line_sep, "Mean")
  } else {
    paste0("Latent", line_sep, "Means")
  }
  out <- stats::setNames(steps, steps)
  out[steps == "configural"] <- "Configural"
  out[steps == "metric"] <- "Metric"
  out[steps == "scalar"] <- "Scalar"
  out[steps == "strict"] <- "Strict"
  out[steps == "lv.variances"] <- latent_variance
  out[steps == "means"] <- latent_mean
  out
}

.mgcfa_measure_label <- function(measure) {
  key <- tolower(as.character(measure))
  switch(
    key,
    "aic" = "AIC",
    "bic" = "BIC",
    "cfi" = "CFI",
    "tli" = "TLI",
    "rmsea" = "RMSEA",
    "rmsea.ci.lower" = "RMSEA CI Lower",
    "rmsea.ci.upper" = "RMSEA CI Upper",
    "srmr" = "SRMR",
    "chisq" = "\u03C7\u00B2",
    "df" = "DF",
    "npar" = "N Parameters",
    gsub("\\.", " ", toupper(key))
  )
}

.mgcfa_default_partial_threshold <- function(criterion) {
  criterion <- match.arg(criterion, choices = c("chisq_pvalue", "delta_cfi", "aic_bic_weight", "measure_change"))
  if (identical(criterion, "chisq_pvalue")) {
    return(0.05)
  }
  if (identical(criterion, "delta_cfi")) {
    return(0.01)
  }
  if (identical(criterion, "aic_bic_weight")) {
    return(0.5)
  }
  0
}

.mgcfa_criterion_higher_is_better <- function(criterion) {
  criterion <- match.arg(criterion, choices = c("chisq_pvalue", "delta_cfi", "aic_bic_weight", "measure_change"))
  criterion %in% c("chisq_pvalue", "aic_bic_weight", "measure_change")
}

.mgcfa_pair_ic_weight <- function(previous_ic, candidate_ic) {
  vals <- c(as.numeric(previous_ic), as.numeric(candidate_ic))
  if (length(vals) != 2L || any(!is.finite(vals))) {
    return(NA_real_)
  }
  delta <- vals - min(vals)
  w <- exp(-0.5 * delta)
  denom <- sum(w)
  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }
  as.numeric(w[[2L]] / denom)
}

.mgcfa_normalize_terms <- function(x) {
  if (is.null(x)) {
    return(character())
  }
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\s+", "", x)
  x <- x[nzchar(x)]
  unique(x)
}

.mgcfa_added_constraints <- function(step_equal, step, prev_step = NULL) {
  cur <- step_equal[[step]] %||% character()
  prev <- if (is.null(prev_step)) character() else step_equal[[prev_step]] %||% character()
  setdiff(cur, prev)
}

.mgcfa_should_run_partial_search <- function(mode, reason) {
  mode <- match.arg(mode, choices = c("prompt", "never", "always"))
  if (identical(mode, "always")) {
    return(TRUE)
  }
  if (identical(mode, "never")) {
    return(FALSE)
  }
  if (!interactive()) {
    return(FALSE)
  }

  ans <- readline(
    prompt = paste0(reason, " Search for a partial invariance model automatically? [y/N]: ")
  )
  tolower(trimws(ans)) %in% c("y", "yes")
}

.mgcfa_resolve_candidate_source <- function(step, source = c("auto", "score", "all"), exhaustive_steps = c("lv.variances", "means")) {
  source <- match.arg(source)
  exhaustive_steps <- unique(as.character(exhaustive_steps %||% character()))
  if (!identical(source, "auto")) {
    return(source)
  }
  if (step %in% exhaustive_steps) {
    return("all")
  }
  "score"
}

.mgcfa_step_eval <- function(
  previous_fit,
  candidate_fit,
  criterion,
  threshold,
  measure = "aic",
  direction = c("decrease", "increase"),
  ic_bic_weight = 0.5
) {
  criterion <- match.arg(criterion, choices = c("chisq_pvalue", "delta_cfi", "aic_bic_weight", "measure_change"))
  direction <- match.arg(direction)
  if (is.null(threshold)) {
    threshold <- .mgcfa_default_partial_threshold(criterion)
  }

  if (identical(criterion, "chisq_pvalue")) {
    chisq_metric <- as.numeric(lavaan::fitMeasures(previous_fit, "chisq"))
    chisq_candidate <- as.numeric(lavaan::fitMeasures(candidate_fit, "chisq"))
    df_metric <- as.numeric(lavaan::fitMeasures(previous_fit, "df"))
    df_candidate <- as.numeric(lavaan::fitMeasures(candidate_fit, "df"))
    d_chisq <- chisq_candidate - chisq_metric
    d_df <- df_candidate - df_metric
    p_value <- if (is.finite(d_chisq) && is.finite(d_df) && d_df > 0) {
      stats::pchisq(d_chisq, df = d_df, lower.tail = FALSE)
    } else {
      NA_real_
    }
    gap <- p_value - threshold
    pass <- is.finite(gap) && gap >= 0
    return(list(
      pass = pass,
      value = p_value,
      threshold = threshold,
      gap = gap,
      details = list(
        chisq_diff = d_chisq,
        df_diff = d_df
      )
    ))
  }

  if (identical(criterion, "delta_cfi")) {
    cfi_metric <- as.numeric(lavaan::fitMeasures(previous_fit, "cfi"))
    cfi_candidate <- as.numeric(lavaan::fitMeasures(candidate_fit, "cfi"))
    delta_cfi <- cfi_metric - cfi_candidate
    gap <- threshold - delta_cfi
    pass <- is.finite(gap) && gap >= 0
    return(list(
      pass = pass,
      value = delta_cfi,
      threshold = threshold,
      gap = gap,
      details = list(
        previous_cfi = cfi_metric,
        candidate_cfi = cfi_candidate
      )
    ))
  }

  if (identical(criterion, "aic_bic_weight")) {
    prev_aic <- as.numeric(lavaan::fitMeasures(previous_fit, "aic"))
    cand_aic <- as.numeric(lavaan::fitMeasures(candidate_fit, "aic"))
    prev_bic <- as.numeric(lavaan::fitMeasures(previous_fit, "bic"))
    cand_bic <- as.numeric(lavaan::fitMeasures(candidate_fit, "bic"))
    aic_weight <- .mgcfa_pair_ic_weight(previous_ic = prev_aic, candidate_ic = cand_aic)
    bic_weight <- .mgcfa_pair_ic_weight(previous_ic = prev_bic, candidate_ic = cand_bic)
    combined_weight <- (1 - ic_bic_weight) * aic_weight + ic_bic_weight * bic_weight
    gap <- combined_weight - threshold
    pass <- is.finite(gap) && gap >= 0
    return(list(
      pass = pass,
      value = combined_weight,
      threshold = threshold,
      gap = gap,
      details = list(
        previous_aic = prev_aic,
        candidate_aic = cand_aic,
        previous_bic = prev_bic,
        candidate_bic = cand_bic,
        aic_weight = aic_weight,
        bic_weight = bic_weight,
        ic_bic_weight = ic_bic_weight
      )
    ))
  }

  measure <- as.character(measure)[[1L]]
  if (!nzchar(measure)) {
    stop("`measure` must be a non-empty fit measure name for `measure_change`.", call. = FALSE)
  }
  previous_value <- as.numeric(lavaan::fitMeasures(previous_fit, measure))
  candidate_value <- as.numeric(lavaan::fitMeasures(candidate_fit, measure))
  change_signed <- if (identical(direction, "decrease")) {
    previous_value - candidate_value
  } else {
    candidate_value - previous_value
  }
  gap <- change_signed - threshold
  pass <- is.finite(gap) && gap >= 0
  list(
    pass = pass,
    value = change_signed,
    threshold = threshold,
    gap = gap,
    details = list(
      measure = measure,
      direction = direction,
      previous_value = previous_value,
      candidate_value = candidate_value,
      raw_change = candidate_value - previous_value
    )
  )
}

.mgcfa_search_partial_step <- function(
  step,
  group_equal,
  added_constraints,
  model_syntax,
  fit_args_base,
  estimator,
  std_lv,
  orthogonal,
  previous_fit,
  base_partial,
  criterion,
  threshold,
  eval_measure,
  eval_direction,
  ic_bic_weight,
  max_free,
  top_n,
  stop_on_accept,
  rank,
  use_best_if_no_pass,
  candidate_source,
  max_models,
  fit_measures
) {
  criterion <- match.arg(criterion, choices = c("chisq_pvalue", "delta_cfi", "aic_bic_weight", "measure_change"))
  eval_direction <- match.arg(eval_direction, choices = c("decrease", "increase"))
  rank <- match.arg(rank, choices = c("closest", "best"))
  candidate_source <- match.arg(candidate_source, choices = c("score", "all"))
  max_models <- as.integer(max_models)
  if (is.na(max_models) || max_models < 1L) {
    max_models <- 5000L
  }
  if (is.null(threshold)) {
    threshold <- .mgcfa_default_partial_threshold(criterion)
  }

  base_partial <- .mgcfa_normalize_terms(base_partial %||% character())
  added_constraints <- unique(as.character(added_constraints %||% character()))
  min_remaining <- 1L
  truncated <- FALSE

  args0 <- c(
    list(
      model = model_syntax,
      std.lv = isTRUE(std_lv),
      orthogonal = isTRUE(orthogonal)
    ),
    fit_args_base
  )
  if (!is.null(estimator)) {
    args0$estimator <- estimator
  }
  if (!is.null(group_equal)) {
    args0$group.equal <- group_equal
  }
  if (length(base_partial) > 0L) {
    args0$group.partial <- base_partial
  }
  fit0 <- tryCatch(
    do.call("sem", args0, envir = asNamespace("lavaan")),
    error = function(e) NULL
  )

  score_candidates <- data.frame()
  all_candidates <- data.frame()
  if (!is.null(fit0)) {
    score_candidates <- .mgcfa_score_release_candidates(
      fit = fit0,
      target_constraints = added_constraints,
      already_freed = base_partial
    )
    all_candidates <- .mgcfa_all_release_candidates(
      fit = fit0,
      target_constraints = added_constraints,
      already_freed = base_partial
    )
  }

  start_candidates <- if (identical(candidate_source, "all")) all_candidates else score_candidates
  releasable_terms <- unique(as.character(start_candidates$term))
  n_releasable <- length(releasable_terms)
  score_lookup <- if (nrow(score_candidates) > 0L) {
    tapply(score_candidates$x2, score_candidates$term, max, na.rm = TRUE)
  } else {
    numeric()
  }
  releasable_scores <- if (n_releasable > 0L) {
    out <- as.numeric(score_lookup[releasable_terms])
    out[is.na(out)] <- 0
    out
  } else {
    numeric()
  }

  max_allow_by_remaining <- max(0L, as.integer(n_releasable) - min_remaining)
  if (is.null(max_free)) {
    max_free <- max_allow_by_remaining
  } else {
    max_free <- min(as.integer(max_free), max_allow_by_remaining)
  }
  max_free <- max(0L, as.integer(max_free))

  candidate_rows <- list()
  candidate_fits <- list()
  candidate_partials <- list()
  candidate_added <- list()
  iter_id <- -1L

  for (k in 0:max_free) {
    if (k == 0L) {
      combos_k <- list(character())
    } else if (k > n_releasable) {
      combos_k <- list()
    } else {
      combos_idx <- utils::combn(seq_len(n_releasable), k, simplify = FALSE)
      if (length(combos_idx) > 1L) {
        score_sum <- vapply(combos_idx, function(ix) sum(releasable_scores[ix]), numeric(1L))
        combos_idx <- combos_idx[order(score_sum, decreasing = TRUE)]
      }
      combos_k <- lapply(combos_idx, function(ix) releasable_terms[ix])
    }

    if (length(combos_k) == 0L) {
      next
    }

    pass_found_k <- FALSE
    for (added_terms in combos_k) {
      if ((iter_id + 1L) >= max_models) {
        truncated <- TRUE
        break
      }
      iter_id <- iter_id + 1L

      partial_now <- unique(c(base_partial, added_terms))
      args <- c(
        list(
          model = model_syntax,
          std.lv = isTRUE(std_lv),
          orthogonal = isTRUE(orthogonal)
        ),
        fit_args_base
      )

      if (!is.null(estimator)) {
        args$estimator <- estimator
      }
      if (!is.null(group_equal)) {
        args$group.equal <- group_equal
      }
      if (length(partial_now) > 0L) {
        args$group.partial <- partial_now
      }

      fit_or_error <- tryCatch(
        do.call("sem", args, envir = asNamespace("lavaan")),
        error = function(e) e
      )

      is_ok <- !inherits(fit_or_error, "error")
      step_eval <- list(pass = FALSE, value = NA_real_, gap = NA_real_)
      fit_stats <- stats::setNames(rep(NA_real_, 5L), c("chisq", "df", "cfi", "aic", "bic"))

      if (is_ok) {
        step_eval <- .mgcfa_step_eval(
          previous_fit = previous_fit,
          candidate_fit = fit_or_error,
          criterion = criterion,
          threshold = threshold,
          measure = eval_measure,
          direction = eval_direction,
          ic_bic_weight = ic_bic_weight
        )
        fit_stats <- lavaan::fitMeasures(fit_or_error, c("chisq", "df", "cfi", "aic", "bic"))
      }

      candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(
        iteration = iter_id,
        step = step,
        added_constraints = paste(added_terms, collapse = ", "),
        constraint_types = paste(added_constraints, collapse = ", "),
        n_added = length(added_terms),
        n_partial_total = length(partial_now),
        criterion = criterion,
        criterion_measure = if (identical(criterion, "measure_change")) eval_measure else "",
        criterion_direction = if (identical(criterion, "measure_change")) eval_direction else "",
        criterion_ic_bic_weight = if (identical(criterion, "aic_bic_weight")) ic_bic_weight else NA_real_,
        criterion_value = as.numeric(step_eval$value),
        threshold = threshold,
        criterion_gap = as.numeric(step_eval$gap),
        pass = isTRUE(step_eval$pass),
        chisq = as.numeric(fit_stats[["chisq"]]),
        df = as.numeric(fit_stats[["df"]]),
        cfi = as.numeric(fit_stats[["cfi"]]),
        aic = as.numeric(fit_stats[["aic"]]),
        bic = as.numeric(fit_stats[["bic"]]),
        status = if (is_ok) "ok" else "error",
        error_message = if (is_ok) NA_character_ else fit_or_error$message,
        stringsAsFactors = FALSE
      )
      candidate_fits[[length(candidate_fits) + 1L]] <- if (is_ok) fit_or_error else NULL
      candidate_partials[[length(candidate_partials) + 1L]] <- partial_now
      candidate_added[[length(candidate_added) + 1L]] <- added_terms

      if (is_ok && isTRUE(step_eval$pass)) {
        pass_found_k <- TRUE
      }
    }

    if (truncated) {
      break
    }
    if (isTRUE(stop_on_accept) && pass_found_k) {
      break
    }
  }

  if (truncated) {
    warning(
      sprintf(
        "Partial search at step `%s` stopped after %d candidate models (limit reached).",
        step, max_models
      ),
      call. = FALSE
    )
  }

  candidates <- if (length(candidate_rows) > 0L) {
    do.call(rbind, candidate_rows)
  } else {
    data.frame()
  }

  base_ok_idx <- which(candidates$iteration == 0L & candidates$status == "ok")
  criterion_higher_is_better <- .mgcfa_criterion_higher_is_better(criterion)
  if (length(base_ok_idx) > 0L) {
    b <- base_ok_idx[[1L]]
    base_chisq <- as.numeric(candidates$chisq[[b]])
    base_cfi <- as.numeric(candidates$cfi[[b]])
    base_criterion <- as.numeric(candidates$criterion_value[[b]])
    candidates$chisq_improve_from_initial <- ifelse(
      candidates$status == "ok" & is.finite(base_chisq),
      base_chisq - as.numeric(candidates$chisq),
      NA_real_
    )
    candidates$cfi_improve_from_initial <- ifelse(
      candidates$status == "ok" & is.finite(base_cfi),
      as.numeric(candidates$cfi) - base_cfi,
      NA_real_
    )
    candidates$criterion_improve_from_initial <- ifelse(
      candidates$status == "ok" & is.finite(base_criterion),
      if (isTRUE(criterion_higher_is_better)) {
        as.numeric(candidates$criterion_value) - base_criterion
      } else {
        base_criterion - as.numeric(candidates$criterion_value)
      },
      NA_real_
    )
  } else {
    candidates$chisq_improve_from_initial <- NA_real_
    candidates$cfi_improve_from_initial <- NA_real_
    candidates$criterion_improve_from_initial <- NA_real_
  }

  ranked_idx <- .mgcfa_rank_partial_candidates(candidates, rank = rank)
  closest_idx <- utils::head(ranked_idx, n = min(length(ranked_idx), top_n))
  closest_models <- if (length(closest_idx) > 0L) candidates[closest_idx, , drop = FALSE] else candidates[0, , drop = FALSE]

  acceptable_ranked_idx <- .mgcfa_rank_acceptable_partial_candidates(candidates)
  top_idx <- if (length(acceptable_ranked_idx) > 0L) {
    utils::head(acceptable_ranked_idx, n = min(length(acceptable_ranked_idx), top_n))
  } else {
    closest_idx
  }
  top_models <- if (length(top_idx) > 0L) candidates[top_idx, , drop = FALSE] else candidates[0, , drop = FALSE]
  acceptable_models <- if (length(acceptable_ranked_idx) > 0L) {
    candidates[acceptable_ranked_idx, , drop = FALSE]
  } else {
    candidates[0, , drop = FALSE]
  }

  pass_idx <- which(candidates$status == "ok" & candidates$pass)
  selected_idx <- NA_integer_
  if (length(pass_idx) > 0L) {
    selected_idx <- if (isTRUE(stop_on_accept)) {
      min_added <- min(candidates$n_added[pass_idx], na.rm = TRUE)
      pool <- pass_idx[candidates$n_added[pass_idx] == min_added]
      ranked_pass <- acceptable_ranked_idx[acceptable_ranked_idx %in% pool]
      if (length(ranked_pass) > 0L) ranked_pass[[1L]] else pool[[1L]]
    } else {
      ranked_pass <- acceptable_ranked_idx[acceptable_ranked_idx %in% pass_idx]
      if (length(ranked_pass) > 0L) ranked_pass[[1L]] else pass_idx[[1L]]
    }
  } else if (isTRUE(use_best_if_no_pass) && length(ranked_idx) > 0L) {
    selected_idx <- ranked_idx[[1L]]
  }

  selected_fit <- if (!is.na(selected_idx)) candidate_fits[[selected_idx]] else NULL
  selected_partial <- if (!is.na(selected_idx)) candidate_partials[[selected_idx]] else character()
  selected_added_terms <- if (!is.na(selected_idx)) candidate_added[[selected_idx]] else character()

  list(
    step = step,
    added_constraints = added_constraints,
    candidate_source = candidate_source,
    total_releasable = as.integer(n_releasable),
    max_free_allowed = as.integer(max_free),
    max_models = max_models,
    truncated = isTRUE(truncated),
    evaluated_models = as.integer(nrow(candidates)),
    criterion = criterion,
    criterion_measure = if (identical(criterion, "measure_change")) eval_measure else NULL,
    criterion_direction = if (identical(criterion, "measure_change")) eval_direction else NULL,
    criterion_ic_bic_weight = if (identical(criterion, "aic_bic_weight")) ic_bic_weight else NULL,
    threshold = threshold,
    candidates = candidates,
    top_models = top_models,
    acceptable_models = acceptable_models,
    closest_models = closest_models,
    selected_index = selected_idx,
    selected_fit = selected_fit,
    selected_partial = selected_partial,
    selected_added_terms = selected_added_terms,
    selected_is_acceptable = if (!is.na(selected_idx)) isTRUE(candidates$pass[[selected_idx]]) else FALSE,
    top_fits = if (length(top_idx) > 0L) candidate_fits[top_idx] else list(),
    closest_fits = if (length(closest_idx) > 0L) candidate_fits[closest_idx] else list()
  )
}

.mgcfa_constraint_class <- function(op, lhs, rhs, ov_names, lv_names) {
  if (identical(op, "=~")) {
    return("loadings")
  }
  if (identical(op, "~1")) {
    if (lhs %in% ov_names) {
      return("intercepts")
    }
    if (lhs %in% lv_names) {
      return("means")
    }
    return(NA_character_)
  }
  if (identical(op, "~~") && identical(lhs, rhs)) {
    if (lhs %in% ov_names) {
      return("residuals")
    }
    if (lhs %in% lv_names) {
      return("lv.variances")
    }
  }
  NA_character_
}

.mgcfa_term_from_parts <- function(lhs, op, rhs) {
  lhs <- as.character(lhs)
  op <- as.character(op)
  rhs <- as.character(rhs)
  if (identical(op, "~1")) {
    return(.mgcfa_normalize_terms(paste(lhs, "~ 1"))[[1L]])
  }
  if (identical(rhs, "")) {
    rhs <- "1"
  }
  .mgcfa_normalize_terms(paste(lhs, op, rhs))[[1L]]
}

.mgcfa_all_release_candidates <- function(fit, target_constraints, already_freed = character()) {
  target_constraints <- unique(as.character(target_constraints %||% character()))
  already_freed <- .mgcfa_normalize_terms(already_freed)
  empty <- data.frame(term = character(), constraint = character(), x2 = numeric(), stringsAsFactors = FALSE)
  if (length(target_constraints) == 0L) {
    return(empty)
  }

  pt <- tryCatch(
    lavaan::parTable(fit),
    error = function(e) NULL
  )
  if (is.null(pt) || nrow(pt) == 0L) {
    return(empty)
  }

  ov_names <- lavaan::lavNames(fit, type = "ov")
  lv_names <- lavaan::lavNames(fit, type = "lv")
  key <- paste(pt$lhs, pt$op, pt$rhs, pt$group, sep = "\r")
  pt <- pt[!duplicated(key), , drop = FALSE]

  rows <- vector("list", nrow(pt))
  out_i <- 0L
  for (i in seq_len(nrow(pt))) {
    cls <- .mgcfa_constraint_class(
      op = as.character(pt$op[[i]]),
      lhs = as.character(pt$lhs[[i]]),
      rhs = as.character(pt$rhs[[i]]),
      ov_names = ov_names,
      lv_names = lv_names
    )
    if (is.na(cls) || !(cls %in% target_constraints)) {
      next
    }
    term <- .mgcfa_term_from_parts(
      lhs = as.character(pt$lhs[[i]]),
      op = as.character(pt$op[[i]]),
      rhs = as.character(pt$rhs[[i]])
    )
    if (term %in% already_freed) {
      next
    }

    out_i <- out_i + 1L
    rows[[out_i]] <- data.frame(
      term = term,
      constraint = cls,
      group = as.integer(pt$group[[i]]),
      label = as.character(pt$label[[i]]),
      plabel = as.character(pt$plabel[[i]]),
      stringsAsFactors = FALSE
    )
  }

  if (out_i == 0L) {
    return(empty)
  }

  out <- do.call(rbind, rows[seq_len(out_i)])
  n_group_by_term <- tapply(out$group, out$term, function(g) length(unique(g)))
  terms_multi_group <- names(n_group_by_term)[n_group_by_term >= 2L]
  if (length(terms_multi_group) == 0L) {
    return(empty)
  }
  out <- out[out$term %in% terms_multi_group, , drop = FALSE]
  if (nrow(out) == 0L) {
    return(empty)
  }

  # Keep terms that still look equality-constrained across groups.
  eq_by_term <- tapply(seq_len(nrow(out)), out$term, function(ix) {
    labs <- out$label[ix]
    labs <- labs[nzchar(labs)]
    plabs <- out$plabel[ix]
    plabs <- plabs[nzchar(plabs)]
    any(duplicated(labs)) || any(duplicated(plabs))
  })
  terms_eq <- names(eq_by_term)[eq_by_term]
  if (length(terms_eq) == 0L) {
    terms_eq <- terms_multi_group
  }

  out <- out[out$term %in% terms_eq, c("term", "constraint"), drop = FALSE]
  out <- unique(out)
  out$x2 <- 0
  rownames(out) <- NULL
  out
}

.mgcfa_score_release_candidates <- function(fit, target_constraints, already_freed = character()) {
  target_constraints <- unique(as.character(target_constraints %||% character()))
  already_freed <- .mgcfa_normalize_terms(already_freed)
  empty <- data.frame(term = character(), constraint = character(), x2 = numeric(), stringsAsFactors = FALSE)
  if (length(target_constraints) == 0L) {
    return(empty)
  }

  score <- tryCatch(
    lavaan::lavTestScore(fit),
    error = function(e) NULL
  )
  if (is.null(score) || is.null(score$uni) || nrow(score$uni) == 0L) {
    return(empty)
  }

  uni <- score$uni
  pt <- lavaan::parTable(fit)
  x2_col <- which(tolower(names(uni)) == "x2")
  if (length(x2_col) != 1L) {
    return(empty)
  }
  ov_names <- lavaan::lavNames(fit, type = "ov")
  lv_names <- lavaan::lavNames(fit, type = "lv")

  labels_lhs <- as.character(uni$lhs)
  labels_rhs <- as.character(uni$rhs)
  x2_vals <- as.numeric(uni[[x2_col]])

  score_rows <- vector("list", length(x2_vals))
  out_i <- 0L
  for (i in seq_along(x2_vals)) {
    labels <- c(labels_lhs[[i]], labels_rhs[[i]])
    hit <- pt$plabel %in% labels
    if (!any(hit)) {
      next
    }
    pt_i <- pt[hit, , drop = FALSE]
    key <- paste(pt_i$lhs, pt_i$op, pt_i$rhs, sep = "\r")
    pt_i <- pt_i[!duplicated(key), , drop = FALSE]
    for (j in seq_len(nrow(pt_i))) {
      cls <- .mgcfa_constraint_class(
        op = as.character(pt_i$op[[j]]),
        lhs = as.character(pt_i$lhs[[j]]),
        rhs = as.character(pt_i$rhs[[j]]),
        ov_names = ov_names,
        lv_names = lv_names
      )
      if (is.na(cls) || !(cls %in% target_constraints)) {
        next
      }
      term <- .mgcfa_term_from_parts(
        lhs = as.character(pt_i$lhs[[j]]),
        op = as.character(pt_i$op[[j]]),
        rhs = as.character(pt_i$rhs[[j]])
      )
      if (term %in% already_freed) {
        next
      }
      out_i <- out_i + 1L
      score_rows[[out_i]] <- data.frame(
        term = term,
        constraint = cls,
        x2 = x2_vals[[i]],
        stringsAsFactors = FALSE
      )
    }
  }

  if (out_i == 0L) {
    return(empty)
  }

  out <- do.call(rbind, score_rows[seq_len(out_i)])
  out <- stats::aggregate(x2 ~ term + constraint, data = out, FUN = max)
  ord <- order(out$x2, decreasing = TRUE, na.last = NA)
  if (length(ord) == 0L) {
    return(empty)
  }
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}

.mgcfa_rank_partial_candidates <- function(candidates, rank = c("closest", "best")) {
  rank <- match.arg(rank)
  if (is.null(candidates) || nrow(candidates) == 0L) {
    return(integer())
  }
  ok_idx <- which(candidates$status == "ok")
  if (length(ok_idx) == 0L) {
    return(integer())
  }

  gap <- as.numeric(candidates$criterion_gap[ok_idx])
  gap_abs <- abs(gap)
  gap_abs[!is.finite(gap_abs)] <- Inf
  gap_best <- gap
  gap_best[!is.finite(gap_best)] <- -Inf
  n_added <- as.numeric(candidates$n_added[ok_idx])
  pass_penalty <- ifelse(candidates$pass[ok_idx], 0L, 1L)

  ord <- if (identical(rank, "closest")) {
    order(pass_penalty, gap_abs, n_added, na.last = TRUE)
  } else {
    order(pass_penalty, -gap_best, n_added, na.last = TRUE)
  }
  ok_idx[ord]
}

.mgcfa_rank_acceptable_partial_candidates <- function(candidates) {
  if (is.null(candidates) || nrow(candidates) == 0L) {
    return(integer())
  }
  ok_pass_idx <- which(candidates$status == "ok" & candidates$pass)
  if (length(ok_pass_idx) == 0L) {
    return(integer())
  }

  n_added <- as.numeric(candidates$n_added[ok_pass_idx])
  gap <- as.numeric(candidates$criterion_gap[ok_pass_idx])
  gap[!is.finite(gap)] <- -Inf
  crit_improve <- as.numeric(candidates$criterion_improve_from_initial[ok_pass_idx])
  crit_improve[!is.finite(crit_improve)] <- -Inf
  chisq_improve <- as.numeric(candidates$chisq_improve_from_initial[ok_pass_idx])
  chisq_improve[!is.finite(chisq_improve)] <- -Inf
  cfi_improve <- as.numeric(candidates$cfi_improve_from_initial[ok_pass_idx])
  cfi_improve[!is.finite(cfi_improve)] <- -Inf

  ord <- order(n_added, -gap, -crit_improve, -chisq_improve, -cfi_improve, na.last = TRUE)
  ok_pass_idx[ord]
}

.mgcfa_prepare_matrix_input <- function(
  sample_cov,
  sample_mean,
  sample_nobs,
  group_labels,
  matrices_are_cor,
  sample_sd
) {
  cov_list <- if (is.matrix(sample_cov)) list(sample_cov) else sample_cov
  if (!is.list(cov_list) || length(cov_list) < 1L) {
    stop("`sample_cov` must be a matrix or a non-empty list of matrices.", call. = FALSE)
  }

  cov_list <- lapply(cov_list, function(m) {
    if (!is.matrix(m) || nrow(m) != ncol(m)) {
      stop("Each element of `sample_cov` must be a square matrix.", call. = FALSE)
    }
    storage.mode(m) <- "double"
    m
  })

  k <- length(cov_list)
  p <- ncol(cov_list[[1L]])

  if (!all(vapply(cov_list, function(m) ncol(m) == p, logical(1L)))) {
    stop("All `sample_cov` matrices must have identical dimensions.", call. = FALSE)
  }

  if (isTRUE(matrices_are_cor)) {
    if (is.null(sample_sd)) {
      stop("`sample_sd` is required when `matrices_are_cor = TRUE`.", call. = FALSE)
    }
    sd_list <- .mgcfa_expand_group_input(sample_sd, k, "sample_sd")
    sd_list <- lapply(sd_list, as.numeric)
    if (!all(vapply(sd_list, length, integer(1L)) == p)) {
      stop("Each `sample_sd` vector must match matrix dimension.", call. = FALSE)
    }
    cov_list <- Map(lavaan::cor2cov, cov_list, sd_list)
  }

  ref_names <- colnames(cov_list[[1L]])
  if (is.null(ref_names)) {
    stop("`sample_cov` matrices must have column names.", call. = FALSE)
  }

  cov_list <- lapply(cov_list, function(m) {
    if (is.null(colnames(m))) {
      colnames(m) <- ref_names
      rownames(m) <- ref_names
    }
    if (!identical(colnames(m), ref_names) || !identical(rownames(m), ref_names)) {
      stop("All `sample_cov` matrices must share identical row/column names and order.", call. = FALSE)
    }
    m
  })

  nobs_list <- .mgcfa_expand_group_input(sample_nobs, k, "sample_nobs")
  nobs_list <- lapply(nobs_list, function(x) {
    x <- as.numeric(x)
    if (length(x) != 1L || is.na(x) || x <= 0) {
      stop("Each `sample_nobs` value must be a positive number.", call. = FALSE)
    }
    x
  })

  mean_list <- NULL
  if (!is.null(sample_mean)) {
    mean_list <- .mgcfa_expand_group_input(sample_mean, k, "sample_mean")
    mean_list <- lapply(mean_list, as.numeric)
    if (!all(vapply(mean_list, length, integer(1L)) == p)) {
      stop("Each `sample_mean` vector must match matrix dimension.", call. = FALSE)
    }
    mean_list <- lapply(mean_list, function(v) {
      names(v) <- ref_names
      v
    })
  }

  if (!is.null(group_labels)) {
    if (length(group_labels) != k) {
      stop("`group_labels` length must match the number of groups.", call. = FALSE)
    }
    group_labels <- as.character(group_labels)
  }

  list(
    sample_cov = cov_list,
    sample_mean = mean_list,
    sample_nobs = nobs_list,
    group_labels = group_labels
  )
}

.mgcfa_expand_group_input <- function(x, k, name) {
  if (is.null(x)) {
    stop(sprintf("`%s` must be supplied in matrix mode.", name), call. = FALSE)
  }
  if (is.list(x)) {
    if (length(x) != k) {
      stop(sprintf("`%s` must have one element per group.", name), call. = FALSE)
    }
    return(x)
  }
  if (length(x) == 1L) {
    return(rep(list(x), k))
  }
  if (length(x) == k) {
    return(as.list(x))
  }
  stop(sprintf("`%s` must be length 1, length K, or a list of length K.", name), call. = FALSE)
}

.mgcfa_exploratory_loadings <- function(
  model_type,
  variables,
  data,
  input_matrix,
  matrix_nobs,
  n_factors,
  rotation
) {
  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      stop("`data` must be a data.frame for EFA/PCA model generation.", call. = FALSE)
    }
    if (!all(variables %in% names(data))) {
      stop("All `variables` must exist in `data`.", call. = FALSE)
    }
    x <- data[, variables, drop = FALSE]
    if (!all(vapply(x, is.numeric, logical(1L)))) {
      stop("All selected `variables` must be numeric for EFA/PCA.", call. = FALSE)
    }
    x <- stats::na.omit(x)
    if (nrow(x) < (length(variables) + 2L)) {
      stop("Insufficient complete rows for EFA/PCA model generation.", call. = FALSE)
    }
  } else if (!is.null(input_matrix)) {
    if (!is.matrix(input_matrix) || nrow(input_matrix) != ncol(input_matrix)) {
      stop("`input_matrix` must be a square matrix for EFA/PCA model generation.", call. = FALSE)
    }
    if (is.null(colnames(input_matrix))) {
      colnames(input_matrix) <- variables
      rownames(input_matrix) <- variables
    }
    x <- stats::cov2cor(input_matrix)
  } else {
    stop("EFA/PCA model generation requires `data` or `input_matrix`.", call. = FALSE)
  }

  if (model_type == "efa") {
    efa_fit <- if (is.data.frame(x)) {
      psych::fa(x, nfactors = n_factors, rotate = rotation, fm = "ml")
    } else {
      psych::fa(r = x, nfactors = n_factors, n.obs = matrix_nobs %||% 1000, rotate = rotation, fm = "ml")
    }
    loadings <- as.matrix(efa_fit$loadings)
  } else {
    pca_fit <- if (is.data.frame(x)) {
      psych::principal(x, nfactors = n_factors, rotate = rotation)
    } else {
      psych::principal(r = x, nfactors = n_factors, n.obs = matrix_nobs %||% 1000, rotate = rotation)
    }
    loadings <- as.matrix(pca_fit$loadings)
  }

  if (is.null(rownames(loadings))) {
    rownames(loadings) <- variables
  }
  loadings
}

.mgcfa_syntax_from_loadings <- function(loadings, loading_threshold, allow_cross_loadings) {
  if (!is.matrix(loadings) || nrow(loadings) < 2L || ncol(loadings) < 1L) {
    stop("Exploratory loadings are malformed.", call. = FALSE)
  }

  abs_l <- abs(loadings)
  keep <- abs_l >= loading_threshold

  if (!allow_cross_loadings) {
    keep[,] <- FALSE
    strongest <- max.col(abs_l, ties.method = "first")
    keep[cbind(seq_len(nrow(abs_l)), strongest)] <- TRUE
  } else {
    none <- rowSums(keep) == 0L
    if (any(none)) {
      strongest <- max.col(abs_l[none, , drop = FALSE], ties.method = "first")
      keep[cbind(which(none), strongest)] <- TRUE
    }
  }

  item_names <- rownames(loadings)
  factor_names <- paste0("F", seq_len(ncol(loadings)))
  lines <- character()
  for (j in seq_len(ncol(keep))) {
    items <- item_names[keep[, j]]
    if (length(items) == 0L) {
      next
    }
    if (length(items) < 2L) {
      warning(sprintf("Factor `%s` has fewer than two indicators after thresholding.", factor_names[j]))
    }
    lines <- c(lines, paste0(factor_names[j], " =~ ", paste(items, collapse = " + ")))
  }

  if (length(lines) == 0L) {
    stop("No loadings survived thresholding; lower `loading_threshold` or inspect input.", call. = FALSE)
  }
  paste(lines, collapse = "\n")
}

.mgcfa_nested_lrt <- function(fits) {
  if (length(fits) < 2L) {
    return(data.frame())
  }

  step_names <- names(fits)
  out <- vector("list", length(fits) - 1L)
  for (i in 2:length(fits)) {
    prev <- fits[[i - 1L]]
    cur <- fits[[i]]
    c_prev <- lavaan::fitMeasures(prev, "chisq")
    c_cur <- lavaan::fitMeasures(cur, "chisq")
    d_prev <- lavaan::fitMeasures(prev, "df")
    d_cur <- lavaan::fitMeasures(cur, "df")
    d_chisq <- as.numeric(c_cur - c_prev)
    d_df <- as.numeric(d_cur - d_prev)
    p_val <- if (is.finite(d_chisq) && is.finite(d_df) && d_df > 0) {
      stats::pchisq(d_chisq, df = d_df, lower.tail = FALSE)
    } else {
      NA_real_
    }
    out[[i - 1L]] <- data.frame(
      from = step_names[i - 1L],
      to = step_names[i],
      chisq_diff = d_chisq,
      df_diff = d_df,
      p_value = p_val,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("step", "value_plot", "series", "variant"))
}
