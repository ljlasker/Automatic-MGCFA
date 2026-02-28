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
#' @param partial_failure_rules Optional list of rule objects for multi-metric
#'   failure decisions. Each rule is a named list with at least
#'   \code{criterion}, and optionally \code{threshold}, \code{measure},
#'   \code{direction}, and \code{ic_bic_weight}.
#' @param partial_failure_rule_policy Aggregation policy for multi-rule failure
#'   decisions: \code{"all"}, \code{"majority"}, \code{"any"}, or
#'   \code{"at_least"}.
#' @param partial_failure_rule_min Minimum number of passing rules required when
#'   \code{partial_failure_rule_policy = "at_least"}.
#' @param partial_failure_rules_by_step Optional named list of step-specific
#'   multi-rule definitions overriding \code{partial_failure_rules}.
#' @param partial_failure_rule_policy_by_step Optional named vector/list of
#'   step-specific aggregation policies for failure decisions.
#' @param partial_failure_rule_min_by_step Optional named vector/list of
#'   step-specific minimum pass counts for failure decisions.
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
#'   free during automatic partial search at the failed step.
#' @param partial_search_top_n Number of closest candidate models returned from
#'   partial search.
#' @param partial_search_stop_on_accept Logical; stop search at first acceptable
#'   candidate when \code{TRUE}.
#' @param partial_search_rank Ranking for candidate output:
#'   \code{"closest"} or \code{"best"}.
#' @param partial_search_use_best_if_no_pass Logical; if no candidate meets the
#'   threshold, continue invariance testing using the best ranked candidate.
#' @param partial_search_rules Optional list of rule objects for multi-metric
#'   partial-search acceptance decisions. Rule format matches
#'   \code{partial_failure_rules}.
#' @param partial_search_rule_policy Aggregation policy for multi-rule search
#'   decisions: \code{"all"}, \code{"majority"}, \code{"any"}, or
#'   \code{"at_least"}.
#' @param partial_search_rule_min Minimum number of passing rules required when
#'   \code{partial_search_rule_policy = "at_least"}.
#' @param partial_search_rules_by_step Optional named list of step-specific
#'   multi-rule definitions overriding \code{partial_search_rules}.
#' @param partial_search_rule_policy_by_step Optional named vector/list of
#'   step-specific aggregation policies for partial-search decisions.
#' @param partial_search_rule_min_by_step Optional named vector/list of
#'   step-specific minimum pass counts for partial-search decisions.
#' @param partial_search_candidate_source Candidate release-term source.
#'   \code{"score"} uses only score-test releasable constraints,
#'   \code{"all"} uses all releasable equality terms detected for the step,
#'   and \code{"auto"} uses \code{"all"} for
#'   \code{partial_search_exhaustive_steps} and \code{"score"} otherwise.
#' @param partial_search_exhaustive_steps Steps that should default to exhaustive
#'   release-term enumeration when \code{partial_search_candidate_source = "auto"}.
#'   Defaults to \code{c("lv.variances", "lv.covariances",
#'   "residual.covariances", "regressions", "means")}.
#' @param partial_search_max_models Maximum number of candidate models evaluated
#'   during automatic partial search for a step.
#' @param partial_search_parallel Logical; if \code{TRUE}, evaluate partial
#'   candidates in parallel where supported.
#' @param partial_search_n_cores Number of worker cores used when
#'   \code{partial_search_parallel = TRUE}.
#' @param partial_search_allow_full_release Logical; if \code{TRUE}, allow
#'   candidate models that free all step-specific equality constraints. These
#'   candidates are flagged as \code{stage_reached = FALSE}. Whether they are
#'   treated as acceptable fallbacks is controlled by
#'   \code{partial_search_full_release_action}. For
#'   \code{"lv.variances"}, \code{"lv.covariances"},
#'   \code{"residual.covariances"}, \code{"regressions"}, and
#'   \code{"means"} with exactly one releasable term, the fully freed
#'   candidate is evaluated automatically as an exploratory comparison even
#'   when this argument is \code{FALSE}.
#' @param partial_search_full_release_action How fully freed candidates are
#'   treated when included. \code{"exploratory"} keeps them non-acceptable
#'   (stage-not-reached), while \code{"eligible"} allows them to be selected as
#'   acceptable fallbacks when decision rules pass.
#' @param means_constrain_lv_variances Logical; if \code{TRUE} (default), the
#'   \code{"means"} stage also constrains latent variances. If \code{FALSE},
#'   latent variances are left unconstrained while testing mean invariance.
#' @param stop_at_first_unacceptable Logical; if \code{TRUE}, stop fitting
#'   higher invariance stages after the first constrained stage that remains
#'   unacceptable relative to the previous fitted stage.
#' @param not_applicable_action Action when a requested stage has no releasable
#'   constraints (\code{"skip"}, \code{"warn"}, or \code{"error"}).
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
  include_steps = c("configural", "metric", "scalar", "strict", "lv.variances", "lv.covariances", "residual.covariances", "means"),
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
  partial_failure_rules = NULL,
  partial_failure_rule_policy = c("all", "majority", "any", "at_least"),
  partial_failure_rule_min = NULL,
  partial_failure_rules_by_step = NULL,
  partial_failure_rule_policy_by_step = NULL,
  partial_failure_rule_min_by_step = NULL,
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
  partial_search_rules = NULL,
  partial_search_rule_policy = c("all", "majority", "any", "at_least"),
  partial_search_rule_min = NULL,
  partial_search_rules_by_step = NULL,
  partial_search_rule_policy_by_step = NULL,
  partial_search_rule_min_by_step = NULL,
  partial_search_candidate_source = c("auto", "score", "all"),
  partial_search_exhaustive_steps = c("lv.variances", "lv.covariances", "residual.covariances", "regressions", "means"),
  partial_search_max_models = 5000L,
  partial_search_parallel = FALSE,
  partial_search_n_cores = 2L,
  partial_search_allow_full_release = FALSE,
  partial_search_full_release_action = c("exploratory", "eligible"),
  means_constrain_lv_variances = TRUE,
  stop_at_first_unacceptable = TRUE,
  not_applicable_action = c("skip", "warn", "error")
) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package `lavaan` is required for MGCFA.", call. = FALSE)
  }

  model_type <- match.arg(model_type)
  available_steps <- c("configural", "metric", "scalar", "strict", "lv.variances", "lv.covariances", "residual.covariances", "regressions", "means")
  include_steps <- unique(match.arg(include_steps, choices = available_steps, several.ok = TRUE))
  include_steps <- available_steps[available_steps %in% include_steps]
  partial_failure_criterion <- match.arg(
    partial_failure_criterion,
    choices = c("chisq_pvalue", "delta_cfi", "aic_bic_weight", "measure_change", "none")
  )
  partial_failure_direction <- match.arg(partial_failure_direction)
  partial_failure_rule_policy <- match.arg(partial_failure_rule_policy, choices = c("all", "majority", "any", "at_least"))
  partial_auto_search <- match.arg(partial_auto_search, choices = c("prompt", "never", "always"))
  partial_search_rank <- match.arg(partial_search_rank)
  partial_search_rule_policy <- match.arg(partial_search_rule_policy, choices = c("all", "majority", "any", "at_least"))
  partial_search_candidate_source <- match.arg(partial_search_candidate_source, choices = c("auto", "score", "all"))
  partial_search_full_release_action <- match.arg(partial_search_full_release_action, choices = c("exploratory", "eligible"))
  not_applicable_action <- match.arg(not_applicable_action, choices = c("skip", "warn", "error"))

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
  if (!is.logical(partial_search_allow_full_release) ||
      length(partial_search_allow_full_release) != 1L ||
      is.na(partial_search_allow_full_release)) {
    stop("`partial_search_allow_full_release` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(partial_failure_rule_min)) {
    partial_failure_rule_min <- as.integer(partial_failure_rule_min)
    if (is.na(partial_failure_rule_min) || partial_failure_rule_min < 1L) {
      stop("`partial_failure_rule_min` must be NULL or a positive integer.", call. = FALSE)
    }
  }
  partial_failure_rules_by_step <- .mgcfa_as_named_list(partial_failure_rules_by_step, "partial_failure_rules_by_step")
  partial_failure_rule_policy_by_step <- .mgcfa_as_named_list(partial_failure_rule_policy_by_step, "partial_failure_rule_policy_by_step")
  partial_failure_rule_min_by_step <- .mgcfa_as_named_list(partial_failure_rule_min_by_step, "partial_failure_rule_min_by_step")
  if (!is.null(partial_search_rule_min)) {
    partial_search_rule_min <- as.integer(partial_search_rule_min)
    if (is.na(partial_search_rule_min) || partial_search_rule_min < 1L) {
      stop("`partial_search_rule_min` must be NULL or a positive integer.", call. = FALSE)
    }
  }
  partial_search_rules_by_step <- .mgcfa_as_named_list(partial_search_rules_by_step, "partial_search_rules_by_step")
  partial_search_rule_policy_by_step <- .mgcfa_as_named_list(partial_search_rule_policy_by_step, "partial_search_rule_policy_by_step")
  partial_search_rule_min_by_step <- .mgcfa_as_named_list(partial_search_rule_min_by_step, "partial_search_rule_min_by_step")
  if (!is.logical(partial_search_parallel) || length(partial_search_parallel) != 1L || is.na(partial_search_parallel)) {
    stop("`partial_search_parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  partial_search_n_cores <- as.integer(partial_search_n_cores)
  if (is.na(partial_search_n_cores) || partial_search_n_cores < 1L) {
    stop("`partial_search_n_cores` must be a positive integer.", call. = FALSE)
  }
  if (!is.logical(means_constrain_lv_variances) ||
      length(means_constrain_lv_variances) != 1L ||
      is.na(means_constrain_lv_variances)) {
    stop("`means_constrain_lv_variances` must be TRUE or FALSE.", call. = FALSE)
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
      no_mean_steps <- c("scalar", "strict", "lv.variances", "lv.covariances", "residual.covariances", "regressions", "means")
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
  post_strict_steps <- c("lv.variances", "lv.covariances", "residual.covariances", "regressions", "means")
  step_geq_used <- list()

  failure_eval_enabled <- !is.null(partial_failure_rules) || partial_failure_criterion != "none" ||
    !is.null(partial_failure_rules_by_step)
  default_failure_rule_set <- .mgcfa_resolve_rule_set(
    rules = partial_failure_rules,
    criterion = if (partial_failure_criterion == "none") "chisq_pvalue" else partial_failure_criterion,
    threshold = partial_failure_threshold,
    measure = partial_failure_measure,
    direction = partial_failure_direction,
    ic_bic_weight = partial_ic_bic_weight
  )
  search_has_single_overrides <- !is.null(partial_search_criterion) ||
    !is.null(partial_search_threshold) ||
    !is.null(partial_search_measure) ||
    !is.null(partial_search_direction)
  default_search_rules <- if (!is.null(partial_search_rules)) {
    .mgcfa_resolve_rule_set(
      rules = partial_search_rules,
      criterion = partial_search_criterion %||% "chisq_pvalue",
      threshold = partial_search_threshold,
      measure = partial_search_measure %||% partial_failure_measure,
      direction = partial_search_direction %||% partial_failure_direction,
      ic_bic_weight = partial_ic_bic_weight
    )
  } else if (!is.null(partial_failure_rules) && !isTRUE(search_has_single_overrides)) {
    default_failure_rule_set
  } else {
    .mgcfa_resolve_rule_set(
      rules = NULL,
      criterion = partial_search_criterion %||% "chisq_pvalue",
      threshold = partial_search_threshold,
      measure = partial_search_measure %||% partial_failure_measure,
      direction = partial_search_direction %||% partial_failure_direction,
      ic_bic_weight = partial_ic_bic_weight
    )
  }
  fits <- list()
  step_failures <- list()
  not_applicable_steps <- list()
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
    step_failure_cfg <- .mgcfa_step_rule_bundle(
      step = step,
      default_rules = default_failure_rule_set,
      default_policy = partial_failure_rule_policy,
      default_min = partial_failure_rule_min,
      rules_by_step = partial_failure_rules_by_step,
      policy_by_step = partial_failure_rule_policy_by_step,
      min_by_step = partial_failure_rule_min_by_step
    )
    failure_rule_set <- step_failure_cfg$rules
    failure_rule_policy_step <- step_failure_cfg$policy
    failure_rule_min_step <- step_failure_cfg$min_pass
    failure_criterion_label <- if (length(failure_rule_set) > 1L) "multi_rule" else failure_rule_set[[1L]]$criterion
    failure_threshold_display <- if (length(failure_rule_set) > 1L) {
      as.numeric(.mgcfa_required_passes(length(failure_rule_set), policy = failure_rule_policy_step, min_pass = failure_rule_min_step) / length(failure_rule_set))
    } else {
      as.numeric(failure_rule_set[[1L]]$threshold)
    }
    step_search_cfg <- .mgcfa_step_rule_bundle(
      step = step,
      default_rules = default_search_rules,
      default_policy = partial_search_rule_policy,
      default_min = partial_search_rule_min,
      rules_by_step = partial_search_rules_by_step,
      policy_by_step = partial_search_rule_policy_by_step,
      min_by_step = partial_search_rule_min_by_step
    )
    search_rules_step <- step_search_cfg$rules
    search_rule_policy_step <- step_search_cfg$policy
    search_rule_min_step <- step_search_cfg$min_pass
    search_criterion_label <- if (length(search_rules_step) > 1L) "multi_rule" else search_rules_step[[1L]]$criterion
    search_threshold_display <- if (length(search_rules_step) > 1L) {
      as.numeric(.mgcfa_required_passes(length(search_rules_step), policy = search_rule_policy_step, min_pass = search_rule_min_step) / length(search_rules_step))
    } else {
      as.numeric(search_rules_step[[1L]]$threshold)
    }
    eval_prev_step <- prev_step
    step_anchor <- "strict"
    if (step %in% post_strict_steps) {
      strict_ok <- .mgcfa_step_is_acceptable(
        step = "strict",
        fits = fits,
        step_failures = step_failures,
        recovered_steps = recovered_steps
      )
      scalar_ok <- .mgcfa_step_is_acceptable(
        step = "scalar",
        fits = fits,
        step_failures = step_failures,
        recovered_steps = recovered_steps
      )
      if (isTRUE(strict_ok)) {
        eval_prev_step <- "strict"
        step_anchor <- "strict"
      } else if (isTRUE(scalar_ok)) {
        eval_prev_step <- "scalar"
        step_anchor <- "scalar"
      } else {
        fail_info <- list(
          failed = TRUE,
          reason = sprintf(
            "Step `%s` requires an acceptable `%s` or `%s` baseline model.",
            step, "strict", "scalar"
          ),
          criterion = NA_character_,
          threshold = NA_real_,
          value = NA_real_,
          gap = NA_real_,
          from_step = NA_character_,
          to_step = step,
          failed_fit = NULL,
          failed_fit_measures = NULL,
          details = list(
            anchor_required = c("strict", "scalar")
          )
        )
        step_failures[[step]] <- fail_info
        if (is.null(partial_failure)) {
          partial_failure <- fail_info
        }
        if (isTRUE(stop_at_first_unacceptable)) {
          stopped_early <- TRUE
          stopped_after_step <- step
          stopped_reason <- fail_info$reason
          break
        }
        next
      }
    }

    prev_fit <- if (!is.null(eval_prev_step)) fits[[eval_prev_step]] else NULL
    prior_requested <- if (i > 1L) include_steps[seq_len(i - 1L)] else character()
    include_regressions_for_means <- identical(step, "means") &&
      ("regressions" %in% prior_requested) &&
      ("regressions" %in% names(fits))
    geq <- .mgcfa_step_group_equal(
      step = step,
      anchor = step_anchor,
      means_constrain_lv_variances = isTRUE(means_constrain_lv_variances),
      include_regressions_for_means = include_regressions_for_means
    )
    prev_geq <- if (!is.null(prev_step)) step_geq_used[[prev_step]] %||% character() else character()
    added_constraints <- setdiff(geq %||% character(), prev_geq %||% character())
    step_geq_used[[step]] <- geq %||% character()

    if (!identical(step, "configural") &&
        !is.null(prev_fit) &&
        length(added_constraints) > 0L &&
        !isTRUE(.mgcfa_stage_has_releasable_terms(prev_fit, added_constraints))) {
      fits[[step]] <- prev_fit
      fail_info <- list(
        failed = FALSE,
        not_applicable = TRUE,
        reason = sprintf(
          "Step `%s` is not applicable for this model (no releasable constraints in: %s).",
          step,
          paste(added_constraints, collapse = ", ")
        ),
        criterion = NA_character_,
        threshold = NA_real_,
        value = NA_real_,
        gap = NA_real_,
        from_step = eval_prev_step %||% NA_character_,
        to_step = step,
        failed_fit = NULL,
        failed_fit_measures = NULL,
        details = list(
          added_constraints = added_constraints
        )
      )
      step_failures[[step]] <- fail_info
      not_applicable_steps[[step]] <- fail_info
      if (identical(not_applicable_action, "warn")) {
        warning(fail_info$reason, call. = FALSE)
      }
      if (identical(not_applicable_action, "error")) {
        stop(fail_info$reason, call. = FALSE)
      }
      next
    }

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
      can_evaluate <- !is.null(prev_fit) && !is.null(geq) && isTRUE(failure_eval_enabled)
      if (can_evaluate) {
        fail_info <- list(
          failed = TRUE,
          reason = fit_or_error$message,
          criterion = failure_criterion_label,
          threshold = failure_threshold_display,
          value = NA_real_,
          gap = NA_real_,
          from_step = eval_prev_step,
          to_step = step,
          failed_fit = NULL,
          failed_fit_measures = NULL,
          details = list(
            policy = failure_rule_policy_step,
            rules = failure_rule_set
          )
        )
        step_failures[[step]] <- fail_info

        do_search <- .mgcfa_should_run_partial_search(
          mode = partial_auto_search,
          reason = sprintf("Step `%s` estimation failed (%s).", step, fit_or_error$message)
        )

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
            rule_set = search_rules_step,
            rule_policy = search_rule_policy_step,
            rule_min = search_rule_min_step,
            max_free = partial_search_max_free,
            top_n = partial_search_top_n,
            stop_on_accept = step_stop_on_accept,
            rank = partial_search_rank,
            use_best_if_no_pass = isTRUE(partial_search_use_best_if_no_pass),
            candidate_source = step_candidate_source,
            max_models = partial_search_max_models,
            use_parallel = isTRUE(partial_search_parallel),
            n_cores = partial_search_n_cores,
            allow_full_release = isTRUE(partial_search_allow_full_release),
            full_release_action = partial_search_full_release_action,
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
            criterion = search_criterion_label,
            threshold = search_threshold_display,
            rule_policy = search_rule_policy_step,
            rules = search_rules_step,
            candidates = data.frame(),
            top_models = data.frame(),
            selected_partial = character(),
            selected_index = NA_integer_,
            selected_fit = NULL,
            selected_added_terms = character(),
            selected_stage_reached = NA
          )
          partial_searches[[step]] <- search_out
        }

        if (is.null(partial_failure)) {
          partial_failure <- fail_info
        }
        if (is.null(partial_search)) {
          partial_search <- partial_searches[[step]]
        }
      } else {
        fail_info <- list(
          failed = TRUE,
          reason = fit_or_error$message,
          criterion = if (isTRUE(failure_eval_enabled)) failure_criterion_label else "none",
          threshold = if (isTRUE(failure_eval_enabled)) failure_threshold_display else NA_real_,
          value = NA_real_,
          gap = NA_real_,
          from_step = eval_prev_step %||% NA_character_,
          to_step = step,
          failed_fit = NULL,
          failed_fit_measures = NULL,
          details = NULL
        )
        step_failures[[step]] <- fail_info
        if (is.null(partial_failure)) {
          partial_failure <- fail_info
        }
      }
      remaining_steps <- if (i < length(include_steps)) include_steps[(i + 1L):length(include_steps)] else character()
      keep_testing_post_strict <- identical(step, "strict") &&
        any(remaining_steps %in% post_strict_steps) &&
        .mgcfa_step_is_acceptable(
          step = "scalar",
          fits = fits,
          step_failures = step_failures,
          recovered_steps = recovered_steps
        )
      if (isTRUE(keep_testing_post_strict)) {
        next
      }

      stop(sprintf("MGCFA failed at step `%s`: %s", step, fit_or_error$message), call. = FALSE)
    }

    fits[[step]] <- fit_or_error

    if (!is.null(prev_fit) && !is.null(geq) && isTRUE(failure_eval_enabled)) {
      step_eval <- .mgcfa_step_eval_rules(
        previous_fit = prev_fit,
        candidate_fit = fit_or_error,
        rules = failure_rule_set,
        policy = failure_rule_policy_step,
        min_pass = failure_rule_min_step
      )

      fail_info <- list(
        failed = !isTRUE(step_eval$pass),
        reason = if (isTRUE(step_eval$pass)) NULL else sprintf("Step `%s` failed `%s` criterion.", step, failure_criterion_label),
        criterion = failure_criterion_label,
        threshold = step_eval$threshold,
        value = step_eval$value,
        gap = step_eval$gap,
        from_step = eval_prev_step,
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
          reason = sprintf("Step `%s` failed `%s` criterion.", step, failure_criterion_label)
        )

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
            rule_set = search_rules_step,
            rule_policy = search_rule_policy_step,
            rule_min = search_rule_min_step,
            max_free = partial_search_max_free,
            top_n = partial_search_top_n,
            stop_on_accept = step_stop_on_accept,
            rank = partial_search_rank,
            use_best_if_no_pass = isTRUE(partial_search_use_best_if_no_pass),
            candidate_source = step_candidate_source,
            max_models = partial_search_max_models,
            use_parallel = isTRUE(partial_search_parallel),
            n_cores = partial_search_n_cores,
            allow_full_release = isTRUE(partial_search_allow_full_release),
            full_release_action = partial_search_full_release_action,
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
            criterion = search_criterion_label,
            threshold = search_threshold_display,
            rule_policy = search_rule_policy_step,
            rules = search_rules_step,
            candidates = data.frame(),
            top_models = data.frame(),
            selected_partial = character(),
            selected_index = NA_integer_,
            selected_fit = NULL,
            selected_added_terms = character(),
            selected_stage_reached = NA
          )
          partial_searches[[step]] <- search_out
          if (is.null(partial_search)) {
            partial_search <- search_out
          }
        }

        remaining_steps <- if (i < length(include_steps)) include_steps[(i + 1L):length(include_steps)] else character()
        keep_testing_means <- identical(step, "lv.variances") &&
          !isTRUE(means_constrain_lv_variances) &&
          ("means" %in% remaining_steps)
        keep_testing_post_strict <- identical(step, "strict") &&
          any(remaining_steps %in% post_strict_steps) &&
          .mgcfa_step_is_acceptable(
            step = "scalar",
            fits = fits,
            step_failures = step_failures,
            recovered_steps = recovered_steps
          )

        if (isTRUE(stop_at_first_unacceptable) &&
            !isTRUE(step_acceptable) &&
            !isTRUE(keep_testing_means) &&
            !isTRUE(keep_testing_post_strict)) {
          stopped_early <- TRUE
          stopped_after_step <- step
          stopped_reason <- sprintf(
            "Stopped after `%s` because invariance remained unacceptable under `%s` (threshold = %s).",
            step,
            failure_criterion_label,
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
  practical_change_table <- .mgcfa_practical_change_table(fits)
  decision_trace <- .mgcfa_build_decision_trace(
    fits = fits,
    step_failures = step_failures,
    partial_searches = partial_searches,
    recovered_steps = recovered_steps,
    not_applicable_steps = not_applicable_steps
  )
  metadata <- .mgcfa_run_metadata(
    call = match.call(),
    include_steps = include_steps,
    fit_measures = fit_measures
  )

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
    not_applicable_steps = not_applicable_steps,
    decision_trace = decision_trace,
    practical_change_table = practical_change_table,
    metadata = metadata,
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

#' Compute Bias Effect Sizes from Partial Invariance Recovery
#'
#' Produces signed bias-effect contrasts (dMACS-style deltas) by comparing the
#' failed non-partial model at each recovered step to the selected bias-adjusted
#' model. Supports observed item/composite outcomes and latent mean/SD
#' comparisons, with optional decomposition by freed parameter class. Also
#' returns item/composite bias effect-size families used in JLLJ RPubs and the
#' cited literature, including \code{dMACS}, \code{dMACS_Signed},
#' \code{dMACS_True}, \code{SDI2}, \code{UDI2}, and \code{SUDI2}.
#'
#' @param x An object returned by \code{mgcfa_auto()}.
#' @param steps Optional character vector of recovered steps to summarize.
#'   Defaults to all steps with both failed non-partial and selected recovered
#'   models.
#' @param groups Optional length-2 vector indicating reference and focal groups
#'   (labels or indices). Defaults to the first two groups.
#' @param composites Optional named list defining observed composites. Each
#'   element may be (a) a character vector of observed variable names (equal
#'   weights) or (b) a named numeric vector of weights.
#' @param include_items Logical; if \code{TRUE}, include one-item composites for
#'   each observed indicator.
#' @param include_class_decomposition Logical; if \code{TRUE}, also estimates
#'   class-specific bias effects by refitting the failed model with only the
#'   freed terms from each parameter class.
#' @param digits Number of digits used for numeric formatting.
#' @param rounding Numeric rounding mode. One of \code{"signif"},
#'   \code{"round"}, or \code{"none"}.
#'
#' @return An object of class \code{"mgcfa_bias_effects"} containing summary
#'   tables for observed and latent bias effects.
#' @export
mgcfa_bias_effects <- function(
  x,
  steps = NULL,
  groups = NULL,
  composites = NULL,
  include_items = TRUE,
  include_class_decomposition = TRUE,
  digits = 3L,
  rounding = c("signif", "round", "none")
) {
  x <- .mgcfa_as_result(x)
  rounding <- match.arg(rounding)
  if (!is.logical(include_items) || length(include_items) != 1L || is.na(include_items)) {
    stop("`include_items` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(include_class_decomposition) || length(include_class_decomposition) != 1L || is.na(include_class_decomposition)) {
    stop("`include_class_decomposition` must be TRUE or FALSE.", call. = FALSE)
  }

  fit_steps <- names(x$fits)
  candidate_steps <- Filter(function(step) {
    sf <- (x$step_failures %||% list())[[step]] %||% NULL
    ps <- (x$partial_searches %||% list())[[step]] %||% NULL
    has_failed <- !is.null(sf) && !is.null(sf$failed_fit)
    has_selected <- !is.null((x$fits %||% list())[[step]])
    has_freed <- length(unique(c(
      .mgcfa_normalize_terms(ps$selected_added_terms %||% character()),
      .mgcfa_normalize_terms((x$freed_parameters %||% list())[[step]] %||% character())
    ))) > 0L
    isTRUE(has_failed && has_selected && has_freed)
  }, fit_steps)
  candidate_steps <- as.character(candidate_steps)
  if (length(candidate_steps) == 0L) {
    stop(
      "No recovered partial-invariance steps with both failed and adjusted models were found in `x`.",
      call. = FALSE
    )
  }

  if (is.null(steps)) {
    steps <- candidate_steps
  } else {
    steps <- unique(as.character(steps))
    bad <- setdiff(steps, candidate_steps)
    if (length(bad) > 0L) {
      stop(
        "Requested `steps` are not available for bias effect computation: ",
        paste(bad, collapse = ", "),
        call. = FALSE
      )
    }
  }

  observed_tables <- list()
  observed_stats_tables <- list()
  metric_tables <- list()
  latent_tables <- list()
  latent_stats_tables <- list()
  class_tables <- list()
  freed_rows <- list()

  for (step in steps) {
    sf <- x$step_failures[[step]]
    ps <- (x$partial_searches %||% list())[[step]] %||% list()
    failed_fit <- sf$failed_fit
    adjusted_fit <- x$fits[[step]]
    if (is.null(failed_fit) || is.null(adjusted_fit)) {
      next
    }

    group_pair <- .mgcfa_resolve_group_pair(failed_fit, groups = groups)
    comp_weights <- .mgcfa_prepare_composites(
      fit = failed_fit,
      composites = composites,
      include_items = isTRUE(include_items)
    )

    obs_biased <- .mgcfa_model_observed_stats(
      fit = failed_fit,
      weights = comp_weights,
      model_label = "Biased (Non-Partial)"
    )
    obs_adjusted <- .mgcfa_model_observed_stats(
      fit = adjusted_fit,
      weights = comp_weights,
      model_label = "Bias-Adjusted (Selected)"
    )
    obs_cmp <- .mgcfa_observed_bias_delta(
      biased_stats = obs_biased,
      adjusted_stats = obs_adjusted,
      ref_group = group_pair$reference,
      focal_group = group_pair$focal,
      step = step
    )
    observed_tables[[step]] <- obs_cmp
    metrics_biased <- .mgcfa_effect_size_table(
      fit = failed_fit,
      weights = comp_weights,
      ref_group = group_pair$reference,
      focal_group = group_pair$focal,
      step = step,
      model_label = "Biased (Non-Partial)"
    )
    metrics_adjusted <- .mgcfa_effect_size_table(
      fit = adjusted_fit,
      weights = comp_weights,
      ref_group = group_pair$reference,
      focal_group = group_pair$focal,
      step = step,
      model_label = "Bias-Adjusted (Selected)"
    )
    metric_tables[[step]] <- rbind(metrics_biased, metrics_adjusted)
    observed_stats_tables[[step]] <- rbind(
      transform(obs_biased, step = step, stringsAsFactors = FALSE),
      transform(obs_adjusted, step = step, stringsAsFactors = FALSE)
    )

    latent_biased <- .mgcfa_model_latent_stats(
      fit = failed_fit,
      model_label = "Biased (Non-Partial)"
    )
    latent_adjusted <- .mgcfa_model_latent_stats(
      fit = adjusted_fit,
      model_label = "Bias-Adjusted (Selected)"
    )
    latent_cmp <- .mgcfa_latent_bias_delta(
      biased_stats = latent_biased,
      adjusted_stats = latent_adjusted,
      ref_group = group_pair$reference,
      focal_group = group_pair$focal,
      step = step
    )
    latent_tables[[step]] <- latent_cmp
    latent_stats_tables[[step]] <- rbind(
      transform(latent_biased, step = step, stringsAsFactors = FALSE),
      transform(latent_adjusted, step = step, stringsAsFactors = FALSE)
    )

    selected_added <- .mgcfa_normalize_terms(
      ps$selected_added_terms %||% ((x$freed_parameters %||% list())[[step]] %||% character())
    )
    selected_partial <- .mgcfa_normalize_terms(ps$selected_partial %||% character())
    base_partial <- setdiff(selected_partial, selected_added)
    term_classes <- .mgcfa_term_classes(adjusted_fit, selected_added)

    if (length(selected_added) > 0L) {
      freed_rows[[step]] <- data.frame(
        step = step,
        term = selected_added,
        class = unname(term_classes[selected_added]),
        stringsAsFactors = FALSE
      )
    }

    if (isTRUE(include_class_decomposition) && length(selected_added) > 0L) {
      cls <- unique(unname(term_classes[selected_added]))
      cls <- cls[nzchar(cls)]
      step_class_rows <- list()
      for (cls_i in cls) {
        cls_terms <- names(term_classes)[term_classes == cls_i]
        cls_terms <- intersect(cls_terms, selected_added)
        if (length(cls_terms) == 0L) {
          next
        }
        fit_cls <- .mgcfa_refit_with_group_partial(
          base_fit = failed_fit,
          group_partial_terms = unique(c(base_partial, cls_terms))
        )
        if (is.null(fit_cls)) {
          next
        }
        obs_cls <- .mgcfa_model_observed_stats(
          fit = fit_cls,
          weights = comp_weights,
          model_label = paste0("Class Partial: ", cls_i)
        )
        cls_cmp <- .mgcfa_observed_delta_against_biased(
          biased_stats = obs_biased,
          compare_stats = obs_cls,
          ref_group = group_pair$reference,
          focal_group = group_pair$focal,
          step = step,
          class_label = cls_i
        )
        if (!is.null(cls_cmp) && nrow(cls_cmp) > 0L) {
          cls_cmp$freed_terms <- paste(cls_terms, collapse = ", ")
          step_class_rows[[length(step_class_rows) + 1L]] <- cls_cmp
        }
      }
      if (length(step_class_rows) > 0L) {
        class_tables[[step]] <- do.call(rbind, step_class_rows)
      }
    }
  }

  observed_effects <- if (length(observed_tables) > 0L) do.call(rbind, observed_tables) else data.frame()
  observed_group_stats <- if (length(observed_stats_tables) > 0L) do.call(rbind, observed_stats_tables) else data.frame()
  effect_size_metrics <- if (length(metric_tables) > 0L) do.call(rbind, metric_tables) else data.frame()
  effect_size_recovery <- .mgcfa_effect_size_recovery(effect_size_metrics)
  latent_effects <- if (length(latent_tables) > 0L) do.call(rbind, latent_tables) else data.frame()
  latent_group_stats <- if (length(latent_stats_tables) > 0L) do.call(rbind, latent_stats_tables) else data.frame()
  class_effects <- if (length(class_tables) > 0L) do.call(rbind, class_tables) else data.frame()
  freed_parameters <- if (length(freed_rows) > 0L) do.call(rbind, freed_rows) else data.frame()

  if (nrow(observed_effects) > 0L) {
    observed_effects <- .mgcfa_format_numeric_df(observed_effects, digits = digits, rounding = rounding)
  }
  if (nrow(observed_group_stats) > 0L) {
    observed_group_stats <- .mgcfa_format_numeric_df(observed_group_stats, digits = digits, rounding = rounding)
  }
  if (nrow(effect_size_metrics) > 0L) {
    effect_size_metrics <- .mgcfa_format_numeric_df(effect_size_metrics, digits = digits, rounding = rounding)
  }
  if (nrow(effect_size_recovery) > 0L) {
    effect_size_recovery <- .mgcfa_format_numeric_df(effect_size_recovery, digits = digits, rounding = rounding)
  }
  if (nrow(latent_effects) > 0L) {
    latent_effects <- .mgcfa_format_numeric_df(latent_effects, digits = digits, rounding = rounding)
  }
  if (nrow(latent_group_stats) > 0L) {
    latent_group_stats <- .mgcfa_format_numeric_df(latent_group_stats, digits = digits, rounding = rounding)
  }
  if (nrow(class_effects) > 0L) {
    class_effects <- .mgcfa_format_numeric_df(class_effects, digits = digits, rounding = rounding)
  }

  out <- list(
    steps = steps,
    groups = if (exists("group_pair")) c(group_pair$reference, group_pair$focal) else NULL,
    observed_effects = observed_effects,
    observed_group_stats = observed_group_stats,
    effect_size_metrics = effect_size_metrics,
    effect_size_recovery = effect_size_recovery,
    latent_effects = latent_effects,
    latent_group_stats = latent_group_stats,
    class_effects = class_effects,
    freed_parameters = freed_parameters
  )
  class(out) <- "mgcfa_bias_effects"
  out
}

#' Print Method for Bias Effect Summaries
#'
#' @param x An object returned by \code{mgcfa_bias_effects()}.
#' @param ... Unused.
#'
#' @return The input object invisibly.
#' @export
print.mgcfa_bias_effects <- function(x, ...) {
  if (!inherits(x, "mgcfa_bias_effects")) {
    stop("`x` must be an object returned by `mgcfa_bias_effects()`.", call. = FALSE)
  }
  cat("Steps:", paste(x$steps, collapse = ", "), "\n")
  if (!is.null(x$groups) && length(x$groups) == 2L) {
    cat("Comparison:", x$groups[[2L]], "vs", x$groups[[1L]], "\n")
  }
  if (!is.null(x$observed_effects) && nrow(x$observed_effects) > 0L) {
    cat("\nObserved bias effects (signed dMACS-style deltas)\n")
    print(utils::head(x$observed_effects, n = 12L), row.names = FALSE)
  }
  if (!is.null(x$effect_size_metrics) && nrow(x$effect_size_metrics) > 0L) {
    cat("\nBias effect-size metrics (dMACS, dMACS_Signed, SDI2/UDI2)\n")
    print(utils::head(x$effect_size_metrics, n = 12L), row.names = FALSE)
  }
  if (!is.null(x$effect_size_recovery) && nrow(x$effect_size_recovery) > 0L) {
    cat("\nEffect-size recovery (Adjusted - Biased)\n")
    print(utils::head(x$effect_size_recovery, n = 12L), row.names = FALSE)
  }
  if (!is.null(x$latent_effects) && nrow(x$latent_effects) > 0L) {
    cat("\nLatent mean/SD comparisons (biased vs adjusted)\n")
    print(utils::head(x$latent_effects, n = 12L), row.names = FALSE)
  }
  if (!is.null(x$latent_group_stats) && nrow(x$latent_group_stats) > 0L) {
    cat("\nLatent group and pooled means/SDs\n")
    print(utils::head(x$latent_group_stats, n = 12L), row.names = FALSE)
  }
  if (!is.null(x$class_effects) && nrow(x$class_effects) > 0L) {
    cat("\nClass decomposition of observed bias effects\n")
    print(utils::head(x$class_effects, n = 12L), row.names = FALSE)
  }
  invisible(x)
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
    status <- if (!is.null((x$not_applicable_steps %||% list())[[step]])) {
      "not_applicable"
    } else if (is.null(fail_rec) || !isTRUE(fail_rec$failed)) {
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

.mgcfa_resolve_group_pair <- function(fit, groups = NULL) {
  gl <- tryCatch(lavaan::lavInspect(fit, "group.label"), error = function(e) NULL)
  if (is.null(gl) || length(gl) < 2L) {
    stop("Bias effects require at least two groups.", call. = FALSE)
  }
  gl <- as.character(gl)
  if (is.null(groups)) {
    return(list(reference = gl[[1L]], focal = gl[[2L]]))
  }
  groups <- unlist(groups, use.names = FALSE)
  if (length(groups) != 2L) {
    stop("`groups` must be NULL or length 2 (reference, focal).", call. = FALSE)
  }
  if (is.numeric(groups)) {
    idx <- as.integer(groups)
    if (any(is.na(idx)) || any(idx < 1L) || any(idx > length(gl))) {
      stop("Numeric `groups` indices are out of range.", call. = FALSE)
    }
    return(list(reference = gl[[idx[[1L]]]], focal = gl[[idx[[2L]]]]))
  }
  groups <- as.character(groups)
  if (!all(groups %in% gl)) {
    stop("Character `groups` must match lavaan group labels.", call. = FALSE)
  }
  list(reference = groups[[1L]], focal = groups[[2L]])
}

.mgcfa_prepare_composites <- function(fit, composites = NULL, include_items = TRUE) {
  ov <- tryCatch(lavaan::lavNames(fit, "ov"), error = function(e) character())
  if (length(ov) == 0L) {
    stop("No observed variables were found in the fitted model.", call. = FALSE)
  }
  empty_w <- stats::setNames(rep(0, length(ov)), ov)
  out <- list()

  if (is.null(composites)) {
    out[["Observed Composite"]] <- stats::setNames(rep(1 / length(ov), length(ov)), ov)
  } else {
    if (!is.list(composites) || length(composites) == 0L) {
      stop("`composites` must be NULL or a non-empty named list.", call. = FALSE)
    }
    if (is.null(names(composites)) || any(!nzchar(names(composites)))) {
      stop("`composites` must be a named list.", call. = FALSE)
    }
    for (nm in names(composites)) {
      def <- composites[[nm]]
      w <- empty_w
      if (is.character(def)) {
        vars <- unique(as.character(def))
        if (!all(vars %in% ov)) {
          stop(sprintf("Composite `%s` includes variables not in the model.", nm), call. = FALSE)
        }
        w[vars] <- 1 / length(vars)
      } else if (is.numeric(def) && !is.null(names(def))) {
        vars <- names(def)
        if (!all(vars %in% ov)) {
          stop(sprintf("Composite `%s` includes weight names not in the model.", nm), call. = FALSE)
        }
        w[vars] <- as.numeric(def)
      } else {
        stop(sprintf(
          "Composite `%s` must be a character vector of variable names or a named numeric weight vector.",
          nm
        ), call. = FALSE)
      }
      out[[nm]] <- w
    }
  }

  if (isTRUE(include_items)) {
    for (v in ov) {
      w <- empty_w
      w[[v]] <- 1
      out[[paste0("Item: ", v)]] <- w
    }
  }
  out
}

.mgcfa_model_observed_stats <- function(fit, weights, model_label) {
  implied <- tryCatch(lavaan::lavInspect(fit, "implied"), error = function(e) NULL)
  gl <- tryCatch(lavaan::lavInspect(fit, "group.label"), error = function(e) NULL)
  nobs <- tryCatch(as.numeric(lavaan::lavInspect(fit, "nobs")), error = function(e) NULL)
  if (is.null(implied) || is.null(gl) || length(gl) < 1L) {
    return(data.frame())
  }
  if (!is.list(implied)) {
    return(data.frame())
  }
  gl <- as.character(gl)
  n_g <- length(gl)
  if (is.null(nobs) || length(nobs) != n_g || any(!is.finite(nobs))) {
    nobs <- rep(1, n_g)
  }

  out <- list()
  out_i <- 0L
  for (comp_name in names(weights)) {
    w <- as.numeric(weights[[comp_name]])
    names(w) <- names(weights[[comp_name]])
    mean_vec <- numeric(n_g)
    sd_vec <- numeric(n_g)
    for (g in seq_len(n_g)) {
      mu_g <- as.numeric(implied[[g]]$mean)
      s_g <- implied[[g]]$cov
      if (is.null(dim(s_g))) {
        s_g <- matrix(s_g, nrow = length(w), ncol = length(w))
      }
      nm <- colnames(s_g) %||% names(w)
      names(mu_g) <- names(mu_g) %||% nm
      w_use <- w[nm]
      w_use[is.na(w_use)] <- 0
      mu_use <- mu_g[nm]
      mu_use[is.na(mu_use)] <- 0
      var_g <- as.numeric(t(w_use) %*% s_g[nm, nm, drop = FALSE] %*% w_use)
      sd_g <- sqrt(max(var_g, 0))
      mean_g <- sum(w_use * mu_use)

      mean_vec[[g]] <- mean_g
      sd_vec[[g]] <- sd_g
      out_i <- out_i + 1L
      out[[out_i]] <- data.frame(
        model = model_label,
        composite = comp_name,
        group = gl[[g]],
        mean = mean_g,
        sd = sd_g,
        n = nobs[[g]],
        stringsAsFactors = FALSE
      )
    }

    n_tot <- sum(nobs)
    w_n <- nobs / n_tot
    mean_pool <- sum(w_n * mean_vec)
    var_pool <- if (n_tot > 1) {
      sum((nobs - 1) * (sd_vec^2) + nobs * (mean_vec - mean_pool)^2) / (n_tot - 1)
    } else {
      sum(w_n * sd_vec^2)
    }
    out_i <- out_i + 1L
    out[[out_i]] <- data.frame(
      model = model_label,
      composite = comp_name,
      group = "Pooled",
      mean = mean_pool,
      sd = sqrt(max(var_pool, 0)),
      n = n_tot,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

.mgcfa_effect_size_table <- function(
  fit,
  weights,
  ref_group,
  focal_group,
  step,
  model_label
) {
  est <- tryCatch(lavaan::lavInspect(fit, "est"), error = function(e) NULL)
  gl <- tryCatch(as.character(lavaan::lavInspect(fit, "group.label")), error = function(e) NULL)
  nobs <- tryCatch(as.numeric(lavaan::lavInspect(fit, "nobs")), error = function(e) NULL)
  if (is.null(est) || is.null(gl) || !is.list(est) || length(gl) < 2L) {
    return(data.frame())
  }

  g_ref <- match(ref_group, gl)
  g_foc <- match(focal_group, gl)
  if (is.na(g_ref) || is.na(g_foc)) {
    return(data.frame())
  }
  est_ref <- est[[g_ref]]
  est_foc <- est[[g_foc]]
  if (is.null(est_ref$lambda) || is.null(est_foc$lambda) ||
      is.null(est_ref$nu) || is.null(est_foc$nu)) {
    return(data.frame())
  }

  lambda_ref <- as.matrix(est_ref$lambda)
  lambda_foc <- as.matrix(est_foc$lambda)
  nu_ref <- as.numeric(est_ref$nu)
  nu_foc <- as.numeric(est_foc$nu)
  ov_names <- rownames(lambda_ref) %||% names(nu_ref)
  if (is.null(ov_names)) {
    ov_names <- paste0("V", seq_len(nrow(lambda_ref)))
  }
  rownames(lambda_ref) <- ov_names
  rownames(lambda_foc) <- ov_names
  names(nu_ref) <- ov_names
  names(nu_foc) <- ov_names

  lv_ref <- colnames(lambda_ref) %||% paste0("LV", seq_len(ncol(lambda_ref)))
  lv_foc <- colnames(lambda_foc) %||% paste0("LV", seq_len(ncol(lambda_foc)))
  common_lv <- intersect(lv_ref, lv_foc)
  if (length(common_lv) == 0L) {
    return(data.frame())
  }
  colnames(lambda_ref) <- lv_ref
  colnames(lambda_foc) <- lv_foc
  lambda_ref <- lambda_ref[, common_lv, drop = FALSE]
  lambda_foc <- lambda_foc[, common_lv, drop = FALSE]

  alpha_raw <- est_foc$alpha %||% rep(0, ncol(lambda_foc))
  alpha_foc <- as.numeric(alpha_raw)
  names(alpha_foc) <- names(alpha_raw) %||% common_lv
  alpha_foc <- alpha_foc[common_lv]
  alpha_foc[is.na(alpha_foc)] <- 0

  psi_foc <- as.matrix(est_foc$psi %||% diag(length(common_lv)))
  rownames(psi_foc) <- rownames(psi_foc) %||% common_lv
  colnames(psi_foc) <- colnames(psi_foc) %||% common_lv
  if (!all(common_lv %in% rownames(psi_foc)) || !all(common_lv %in% colnames(psi_foc))) {
    psi_foc <- diag(length(common_lv))
    rownames(psi_foc) <- common_lv
    colnames(psi_foc) <- common_lv
  } else {
    psi_foc <- psi_foc[common_lv, common_lv, drop = FALSE]
  }

  # Item/composite SDs from implied observed moments.
  obs_stats <- .mgcfa_model_observed_stats(fit, weights = weights, model_label = model_label)
  if (nrow(obs_stats) == 0L) {
    return(data.frame())
  }

  if (is.null(nobs) || length(nobs) != length(gl) || any(!is.finite(nobs))) {
    nobs <- rep(1, length(gl))
  }
  n_ref <- nobs[[g_ref]]
  n_foc <- nobs[[g_foc]]

  out <- list()
  out_i <- 0L
  for (nm in names(weights)) {
    w <- as.numeric(weights[[nm]])
    names(w) <- names(weights[[nm]])
    w <- w[ov_names]
    w[is.na(w)] <- 0

    delta_nu <- nu_ref - nu_foc
    delta_lambda <- lambda_ref - lambda_foc
    a <- sum(w * delta_nu)
    b <- as.numeric(colSums(w * delta_lambda))

    mu_delta <- a + sum(b * alpha_foc)
    var_delta <- as.numeric(t(b) %*% psi_foc %*% b)
    if (!is.finite(var_delta) || var_delta < 0) {
      var_delta <- 0
    }
    sd_delta <- sqrt(var_delta)
    abs_expect <- .mgcfa_abs_normal_expectation(mu = mu_delta, sd = sd_delta)

    sd_ref <- obs_stats$sd[obs_stats$composite == nm & obs_stats$group == ref_group][1]
    sd_foc <- obs_stats$sd[obs_stats$composite == nm & obs_stats$group == focal_group][1]
    sd_pool <- .mgcfa_two_group_pooled_sd(
      sd1 = sd_ref,
      sd2 = sd_foc,
      n1 = n_ref,
      n2 = n_foc
    )

    dmacs <- if (is.finite(sd_pool) && sd_pool > 0) sqrt(mu_delta^2 + var_delta) / sd_pool else NA_real_
    dmacs_signed <- if (is.finite(sd_pool) && sd_pool > 0) mu_delta / sd_pool else NA_real_
    dmacs_true <- if (is.finite(dmacs_signed)) sign(dmacs_signed) * dmacs else NA_real_

    sdi2 <- if (is.finite(sd_foc) && sd_foc > 0) mu_delta / sd_foc else NA_real_
    udi2 <- if (is.finite(sd_foc) && sd_foc > 0) abs_expect / sd_foc else NA_real_
    sudi2 <- if (is.finite(sdi2)) sign(sdi2) * udi2 else NA_real_

    out_i <- out_i + 1L
    out[[out_i]] <- data.frame(
      step = step,
      model = model_label,
      comparison = paste(ref_group, "-", focal_group),
      outcome = nm,
      dmacs = as.numeric(dmacs),
      dmacs_signed = as.numeric(dmacs_signed),
      dmacs_true = as.numeric(dmacs_true),
      sdi2 = as.numeric(sdi2),
      udi2 = as.numeric(udi2),
      sudi2 = as.numeric(sudi2),
      delta_mean = as.numeric(mu_delta),
      delta_sd = as.numeric(sd_delta),
      denom_pooled_sd = as.numeric(sd_pool),
      denom_focal_sd = as.numeric(sd_foc),
      stringsAsFactors = FALSE
    )
  }

  if (out_i == 0L) {
    return(data.frame())
  }
  out_df <- do.call(rbind, out)

  agg <- data.frame(
    step = step,
    model = model_label,
    comparison = paste(ref_group, "-", focal_group),
    outcome = "Average Across Outcomes",
    dmacs = mean(out_df$dmacs, na.rm = TRUE),
    dmacs_signed = mean(out_df$dmacs_signed, na.rm = TRUE),
    dmacs_true = mean(out_df$dmacs_true, na.rm = TRUE),
    sdi2 = mean(out_df$sdi2, na.rm = TRUE),
    udi2 = mean(out_df$udi2, na.rm = TRUE),
    sudi2 = mean(out_df$sudi2, na.rm = TRUE),
    delta_mean = mean(out_df$delta_mean, na.rm = TRUE),
    delta_sd = mean(out_df$delta_sd, na.rm = TRUE),
    denom_pooled_sd = mean(out_df$denom_pooled_sd, na.rm = TRUE),
    denom_focal_sd = mean(out_df$denom_focal_sd, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  rbind(out_df, agg)
}

.mgcfa_two_group_pooled_sd <- function(sd1, sd2, n1, n2) {
  sd1 <- as.numeric(sd1)
  sd2 <- as.numeric(sd2)
  n1 <- as.numeric(n1)
  n2 <- as.numeric(n2)
  if (!is.finite(sd1) || !is.finite(sd2) || !is.finite(n1) || !is.finite(n2)) {
    return(NA_real_)
  }
  if (n1 <= 1 || n2 <= 1) {
    return(sqrt((sd1^2 + sd2^2) / 2))
  }
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

.mgcfa_abs_normal_expectation <- function(mu, sd) {
  mu <- as.numeric(mu)
  sd <- as.numeric(sd)
  if (!is.finite(mu) || !is.finite(sd) || sd < 0) {
    return(NA_real_)
  }
  if (sd == 0) {
    return(abs(mu))
  }
  z <- mu / sd
  sd * sqrt(2 / pi) * exp(-0.5 * z^2) + mu * (2 * stats::pnorm(z) - 1)
}

.mgcfa_effect_size_recovery <- function(effect_size_metrics) {
  if (is.null(effect_size_metrics) || nrow(effect_size_metrics) == 0L) {
    return(data.frame())
  }
  need_cols <- c("step", "comparison", "outcome", "model", "dmacs", "dmacs_signed", "dmacs_true", "sdi2", "udi2", "sudi2")
  if (!all(need_cols %in% names(effect_size_metrics))) {
    return(data.frame())
  }
  b <- effect_size_metrics[effect_size_metrics$model == "Biased (Non-Partial)", , drop = FALSE]
  a <- effect_size_metrics[effect_size_metrics$model == "Bias-Adjusted (Selected)", , drop = FALSE]
  if (nrow(b) == 0L || nrow(a) == 0L) {
    return(data.frame())
  }
  key <- c("step", "comparison", "outcome")
  m <- merge(
    b,
    a,
    by = key,
    suffixes = c("_biased", "_adjusted")
  )
  data.frame(
    step = m$step,
    comparison = m$comparison,
    outcome = m$outcome,
    dmacs_change = m$dmacs_adjusted - m$dmacs_biased,
    dmacs_signed_change = m$dmacs_signed_adjusted - m$dmacs_signed_biased,
    dmacs_true_change = m$dmacs_true_adjusted - m$dmacs_true_biased,
    sdi2_change = m$sdi2_adjusted - m$sdi2_biased,
    udi2_change = m$udi2_adjusted - m$udi2_biased,
    sudi2_change = m$sudi2_adjusted - m$sudi2_biased,
    stringsAsFactors = FALSE
  )
}

.mgcfa_group_d_from_stats <- function(stats_df, ref_group, focal_group) {
  ref <- stats_df[stats_df$group == ref_group, , drop = FALSE]
  focal <- stats_df[stats_df$group == focal_group, , drop = FALSE]
  if (nrow(ref) < 1L || nrow(focal) < 1L) {
    return(data.frame())
  }
  merge_df <- merge(
    ref[, c("composite", "mean", "sd"), drop = FALSE],
    focal[, c("composite", "mean", "sd"), drop = FALSE],
    by = "composite",
    suffixes = c("_ref", "_focal")
  )
  pooled_sd <- sqrt((merge_df$sd_ref^2 + merge_df$sd_focal^2) / 2)
  d <- (merge_df$mean_focal - merge_df$mean_ref) / pooled_sd
  data.frame(
    composite = merge_df$composite,
    mean_diff = merge_df$mean_focal - merge_df$mean_ref,
    d = d,
    stringsAsFactors = FALSE
  )
}

.mgcfa_observed_bias_delta <- function(
  biased_stats,
  adjusted_stats,
  ref_group,
  focal_group,
  step
) {
  d_biased <- .mgcfa_group_d_from_stats(biased_stats, ref_group = ref_group, focal_group = focal_group)
  d_adjusted <- .mgcfa_group_d_from_stats(adjusted_stats, ref_group = ref_group, focal_group = focal_group)
  if (nrow(d_biased) == 0L || nrow(d_adjusted) == 0L) {
    return(data.frame())
  }
  m <- merge(
    d_biased,
    d_adjusted,
    by = "composite",
    suffixes = c("_biased", "_adjusted")
  )
  data.frame(
    step = step,
    comparison = paste(focal_group, "vs", ref_group),
    composite = m$composite,
    d_biased = as.numeric(m$d_biased),
    d_adjusted = as.numeric(m$d_adjusted),
    dmacs_signed = as.numeric(m$d_adjusted - m$d_biased),
    mean_diff_biased = as.numeric(m$mean_diff_biased),
    mean_diff_adjusted = as.numeric(m$mean_diff_adjusted),
    stringsAsFactors = FALSE
  )
}

.mgcfa_observed_delta_against_biased <- function(
  biased_stats,
  compare_stats,
  ref_group,
  focal_group,
  step,
  class_label
) {
  d_biased <- .mgcfa_group_d_from_stats(biased_stats, ref_group = ref_group, focal_group = focal_group)
  d_cmp <- .mgcfa_group_d_from_stats(compare_stats, ref_group = ref_group, focal_group = focal_group)
  if (nrow(d_biased) == 0L || nrow(d_cmp) == 0L) {
    return(data.frame())
  }
  m <- merge(d_biased, d_cmp, by = "composite", suffixes = c("_biased", "_class"))
  data.frame(
    step = step,
    class = class_label,
    comparison = paste(focal_group, "vs", ref_group),
    composite = m$composite,
    d_biased = as.numeric(m$d_biased),
    d_class_partial = as.numeric(m$d_class),
    dmacs_signed = as.numeric(m$d_class - m$d_biased),
    stringsAsFactors = FALSE
  )
}

.mgcfa_model_latent_stats <- function(fit, model_label) {
  mean_lv <- tryCatch(lavaan::lavInspect(fit, "mean.lv"), error = function(e) NULL)
  cov_lv <- tryCatch(lavaan::lavInspect(fit, "cov.lv"), error = function(e) NULL)
  gl <- tryCatch(as.character(lavaan::lavInspect(fit, "group.label")), error = function(e) NULL)
  nobs <- tryCatch(as.numeric(lavaan::lavInspect(fit, "nobs")), error = function(e) NULL)
  if (is.null(mean_lv) || is.null(cov_lv) || is.null(gl)) {
    return(data.frame())
  }
  if (!is.list(mean_lv)) {
    mean_lv <- list(mean_lv)
  }
  if (!is.list(cov_lv)) {
    cov_lv <- list(cov_lv)
  }
  if (length(mean_lv) != length(gl) || length(cov_lv) != length(gl)) {
    return(data.frame())
  }
  if (is.null(nobs) || length(nobs) != length(gl) || any(!is.finite(nobs))) {
    nobs <- rep(1, length(gl))
  }

  out <- list()
  out_i <- 0L
  lv_names <- colnames(mean_lv[[1L]]) %||% names(mean_lv[[1L]])
  for (g in seq_along(gl)) {
    mu <- as.numeric(mean_lv[[g]])
    names(mu) <- lv_names
    sigma <- cov_lv[[g]]
    sdv <- sqrt(pmax(diag(sigma), 0))
    names(sdv) <- colnames(sigma)
    common <- intersect(names(mu), names(sdv))
    for (lv in common) {
      out_i <- out_i + 1L
      out[[out_i]] <- data.frame(
        model = model_label,
        latent = lv,
        group = gl[[g]],
        mean = mu[[lv]],
        sd = sdv[[lv]],
        n = nobs[[g]],
        stringsAsFactors = FALSE
      )
    }
  }
  if (out_i == 0L) {
    return(data.frame())
  }
  df <- do.call(rbind, out)

  pooled <- do.call(rbind, lapply(split(df, df$latent), function(dfi) {
    n_tot <- sum(dfi$n)
    w <- dfi$n / n_tot
    m <- sum(w * dfi$mean)
    v <- if (n_tot > 1) {
      sum((dfi$n - 1) * (dfi$sd^2) + dfi$n * (dfi$mean - m)^2) / (n_tot - 1)
    } else {
      sum(w * dfi$sd^2)
    }
    data.frame(
      model = unique(dfi$model)[[1L]],
      latent = unique(dfi$latent)[[1L]],
      group = "Pooled",
      mean = m,
      sd = sqrt(max(v, 0)),
      n = n_tot,
      stringsAsFactors = FALSE
    )
  }))
  rbind(df, pooled)
}

.mgcfa_latent_group_d <- function(stats_df, ref_group, focal_group) {
  ref <- stats_df[stats_df$group == ref_group, , drop = FALSE]
  focal <- stats_df[stats_df$group == focal_group, , drop = FALSE]
  if (nrow(ref) == 0L || nrow(focal) == 0L) {
    return(data.frame())
  }
  m <- merge(
    ref[, c("latent", "mean", "sd"), drop = FALSE],
    focal[, c("latent", "mean", "sd"), drop = FALSE],
    by = "latent",
    suffixes = c("_ref", "_focal")
  )
  pooled_sd <- sqrt((m$sd_ref^2 + m$sd_focal^2) / 2)
  data.frame(
    latent = m$latent,
    mean_diff = m$mean_focal - m$mean_ref,
    d = (m$mean_focal - m$mean_ref) / pooled_sd,
    stringsAsFactors = FALSE
  )
}

.mgcfa_latent_bias_delta <- function(biased_stats, adjusted_stats, ref_group, focal_group, step) {
  if (nrow(biased_stats) == 0L || nrow(adjusted_stats) == 0L) {
    return(data.frame())
  }
  d_b <- .mgcfa_latent_group_d(biased_stats, ref_group = ref_group, focal_group = focal_group)
  d_a <- .mgcfa_latent_group_d(adjusted_stats, ref_group = ref_group, focal_group = focal_group)
  if (nrow(d_b) == 0L || nrow(d_a) == 0L) {
    return(data.frame())
  }
  m <- merge(d_b, d_a, by = "latent", suffixes = c("_biased", "_adjusted"))
  data.frame(
    step = step,
    comparison = paste(focal_group, "vs", ref_group),
    latent = m$latent,
    d_biased = as.numeric(m$d_biased),
    d_adjusted = as.numeric(m$d_adjusted),
    dmacs_signed = as.numeric(m$d_adjusted - m$d_biased),
    mean_diff_biased = as.numeric(m$mean_diff_biased),
    mean_diff_adjusted = as.numeric(m$mean_diff_adjusted),
    stringsAsFactors = FALSE
  )
}

.mgcfa_term_classes <- function(fit, terms) {
  terms <- .mgcfa_normalize_terms(terms)
  if (length(terms) == 0L) {
    return(stats::setNames(character(), character()))
  }
  pt <- tryCatch(lavaan::parTable(fit), error = function(e) NULL)
  if (is.null(pt) || nrow(pt) == 0L) {
    return(stats::setNames(rep("unknown", length(terms)), terms))
  }
  ov <- tryCatch(lavaan::lavNames(fit, type = "ov"), error = function(e) character())
  lv <- tryCatch(lavaan::lavNames(fit, type = "lv"), error = function(e) character())
  rows <- vector("list", nrow(pt))
  out_i <- 0L
  for (i in seq_len(nrow(pt))) {
    term <- .mgcfa_term_from_parts(
      lhs = as.character(pt$lhs[[i]]),
      op = as.character(pt$op[[i]]),
      rhs = as.character(pt$rhs[[i]])
    )
    cls <- .mgcfa_constraint_class(
      op = as.character(pt$op[[i]]),
      lhs = as.character(pt$lhs[[i]]),
      rhs = as.character(pt$rhs[[i]]),
      ov_names = ov,
      lv_names = lv
    )
    if (is.na(cls)) {
      cls <- "unknown"
    }
    out_i <- out_i + 1L
    rows[[out_i]] <- data.frame(term = term, class = cls, stringsAsFactors = FALSE)
  }
  map <- do.call(rbind, rows[seq_len(out_i)])
  map <- map[!duplicated(map$term), , drop = FALSE]
  cls <- map$class[match(terms, map$term)]
  cls[is.na(cls) | !nzchar(cls)] <- "unknown"
  stats::setNames(cls, terms)
}

.mgcfa_refit_with_group_partial <- function(base_fit, group_partial_terms) {
  gp <- .mgcfa_normalize_terms(group_partial_terms)
  tryCatch(
    stats::update(base_fit, group.partial = gp),
    error = function(e) NULL
  )
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
  out[steps == "lv.covariances"] <- paste0("Latent", line_sep, "Covariances")
  out[steps == "residual.covariances"] <- paste0("Residual", line_sep, "Covariances")
  out[steps == "regressions"] <- "Regressions"
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

.mgcfa_as_named_list <- function(x, arg_name) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.list(x)) {
    x <- as.list(x)
  }
  nm <- names(x)
  if (is.null(nm) || any(!nzchar(nm))) {
    stop(sprintf("`%s` must be NULL or a named list/vector.", arg_name), call. = FALSE)
  }
  x
}

.mgcfa_step_rule_bundle <- function(
  step,
  default_rules,
  default_policy,
  default_min,
  rules_by_step = NULL,
  policy_by_step = NULL,
  min_by_step = NULL
) {
  step <- as.character(step)[[1L]]
  rules_step <- rules_by_step[[step]] %||% NULL
  rules <- if (is.null(rules_step)) {
    default_rules
  } else {
    .mgcfa_resolve_rule_set(rules = rules_step)
  }

  policy <- as.character(policy_by_step[[step]] %||% default_policy)[[1L]]
  policy <- match.arg(policy, choices = c("all", "majority", "any", "at_least"))
  min_step <- min_by_step[[step]] %||% default_min
  if (!is.null(min_step)) {
    min_step <- as.integer(min_step)
    if (is.na(min_step) || min_step < 1L) {
      stop(sprintf("Step-specific minimum pass count for `%s` must be a positive integer.", step), call. = FALSE)
    }
  }
  list(
    rules = rules,
    policy = policy,
    min_pass = min_step
  )
}

.mgcfa_practical_change_table <- function(fits, measures = c("cfi", "rmsea", "srmr", "aic", "bic")) {
  if (is.null(fits) || length(fits) < 2L) {
    return(data.frame())
  }
  steps <- names(fits)
  out <- vector("list", length(fits) - 1L)
  for (i in 2:length(fits)) {
    prev <- fits[[i - 1L]]
    cur <- fits[[i]]
    prev_vals <- suppressWarnings(lavaan::fitMeasures(prev, measures))
    cur_vals <- suppressWarnings(lavaan::fitMeasures(cur, measures))
    out[[i - 1L]] <- data.frame(
      from_step = steps[[i - 1L]],
      to_step = steps[[i]],
      delta_cfi = as.numeric(cur_vals[["cfi"]] - prev_vals[["cfi"]]),
      delta_rmsea = as.numeric(cur_vals[["rmsea"]] - prev_vals[["rmsea"]]),
      delta_srmr = as.numeric(cur_vals[["srmr"]] - prev_vals[["srmr"]]),
      delta_aic = as.numeric(cur_vals[["aic"]] - prev_vals[["aic"]]),
      delta_bic = as.numeric(cur_vals[["bic"]] - prev_vals[["bic"]]),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

.mgcfa_build_decision_trace <- function(
  fits,
  step_failures,
  partial_searches,
  recovered_steps,
  not_applicable_steps
) {
  if (is.null(fits) || length(fits) == 0L) {
    return(data.frame())
  }
  steps <- names(fits)
  out <- vector("list", length(steps))
  for (i in seq_along(steps)) {
    step <- steps[[i]]
    sf <- step_failures[[step]] %||% list()
    ps <- partial_searches[[step]] %||% list()
    status <- if (!is.null(not_applicable_steps[[step]])) {
      "not_applicable"
    } else if (isTRUE(sf$failed) && !(step %in% (recovered_steps %||% character()))) {
      "failed"
    } else if (step %in% (recovered_steps %||% character())) {
      "recovered_partial"
    } else {
      "ok"
    }
    out[[i]] <- data.frame(
      step = step,
      status = status,
      failed_non_partial = isTRUE(sf$failed),
      recovered = step %in% (recovered_steps %||% character()),
      search_triggered = isTRUE(ps$triggered),
      search_decision = as.character(ps$decision %||% NA_character_),
      selected_terms = paste(.mgcfa_normalize_terms(ps$selected_added_terms %||% character()), collapse = ", "),
      criterion = as.character(sf$criterion %||% NA_character_),
      criterion_value = as.numeric(sf$value %||% NA_real_),
      criterion_threshold = as.numeric(sf$threshold %||% NA_real_),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

.mgcfa_run_metadata <- function(call, include_steps, fit_measures) {
  list(
    timestamp_utc = format(as.POSIXct(Sys.time(), tz = "UTC"), tz = "UTC", usetz = TRUE),
    r_version = R.version.string,
    platform = paste(R.version$platform, R.version$arch, sep = " / "),
    include_steps = as.character(include_steps),
    fit_measures = as.character(fit_measures),
    call = paste(deparse(call), collapse = "")
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
  criterion <- match.arg(criterion, choices = c("chisq_pvalue", "delta_cfi", "aic_bic_weight", "measure_change", "multi_rule"))
  if (identical(criterion, "multi_rule")) {
    return(TRUE)
  }
  criterion %in% c("chisq_pvalue", "aic_bic_weight", "measure_change")
}

.mgcfa_required_passes <- function(n_rules, policy = c("all", "majority", "any", "at_least"), min_pass = NULL) {
  policy <- match.arg(policy)
  n_rules <- as.integer(n_rules)
  if (is.na(n_rules) || n_rules < 1L) {
    stop("`n_rules` must be a positive integer.", call. = FALSE)
  }
  if (identical(policy, "all")) {
    return(n_rules)
  }
  if (identical(policy, "any")) {
    return(1L)
  }
  if (identical(policy, "majority")) {
    return(floor(n_rules / 2) + 1L)
  }
  min_pass <- as.integer(min_pass %||% ceiling(n_rules / 2))
  if (is.na(min_pass) || min_pass < 1L) {
    min_pass <- 1L
  }
  min(min_pass, n_rules)
}

.mgcfa_resolve_rule_set <- function(
  rules = NULL,
  criterion = "chisq_pvalue",
  threshold = NULL,
  measure = "aic",
  direction = "decrease",
  ic_bic_weight = 0.5
) {
  normalize_one <- function(r) {
    if (!is.list(r)) {
      stop("Each rule in `rules` must be a named list.", call. = FALSE)
    }
    crit <- as.character(r$criterion %||% criterion)[[1L]]
    if (!nzchar(crit)) {
      stop("Each rule must define a non-empty `criterion`.", call. = FALSE)
    }
    crit <- match.arg(crit, choices = c("chisq_pvalue", "delta_cfi", "aic_bic_weight", "measure_change"))
    thr <- r$threshold %||% if (!is.null(threshold)) threshold else .mgcfa_default_partial_threshold(crit)
    if (!is.numeric(thr) || length(thr) != 1L || !is.finite(thr)) {
      stop("Each rule `threshold` must be a finite numeric scalar.", call. = FALSE)
    }

    out <- list(
      criterion = crit,
      threshold = as.numeric(thr),
      measure = NULL,
      direction = NULL,
      ic_bic_weight = NULL,
      label = as.character(r$label %||% crit)[[1L]]
    )

    if (identical(crit, "measure_change")) {
      m <- as.character(r$measure %||% measure)[[1L]]
      if (!nzchar(m)) {
        stop("`measure_change` rules require a non-empty `measure`.", call. = FALSE)
      }
      dir <- match.arg(as.character(r$direction %||% direction)[[1L]], choices = c("decrease", "increase"))
      out$measure <- m
      out$direction <- dir
    }
    if (identical(crit, "aic_bic_weight")) {
      w <- as.numeric(r$ic_bic_weight %||% ic_bic_weight)
      if (length(w) != 1L || !is.finite(w) || w < 0 || w > 1) {
        stop("`ic_bic_weight` in each rule must be in [0, 1].", call. = FALSE)
      }
      out$ic_bic_weight <- w
    }
    out
  }

  if (is.null(rules)) {
    return(list(normalize_one(list(
      criterion = criterion,
      threshold = threshold,
      measure = measure,
      direction = direction,
      ic_bic_weight = ic_bic_weight
    ))))
  }
  if (!is.list(rules) || length(rules) < 1L) {
    stop("`rules` must be NULL or a non-empty list of rule definitions.", call. = FALSE)
  }
  out <- lapply(rules, normalize_one)
  if (length(out) < 1L) {
    stop("At least one decision rule is required.", call. = FALSE)
  }
  out
}

.mgcfa_step_eval_rules <- function(previous_fit, candidate_fit, rules, policy = c("all", "majority", "any", "at_least"), min_pass = NULL) {
  policy <- match.arg(policy)
  if (!is.list(rules) || length(rules) < 1L) {
    stop("`rules` must be a non-empty list.", call. = FALSE)
  }

  evals <- lapply(rules, function(r) {
    .mgcfa_step_eval(
      previous_fit = previous_fit,
      candidate_fit = candidate_fit,
      criterion = r$criterion,
      threshold = r$threshold,
      measure = r$measure %||% "aic",
      direction = r$direction %||% "decrease",
      ic_bic_weight = r$ic_bic_weight %||% 0.5
    )
  })

  n_rules <- length(evals)
  pass_vec <- vapply(evals, function(x) isTRUE(x$pass), logical(1L))
  n_pass <- sum(pass_vec)
  required_pass <- .mgcfa_required_passes(n_rules = n_rules, policy = policy, min_pass = min_pass)
  pass_all <- n_pass >= required_pass
  single_rule <- (n_rules == 1L)

  details_tbl <- data.frame(
    rule = seq_len(n_rules),
    criterion = vapply(rules, function(r) as.character(r$criterion), character(1L)),
    pass = pass_vec,
    value = vapply(evals, function(x) as.numeric(x$value), numeric(1L)),
    threshold = vapply(evals, function(x) as.numeric(x$threshold), numeric(1L)),
    gap = vapply(evals, function(x) as.numeric(x$gap), numeric(1L)),
    measure = vapply(rules, function(r) as.character(r$measure %||% ""), character(1L)),
    direction = vapply(rules, function(r) as.character(r$direction %||% ""), character(1L)),
    stringsAsFactors = FALSE
  )

  list(
    pass = pass_all,
    value = if (single_rule) as.numeric(evals[[1L]]$value) else as.numeric(n_pass / n_rules),
    threshold = if (single_rule) as.numeric(evals[[1L]]$threshold) else as.numeric(required_pass / n_rules),
    gap = if (single_rule) as.numeric(evals[[1L]]$gap) else as.numeric(n_pass - required_pass),
    n_pass = as.integer(n_pass),
    n_rules = as.integer(n_rules),
    required_pass = as.integer(required_pass),
    details = list(
      policy = policy,
      required_pass = required_pass,
      n_pass = n_pass,
      n_rules = n_rules,
      rule_results = details_tbl,
      per_rule = evals
    )
  )
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

.mgcfa_step_is_acceptable <- function(step, fits, step_failures, recovered_steps = character()) {
  step <- as.character(step)[[1L]]
  if (!nzchar(step) || is.null(fits[[step]])) {
    return(FALSE)
  }
  fail_rec <- step_failures[[step]] %||% NULL
  if (is.null(fail_rec) || !isTRUE(fail_rec$failed)) {
    return(TRUE)
  }
  step %in% as.character(recovered_steps %||% character())
}

.mgcfa_stage_has_releasable_terms <- function(fit, target_constraints) {
  target_constraints <- unique(as.character(target_constraints %||% character()))
  if (is.null(fit) || length(target_constraints) == 0L) {
    return(FALSE)
  }

  pt <- tryCatch(
    lavaan::parTable(fit),
    error = function(e) NULL
  )
  if (is.null(pt) || nrow(pt) == 0L) {
    return(FALSE)
  }

  ov_names <- tryCatch(lavaan::lavNames(fit, type = "ov"), error = function(e) character())
  lv_names <- tryCatch(lavaan::lavNames(fit, type = "lv"), error = function(e) character())

  rows <- vector("list", nrow(pt))
  out_i <- 0L
  for (i in seq_len(nrow(pt))) {
    free_i <- suppressWarnings(as.numeric(pt$free[[i]]))
    if (!is.finite(free_i) || free_i <= 0) {
      next
    }
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
    out_i <- out_i + 1L
    rows[[out_i]] <- data.frame(
      term = .mgcfa_term_from_parts(
        lhs = as.character(pt$lhs[[i]]),
        op = as.character(pt$op[[i]]),
        rhs = as.character(pt$rhs[[i]])
      ),
      group = as.integer(pt$group[[i]]),
      stringsAsFactors = FALSE
    )
  }

  if (out_i == 0L) {
    return(FALSE)
  }
  out <- do.call(rbind, rows[seq_len(out_i)])
  out <- unique(out)
  n_group_by_term <- tapply(out$group, out$term, function(g) length(unique(g)))
  any(n_group_by_term >= 2L)
}

.mgcfa_step_group_equal <- function(
  step,
  anchor = c("strict", "scalar"),
  means_constrain_lv_variances = TRUE,
  include_regressions_for_means = FALSE
) {
  step <- as.character(step)[[1L]]
  anchor <- match.arg(anchor)

  if (identical(step, "configural")) {
    return(NULL)
  }
  if (identical(step, "metric")) {
    return(c("loadings"))
  }
  if (identical(step, "scalar")) {
    return(c("loadings", "intercepts"))
  }
  if (identical(step, "strict")) {
    return(c("loadings", "intercepts", "residuals"))
  }

  base <- if (identical(anchor, "strict")) {
    c("loadings", "intercepts", "residuals")
  } else {
    c("loadings", "intercepts")
  }

  if (identical(step, "lv.variances")) {
    return(unique(c(base, "lv.variances")))
  }
  if (identical(step, "lv.covariances")) {
    return(unique(c(base, "lv.variances", "lv.covariances")))
  }
  if (identical(step, "residual.covariances")) {
    return(unique(c(base, "lv.variances", "lv.covariances", "residual.covariances")))
  }
  if (identical(step, "regressions")) {
    return(unique(c(base, "lv.variances", "lv.covariances", "residual.covariances", "regressions")))
  }
  if (identical(step, "means")) {
    eq <- unique(c(base, "lv.variances", "lv.covariances", "residual.covariances", "means"))
    if (!isTRUE(means_constrain_lv_variances)) {
      eq <- setdiff(eq, "lv.variances")
    }
    if (isTRUE(include_regressions_for_means)) {
      eq <- unique(c(eq, "regressions"))
    }
    return(eq)
  }

  stop("Unsupported invariance step: `", step, "`.", call. = FALSE)
}

.mgcfa_nested_lrt_compare <- function(previous_fit, candidate_fit) {
  fallback <- function() {
    chisq_prev <- as.numeric(lavaan::fitMeasures(previous_fit, "chisq"))
    chisq_cand <- as.numeric(lavaan::fitMeasures(candidate_fit, "chisq"))
    df_prev <- as.numeric(lavaan::fitMeasures(previous_fit, "df"))
    df_cand <- as.numeric(lavaan::fitMeasures(candidate_fit, "df"))
    d_chisq <- chisq_cand - chisq_prev
    d_df <- df_cand - df_prev
    p_val <- if (is.finite(d_chisq) && is.finite(d_df) && d_df > 0) {
      stats::pchisq(d_chisq, df = d_df, lower.tail = FALSE)
    } else if (is.finite(d_df) && d_df <= 0) {
      1
    } else {
      NA_real_
    }
    list(
      chisq_diff = as.numeric(d_chisq),
      df_diff = as.numeric(d_df),
      p_value = as.numeric(p_val),
      method = "manual_diff"
    )
  }

  lrt <- tryCatch(
    lavaan::lavTestLRT(previous_fit, candidate_fit),
    error = function(e) NULL
  )
  if (is.null(lrt)) {
    return(fallback())
  }

  df_lrt <- tryCatch(as.data.frame(lrt), error = function(e) NULL)
  if (is.null(df_lrt) || nrow(df_lrt) < 2L) {
    return(fallback())
  }

  norm_names <- gsub("[^a-z0-9]+", "", tolower(names(df_lrt)))
  get_col <- function(candidates) {
    idx <- match(candidates, norm_names, nomatch = 0L)
    idx <- idx[idx > 0L]
    if (length(idx) < 1L) return(integer())
    idx[[1L]]
  }

  chisq_col <- get_col(c("chisqdiff", "scaledchisqdiff"))
  df_col <- get_col(c("dfdiff"))
  p_col <- get_col(c("prchisq", "pvalue"))

  row_idx <- nrow(df_lrt)
  d_chisq <- if (length(chisq_col) == 1L) as.numeric(df_lrt[[chisq_col]][row_idx]) else NA_real_
  d_df <- if (length(df_col) == 1L) as.numeric(df_lrt[[df_col]][row_idx]) else NA_real_
  p_val <- if (length(p_col) == 1L) as.numeric(df_lrt[[p_col]][row_idx]) else NA_real_

  if (!is.finite(p_val)) {
    p_val <- if (is.finite(d_chisq) && is.finite(d_df) && d_df > 0) {
      stats::pchisq(d_chisq, df = d_df, lower.tail = FALSE)
    } else if (is.finite(d_df) && d_df <= 0) {
      1
    } else {
      NA_real_
    }
  }
  if (!is.finite(d_chisq) || !is.finite(d_df)) {
    return(fallback())
  }

  list(
    chisq_diff = as.numeric(d_chisq),
    df_diff = as.numeric(d_df),
    p_value = as.numeric(p_val),
    method = "lavTestLRT"
  )
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

.mgcfa_resolve_candidate_source <- function(
  step,
  source = c("auto", "score", "all"),
  exhaustive_steps = c("lv.variances", "lv.covariances", "residual.covariances", "regressions", "means")
) {
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
    nested <- .mgcfa_nested_lrt_compare(
      previous_fit = previous_fit,
      candidate_fit = candidate_fit
    )
    d_chisq <- as.numeric(nested$chisq_diff)
    d_df <- as.numeric(nested$df_diff)
    p_value <- as.numeric(nested$p_value)
    gap <- p_value - threshold
    pass <- is.finite(gap) && gap >= 0
    return(list(
      pass = pass,
      value = p_value,
      threshold = threshold,
      gap = gap,
      details = list(
        chisq_diff = d_chisq,
        df_diff = d_df,
        method = as.character(nested$method %||% NA_character_)
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
  rule_set,
  rule_policy,
  rule_min,
  max_free,
  top_n,
  stop_on_accept,
  rank,
  use_best_if_no_pass,
  candidate_source,
  max_models,
  use_parallel,
  n_cores,
  allow_full_release,
  full_release_action,
  fit_measures
) {
  rule_policy <- match.arg(rule_policy, choices = c("all", "majority", "any", "at_least"))
  rank <- match.arg(rank, choices = c("closest", "best"))
  candidate_source <- match.arg(candidate_source, choices = c("score", "all"))
  full_release_action <- match.arg(full_release_action, choices = c("exploratory", "eligible"))
  max_models <- as.integer(max_models)
  if (is.na(max_models) || max_models < 1L) {
    max_models <- 5000L
  }
  use_parallel <- isTRUE(use_parallel)
  n_cores <- as.integer(n_cores %||% 1L)
  if (is.na(n_cores) || n_cores < 1L) {
    n_cores <- 1L
  }
  parallel_enabled <- isTRUE(use_parallel) && n_cores > 1L &&
    !identical(.Platform$OS.type, "windows") &&
    !isTRUE(stop_on_accept) &&
    requireNamespace("parallel", quietly = TRUE)
  rule_set <- .mgcfa_resolve_rule_set(rules = rule_set)
  criterion <- if (length(rule_set) > 1L) "multi_rule" else rule_set[[1L]]$criterion
  threshold <- if (length(rule_set) > 1L) {
    as.numeric(.mgcfa_required_passes(length(rule_set), policy = rule_policy, min_pass = rule_min) / length(rule_set))
  } else {
    as.numeric(rule_set[[1L]]$threshold)
  }
  allow_full_release <- isTRUE(allow_full_release)

  base_partial <- .mgcfa_normalize_terms(base_partial %||% character())
  added_constraints <- unique(as.character(added_constraints %||% character()))
  force_single_release_test <- FALSE
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
  force_single_release_test <- !isTRUE(allow_full_release) &&
    n_releasable == 1L &&
    (step %in% c("lv.variances", "lv.covariances", "residual.covariances", "regressions", "means"))
  min_remaining <- if (isTRUE(allow_full_release) || isTRUE(force_single_release_test)) 0L else 1L
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
      step_eval <- list(
        pass = FALSE,
        value = NA_real_,
        gap = NA_real_,
        threshold = threshold,
        n_pass = NA_integer_,
        n_rules = length(rule_set),
        required_pass = .mgcfa_required_passes(length(rule_set), policy = rule_policy, min_pass = rule_min),
        details = NULL
      )
      fit_stats <- stats::setNames(rep(NA_real_, 5L), c("chisq", "df", "cfi", "aic", "bic"))

      if (is_ok) {
        step_eval <- .mgcfa_step_eval_rules(
          previous_fit = previous_fit,
          candidate_fit = fit_or_error,
          rules = rule_set,
          policy = rule_policy,
          min_pass = rule_min
        )
        fit_stats <- lavaan::fitMeasures(fit_or_error, c("chisq", "df", "cfi", "aic", "bic"))
      }

      stage_reached_i <- if (n_releasable <= 0L) TRUE else (length(added_terms) < n_releasable)
      full_release_eligible <- identical(full_release_action, "eligible") && !isTRUE(stage_reached_i)
      stage_acceptable_i <- is_ok && isTRUE(step_eval$pass) &&
        (isTRUE(stage_reached_i) || isTRUE(full_release_eligible))

      candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(
        iteration = iter_id,
        step = step,
        added_constraints = paste(added_terms, collapse = ", "),
        constraint_types = paste(added_constraints, collapse = ", "),
        n_added = length(added_terms),
        n_partial_total = length(partial_now),
        criterion = criterion,
        criterion_measure = if (length(rule_set) == 1L && identical(rule_set[[1L]]$criterion, "measure_change")) rule_set[[1L]]$measure else "",
        criterion_direction = if (length(rule_set) == 1L && identical(rule_set[[1L]]$criterion, "measure_change")) rule_set[[1L]]$direction else "",
        criterion_ic_bic_weight = if (length(rule_set) == 1L && identical(rule_set[[1L]]$criterion, "aic_bic_weight")) rule_set[[1L]]$ic_bic_weight else NA_real_,
        criterion_value = as.numeric(step_eval$value),
        threshold = as.numeric(step_eval$threshold),
        criterion_gap = as.numeric(step_eval$gap),
        pass = isTRUE(step_eval$pass),
        pass_count = as.integer(step_eval$n_pass),
        rule_count = as.integer(step_eval$n_rules),
        required_pass = as.integer(step_eval$required_pass),
        stage_reached = isTRUE(stage_reached_i),
        stage_not_reached = !isTRUE(stage_reached_i),
        full_release_selected = !isTRUE(stage_reached_i),
        stage_acceptable = isTRUE(stage_acceptable_i),
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

      if (isTRUE(stage_acceptable_i)) {
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

  stage_pass <- if ("stage_acceptable" %in% names(candidates)) {
    as.logical(candidates$stage_acceptable)
  } else {
    candidates$status == "ok" & candidates$pass
  }
  stage_pass[is.na(stage_pass)] <- FALSE
  pass_idx <- which(stage_pass)
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
    rules = rule_set,
    rule_policy = rule_policy,
    rule_min = rule_min,
    candidate_source = candidate_source,
    allow_full_release = isTRUE(allow_full_release),
    full_release_action = full_release_action,
    forced_single_release_test = isTRUE(force_single_release_test),
    total_releasable = as.integer(n_releasable),
    min_active_constraints = as.integer(min_remaining),
    max_free_allowed = as.integer(max_free),
    max_models = max_models,
    parallel_requested = use_parallel,
    parallel_enabled = parallel_enabled,
    parallel_n_cores = n_cores,
    truncated = isTRUE(truncated),
    evaluated_models = as.integer(nrow(candidates)),
    criterion = criterion,
    criterion_measure = if (length(rule_set) == 1L && identical(rule_set[[1L]]$criterion, "measure_change")) rule_set[[1L]]$measure else NULL,
    criterion_direction = if (length(rule_set) == 1L && identical(rule_set[[1L]]$criterion, "measure_change")) rule_set[[1L]]$direction else NULL,
    criterion_ic_bic_weight = if (length(rule_set) == 1L && identical(rule_set[[1L]]$criterion, "aic_bic_weight")) rule_set[[1L]]$ic_bic_weight else NULL,
    threshold = threshold,
    candidates = candidates,
    top_models = top_models,
    acceptable_models = acceptable_models,
    closest_models = closest_models,
    selected_index = selected_idx,
    selected_fit = selected_fit,
    selected_partial = selected_partial,
    selected_added_terms = selected_added_terms,
    selected_stage_reached = if (!is.na(selected_idx)) isTRUE(candidates$stage_reached[[selected_idx]]) else NA,
    selected_is_acceptable = if (!is.na(selected_idx)) isTRUE(stage_pass[[selected_idx]]) else FALSE,
    top_fits = if (length(top_idx) > 0L) candidate_fits[top_idx] else list(),
    closest_fits = if (length(closest_idx) > 0L) candidate_fits[closest_idx] else list()
  )
}

.mgcfa_constraint_class <- function(op, lhs, rhs, ov_names, lv_names) {
  if (identical(op, "=~")) {
    return("loadings")
  }
  if (identical(op, "~")) {
    return("regressions")
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
  if (identical(op, "~~") && !identical(lhs, rhs)) {
    if ((lhs %in% ov_names) && (rhs %in% ov_names)) {
      return("residual.covariances")
    }
    if ((lhs %in% lv_names) && (rhs %in% lv_names)) {
      return("lv.covariances")
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
    free_i <- suppressWarnings(as.numeric(pt$free[[i]]))
    if (!is.finite(free_i) || free_i <= 0) {
      next
    }
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
  pass_ok <- if ("stage_acceptable" %in% names(candidates)) {
    as.logical(candidates$stage_acceptable[ok_idx])
  } else {
    as.logical(candidates$pass[ok_idx])
  }
  pass_ok[is.na(pass_ok)] <- FALSE
  pass_penalty <- ifelse(pass_ok, 0L, 1L)

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
  if ("stage_acceptable" %in% names(candidates)) {
    ok_pass_idx <- which(candidates$status == "ok" & as.logical(candidates$stage_acceptable))
  } else {
    ok_pass_idx <- which(candidates$status == "ok" & candidates$pass)
  }
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
    nested <- .mgcfa_nested_lrt_compare(previous_fit = prev, candidate_fit = cur)
    out[[i - 1L]] <- data.frame(
      from = step_names[i - 1L],
      to = step_names[i],
      chisq_diff = as.numeric(nested$chisq_diff),
      df_diff = as.numeric(nested$df_diff),
      p_value = as.numeric(nested$p_value),
      method = as.character(nested$method %||% NA_character_),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("step", "value_plot", "series", "variant", "series_plot", "measure_label", "variant_plot"))
}
