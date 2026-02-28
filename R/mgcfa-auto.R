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
#' @param input Optional object returned by \code{mgcfa_prepare_input()}.
#'   When supplied, raw/summary input arguments are taken from this object.
#' @param ordered Optional ordered/categorical indicator specification for raw
#'   data mode. Use \code{TRUE} to treat all selected indicators as ordered, or
#'   provide a character vector of ordered variable names.
#' @param parameterization Parameterization passed to lavaan when
#'   \code{ordered} is used (\code{"delta"} or \code{"theta"}).
#' @param data Raw-data input (data.frame) for multi-group SEM.
#' @param group Grouping variable name in \code{data}.
#' @param sample_cov Matrix or list of group covariance/correlation matrices.
#' @param sample_cor Optional alias for correlation-matrix input. Equivalent to
#'   supplying \code{sample_cov = sample_cor} with
#'   \code{matrices_are_cor = TRUE}.
#' @param sample_mean Optional vector/list of group means (required for scalar
#'   and stricter steps in matrix mode).
#' @param sample_nobs Scalar, vector, or list of group sample sizes in matrix
#'   mode.
#' @param group_labels Optional explicit group labels.
#' @param matrices_are_cor Logical. If \code{TRUE}, \code{sample_cov} is treated
#'   as correlation matrices and converted using \code{sample_sd}.
#' @param sample_sd Group SD vectors used when \code{matrices_are_cor = TRUE}.
#' @param summary_check_action Validation action for summary matrices:
#'   \code{"warn"}, \code{"error"}, or \code{"none"}.
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
  input = NULL,
  ordered = NULL,
  parameterization = c("delta", "theta"),
  data = NULL,
  group = NULL,
  sample_cov = NULL,
  sample_cor = NULL,
  sample_mean = NULL,
  sample_nobs = NULL,
  group_labels = NULL,
  matrices_are_cor = FALSE,
  sample_sd = NULL,
  summary_check_action = c("warn", "error", "none"),
  include_steps = c("configural", "metric", "scalar", "strict", "lv.variances", "lv.covariances", "residual.covariances", "regressions", "means"),
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

  parameterization <- match.arg(parameterization, choices = c("delta", "theta"))
  summary_check_action <- match.arg(summary_check_action, choices = c("warn", "error", "none"))
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

  if (!is.null(sample_cor)) {
    if (!is.null(sample_cov)) {
      stop("Provide `sample_cov` or `sample_cor`, not both.", call. = FALSE)
    }
    sample_cov <- sample_cor
    matrices_are_cor <- TRUE
  }

  input_summary <- NULL
  input_notes <- character()
  if (!is.null(input)) {
    has_explicit <- !is.null(data) || !is.null(group) || !is.null(sample_cov) ||
      !is.null(sample_mean) || !is.null(sample_nobs) || !is.null(group_labels) || !is.null(sample_sd)
    if (has_explicit) {
      stop(
        "When `input` is supplied, do not also pass raw/summary input arguments directly.",
        call. = FALSE
      )
    }
    input_obj <- if (inherits(input, "mgcfa_input")) {
      input
    } else if (is.list(input) && !is.null(input$mgcfa_args)) {
      input
    } else if (is.list(input)) {
      list(
        mode = if (!is.null(input$data)) "raw_data" else "summary_matrices",
        mgcfa_args = input,
        summary = NULL,
        notes = character()
      )
    } else {
      stop("`input` must be an `mgcfa_input` object or a list of input arguments.", call. = FALSE)
    }
    in_args <- input_obj$mgcfa_args %||% list()
    data <- in_args$data %||% NULL
    group <- in_args$group %||% NULL
    sample_cov <- in_args$sample_cov %||% NULL
    sample_mean <- in_args$sample_mean %||% NULL
    sample_nobs <- in_args$sample_nobs %||% NULL
    group_labels <- in_args$group_labels %||% NULL
    sample_sd <- in_args$sample_sd %||% NULL
    if (!is.null(in_args$matrices_are_cor)) {
      matrices_are_cor <- isTRUE(in_args$matrices_are_cor)
    }
    input_summary <- input_obj$summary %||% NULL
    input_notes <- as.character(input_obj$notes %||% character())
  }

  using_raw <- !is.null(data)
  using_matrices <- !is.null(sample_cov)
  if (identical(using_raw, using_matrices)) {
    stop("Provide either raw `data` (with `group`) or matrix input via `sample_cov`.", call. = FALSE)
  }
  ordered_vars <- character()

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
    ordered_vars <- .mgcfa_resolve_ordered_vars(
      ordered = ordered,
      data = data,
      group = group,
      variables = variables
    )
    if (length(ordered_vars) > 0L && is.null(estimator)) {
      estimator <- "WLSMV"
      input_notes <- unique(c(input_notes, "Set estimator to WLSMV because ordered indicators were supplied."))
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
    if (length(ordered_vars) > 0L) {
      fit_args_base$ordered <- ordered_vars
      fit_args_base$parameterization <- parameterization
    }
    if (!is.null(group_labels)) {
      fit_args_base$group.label <- as.character(group_labels)
    }
    if (is.null(input_summary)) {
      g_chr <- as.character(data[[group]])
      g_chr <- g_chr[!is.na(g_chr)]
      g_tab <- table(g_chr)
      input_summary <- list(
        mode = "raw_data",
        n_rows = nrow(data),
        n_columns = ncol(data),
        group = group,
        n_groups = length(g_tab),
        group_labels = as.character(names(g_tab)),
        nobs_by_group = as.integer(g_tab),
        ordered = length(ordered_vars) > 0L,
        ordered_variables = ordered_vars
      )
    }
  } else {
    if (!is.null(ordered)) {
      stop("`ordered` is currently supported in raw-data mode only.", call. = FALSE)
    }
    prepared <- .mgcfa_prepare_matrix_input(
      sample_cov = sample_cov,
      sample_mean = sample_mean,
      sample_nobs = sample_nobs,
      group_labels = group_labels,
      matrices_are_cor = matrices_are_cor,
      sample_sd = sample_sd
    )
    .mgcfa_validate_summary_matrices(
      sample_cov = prepared$sample_cov,
      sample_nobs = prepared$sample_nobs,
      sample_mean = prepared$sample_mean,
      action = summary_check_action
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
        input_notes <- unique(c(
          input_notes,
          paste0("Dropped mean-dependent stages in matrix mode: ", paste(dropped, collapse = ", "))
        ))
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
    if (is.null(input_summary)) {
      input_summary <- list(
        mode = "summary_matrices",
        n_groups = length(prepared$sample_cov),
        n_variables = ncol(prepared$sample_cov[[1L]]),
        matrix_type = if (isTRUE(matrices_are_cor)) "correlation" else "covariance",
        has_means = !is.null(prepared$sample_mean),
        has_sds = !is.null(sample_sd),
        group_labels = as.character(prepared$group_labels %||% paste0("g", seq_along(prepared$sample_cov))),
        nobs_by_group = as.numeric(unlist(prepared$sample_nobs))
      )
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
      use_thresholds = length(ordered_vars) > 0L,
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
    input_summary = input_summary,
    input_notes = unique(input_notes),
    stopped_early = stopped_early,
    stopped_after_step = stopped_after_step,
    stopped_reason = stopped_reason
  )
  class(out) <- "mgcfa_result"
  out
}
