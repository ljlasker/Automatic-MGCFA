\name{mgcfa_auto}
\alias{mgcfa_auto}
\title{Automatic Multi-Group CFA/SEM Pipeline}
\usage{
mgcfa_auto(
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
  include_steps = c("configural", "metric", "scalar", "strict",
    "lv.variances", "means"),
  partial = NULL,
  estimator = NULL,
  std.lv_configural = TRUE,
  std.lv_constrained = FALSE,
  orthogonal = FALSE,
  fit_measures = c("chisq", "df", "npar", "cfi", "tli", "rmsea",
    "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic", "bic"),
  partial_failure_criterion = c("chisq_pvalue", "delta_cfi", "aic_bic_weight",
    "measure_change", "none"),
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
  partial_search_allow_full_release = FALSE,
  means_constrain_lv_variances = TRUE,
  stop_at_first_unacceptable = TRUE
)
}
\arguments{
\item{model_type}{One of \code{"custom"}, \code{"single_factor"},
\code{"efa"}, or \code{"pca"}.}

\item{model}{Lavaan model syntax when \code{model_type = "custom"}.}

\item{variables}{Observed indicators used in non-custom model generation.}

\item{n_factors}{Number of factors/components for \code{"efa"} or \code{"pca"}.}

\item{loading_threshold}{Absolute loading threshold for path retention in
\code{"efa"} and \code{"pca"} syntax generation.}

\item{rotation}{Rotation passed to \code{psych::fa()} or
\code{psych::principal()}.}

\item{allow_cross_loadings}{If \code{TRUE}, keep all loadings above threshold.
If \code{FALSE}, keep only one loading per indicator.}

\item{data}{Raw-data input (data.frame).}

\item{group}{Grouping column name in \code{data}.}

\item{sample_cov}{Matrix or list of group covariance/correlation matrices
for summary-statistics mode.}

\item{sample_mean}{Optional vector/list of group means in matrix mode.}

\item{sample_nobs}{Sample size input for matrix mode (scalar, vector, or list).}

\item{group_labels}{Optional explicit group labels for matrix mode.}

\item{matrices_are_cor}{Logical; treat \code{sample_cov} as correlation
matrices if \code{TRUE}.}

\item{sample_sd}{Group SD vectors used when \code{matrices_are_cor = TRUE}.}

\item{include_steps}{Invariance steps to run, chosen from
\code{"configural"}, \code{"metric"}, \code{"scalar"}, \code{"strict"},
\code{"lv.variances"}, and \code{"means"}.}

\item{partial}{Optional named list of \code{group.partial} vectors keyed by
step name.}

\item{estimator}{Optional lavaan estimator (for example \code{"ML"} or
\code{"DWLS"}).}

\item{std.lv_configural}{\code{std.lv} setting for the configural step.}

\item{std.lv_constrained}{\code{std.lv} setting for all constrained steps.}

\item{orthogonal}{Passed to lavaan.}

\item{fit_measures}{Fit statistics extracted from each step with
\code{lavaan::fitMeasures()}.}

\item{partial_failure_criterion}{Criterion used to detect failure of a
constrained invariance step relative to the previous step in sequence. One of
\code{"chisq_pvalue"}, \code{"delta_cfi"}, \code{"aic_bic_weight"},
\code{"measure_change"}, or \code{"none"}.}

\item{partial_failure_threshold}{Numeric threshold for
\code{partial_failure_criterion}. Defaults: \code{0.05} for
\code{"chisq_pvalue"}, \code{0.01} for \code{"delta_cfi"},
\code{0.5} for \code{"aic_bic_weight"}, and \code{0} for
\code{"measure_change"}.}

\item{partial_failure_measure}{Fit measure used when
\code{partial_failure_criterion = "measure_change"}.}

\item{partial_failure_direction}{Direction used when
\code{partial_failure_criterion = "measure_change"}.}

\item{partial_auto_search}{Behavior when a constrained step fails:
\code{"prompt"}, \code{"never"}, or \code{"always"}.}

\item{partial_search_criterion}{Criterion used for automatic partial-model
search. One of \code{"chisq_pvalue"}, \code{"delta_cfi"},
\code{"aic_bic_weight"}, or \code{"measure_change"}.}

\item{partial_search_threshold}{Numeric threshold for
\code{partial_search_criterion}.}

\item{partial_search_measure}{Fit measure used when
\code{partial_search_criterion = "measure_change"}. Defaults to
\code{partial_failure_measure}.}

\item{partial_search_direction}{Direction used when
\code{partial_search_criterion = "measure_change"}. Defaults to
\code{partial_failure_direction}.}

\item{partial_ic_bic_weight}{Weight on BIC when
\code{partial_failure_criterion} or \code{partial_search_criterion} is
\code{"aic_bic_weight"}. The AIC weight is
\code{1 - partial_ic_bic_weight}.}

\item{partial_search_max_free}{Maximum number of additional constraints to free
during automatic partial search at a failed step.}

\item{partial_search_top_n}{Number of closest candidate models to return.}

\item{partial_search_stop_on_accept}{If \code{TRUE}, stop searching at the first
acceptable candidate.}

\item{partial_search_rank}{Candidate ranking for output:
\code{"closest"} or \code{"best"}.}

\item{partial_search_use_best_if_no_pass}{If \code{TRUE} and no candidate passes
the threshold, continue with the best-ranked candidate.}

\item{partial_search_candidate_source}{Candidate release-term source.
\code{"score"} uses only score-test releasable constraints,
\code{"all"} uses all releasable equality terms detected for the step, and
\code{"auto"} uses \code{"all"} for
\code{partial_search_exhaustive_steps} and \code{"score"} otherwise.}

\item{partial_search_exhaustive_steps}{Steps that should default to exhaustive
release-term enumeration when
\code{partial_search_candidate_source = "auto"}. Defaults to
\code{c("lv.variances", "means")}.}

\item{partial_search_max_models}{Maximum number of candidate models evaluated
during automatic partial search for a step.}

\item{partial_search_allow_full_release}{If \code{TRUE}, allow candidate models
that free all step-specific equality constraints. These candidates are flagged
as \code{stage_reached = FALSE} and are not counted as acceptable
invariance-stage recoveries. For \code{"lv.variances"} and \code{"means"}
with exactly one releasable term, the fully freed candidate is evaluated
automatically as an exploratory comparison even when this argument is
\code{FALSE}.}

\item{means_constrain_lv_variances}{If \code{TRUE} (default), the
\code{"means"} stage also constrains latent variances. If \code{FALSE},
latent variances are left unconstrained while testing mean invariance.}

\item{stop_at_first_unacceptable}{If \code{TRUE}, stop fitting higher
invariance stages after the first constrained stage that remains
unacceptable relative to the previous fitted stage.}
}
\value{
An object of class \code{mgcfa_result} containing:
\itemize{
\item \code{model_syntax}: the fitted lavaan syntax.
\item \code{fits}: named list of lavaan fit objects by step.
\item \code{fit_table}: matrix of requested fit measures.
\item \code{lrt_table}: nested chi-square differences between adjacent steps.
\item \code{step_failures}: per-step criterion diagnostics for constrained
steps.
\item \code{failed_step_outputs}: failed non-partial model outputs for all
failed constrained steps, including failed fit measures when available.
\item \code{partial_searches}: per-step automatic partial-model discovery
results.
\item \code{freed_parameters}: per-step list of selected freed constraints.
\item \code{recovered_steps}: constrained steps where automatic discovery
recovered a selected model.
\item \code{partial_failure}: diagnostics for the first failed constrained
step encountered in sequence (if any).
\item \code{partial_search}: automatic partial-search results for that first
failed constrained step (if triggered).
\item \code{recovered_partial}: logical flag indicating whether any constrained
step was recovered by partial search.
\item \code{stopped_early}: logical flag indicating whether progression was
stopped after an unacceptable stage.
\item \code{stopped_after_step}: step name where progression stopped, when
applicable.
\item \code{stopped_reason}: brief explanation of why progression stopped.
}
}
\description{
Fits a sequential MGCFA/SEM invariance workflow from either:
\itemize{
\item Raw data (\code{data} + \code{group}), or
\item Group summary matrices (\code{sample_cov}, \code{sample_mean}, \code{sample_nobs}).
}

The model can be user-provided (\code{model_type = "custom"}) or generated
automatically from a single-factor specification, EFA, or PCA with
thresholded loading retention.

If a constrained invariance step fails, \code{mgcfa_auto()} can prompt for or
automatically run a partial-model search for that step. The search frees
eligible equality constraints iteratively, prioritizes acceptable candidates
with fewer freed constraints. By default it retains at least one active
equality constraint for a level; setting
\code{partial_search_allow_full_release = TRUE} includes fully free
exploratory candidates (flagged as stage not reached).
By default, latent-variance and latent-mean stages use exhaustive candidate
enumeration in automatic mode.
}
