# AutomaticMGCFA

`AutomaticMGCFA` automates multi-group CFA/SEM measurement invariance workflows
using `lavaan`.

## Install from GitHub (Public)

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("ljlasker/Automatic-MGCFA")
```

## Core Features

- Fit invariance models from raw data (`data` + `group`).
- Fit invariance models from summary inputs (`sample_cov`, `sample_mean`, `sample_nobs`).
- Build models from custom lavaan syntax, a single-factor template, or thresholded EFA/PCA loadings.
- Run staged invariance testing across `configural`, `metric`, `scalar`, `strict`, `lv.variances`, `lv.covariances`, `residual.covariances`, and `means` by default (with optional `regressions` stage).
- For post-strict stages (`lv.variances`, `lv.covariances`, `residual.covariances`, `regressions`, `means`), automatically use `strict` as the baseline when acceptable, otherwise fall back to `scalar`.
- Detect failed constrained stages using chi-square change p-value or delta CFI.
- Run automatic partial-invariance search at failed stages (prompt, never, or always).
- Automatically evaluate exhaustive higher-stage subsets for `lv.variances`, `lv.covariances`, `residual.covariances`, `regressions`, and `means` (for example, general-only, general+subset, up to nearly unconstrained).
- By default, stop progression after the first constrained stage that remains unacceptable.
- Preserve failed non-partial outputs for failed constrained stages in `out$failed_step_outputs`.
- Return tidy fit tables and plotting-ready outputs.

## Quick Start

```r
library(AutomaticMGCFA)

out <- mgcfa_auto(
  model_type = "custom",
  model = "g =~ x1 + x2 + x3 + x4",
  data = dat,
  group = "grp",
  include_steps = c("configural", "metric", "scalar", "strict")
)

print(out)
print(out, show_freed_parameters = TRUE)
print(out, digits = 4, rounding = "signif")
print(out, verbose = TRUE)

# Continue testing higher levels even after an unacceptable stage:
# out <- mgcfa_auto(..., stop_at_first_unacceptable = FALSE)

# Failed non-partial outputs for failed constrained steps:
# out$failed_step_outputs
```

## Alternative Partial-Search Criteria

```r
# 60% BIC + 40% AIC weighting for partial-model selection
out_weighted <- mgcfa_auto(
  ...,
  partial_search_criterion = "aic_bic_weight",
  partial_search_threshold = 0.55,  # optional (default 0.5)
  partial_ic_bic_weight = 0.60
)

# Require AIC to improve (decrease) by at least 2 units
out_aic_drop <- mgcfa_auto(
  ...,
  partial_search_criterion = "measure_change",
  partial_search_measure = "aic",
  partial_search_direction = "decrease",
  partial_search_threshold = 2
)

# Same idea can be used for stage-failure detection too
out_failure_rule <- mgcfa_auto(
  ...,
  partial_failure_criterion = "measure_change",
  partial_failure_measure = "bic",
  partial_failure_direction = "decrease",
  partial_failure_threshold = 0
)

# Force exhaustive release-term enumeration at all stages (not only latent stages)
out_all_terms <- mgcfa_auto(
  ...,
  partial_search_candidate_source = "all",
  partial_search_stop_on_accept = FALSE,
  partial_search_max_models = 20000L
)

# Include fully free exploratory candidates (flagged as stage_not_reached)
out_with_fully_free <- mgcfa_auto(
  ...,
  partial_search_allow_full_release = TRUE
)

# Allow fully free candidates to be selected as eligible fallback models
out_with_fully_free_fallback <- mgcfa_auto(
  ...,
  partial_search_allow_full_release = TRUE,
  partial_search_full_release_action = "eligible"
)

# Note: for one-term stages (`lv.variances`, `lv.covariances`,
# `residual.covariances`, `regressions`, `means`), the fully freed
# exploratory candidate is evaluated automatically.

# Multi-metric voting (3 rules, require at least 2 to pass)
out_multi_rule <- mgcfa_auto(
  ...,
  partial_failure_rules = list(
    list(criterion = "chisq_pvalue", threshold = 0.05),
    list(criterion = "delta_cfi", threshold = 0.01),
    list(criterion = "measure_change", measure = "aic", direction = "decrease", threshold = 0)
  ),
  partial_failure_rule_policy = "at_least",
  partial_failure_rule_min = 2L,
  partial_search_rules = list(
    list(criterion = "chisq_pvalue", threshold = 0.05),
    list(criterion = "delta_cfi", threshold = 0.01),
    list(criterion = "measure_change", measure = "aic", direction = "decrease", threshold = 0)
  ),
  partial_search_rule_policy = "at_least",
  partial_search_rule_min = 2L
)

# Test mean invariance without constraining latent variances
out_means_free_lvvar <- mgcfa_auto(
  ...,
  means_constrain_lv_variances = FALSE
)
```

## Fit-Improvement Graphics

```r
# Delta AIC/BIC from configural
# solid = selected model, dashed = Non-Partial (Failed)
p1 <- mgcfa_plot_fit(
  out,
  measures = c("aic", "bic"),
  plot_type = "delta",
  include_non_partial = TRUE
)
print(p1)

# Raw CFI/RMSEA trajectories
p2 <- mgcfa_plot_fit(
  out,
  measures = c("cfi", "rmsea"),
  plot_type = "raw"
)
print(p2)
```

Default plotting labels use stage names such as `Configural`, `Metric`,
`Scalar`, `Strict`, `Latent Variance(s)`, and `Latent Mean(s)`.

## Tidy Output for Custom Reporting

```r
df_fit <- mgcfa_tidy_fit(
  out,
  measures = c("aic", "bic", "cfi"),
  add_delta = TRUE,
  digits = 3,
  rounding = "signif"
)
head(df_fit)
```

If you see `x must be an object returned by mgcfa_auto()`, pass the full
`mgcfa_auto()` result object (for example, `out`), not an individual
`lavaan` fit (for example, `out$fits$scalar`).

## Summary-Matrix Workflow

```r
# Build paper-style summary inputs (e.g., r = 0.45, SD = 1.01)
sum_in <- mgcfa_make_summary(
  data = dat,
  group = "grp",
  variables = c("x1", "x2", "x3", "x4"),
  matrix_type = "cor",
  format = "truncate",
  cor_digits = 2,
  sd_digits = 2,
  mean_digits = 2
)

out <- mgcfa_auto(
  model_type = "custom",
  model = "g =~ x1 + x2 + x3 + x4",
  sample_cov = sum_in$sample_cov,
  sample_sd = sum_in$sample_sd,
  matrices_are_cor = sum_in$matrices_are_cor,
  sample_mean = sum_in$sample_mean,
  sample_nobs = sum_in$sample_nobs,
  group_labels = sum_in$group_labels
)
```

## Automatic Model Generation

```r
syntax <- mgcfa_build_model(
  model_type = "pca",
  data = dat[c("x1", "x2", "x3", "x4", "x5")],
  variables = c("x1", "x2", "x3", "x4", "x5"),
  n_factors = 2,
  loading_threshold = 0.30
)
```

## An Open-Source Demonstration

Run the package demo script:

```r
source("inst/examples/MatrixExperiment1Demo.R")
```

This script:

- Runs a single measurement invariance test from an experiment
- Fits through testing the equivalence of latent means
- Prints a concise summary of the results
- Shows select free constraints when partial invariance is achieved
- Generates AIC and BIC delta plots
