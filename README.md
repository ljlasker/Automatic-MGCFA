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
- Run staged invariance testing across `configural`, `metric`, `scalar`, `strict`, `lv.variances`, `lv.covariances`, `residual.covariances`, `regressions`, and `means` by default.
- For post-strict stages (`lv.variances`, `lv.covariances`, `residual.covariances`, `regressions`, `means`), automatically use `strict` as the baseline when acceptable, otherwise fall back to `scalar`.
- Detect failed constrained stages using chi-square change p-value or delta CFI.
- Run automatic partial-invariance search at failed stages (prompt, never, or always).
- Automatically evaluate exhaustive higher-stage subsets for `lv.variances`, `lv.covariances`, `residual.covariances`, `regressions`, and `means` (for example, general-only, general+subset, up to nearly unconstrained).
- Automatically mark stage tests as not-applicable when a constraint class is absent (for example, latent covariances in a one-factor model).
- Support step-specific decision rules and voting policies (for example, scalar uses multi-rule `any`, strict uses `all`).
- By default, stop progression after the first constrained stage that remains unacceptable.
- Preserve failed non-partial outputs for failed constrained stages in `out$failed_step_outputs`.
- Return tidy fit tables, practical step-change summaries (`out$practical_change_table`), and decision traces (`out$decision_trace`).
- Compute signed bias effect sizes (`mgcfa_bias_effects`) for observed composites/items and latent mean/SD comparisons.

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

## Easiest Input Route

Use `mgcfa_prepare_input()` to normalize whatever you have, then pass that
into `mgcfa_auto(input = ...)`.

```r
# 1) Raw data
inp_raw <- mgcfa_prepare_input(data = dat, group = "grp")
print(inp_raw)

out_raw <- mgcfa_auto(
  model_type = "custom",
  model = "g =~ x1 + x2 + x3 + x4",
  input = inp_raw
)

# 2) Correlation matrices (+ Ns; SDs auto-filled if needed)
inp_cor <- mgcfa_prepare_input(
  sample_cor = cor_list,
  sample_nobs = n_list,
  sample_mean = mean_list,
  cor_without_sd_action = "auto_unit"
)

out_cor <- mgcfa_auto(
  model_type = "custom",
  model = "g =~ x1 + x2 + x3 + x4",
  input = inp_cor,
  include_steps = c("configural", "metric")
)

# 3) Ordered indicators (threshold invariance workflow)
out_ord <- mgcfa_run_simple(
  model = "g =~ i1 + i2 + i3 + i4",
  data = dat_ord,
  group = "grp",
  ordered = c("i1", "i2", "i3", "i4")
)

# 4) Quick guidance table
mgcfa_help_inputs()
```

## Simplest End-to-End Runner

```r
out <- mgcfa_run_simple(
  model = "g =~ x1 + x2 + x3 + x4",
  data = dat,
  group = "grp",
  partial_auto_search = "always"
)

rep <- mgcfa_report(out, include_plots = TRUE)
print(rep)
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

# Enable true parallel candidate evaluation for large partial searches
out_parallel <- mgcfa_auto(
  ...,
  partial_search_parallel = TRUE,
  partial_search_n_cores = 4,
  partial_search_stop_on_accept = FALSE
)

# Step-specific rules (example: scalar uses any-of-two rules)
out_step_specific <- mgcfa_auto(
  ...,
  partial_failure_rules_by_step = list(
    scalar = list(
      list(criterion = "chisq_pvalue", threshold = 0.05),
      list(criterion = "delta_cfi", threshold = 0.01)
    )
  ),
  partial_failure_rule_policy_by_step = list(scalar = "any")
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
`Scalar`, `Strict`, `Latent Variance(s)`, `Latent Covariances`,
`Residual Covariances`, `Regressions`, and `Latent Mean(s)`.

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

## Bias Effect Sizes

```r
bias_out <- mgcfa_bias_effects(
  out,
  composites = list(
    Total = c("x1", "x2", "x3", "x4"),
    Subtest_A = c("x1", "x2")
  ),
  include_items = TRUE
)

print(bias_out)

# Signed dMACS-style bias effects on observed outcomes
head(bias_out$observed_effects)

# Published bias effect-size families:
# dMACS, dMACS_Signed, dMACS_True, SDI2, UDI2, SUDI2
head(bias_out$effect_size_metrics)

# Optional bootstrap confidence intervals for effect-size metrics
bias_out_ci <- mgcfa_bias_effects(out, ci = TRUE, n_boot = 200, conf_level = 0.95)
head(bias_out_ci$effect_size_metrics_ci)

# Per-group and pooled latent means/SDs for biased vs adjusted models
head(bias_out$latent_group_stats)
```

If you see `x must be an object returned by mgcfa_auto()`, pass the full
`mgcfa_auto()` result object (for example, `out`), not an individual
`lavaan` fit (for example, `out$fits$scalar`).

## Summary-Matrix Workflow

```r
# Build raw-equivalent summary inputs (default): optimized to reproduce
# raw-data MGCFA decisions as closely as possible.
sum_in <- mgcfa_make_summary(
  data = dat,
  group = "grp",
  variables = c("x1", "x2", "x3", "x4")
)

out_summary <- mgcfa_auto(
  model_type = "custom",
  model = "g =~ x1 + x2 + x3 + x4",
  sample_cov = sum_in$sample_cov,
  sample_sd = sum_in$sample_sd,
  matrices_are_cor = sum_in$matrices_are_cor,
  sample_mean = sum_in$sample_mean,
  sample_nobs = sum_in$sample_nobs,
  group_labels = sum_in$group_labels
)

# Equivalent shortcut:
# out_summary <- do.call(mgcfa_auto, c(list(model_type = "custom", model = "g =~ x1 + x2 + x3 + x4"), sum_in$mgcfa_args))

# Optional paper-style compact summary tables:
sum_paper <- mgcfa_make_summary(
  data = dat,
  group = "grp",
  variables = c("x1", "x2", "x3", "x4"),
  summary_profile = "paper"   # defaults to truncated 2-decimal compact output
)

# Compare raw-data and summary-matrix outputs
out_raw <- mgcfa_auto(
  model_type = "custom",
  model = "g =~ x1 + x2 + x3 + x4",
  data = dat,
  group = "grp"
)
cmp <- mgcfa_compare_results(out_raw, out_summary, tolerance = 1e-4)
head(cmp)
attr(cmp, "overall_within_tolerance")
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
