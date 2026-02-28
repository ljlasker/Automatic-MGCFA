make_bias_data <- function(n = 450L, seed = 2026L) {
  set.seed(seed)
  lam <- c(0.80, 0.75, 0.70, 0.85, 0.78, 0.73)
  mk <- function(eta, int_shift = rep(0, length(lam))) {
    y <- sapply(seq_along(lam), function(j) {
      int_shift[[j]] + lam[[j]] * eta + rnorm(length(eta), 0, sqrt(1 - lam[[j]]^2))
    })
    y <- as.data.frame(y)
    names(y) <- paste0("y", seq_len(ncol(y)))
    y
  }

  g1 <- mk(rnorm(n))
  g2 <- mk(rnorm(n), int_shift = c(0.65, 0.55, 0, 0, 0, 0))
  g1$grp <- "g1"
  g2$grp <- "g2"
  out <- rbind(g1, g2)
  out$grp <- factor(out$grp)
  out
}

test_that("bias effects are produced for recovered scalar stage", {
  dat <- make_bias_data(n = 500L, seed = 90)
  model <- "F =~ y1 + y2 + y3 + y4 + y5 + y6"

  out <- suppressWarnings(
    mgcfa_auto(
      model_type = "custom",
      model = model,
      data = dat,
      group = "grp",
      include_steps = c("configural", "metric", "scalar", "strict", "lv.variances", "means"),
      partial_failure_criterion = "chisq_pvalue",
      partial_auto_search = "always",
      stop_at_first_unacceptable = FALSE
    )
  )

  be <- mgcfa_bias_effects(out, include_items = FALSE)
  expect_s3_class(be, "mgcfa_bias_effects")
  expect_true(nrow(be$observed_effects) >= 1L)
  expect_true("dmacs_signed" %in% names(be$observed_effects))
  expect_true(nrow(be$effect_size_metrics) >= 1L)
  expect_true(all(c("dmacs", "dmacs_signed", "dmacs_true", "sdi2", "udi2", "sudi2") %in% names(be$effect_size_metrics)))
  expect_true(nrow(be$latent_effects) >= 1L)
  expect_true(nrow(be$latent_group_stats) >= 2L)
  expect_true(nrow(be$freed_parameters) >= 1L)

  be_ci <- suppressWarnings(
    mgcfa_bias_effects(out, include_items = FALSE, ci = TRUE, n_boot = 20L, boot_seed = 123)
  )
  expect_true(nrow(be_ci$effect_size_metrics_ci) >= 1L)
  expect_true("dmacs_ci_low" %in% names(be_ci$effect_size_metrics_ci))
  expect_true("dmacs_change_ci_low" %in% names(be_ci$effect_size_recovery_ci))
})

test_that("step-specific decision rules are accepted", {
  dat <- make_bias_data(n = 420L, seed = 902)
  model <- "F =~ y1 + y2 + y3 + y4 + y5 + y6"

  out <- suppressWarnings(
    mgcfa_auto(
      model_type = "custom",
      model = model,
      data = dat,
      group = "grp",
      include_steps = c("configural", "metric", "scalar"),
      partial_failure_criterion = "delta_cfi",
      partial_failure_threshold = 0.01,
      partial_failure_rules_by_step = list(
        scalar = list(
          list(criterion = "chisq_pvalue", threshold = 0.05),
          list(criterion = "delta_cfi", threshold = 0.01)
        )
      ),
      partial_failure_rule_policy_by_step = list(scalar = "any"),
      partial_auto_search = "never",
      stop_at_first_unacceptable = FALSE
    )
  )

  expect_true("decision_trace" %in% names(out))
  expect_true(nrow(out$decision_trace) == length(out$fits))
  expect_true("practical_change_table" %in% names(out))
  expect_true(nrow(out$practical_change_table) >= 1L)
})
