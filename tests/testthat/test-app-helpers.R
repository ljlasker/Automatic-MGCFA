test_that("app summary unpack accepts mgcfa_make_summary outputs", {
  dat <- data.frame(
    x1 = c(rnorm(60), rnorm(60)),
    x2 = c(rnorm(60), rnorm(60)),
    x3 = c(rnorm(60), rnorm(60)),
    grp = factor(rep(c("g1", "g2"), each = 60))
  )

  sum_obj <- mgcfa_make_summary(
    data = dat,
    group = "grp",
    variables = c("x1", "x2", "x3"),
    summary_profile = "raw_equivalent"
  )

  out <- .mgcfa_app_unpack_summary_input(sum_obj)
  expect_true(is.list(out$sample_cov))
  expect_equal(length(out$sample_cov), 2L)
  expect_true(!is.null(out$sample_nobs))
})

test_that("app summary unpack errors on missing required fields", {
  expect_error(
    .mgcfa_app_unpack_summary_input(list(sample_cov = list(diag(2)))),
    "missing required fields"
  )
})

test_that("app raw reader handles csv and rds", {
  dat <- data.frame(a = 1:4, b = c(2, 3, 4, 5), grp = c("g1", "g1", "g2", "g2"))

  csv <- tempfile(fileext = ".csv")
  rds <- tempfile(fileext = ".rds")
  utils::write.csv(dat, csv, row.names = FALSE)
  saveRDS(dat, rds)

  d_csv <- .mgcfa_app_read_raw_data(csv, header = TRUE, sep = ",")
  d_rds <- .mgcfa_app_read_raw_data(rds)

  expect_true(is.data.frame(d_csv))
  expect_true(is.data.frame(d_rds))
  expect_equal(nrow(d_csv), nrow(dat))
  expect_equal(nrow(d_rds), nrow(dat))
})

test_that("app criterion parser supports AIC/BIC aliases and measure_change", {
  aic_cfg <- .mgcfa_app_parse_criterion(
    criterion = "aic_weight",
    threshold = 0.55,
    allow_none = FALSE
  )
  expect_identical(aic_cfg$criterion, "aic_bic_weight")
  expect_equal(aic_cfg$ic_bic_weight, 0)
  expect_equal(aic_cfg$threshold, 0.55)

  bic_cfg <- .mgcfa_app_parse_criterion(
    criterion = "bic_weight",
    threshold = NA_real_,
    allow_none = FALSE
  )
  expect_identical(bic_cfg$criterion, "aic_bic_weight")
  expect_equal(bic_cfg$ic_bic_weight, 1)
  expect_true(is.null(bic_cfg$threshold))

  m_cfg <- .mgcfa_app_parse_criterion(
    criterion = "measure_change",
    threshold = 0.02,
    measure = "rmsea",
    direction = "increase",
    allow_none = FALSE
  )
  expect_identical(m_cfg$criterion, "measure_change")
  expect_identical(m_cfg$measure, "rmsea")
  expect_identical(m_cfg$direction, "increase")
  expect_equal(m_cfg$threshold, 0.02)
})

test_that("app rule builder creates normalized multi-rule lists", {
  rules <- .mgcfa_app_build_rules(
    criteria = c("aic_weight", "bic_weight", "measure_change"),
    thresholds = c(0.6, 0.55, 0.01),
    measures = c("aic", "bic", "rmsea"),
    directions = c("decrease", "decrease", "decrease"),
    ic_bic_weights = c(0.5, 0.5, 0.5)
  )

  expect_equal(length(rules), 3L)
  expect_identical(rules[[1]]$criterion, "aic_bic_weight")
  expect_equal(rules[[1]]$ic_bic_weight, 0)
  expect_identical(rules[[2]]$criterion, "aic_bic_weight")
  expect_equal(rules[[2]]$ic_bic_weight, 1)
  expect_identical(rules[[3]]$criterion, "measure_change")
  expect_identical(rules[[3]]$measure, "rmsea")
  expect_identical(rules[[3]]$direction, "decrease")
})

test_that("app rule builder validates length mismatches", {
  expect_error(
    .mgcfa_app_build_rules(
      criteria = c("chisq_pvalue", "delta_cfi"),
      thresholds = c(0.05, 0.01, 0.2)
    ),
    "length 1 or match"
  )
})

test_that("step-specific override builder returns named step rule bundles", {
  vals <- list(
    failure_step_use_metric = TRUE,
    failure_step_criterion_metric = "measure_change",
    failure_step_threshold_metric = 0.01,
    failure_step_measure_metric = "rmsea",
    failure_step_direction_metric = "decrease",
    failure_step_ic_bic_weight_metric = 0.5,
    failure_step_policy_metric = "at_least",
    failure_step_min_metric = 1,
    failure_step_use_scalar = FALSE
  )
  out <- .mgcfa_app_build_step_rule_overrides(
    values = vals,
    steps = c("metric", "scalar"),
    prefix = "failure"
  )
  expect_true("metric" %in% names(out$rules_by_step))
  expect_false("scalar" %in% names(out$rules_by_step))
  expect_identical(out$rules_by_step$metric[[1]]$criterion, "measure_change")
  expect_identical(out$rules_by_step$metric[[1]]$measure, "rmsea")
  expect_identical(out$policy_by_step$metric, "at_least")
  expect_identical(unname(out$min_by_step$metric), 1L)
})

test_that("run-arg validation rejects unknown measure names when strict", {
  args <- list(
    include_steps = c("configural", "metric"),
    partial_failure_criterion = "measure_change",
    partial_failure_measure = "not_a_real_measure",
    partial_search_criterion = "chisq_pvalue",
    partial_ic_bic_weight = 0.5
  )
  expect_error(
    .mgcfa_app_validate_run_args(args, allow_nonstandard_measures = FALSE),
    "Unknown fit-measure name"
  )
  expect_true(isTRUE(.mgcfa_app_validate_run_args(args, allow_nonstandard_measures = TRUE)))
})

test_that("run-arg validation catches invalid at_least minima", {
  args <- list(
    include_steps = c("configural", "metric"),
    partial_failure_criterion = "chisq_pvalue",
    partial_search_criterion = "chisq_pvalue",
    partial_ic_bic_weight = 0.5,
    partial_failure_rules = list(
      list(criterion = "chisq_pvalue", threshold = 0.05),
      list(criterion = "delta_cfi", threshold = 0.01)
    ),
    partial_failure_rule_policy = "at_least",
    partial_failure_rule_min = 3L
  )
  expect_error(
    .mgcfa_app_validate_run_args(args, allow_nonstandard_measures = FALSE),
    "between 1 and number of rules"
  )
})
