make_easy_input_data <- function(n = 400L, seed = 222L) {
  set.seed(seed)
  eta1 <- rnorm(n)
  eta2 <- rnorm(n)
  lam <- c(0.8, 0.75, 0.7, 0.85)
  mk <- function(eta) {
    y <- sapply(seq_along(lam), function(j) lam[j] * eta + rnorm(length(eta), 0, sqrt(1 - lam[j]^2)))
    y <- as.data.frame(y)
    names(y) <- paste0("x", seq_len(ncol(y)))
    y
  }
  g1 <- mk(eta1)
  g2 <- mk(eta2)
  g1$grp <- "g1"
  g2$grp <- "g2"
  out <- rbind(g1, g2)
  out$grp <- factor(out$grp)
  out
}

test_that("mgcfa_prepare_input supports raw-data workflow via mgcfa_auto(input=...)", {
  dat <- make_easy_input_data(n = 450L, seed = 31)
  prep <- mgcfa_prepare_input(data = dat, group = "grp")
  expect_s3_class(prep, "mgcfa_input")
  expect_identical(prep$mode, "raw_data")

  out <- mgcfa_auto(
    model_type = "custom",
    model = "F =~ x1 + x2 + x3 + x4",
    input = prep,
    include_steps = c("configural", "metric", "scalar"),
    partial_failure_criterion = "none",
    partial_auto_search = "never"
  )

  expect_s3_class(out, "mgcfa_result")
  expect_identical(out$mode, "raw_data")
  expect_true(!is.null(out$input_summary))
  expect_identical(out$input_summary$mode, "raw_data")
})

test_that("summary input helper auto-fills unit SDs for correlation input when requested", {
  dat <- make_easy_input_data(n = 350L, seed = 41)
  sum_in <- mgcfa_make_summary(
    data = dat,
    group = "grp",
    variables = c("x1", "x2", "x3", "x4"),
    summary_profile = "paper"
  )

  prep <- mgcfa_prepare_input(
    sample_cor = sum_in$sample_cov,
    sample_nobs = sum_in$sample_nobs,
    sample_mean = sum_in$sample_mean,
    cor_without_sd_action = "auto_unit"
  )
  expect_s3_class(prep, "mgcfa_input")
  expect_identical(prep$mode, "summary_matrices")
  expect_true(any(grepl("unit SDs", prep$notes, fixed = TRUE)))
  expect_true(!is.null(prep$mgcfa_args$sample_sd))
})

test_that("mgcfa_auto accepts sample_cor alias", {
  dat <- make_easy_input_data(n = 380L, seed = 51)
  sum_in <- mgcfa_make_summary(
    data = dat,
    group = "grp",
    variables = c("x1", "x2", "x3", "x4"),
    summary_profile = "paper"
  )

  out <- mgcfa_auto(
    model_type = "custom",
    model = "F =~ x1 + x2 + x3 + x4",
    sample_cor = sum_in$sample_cov,
    sample_sd = sum_in$sample_sd,
    sample_mean = sum_in$sample_mean,
    sample_nobs = sum_in$sample_nobs,
    group_labels = sum_in$group_labels,
    include_steps = c("configural", "metric"),
    partial_failure_criterion = "none",
    partial_auto_search = "never"
  )

  expect_identical(out$mode, "summary_matrices")
  expect_identical(out$input_summary$matrix_type, "correlation")
})

test_that("mgcfa_help_inputs returns guidance table", {
  g <- suppressMessages(capture.output(tbl <- mgcfa_help_inputs()))
  expect_true(length(g) > 0L)
  expect_true(is.data.frame(tbl))
  expect_true(all(c("scenario", "minimal_requirements") %in% names(tbl)))
})

test_that("mgcfa_run_simple and mgcfa_report produce expected outputs", {
  dat <- make_easy_input_data(n = 300L, seed = 88)
  out <- mgcfa_run_simple(
    model = "F =~ x1 + x2 + x3 + x4",
    data = dat,
    group = "grp",
    include_steps = c("configural", "metric", "scalar"),
    partial_auto_search = "never",
    partial_failure_criterion = "none"
  )
  expect_s3_class(out, "mgcfa_result")

  rep <- mgcfa_report(out, include_plots = FALSE)
  expect_s3_class(rep, "mgcfa_report")
  expect_true(is.data.frame(rep$overview))
  expect_true(is.data.frame(rep$tidy_fit))
  expect_true(is.data.frame(rep$decision_trace))
})

test_that("ordered workflow sets categorical estimator and thresholds in scalar stage", {
  dat <- make_easy_input_data(n = 320L, seed = 99)
  ord_vars <- c("x1", "x2", "x3", "x4")
  for (v in ord_vars) {
    dat[[v]] <- ordered(as.integer(cut(dat[[v]], breaks = 5)))
  }
  out <- mgcfa_run_simple(
    model = "F =~ x1 + x2 + x3 + x4",
    data = dat,
    group = "grp",
    ordered = ord_vars,
    include_steps = c("configural", "metric", "scalar"),
    partial_auto_search = "never",
    partial_failure_criterion = "none",
    stop_at_first_unacceptable = FALSE
  )
  expect_s3_class(out, "mgcfa_result")
  est <- tryCatch(lavaan::lavInspect(out$fits$configural, "options")$estimator, error = function(e) NA_character_)
  expect_true(est %in% c("WLSMV", "DWLS"))
  expect_true(!is.null(out$step_failures$scalar) || "scalar" %in% names(out$fits))
})
