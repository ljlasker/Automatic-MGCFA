make_summary_test_data <- function(n = 500L, seed = 777L) {
  set.seed(seed)
  eta1 <- rnorm(n)
  eta2 <- rnorm(n)
  lam <- c(0.80, 0.75, 0.70, 0.85)

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

test_that("raw_equivalent summary profile reproduces raw-data fits closely", {
  dat <- make_summary_test_data(n = 450L, seed = 18)
  vars <- paste0("x", 1:4)
  model <- "F =~ x1 + x2 + x3 + x4"

  sum_in <- mgcfa_make_summary(
    data = dat,
    group = "grp",
    variables = vars
  )

  expect_identical(sum_in$summary_profile, "raw_equivalent")
  expect_false(isTRUE(sum_in$matrices_are_cor))
  expect_null(sum_in$sample_sd)
  expect_true(is.list(sum_in$mgcfa_args))

  out_raw <- mgcfa_auto(
    model_type = "custom",
    model = model,
    data = dat,
    group = "grp",
    include_steps = c("configural", "metric", "scalar", "strict"),
    partial_failure_criterion = "none",
    partial_auto_search = "never"
  )

  out_sum <- do.call(
    mgcfa_auto,
    c(
      list(
        model_type = "custom",
        model = model,
        include_steps = c("configural", "metric", "scalar", "strict"),
        partial_failure_criterion = "none",
        partial_auto_search = "never"
      ),
      sum_in$mgcfa_args
    )
  )

  cmp <- mgcfa_compare_results(
    out_raw,
    out_sum,
    measures = c("cfi", "rmsea", "srmr", "aic", "bic"),
    tolerance = 1e-4
  )

  expect_true(all(cmp$within_tolerance))
  expect_true(isTRUE(attr(cmp, "overall_within_tolerance")))
})

test_that("paper profile defaults to compact correlation summaries", {
  dat <- make_summary_test_data(n = 300L, seed = 41)
  vars <- paste0("x", 1:4)

  sum_paper <- mgcfa_make_summary(
    data = dat,
    group = "grp",
    variables = vars,
    summary_profile = "paper"
  )

  expect_identical(sum_paper$summary_profile, "paper")
  expect_true(isTRUE(sum_paper$matrices_are_cor))
  expect_true(is.list(sum_paper$sample_sd))
  expect_true(all(vapply(sum_paper$sample_cov, is.matrix, logical(1L))))
})
