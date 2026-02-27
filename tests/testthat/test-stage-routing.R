make_two_group_onefactor <- function(n = 500L, res_scale_group2 = rep(1, 6), seed = 1L) {
  set.seed(seed)
  stopifnot(length(res_scale_group2) == 6L)

  make_group <- function(res_scale) {
    eta <- rnorm(n)
    lam <- c(0.80, 0.75, 0.70, 0.85, 0.78, 0.73)
    err_sd <- sqrt((1 - lam^2) * res_scale)
    y <- sapply(seq_along(lam), function(j) lam[j] * eta + rnorm(n, 0, err_sd[j]))
    y <- as.data.frame(y)
    names(y) <- paste0("y", seq_len(ncol(y)))
    y
  }

  g1 <- make_group(rep(1, 6))
  g2 <- make_group(res_scale_group2)
  g1$grp <- "g1"
  g2$grp <- "g2"
  out <- rbind(g1, g2)
  out$grp <- factor(out$grp)
  out
}

test_that("include_steps are routed in canonical order and post-strict stages fall back to scalar", {
  dat <- make_two_group_onefactor(n = 700L, res_scale_group2 = c(1.9, 1.7, 1.6, 1.8, 1.9, 1.7), seed = 42)
  model <- "F =~ y1 + y2 + y3 + y4 + y5 + y6"

  out <- suppressWarnings(
    mgcfa_auto(
      model_type = "custom",
      model = model,
      data = dat,
      group = "grp",
      include_steps = c("lv.covariances", "configural", "scalar", "means", "metric", "lv.variances", "residual.covariances", "strict"),
      partial_failure_criterion = "chisq_pvalue",
      partial_auto_search = "never",
      stop_at_first_unacceptable = FALSE
    )
  )

  expect_identical(
    names(out$fits),
    c("configural", "metric", "scalar", "strict", "lv.variances", "lv.covariances", "residual.covariances", "means")
  )
  expect_true(isTRUE(out$step_failures$strict$failed))
  expect_identical(out$step_failures$lv.variances$from_step, "scalar")
})

test_that("no-op stages are flagged as not applicable", {
  dat <- make_two_group_onefactor(n = 400L, res_scale_group2 = rep(1, 6), seed = 99)
  model <- "F =~ y1 + y2 + y3 + y4 + y5 + y6"

  out <- suppressWarnings(
    mgcfa_auto(
      model_type = "custom",
      model = model,
      data = dat,
      group = "grp",
      include_steps = c("configural", "metric", "scalar", "strict", "lv.variances", "lv.covariances", "residual.covariances", "regressions", "means"),
      partial_failure_criterion = "none",
      partial_auto_search = "never",
      stop_at_first_unacceptable = FALSE
    )
  )

  expect_true("lv.covariances" %in% names(out$not_applicable_steps))
  expect_true("residual.covariances" %in% names(out$not_applicable_steps))
  expect_true("regressions" %in% names(out$not_applicable_steps))
  expect_true(isTRUE(out$step_failures$lv.covariances$not_applicable))
})
