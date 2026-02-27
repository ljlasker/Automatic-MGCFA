make_two_group_data <- function(n = 350L, seed = 123L) {
  set.seed(seed)
  eta1 <- rnorm(n)
  eta2 <- rnorm(n)
  lam <- c(0.8, 0.75, 0.7, 0.85)

  make_group <- function(eta) {
    y <- sapply(seq_along(lam), function(j) lam[j] * eta + rnorm(n, 0, sqrt(1 - lam[j]^2)))
    y <- as.data.frame(y)
    names(y) <- paste0("x", seq_len(ncol(y)))
    y
  }

  g1 <- make_group(eta1)
  g2 <- make_group(eta2)
  g1$grp <- "g1"
  g2$grp <- "g2"
  out <- rbind(g1, g2)
  out$grp <- factor(out$grp)
  out
}

test_that("chisq-pvalue criterion records lavTestLRT-based method", {
  dat <- make_two_group_data(n = 450L, seed = 77)
  model <- "F =~ x1 + x2 + x3 + x4"

  out <- mgcfa_auto(
    model_type = "custom",
    model = model,
    data = dat,
    group = "grp",
    include_steps = c("configural", "metric", "scalar"),
    estimator = "MLR",
    std.lv_configural = FALSE,
    std.lv_constrained = FALSE,
    partial_failure_criterion = "chisq_pvalue",
    partial_auto_search = "never",
    stop_at_first_unacceptable = FALSE
  )

  method_metric <- out$step_failures$metric$details$per_rule[[1]]$details$method
  expect_identical(method_metric, "lavTestLRT")

  expect_true("method" %in% names(out$lrt_table))
  expect_true(all(out$lrt_table$method %in% c("lavTestLRT", "manual_diff")))
})
