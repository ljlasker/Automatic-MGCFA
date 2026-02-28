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
