# Scalar Partial Invariance Demo
# Data from https://rpubs.com/JLLJ/MatMI

library(lavaan)

load_automaticmgcfa_functions <- function() {
  needed <- c("mgcfa_auto", "mgcfa_build_model", "mgcfa_tidy_fit", "mgcfa_plot_fit")

  if (requireNamespace("AutomaticMGCFA", quietly = TRUE)) {
    suppressPackageStartupMessages(library(AutomaticMGCFA))
  }

  has_needed <- all(vapply(needed, exists, logical(1L), mode = "function"))
  if (isTRUE(has_needed)) {
    return(invisible(TRUE))
  }

  script_file <- NULL
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    of <- frames[[i]]$ofile
    if (!is.null(of) && nzchar(as.character(of))) {
      script_file <- normalizePath(as.character(of), winslash = "/", mustWork = FALSE)
      break
    }
  }

  repo_candidates <- c(
    Sys.getenv("AUTOMATIC_MGCFA_REPO", ""),
    if (!is.null(script_file)) normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = FALSE) else "",
    getwd(),
    normalizePath(file.path(getwd(), "Automatic-MGCFA"), winslash = "/", mustWork = FALSE)
  )
  repo_candidates <- unique(repo_candidates[nzchar(repo_candidates)])

  loaded <- FALSE
  for (repo in repo_candidates) {
    impl <- file.path(repo, "R", "mgcfa.R")
    if (!file.exists(impl)) {
      next
    }
    source(impl, local = .GlobalEnv)
    loaded <- all(vapply(needed, exists, logical(1L), mode = "function"))
    if (isTRUE(loaded)) {
      break
    }
  }

  if (!isTRUE(loaded)) {
    stop(
      "AutomaticMGCFA functions could not be loaded. Install `AutomaticMGCFA` or set `AUTOMATIC_MGCFA_REPO` to the package source directory.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

load_automaticmgcfa_functions()

lowerExp1 <- "
1
0.82    1
0.42    0.42    1
0.5     0.6     0.4     1
0.19    0.22    0.24    0.39    1
0.46    0.51    0.33    0.5     0.26    1"

lowerCon1 <- "
1
0.73    1
0.49    0.38    1
0.37    0.42    0.39    1
0.35    0.44    0.44    0.34    1
0.45    0.37    0.41    0.57    0.36    1"

vars1 <- c("MatCon", "MatDis", "LPS2GK", "LPS2NumSeq", "LPS2MenRot", "LPS2Add")

Exp1Cor <- lavaan::getCov(lowerExp1, names = vars1)
Con1Cor <- lavaan::getCov(lowerCon1, names = vars1)

Exp1SDs <- c(4.39, 4.10, 7.86, 3.65, 6.69, 6.00)
Con1SDs <- c(5.48, 4.83, 9.66, 3.37, 7.54, 5.98)

Exp1Cov <- lavaan::cor2cov(Exp1Cor, Exp1SDs)
Con1Cov <- lavaan::cor2cov(Con1Cor, Con1SDs)

Exp1Means <- c(13.59, 14.16, 34.59, 20.59, 20.98, 14.46)
Con1Means <- c(7.21, 9.07, 35.91, 20.89, 22.95, 14.79)

model_exp1 <- "g =~ MatCon + MatDis + LPS2GK + LPS2NumSeq + LPS2MenRot + LPS2Add"

partial_auto_mode <- if (interactive()) "prompt" else "always"

message("Running Experiment 1 invariance test.")
if (identical(partial_auto_mode, "prompt")) {
  message("If a constrained step fails, you will be prompted to run automatic partial search (y/n).")
  message("Answer 'y' to include both partial and failed non-partial models in plots.")
} else {
  message("Non-interactive session detected: partial search set to 'always'.")
}

fit_exp1 <- mgcfa_auto(
  model_type = "custom",
  model = model_exp1,
  sample_cov = list(Exp1Cov, Con1Cov),
  sample_mean = list(Exp1Means, Con1Means),
  sample_nobs = list(56, 56),
  group_labels = c("Experimental", "Control"),
  include_steps = c("configural", "metric", "scalar", "strict", "lv.variances", "means"),
  estimator = "ML",
  partial_failure_criterion = "chisq_pvalue",
  partial_failure_threshold = 0.05,
  partial_auto_search = partial_auto_mode,
  partial_search_criterion = "chisq_pvalue",
  partial_search_threshold = 0.05,
  partial_search_top_n = 5L,
  partial_search_stop_on_accept = TRUE,
  partial_search_rank = "closest"
)

message("\nExperiment 1 results (3 significant figures):")
print(
  fit_exp1,
  digits = 3,
  rounding = "signif",
  show_freed_parameters = TRUE,
  verbose = FALSE
)

single_partial <- fit_exp1$partial_search
if (!is.null(single_partial) && isTRUE(single_partial$triggered)) {
  message("\nPartial search decision: ", single_partial$decision)
  if (length(single_partial$selected_partial) > 0L) {
    message(
      "Selected freed constraints: ",
      paste(single_partial$selected_partial, collapse = ", ")
    )
  }
}

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p_aic <- mgcfa_plot_fit(
    fit_exp1,
    measures = "aic",
    include_non_partial = TRUE,
    baseline_step = "configural",
    plot_type = "delta",
    color_values = c(selected = "#1B6CA8", non_partial = "#C0392B"),
    linetype_values = c(selected = "solid", non_partial = "dashed")
  )

  p_bic <- mgcfa_plot_fit(
    fit_exp1,
    measures = "bic",
    include_non_partial = TRUE,
    baseline_step = "configural",
    plot_type = "delta",
    color_values = c(selected = "#1B6CA8", non_partial = "#C0392B"),
    linetype_values = c(selected = "solid", non_partial = "dashed")
  )

  print(p_aic)
  print(p_bic)
} else {
  message("Install ggplot2 to generate AIC/BIC plots.")
}
