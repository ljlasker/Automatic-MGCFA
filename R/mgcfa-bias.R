#' Compute Bias Effect Sizes from Partial Invariance Recovery
#'
#' Produces signed bias-effect contrasts (dMACS-style deltas) by comparing the
#' failed non-partial model at each recovered step to the selected bias-adjusted
#' model. Supports observed item/composite outcomes and latent mean/SD
#' comparisons, with optional decomposition by freed parameter class. Also
#' returns item/composite bias effect-size families used in JLLJ RPubs and the
#' cited literature, including \code{dMACS}, \code{dMACS_Signed},
#' \code{dMACS_True}, \code{SDI2}, \code{UDI2}, and \code{SUDI2}.
#'
#' @param x An object returned by \code{mgcfa_auto()}.
#' @param steps Optional character vector of recovered steps to summarize.
#'   Defaults to all steps with both failed non-partial and selected recovered
#'   models.
#' @param groups Optional length-2 vector indicating reference and focal groups
#'   (labels or indices). Defaults to the first two groups.
#' @param composites Optional named list defining observed composites. Each
#'   element may be (a) a character vector of observed variable names (equal
#'   weights) or (b) a named numeric vector of weights.
#' @param include_items Logical; if \code{TRUE}, include one-item composites for
#'   each observed indicator.
#' @param include_class_decomposition Logical; if \code{TRUE}, also estimates
#'   class-specific bias effects by refitting the failed model with only the
#'   freed terms from each parameter class.
#' @param ci Logical; if \code{TRUE}, computes parametric bootstrap confidence
#'   intervals for effect-size metrics.
#' @param n_boot Number of bootstrap draws when \code{ci = TRUE}.
#' @param conf_level Confidence level for bootstrap intervals.
#' @param boot_seed Optional random seed for bootstrap replicates.
#' @param digits Number of digits used for numeric formatting.
#' @param rounding Numeric rounding mode. One of \code{"signif"},
#'   \code{"round"}, or \code{"none"}.
#'
#' @return An object of class \code{"mgcfa_bias_effects"} containing summary
#'   tables for observed and latent bias effects.
#' @export
mgcfa_bias_effects <- function(
  x,
  steps = NULL,
  groups = NULL,
  composites = NULL,
  include_items = TRUE,
  include_class_decomposition = TRUE,
  ci = FALSE,
  n_boot = 200L,
  conf_level = 0.95,
  boot_seed = NULL,
  digits = 3L,
  rounding = c("signif", "round", "none")
) {
  x <- .mgcfa_as_result(x)
  rounding <- match.arg(rounding)
  if (!is.logical(include_items) || length(include_items) != 1L || is.na(include_items)) {
    stop("`include_items` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(include_class_decomposition) || length(include_class_decomposition) != 1L || is.na(include_class_decomposition)) {
    stop("`include_class_decomposition` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(ci) || length(ci) != 1L || is.na(ci)) {
    stop("`ci` must be TRUE or FALSE.", call. = FALSE)
  }
  n_boot <- as.integer(n_boot)
  if (is.na(n_boot) || n_boot < 20L) {
    n_boot <- 200L
  }
  conf_level <- as.numeric(conf_level)[[1L]]
  if (!is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be in (0, 1).", call. = FALSE)
  }
  if (!is.null(boot_seed)) {
    boot_seed <- as.integer(boot_seed)[[1L]]
    if (is.na(boot_seed)) {
      stop("`boot_seed` must be NULL or an integer.", call. = FALSE)
    }
  }

  fit_steps <- names(x$fits)
  candidate_steps <- Filter(function(step) {
    sf <- (x$step_failures %||% list())[[step]] %||% NULL
    ps <- (x$partial_searches %||% list())[[step]] %||% NULL
    has_failed <- !is.null(sf) && !is.null(sf$failed_fit)
    has_selected <- !is.null((x$fits %||% list())[[step]])
    has_freed <- length(unique(c(
      .mgcfa_normalize_terms(ps$selected_added_terms %||% character()),
      .mgcfa_normalize_terms((x$freed_parameters %||% list())[[step]] %||% character())
    ))) > 0L
    isTRUE(has_failed && has_selected && has_freed)
  }, fit_steps)
  candidate_steps <- as.character(candidate_steps)
  if (length(candidate_steps) == 0L) {
    stop(
      "No recovered partial-invariance steps with both failed and adjusted models were found in `x`.",
      call. = FALSE
    )
  }

  if (is.null(steps)) {
    steps <- candidate_steps
  } else {
    steps <- unique(as.character(steps))
    bad <- setdiff(steps, candidate_steps)
    if (length(bad) > 0L) {
      stop(
        "Requested `steps` are not available for bias effect computation: ",
        paste(bad, collapse = ", "),
        call. = FALSE
      )
    }
  }

  observed_tables <- list()
  observed_stats_tables <- list()
  metric_tables <- list()
  metric_ci_tables <- list()
  recovery_ci_tables <- list()
  latent_tables <- list()
  latent_stats_tables <- list()
  class_tables <- list()
  freed_rows <- list()

  for (step in steps) {
    sf <- x$step_failures[[step]]
    ps <- (x$partial_searches %||% list())[[step]] %||% list()
    failed_fit <- sf$failed_fit
    adjusted_fit <- x$fits[[step]]
    if (is.null(failed_fit) || is.null(adjusted_fit)) {
      next
    }

    group_pair <- .mgcfa_resolve_group_pair(failed_fit, groups = groups)
    comp_weights <- .mgcfa_prepare_composites(
      fit = failed_fit,
      composites = composites,
      include_items = isTRUE(include_items)
    )

    obs_biased <- .mgcfa_model_observed_stats(
      fit = failed_fit,
      weights = comp_weights,
      model_label = "Biased (Non-Partial)"
    )
    obs_adjusted <- .mgcfa_model_observed_stats(
      fit = adjusted_fit,
      weights = comp_weights,
      model_label = "Bias-Adjusted (Selected)"
    )
    obs_cmp <- .mgcfa_observed_bias_delta(
      biased_stats = obs_biased,
      adjusted_stats = obs_adjusted,
      ref_group = group_pair$reference,
      focal_group = group_pair$focal,
      step = step
    )
    observed_tables[[step]] <- obs_cmp
    metrics_biased <- .mgcfa_effect_size_table(
      fit = failed_fit,
      weights = comp_weights,
      ref_group = group_pair$reference,
      focal_group = group_pair$focal,
      step = step,
      model_label = "Biased (Non-Partial)"
    )
    metrics_adjusted <- .mgcfa_effect_size_table(
      fit = adjusted_fit,
      weights = comp_weights,
      ref_group = group_pair$reference,
      focal_group = group_pair$focal,
      step = step,
      model_label = "Bias-Adjusted (Selected)"
    )
    metric_tables[[step]] <- rbind(metrics_biased, metrics_adjusted)
    if (isTRUE(ci)) {
      seed_step <- if (is.null(boot_seed)) NULL else (boot_seed + which(steps == step) - 1L)
      boot_b <- .mgcfa_boot_effect_sizes(
        fit = failed_fit,
        weights = comp_weights,
        ref_group = group_pair$reference,
        focal_group = group_pair$focal,
        step = step,
        model_label = "Biased (Non-Partial)",
        n_boot = n_boot,
        conf_level = conf_level,
        seed = seed_step
      )
      boot_a <- .mgcfa_boot_effect_sizes(
        fit = adjusted_fit,
        weights = comp_weights,
        ref_group = group_pair$reference,
        focal_group = group_pair$focal,
        step = step,
        model_label = "Bias-Adjusted (Selected)",
        n_boot = n_boot,
        conf_level = conf_level,
        seed = if (is.null(seed_step)) NULL else seed_step + 10000L
      )
      if (nrow(boot_b$summary) > 0L || nrow(boot_a$summary) > 0L) {
        metric_ci_tables[[step]] <- rbind(boot_b$summary, boot_a$summary)
      }
      rec_ci <- .mgcfa_boot_recovery_ci(
        boot_b = boot_b$draws,
        boot_a = boot_a$draws,
        conf_level = conf_level
      )
      if (nrow(rec_ci) > 0L) {
        recovery_ci_tables[[step]] <- rec_ci
      }
    }
    observed_stats_tables[[step]] <- rbind(
      transform(obs_biased, step = step, stringsAsFactors = FALSE),
      transform(obs_adjusted, step = step, stringsAsFactors = FALSE)
    )

    latent_biased <- .mgcfa_model_latent_stats(
      fit = failed_fit,
      model_label = "Biased (Non-Partial)"
    )
    latent_adjusted <- .mgcfa_model_latent_stats(
      fit = adjusted_fit,
      model_label = "Bias-Adjusted (Selected)"
    )
    latent_cmp <- .mgcfa_latent_bias_delta(
      biased_stats = latent_biased,
      adjusted_stats = latent_adjusted,
      ref_group = group_pair$reference,
      focal_group = group_pair$focal,
      step = step
    )
    latent_tables[[step]] <- latent_cmp
    latent_stats_tables[[step]] <- rbind(
      transform(latent_biased, step = step, stringsAsFactors = FALSE),
      transform(latent_adjusted, step = step, stringsAsFactors = FALSE)
    )

    selected_added <- .mgcfa_normalize_terms(
      ps$selected_added_terms %||% ((x$freed_parameters %||% list())[[step]] %||% character())
    )
    selected_partial <- .mgcfa_normalize_terms(ps$selected_partial %||% character())
    base_partial <- setdiff(selected_partial, selected_added)
    term_classes <- .mgcfa_term_classes(adjusted_fit, selected_added)

    if (length(selected_added) > 0L) {
      freed_rows[[step]] <- data.frame(
        step = step,
        term = selected_added,
        class = unname(term_classes[selected_added]),
        stringsAsFactors = FALSE
      )
    }

    if (isTRUE(include_class_decomposition) && length(selected_added) > 0L) {
      cls <- unique(unname(term_classes[selected_added]))
      cls <- cls[nzchar(cls)]
      step_class_rows <- list()
      for (cls_i in cls) {
        cls_terms <- names(term_classes)[term_classes == cls_i]
        cls_terms <- intersect(cls_terms, selected_added)
        if (length(cls_terms) == 0L) {
          next
        }
        fit_cls <- .mgcfa_refit_with_group_partial(
          base_fit = failed_fit,
          group_partial_terms = unique(c(base_partial, cls_terms))
        )
        if (is.null(fit_cls)) {
          next
        }
        obs_cls <- .mgcfa_model_observed_stats(
          fit = fit_cls,
          weights = comp_weights,
          model_label = paste0("Class Partial: ", cls_i)
        )
        cls_cmp <- .mgcfa_observed_delta_against_biased(
          biased_stats = obs_biased,
          compare_stats = obs_cls,
          ref_group = group_pair$reference,
          focal_group = group_pair$focal,
          step = step,
          class_label = cls_i
        )
        if (!is.null(cls_cmp) && nrow(cls_cmp) > 0L) {
          cls_cmp$freed_terms <- paste(cls_terms, collapse = ", ")
          step_class_rows[[length(step_class_rows) + 1L]] <- cls_cmp
        }
      }
      if (length(step_class_rows) > 0L) {
        class_tables[[step]] <- do.call(rbind, step_class_rows)
      }
    }
  }

  observed_effects <- if (length(observed_tables) > 0L) do.call(rbind, observed_tables) else data.frame()
  observed_group_stats <- if (length(observed_stats_tables) > 0L) do.call(rbind, observed_stats_tables) else data.frame()
  effect_size_metrics <- if (length(metric_tables) > 0L) do.call(rbind, metric_tables) else data.frame()
  effect_size_metrics_ci <- if (length(metric_ci_tables) > 0L) do.call(rbind, metric_ci_tables) else data.frame()
  effect_size_recovery <- .mgcfa_effect_size_recovery(effect_size_metrics)
  effect_size_recovery_ci <- if (length(recovery_ci_tables) > 0L) do.call(rbind, recovery_ci_tables) else data.frame()
  latent_effects <- if (length(latent_tables) > 0L) do.call(rbind, latent_tables) else data.frame()
  latent_group_stats <- if (length(latent_stats_tables) > 0L) do.call(rbind, latent_stats_tables) else data.frame()
  class_effects <- if (length(class_tables) > 0L) do.call(rbind, class_tables) else data.frame()
  freed_parameters <- if (length(freed_rows) > 0L) do.call(rbind, freed_rows) else data.frame()

  if (nrow(observed_effects) > 0L) {
    observed_effects <- .mgcfa_format_numeric_df(observed_effects, digits = digits, rounding = rounding)
  }
  if (nrow(observed_group_stats) > 0L) {
    observed_group_stats <- .mgcfa_format_numeric_df(observed_group_stats, digits = digits, rounding = rounding)
  }
  if (nrow(effect_size_metrics) > 0L) {
    effect_size_metrics <- .mgcfa_format_numeric_df(effect_size_metrics, digits = digits, rounding = rounding)
  }
  if (nrow(effect_size_metrics_ci) > 0L) {
    effect_size_metrics_ci <- .mgcfa_format_numeric_df(effect_size_metrics_ci, digits = digits, rounding = rounding)
  }
  if (nrow(effect_size_recovery) > 0L) {
    effect_size_recovery <- .mgcfa_format_numeric_df(effect_size_recovery, digits = digits, rounding = rounding)
  }
  if (nrow(effect_size_recovery_ci) > 0L) {
    effect_size_recovery_ci <- .mgcfa_format_numeric_df(effect_size_recovery_ci, digits = digits, rounding = rounding)
  }
  if (nrow(latent_effects) > 0L) {
    latent_effects <- .mgcfa_format_numeric_df(latent_effects, digits = digits, rounding = rounding)
  }
  if (nrow(latent_group_stats) > 0L) {
    latent_group_stats <- .mgcfa_format_numeric_df(latent_group_stats, digits = digits, rounding = rounding)
  }
  if (nrow(class_effects) > 0L) {
    class_effects <- .mgcfa_format_numeric_df(class_effects, digits = digits, rounding = rounding)
  }

  out <- list(
    steps = steps,
    groups = if (exists("group_pair")) c(group_pair$reference, group_pair$focal) else NULL,
    observed_effects = observed_effects,
    observed_group_stats = observed_group_stats,
    effect_size_metrics = effect_size_metrics,
    effect_size_metrics_ci = effect_size_metrics_ci,
    effect_size_recovery = effect_size_recovery,
    effect_size_recovery_ci = effect_size_recovery_ci,
    latent_effects = latent_effects,
    latent_group_stats = latent_group_stats,
    class_effects = class_effects,
    freed_parameters = freed_parameters
  )
  class(out) <- "mgcfa_bias_effects"
  out
}

#' Print Method for Bias Effect Summaries
#'
#' @param x An object returned by \code{mgcfa_bias_effects()}.
#' @param ... Unused.
#'
#' @return The input object invisibly.
#' @export
print.mgcfa_bias_effects <- function(x, ...) {
  if (!inherits(x, "mgcfa_bias_effects")) {
    stop("`x` must be an object returned by `mgcfa_bias_effects()`.", call. = FALSE)
  }
  cat("Steps:", paste(x$steps, collapse = ", "), "\n")
  if (!is.null(x$groups) && length(x$groups) == 2L) {
    cat("Comparison:", x$groups[[2L]], "vs", x$groups[[1L]], "\n")
  }
  if (!is.null(x$observed_effects) && nrow(x$observed_effects) > 0L) {
    cat("\nObserved bias effects (signed dMACS-style deltas)\n")
    print(utils::head(x$observed_effects, n = 12L), row.names = FALSE)
  }
  if (!is.null(x$effect_size_metrics) && nrow(x$effect_size_metrics) > 0L) {
    cat("\nBias effect-size metrics (dMACS, dMACS_Signed, SDI2/UDI2)\n")
    print(utils::head(x$effect_size_metrics, n = 12L), row.names = FALSE)
  }
  if (!is.null(x$effect_size_metrics_ci) && nrow(x$effect_size_metrics_ci) > 0L) {
    cat("\nBootstrap CIs for effect-size metrics\n")
    print(utils::head(x$effect_size_metrics_ci, n = 12L), row.names = FALSE)
  }
  if (!is.null(x$effect_size_recovery) && nrow(x$effect_size_recovery) > 0L) {
    cat("\nEffect-size recovery (Adjusted - Biased)\n")
    print(utils::head(x$effect_size_recovery, n = 12L), row.names = FALSE)
  }
  if (!is.null(x$effect_size_recovery_ci) && nrow(x$effect_size_recovery_ci) > 0L) {
    cat("\nBootstrap CIs for recovery deltas\n")
    print(utils::head(x$effect_size_recovery_ci, n = 12L), row.names = FALSE)
  }
  if (!is.null(x$latent_effects) && nrow(x$latent_effects) > 0L) {
    cat("\nLatent mean/SD comparisons (biased vs adjusted)\n")
    print(utils::head(x$latent_effects, n = 12L), row.names = FALSE)
  }
  if (!is.null(x$latent_group_stats) && nrow(x$latent_group_stats) > 0L) {
    cat("\nLatent group and pooled means/SDs\n")
    print(utils::head(x$latent_group_stats, n = 12L), row.names = FALSE)
  }
  if (!is.null(x$class_effects) && nrow(x$class_effects) > 0L) {
    cat("\nClass decomposition of observed bias effects\n")
    print(utils::head(x$class_effects, n = 12L), row.names = FALSE)
  }
  invisible(x)
}

#' Print Method for MGCFA Reports
#'
#' @param x An object returned by \code{mgcfa_report()}.
#' @param ... Unused.
#'
#' @return The input object invisibly.
#' @export
print.mgcfa_report <- function(x, ...) {
  if (!inherits(x, "mgcfa_report")) {
    stop("`x` must be an object returned by `mgcfa_report()`.", call. = FALSE)
  }
  cat("AutomaticMGCFA Report\n")
  if (!is.null(x$overview) && nrow(x$overview) > 0L) {
    cat("\nOverview\n")
    print(x$overview, row.names = FALSE)
  }
  if (!is.null(x$failures) && nrow(x$failures) > 0L) {
    cat("\nFailures\n")
    print(x$failures, row.names = FALSE)
  }
  if (!is.null(x$input_notes) && length(x$input_notes) > 0L) {
    cat("\nInput notes:\n")
    for (msg in x$input_notes) {
      cat("- ", msg, "\n", sep = "")
    }
  }
  invisible(x)
}

.mgcfa_step_overview <- function(x) {
  if (!inherits(x, "mgcfa_result") || is.null(x$fits) || length(x$fits) == 0L) {
    return(data.frame())
  }
  steps <- names(x$fits)
  out <- vector("list", length(steps))
  for (i in seq_along(steps)) {
    step <- steps[[i]]
    fm <- suppressWarnings(lavaan::fitMeasures(x$fits[[step]], c("cfi", "rmsea", "srmr", "aic", "bic")))
    freed <- .mgcfa_normalize_terms(x$freed_parameters[[step]] %||% character())
    fail_rec <- x$step_failures[[step]]
    status <- if (!is.null((x$not_applicable_steps %||% list())[[step]])) {
      "not_applicable"
    } else if (is.null(fail_rec) || !isTRUE(fail_rec$failed)) {
      "ok"
    } else if (step %in% (x$recovered_steps %||% character())) {
      "recovered_partial"
    } else {
      "failed"
    }
    out[[i]] <- data.frame(
      step = step,
      status = status,
      model_used = if (length(freed) > 0L) "partial" else "full",
      n_freed = length(freed),
      cfi = as.numeric(fm[["cfi"]]),
      rmsea = as.numeric(fm[["rmsea"]]),
      srmr = as.numeric(fm[["srmr"]]),
      aic = as.numeric(fm[["aic"]]),
      bic = as.numeric(fm[["bic"]]),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

.mgcfa_failure_summary <- function(x, digits = 3L, rounding = c("signif", "round", "none")) {
  rounding <- match.arg(rounding)
  sf <- x$step_failures %||% list()
  if (length(sf) == 0L) {
    return(data.frame())
  }
  rows <- list()
  out_i <- 0L
  for (step in names(sf)) {
    z <- sf[[step]]
    if (is.null(z) || (!isTRUE(z$failed) && !isTRUE(z$not_applicable))) {
      next
    }
    out_i <- out_i + 1L
    rows[[out_i]] <- data.frame(
      step = step,
      failed = isTRUE(z$failed),
      not_applicable = isTRUE(z$not_applicable),
      criterion = as.character(z$criterion %||% NA_character_),
      threshold = as.numeric(z$threshold %||% NA_real_),
      value = as.numeric(z$value %||% NA_real_),
      reason = as.character(z$reason %||% NA_character_),
      stringsAsFactors = FALSE
    )
  }
  if (out_i == 0L) {
    return(data.frame())
  }
  .mgcfa_format_numeric_df(do.call(rbind, rows), digits = digits, rounding = rounding)
}

.mgcfa_print_result_input_summary <- function(input_summary, input_notes = character()) {
  if (is.null(input_summary) || !is.list(input_summary)) {
    if (length(input_notes) > 0L) {
      cat("Input notes:\n")
      for (msg in unique(as.character(input_notes))) {
        cat("- ", msg, "\n", sep = "")
      }
    }
    return(invisible(NULL))
  }

  if (identical(input_summary$mode %||% "", "raw_data")) {
    cat(
      "Input: raw data | rows =",
      input_summary$n_rows %||% NA_integer_,
      "| groups =",
      input_summary$n_groups %||% NA_integer_,
      "| group var =",
      input_summary$group %||% "",
      "\n"
    )
    if (isTRUE(input_summary$ordered)) {
      ov <- input_summary$ordered_variables %||% character()
      cat("Ordered indicators:", if (length(ov) > 0L) paste(ov, collapse = ", ") else "TRUE", "\n")
    }
  } else if (identical(input_summary$mode %||% "", "summary_matrices")) {
    cat(
      "Input: summary matrices | type =",
      input_summary$matrix_type %||% "",
      "| groups =",
      input_summary$n_groups %||% NA_integer_,
      "| variables =",
      input_summary$n_variables %||% NA_integer_,
      "\n"
    )
    cat(
      "Summary details: means =",
      isTRUE(input_summary$has_means),
      "| sds =",
      isTRUE(input_summary$has_sds),
      "\n"
    )
  }

  if (length(input_notes) > 0L) {
    cat("Input notes:\n")
    for (msg in unique(as.character(input_notes))) {
      cat("- ", msg, "\n", sep = "")
    }
  }
  invisible(NULL)
}

.mgcfa_resolve_group_pair <- function(fit, groups = NULL) {
  gl <- tryCatch(lavaan::lavInspect(fit, "group.label"), error = function(e) NULL)
  if (is.null(gl) || length(gl) < 2L) {
    stop("Bias effects require at least two groups.", call. = FALSE)
  }
  gl <- as.character(gl)
  if (is.null(groups)) {
    return(list(reference = gl[[1L]], focal = gl[[2L]]))
  }
  groups <- unlist(groups, use.names = FALSE)
  if (length(groups) != 2L) {
    stop("`groups` must be NULL or length 2 (reference, focal).", call. = FALSE)
  }
  if (is.numeric(groups)) {
    idx <- as.integer(groups)
    if (any(is.na(idx)) || any(idx < 1L) || any(idx > length(gl))) {
      stop("Numeric `groups` indices are out of range.", call. = FALSE)
    }
    return(list(reference = gl[[idx[[1L]]]], focal = gl[[idx[[2L]]]]))
  }
  groups <- as.character(groups)
  if (!all(groups %in% gl)) {
    stop("Character `groups` must match lavaan group labels.", call. = FALSE)
  }
  list(reference = groups[[1L]], focal = groups[[2L]])
}

.mgcfa_prepare_composites <- function(fit, composites = NULL, include_items = TRUE) {
  ov <- tryCatch(lavaan::lavNames(fit, "ov"), error = function(e) character())
  if (length(ov) == 0L) {
    stop("No observed variables were found in the fitted model.", call. = FALSE)
  }
  empty_w <- stats::setNames(rep(0, length(ov)), ov)
  out <- list()

  if (is.null(composites)) {
    out[["Observed Composite"]] <- stats::setNames(rep(1 / length(ov), length(ov)), ov)
  } else {
    if (!is.list(composites) || length(composites) == 0L) {
      stop("`composites` must be NULL or a non-empty named list.", call. = FALSE)
    }
    if (is.null(names(composites)) || any(!nzchar(names(composites)))) {
      stop("`composites` must be a named list.", call. = FALSE)
    }
    for (nm in names(composites)) {
      def <- composites[[nm]]
      w <- empty_w
      if (is.character(def)) {
        vars <- unique(as.character(def))
        if (!all(vars %in% ov)) {
          stop(sprintf("Composite `%s` includes variables not in the model.", nm), call. = FALSE)
        }
        w[vars] <- 1 / length(vars)
      } else if (is.numeric(def) && !is.null(names(def))) {
        vars <- names(def)
        if (!all(vars %in% ov)) {
          stop(sprintf("Composite `%s` includes weight names not in the model.", nm), call. = FALSE)
        }
        w[vars] <- as.numeric(def)
      } else {
        stop(sprintf(
          "Composite `%s` must be a character vector of variable names or a named numeric weight vector.",
          nm
        ), call. = FALSE)
      }
      out[[nm]] <- w
    }
  }

  if (isTRUE(include_items)) {
    for (v in ov) {
      w <- empty_w
      w[[v]] <- 1
      out[[paste0("Item: ", v)]] <- w
    }
  }
  out
}

.mgcfa_model_observed_stats <- function(fit, weights, model_label) {
  implied <- tryCatch(lavaan::lavInspect(fit, "implied"), error = function(e) NULL)
  gl <- tryCatch(lavaan::lavInspect(fit, "group.label"), error = function(e) NULL)
  nobs <- tryCatch(as.numeric(lavaan::lavInspect(fit, "nobs")), error = function(e) NULL)
  if (is.null(implied) || is.null(gl) || length(gl) < 1L) {
    return(data.frame())
  }
  if (!is.list(implied)) {
    return(data.frame())
  }
  gl <- as.character(gl)
  n_g <- length(gl)
  if (is.null(nobs) || length(nobs) != n_g || any(!is.finite(nobs))) {
    nobs <- rep(1, n_g)
  }

  out <- list()
  out_i <- 0L
  for (comp_name in names(weights)) {
    w <- as.numeric(weights[[comp_name]])
    names(w) <- names(weights[[comp_name]])
    mean_vec <- numeric(n_g)
    sd_vec <- numeric(n_g)
    for (g in seq_len(n_g)) {
      mu_g <- as.numeric(implied[[g]]$mean)
      s_g <- implied[[g]]$cov
      if (is.null(dim(s_g))) {
        s_g <- matrix(s_g, nrow = length(w), ncol = length(w))
      }
      nm <- colnames(s_g) %||% names(w)
      names(mu_g) <- names(mu_g) %||% nm
      w_use <- w[nm]
      w_use[is.na(w_use)] <- 0
      mu_use <- mu_g[nm]
      mu_use[is.na(mu_use)] <- 0
      var_g <- as.numeric(t(w_use) %*% s_g[nm, nm, drop = FALSE] %*% w_use)
      sd_g <- sqrt(max(var_g, 0))
      mean_g <- sum(w_use * mu_use)

      mean_vec[[g]] <- mean_g
      sd_vec[[g]] <- sd_g
      out_i <- out_i + 1L
      out[[out_i]] <- data.frame(
        model = model_label,
        composite = comp_name,
        group = gl[[g]],
        mean = mean_g,
        sd = sd_g,
        n = nobs[[g]],
        stringsAsFactors = FALSE
      )
    }

    n_tot <- sum(nobs)
    w_n <- nobs / n_tot
    mean_pool <- sum(w_n * mean_vec)
    var_pool <- if (n_tot > 1) {
      sum((nobs - 1) * (sd_vec^2) + nobs * (mean_vec - mean_pool)^2) / (n_tot - 1)
    } else {
      sum(w_n * sd_vec^2)
    }
    out_i <- out_i + 1L
    out[[out_i]] <- data.frame(
      model = model_label,
      composite = comp_name,
      group = "Pooled",
      mean = mean_pool,
      sd = sqrt(max(var_pool, 0)),
      n = n_tot,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

.mgcfa_effect_size_table <- function(
  fit,
  weights,
  ref_group,
  focal_group,
  step,
  model_label
) {
  est <- tryCatch(lavaan::lavInspect(fit, "est"), error = function(e) NULL)
  gl <- tryCatch(as.character(lavaan::lavInspect(fit, "group.label")), error = function(e) NULL)
  nobs <- tryCatch(as.numeric(lavaan::lavInspect(fit, "nobs")), error = function(e) NULL)
  if (is.null(est) || is.null(gl) || !is.list(est) || length(gl) < 2L) {
    return(data.frame())
  }

  g_ref <- match(ref_group, gl)
  g_foc <- match(focal_group, gl)
  if (is.na(g_ref) || is.na(g_foc)) {
    return(data.frame())
  }
  est_ref <- est[[g_ref]]
  est_foc <- est[[g_foc]]
  if (is.null(est_ref$lambda) || is.null(est_foc$lambda) ||
      is.null(est_ref$nu) || is.null(est_foc$nu)) {
    return(data.frame())
  }

  lambda_ref <- as.matrix(est_ref$lambda)
  lambda_foc <- as.matrix(est_foc$lambda)
  nu_ref <- as.numeric(est_ref$nu)
  nu_foc <- as.numeric(est_foc$nu)
  ov_names <- rownames(lambda_ref) %||% names(nu_ref)
  if (is.null(ov_names)) {
    ov_names <- paste0("V", seq_len(nrow(lambda_ref)))
  }
  rownames(lambda_ref) <- ov_names
  rownames(lambda_foc) <- ov_names
  names(nu_ref) <- ov_names
  names(nu_foc) <- ov_names

  lv_ref <- colnames(lambda_ref) %||% paste0("LV", seq_len(ncol(lambda_ref)))
  lv_foc <- colnames(lambda_foc) %||% paste0("LV", seq_len(ncol(lambda_foc)))
  common_lv <- intersect(lv_ref, lv_foc)
  if (length(common_lv) == 0L) {
    return(data.frame())
  }
  colnames(lambda_ref) <- lv_ref
  colnames(lambda_foc) <- lv_foc
  lambda_ref <- lambda_ref[, common_lv, drop = FALSE]
  lambda_foc <- lambda_foc[, common_lv, drop = FALSE]

  alpha_raw <- est_foc$alpha %||% rep(0, ncol(lambda_foc))
  alpha_foc <- as.numeric(alpha_raw)
  names(alpha_foc) <- names(alpha_raw) %||% common_lv
  alpha_foc <- alpha_foc[common_lv]
  alpha_foc[is.na(alpha_foc)] <- 0

  psi_foc <- as.matrix(est_foc$psi %||% diag(length(common_lv)))
  rownames(psi_foc) <- rownames(psi_foc) %||% common_lv
  colnames(psi_foc) <- colnames(psi_foc) %||% common_lv
  if (!all(common_lv %in% rownames(psi_foc)) || !all(common_lv %in% colnames(psi_foc))) {
    psi_foc <- diag(length(common_lv))
    rownames(psi_foc) <- common_lv
    colnames(psi_foc) <- common_lv
  } else {
    psi_foc <- psi_foc[common_lv, common_lv, drop = FALSE]
  }

  # Item/composite SDs from implied observed moments.
  obs_stats <- .mgcfa_model_observed_stats(fit, weights = weights, model_label = model_label)
  if (nrow(obs_stats) == 0L) {
    return(data.frame())
  }

  if (is.null(nobs) || length(nobs) != length(gl) || any(!is.finite(nobs))) {
    nobs <- rep(1, length(gl))
  }
  n_ref <- nobs[[g_ref]]
  n_foc <- nobs[[g_foc]]

  out <- list()
  out_i <- 0L
  for (nm in names(weights)) {
    w <- as.numeric(weights[[nm]])
    names(w) <- names(weights[[nm]])
    w <- w[ov_names]
    w[is.na(w)] <- 0

    delta_nu <- nu_ref - nu_foc
    delta_lambda <- lambda_ref - lambda_foc
    a <- sum(w * delta_nu)
    b <- as.numeric(colSums(w * delta_lambda))

    mu_delta <- a + sum(b * alpha_foc)
    var_delta <- as.numeric(t(b) %*% psi_foc %*% b)
    if (!is.finite(var_delta) || var_delta < 0) {
      var_delta <- 0
    }
    sd_delta <- sqrt(var_delta)
    abs_expect <- .mgcfa_abs_normal_expectation(mu = mu_delta, sd = sd_delta)

    sd_ref <- obs_stats$sd[obs_stats$composite == nm & obs_stats$group == ref_group][1]
    sd_foc <- obs_stats$sd[obs_stats$composite == nm & obs_stats$group == focal_group][1]
    sd_pool <- .mgcfa_two_group_pooled_sd(
      sd1 = sd_ref,
      sd2 = sd_foc,
      n1 = n_ref,
      n2 = n_foc
    )

    dmacs <- if (is.finite(sd_pool) && sd_pool > 0) sqrt(mu_delta^2 + var_delta) / sd_pool else NA_real_
    dmacs_signed <- if (is.finite(sd_pool) && sd_pool > 0) mu_delta / sd_pool else NA_real_
    dmacs_true <- if (is.finite(dmacs_signed)) sign(dmacs_signed) * dmacs else NA_real_

    sdi2 <- if (is.finite(sd_foc) && sd_foc > 0) mu_delta / sd_foc else NA_real_
    udi2 <- if (is.finite(sd_foc) && sd_foc > 0) abs_expect / sd_foc else NA_real_
    sudi2 <- if (is.finite(sdi2)) sign(sdi2) * udi2 else NA_real_

    out_i <- out_i + 1L
    out[[out_i]] <- data.frame(
      step = step,
      model = model_label,
      comparison = paste(ref_group, "-", focal_group),
      outcome = nm,
      dmacs = as.numeric(dmacs),
      dmacs_signed = as.numeric(dmacs_signed),
      dmacs_true = as.numeric(dmacs_true),
      sdi2 = as.numeric(sdi2),
      udi2 = as.numeric(udi2),
      sudi2 = as.numeric(sudi2),
      delta_mean = as.numeric(mu_delta),
      delta_sd = as.numeric(sd_delta),
      denom_pooled_sd = as.numeric(sd_pool),
      denom_focal_sd = as.numeric(sd_foc),
      stringsAsFactors = FALSE
    )
  }

  if (out_i == 0L) {
    return(data.frame())
  }
  out_df <- do.call(rbind, out)

  agg <- data.frame(
    step = step,
    model = model_label,
    comparison = paste(ref_group, "-", focal_group),
    outcome = "Average Across Outcomes",
    dmacs = mean(out_df$dmacs, na.rm = TRUE),
    dmacs_signed = mean(out_df$dmacs_signed, na.rm = TRUE),
    dmacs_true = mean(out_df$dmacs_true, na.rm = TRUE),
    sdi2 = mean(out_df$sdi2, na.rm = TRUE),
    udi2 = mean(out_df$udi2, na.rm = TRUE),
    sudi2 = mean(out_df$sudi2, na.rm = TRUE),
    delta_mean = mean(out_df$delta_mean, na.rm = TRUE),
    delta_sd = mean(out_df$delta_sd, na.rm = TRUE),
    denom_pooled_sd = mean(out_df$denom_pooled_sd, na.rm = TRUE),
    denom_focal_sd = mean(out_df$denom_focal_sd, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  rbind(out_df, agg)
}

.mgcfa_two_group_pooled_sd <- function(sd1, sd2, n1, n2) {
  sd1 <- as.numeric(sd1)
  sd2 <- as.numeric(sd2)
  n1 <- as.numeric(n1)
  n2 <- as.numeric(n2)
  if (!is.finite(sd1) || !is.finite(sd2) || !is.finite(n1) || !is.finite(n2)) {
    return(NA_real_)
  }
  if (n1 <= 1 || n2 <= 1) {
    return(sqrt((sd1^2 + sd2^2) / 2))
  }
  sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
}

.mgcfa_abs_normal_expectation <- function(mu, sd) {
  mu <- as.numeric(mu)
  sd <- as.numeric(sd)
  if (!is.finite(mu) || !is.finite(sd) || sd < 0) {
    return(NA_real_)
  }
  if (sd == 0) {
    return(abs(mu))
  }
  z <- mu / sd
  sd * sqrt(2 / pi) * exp(-0.5 * z^2) + mu * (2 * stats::pnorm(z) - 1)
}

.mgcfa_effect_size_recovery <- function(effect_size_metrics) {
  if (is.null(effect_size_metrics) || nrow(effect_size_metrics) == 0L) {
    return(data.frame())
  }
  need_cols <- c("step", "comparison", "outcome", "model", "dmacs", "dmacs_signed", "dmacs_true", "sdi2", "udi2", "sudi2")
  if (!all(need_cols %in% names(effect_size_metrics))) {
    return(data.frame())
  }
  b <- effect_size_metrics[effect_size_metrics$model == "Biased (Non-Partial)", , drop = FALSE]
  a <- effect_size_metrics[effect_size_metrics$model == "Bias-Adjusted (Selected)", , drop = FALSE]
  if (nrow(b) == 0L || nrow(a) == 0L) {
    return(data.frame())
  }
  key <- c("step", "comparison", "outcome")
  m <- merge(
    b,
    a,
    by = key,
    suffixes = c("_biased", "_adjusted")
  )
  data.frame(
    step = m$step,
    comparison = m$comparison,
    outcome = m$outcome,
    dmacs_change = m$dmacs_adjusted - m$dmacs_biased,
    dmacs_signed_change = m$dmacs_signed_adjusted - m$dmacs_signed_biased,
    dmacs_true_change = m$dmacs_true_adjusted - m$dmacs_true_biased,
    sdi2_change = m$sdi2_adjusted - m$sdi2_biased,
    udi2_change = m$udi2_adjusted - m$udi2_biased,
    sudi2_change = m$sudi2_adjusted - m$sudi2_biased,
    stringsAsFactors = FALSE
  )
}

.mgcfa_boot_group_name <- function(fit) {
  grp <- tryCatch(as.character(fit@call$group), error = function(e) NULL)
  if (is.null(grp) || length(grp) < 1L || !nzchar(grp[[1L]])) {
    return(".group")
  }
  grp[[1L]]
}

.mgcfa_sample_mvn <- function(n, mu, sigma) {
  n <- as.integer(n)
  if (is.na(n) || n < 1L) {
    return(matrix(numeric(), nrow = 0L, ncol = length(mu)))
  }
  mu <- as.numeric(mu)
  p <- length(mu)
  sigma <- as.matrix(sigma)
  if (!is.matrix(sigma) || nrow(sigma) != p || ncol(sigma) != p) {
    sigma <- diag(p)
  }
  sigma <- (sigma + t(sigma)) / 2

  root <- tryCatch(chol(sigma), error = function(e) NULL)
  if (is.null(root)) {
    ev <- eigen(sigma, symmetric = TRUE)
    vals <- pmax(ev$values, 1e-8)
    root <- ev$vectors %*% diag(sqrt(vals), nrow = length(vals), ncol = length(vals))
    root <- t(root)
  }

  z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  x <- z %*% root
  sweep(x, 2, mu, "+")
}

.mgcfa_bootstrap_data <- function(fit) {
  gl <- tryCatch(as.character(lavaan::lavInspect(fit, "group.label")), error = function(e) NULL)
  if (is.null(gl) || length(gl) < 2L) {
    return(NULL)
  }
  group_name <- .mgcfa_boot_group_name(fit)

  observed <- tryCatch(lavaan::lavInspect(fit, "data"), error = function(e) NULL)
  if (is.list(observed) && length(observed) == length(gl)) {
    pieces <- vector("list", length(gl))
    for (i in seq_along(gl)) {
      d_i <- as.data.frame(observed[[i]])
      if (nrow(d_i) < 1L) {
        return(NULL)
      }
      idx <- sample.int(nrow(d_i), size = nrow(d_i), replace = TRUE)
      d_i <- d_i[idx, , drop = FALSE]
      d_i[[group_name]] <- gl[[i]]
      pieces[[i]] <- d_i
    }
    dat <- do.call(rbind, pieces)
    dat[[group_name]] <- factor(dat[[group_name]], levels = gl)
    return(list(data = dat, group = group_name))
  }

  implied <- tryCatch(lavaan::lavInspect(fit, "implied"), error = function(e) NULL)
  nobs <- tryCatch(as.numeric(lavaan::lavInspect(fit, "nobs")), error = function(e) NULL)
  if (!is.list(implied) || is.null(nobs) || length(nobs) != length(gl)) {
    return(NULL)
  }
  pieces <- vector("list", length(gl))
  for (i in seq_along(gl)) {
    mu_i <- implied[[i]]$mean %||% NULL
    sigma_i <- implied[[i]]$cov %||% NULL
    if (is.null(mu_i) || is.null(sigma_i)) {
      return(NULL)
    }
    x_i <- .mgcfa_sample_mvn(n = nobs[[i]], mu = mu_i, sigma = sigma_i)
    d_i <- as.data.frame(x_i)
    names(d_i) <- names(mu_i) %||% colnames(sigma_i) %||% paste0("V", seq_len(ncol(d_i)))
    d_i[[group_name]] <- gl[[i]]
    pieces[[i]] <- d_i
  }
  dat <- do.call(rbind, pieces)
  dat[[group_name]] <- factor(dat[[group_name]], levels = gl)
  list(data = dat, group = group_name)
}

.mgcfa_boot_effect_sizes <- function(
  fit,
  weights,
  ref_group,
  focal_group,
  step,
  model_label,
  n_boot = 200L,
  conf_level = 0.95,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  gl <- tryCatch(as.character(lavaan::lavInspect(fit, "group.label")), error = function(e) NULL)
  nobs <- tryCatch(as.numeric(lavaan::lavInspect(fit, "nobs")), error = function(e) NULL)
  if (is.null(gl) || is.null(nobs) || length(gl) < 2L) {
    return(list(summary = data.frame(), draws = data.frame()))
  }

  draws <- list()
  out_i <- 0L
  for (b in seq_len(n_boot)) {
    boot_data <- .mgcfa_bootstrap_data(fit)
    if (is.null(boot_data) || is.null(boot_data$data)) {
      next
    }
    refit <- tryCatch(
      stats::update(fit, data = boot_data$data, group = boot_data$group),
      error = function(e) NULL
    )
    if (is.null(refit)) {
      next
    }
    tbl_b <- .mgcfa_effect_size_table(
      fit = refit,
      weights = weights,
      ref_group = ref_group,
      focal_group = focal_group,
      step = step,
      model_label = model_label
    )
    if (nrow(tbl_b) == 0L) {
      next
    }
    tbl_b$boot <- b
    out_i <- out_i + 1L
    draws[[out_i]] <- tbl_b
  }
  if (out_i == 0L) {
    return(list(summary = data.frame(), draws = data.frame()))
  }
  d <- do.call(rbind, draws)
  alpha <- (1 - conf_level) / 2
  qlow <- alpha
  qhigh <- 1 - alpha
  metrics <- c("dmacs", "dmacs_signed", "dmacs_true", "sdi2", "udi2", "sudi2")
  sum_rows <- list()
  s_i <- 0L
  keys <- unique(d[, c("step", "model", "comparison", "outcome"), drop = FALSE])
  for (k in seq_len(nrow(keys))) {
    key <- keys[k, , drop = FALSE]
    dk <- d[
      d$step == key$step &
        d$model == key$model &
        d$comparison == key$comparison &
        d$outcome == key$outcome,
      ,
      drop = FALSE
    ]
    row <- key
    for (m in metrics) {
      vals <- as.numeric(dk[[m]])
      vals <- vals[is.finite(vals)]
      row[[paste0(m, "_ci_low")]] <- if (length(vals) > 0L) stats::quantile(vals, probs = qlow, na.rm = TRUE, names = FALSE) else NA_real_
      row[[paste0(m, "_ci_high")]] <- if (length(vals) > 0L) stats::quantile(vals, probs = qhigh, na.rm = TRUE, names = FALSE) else NA_real_
    }
    row$boot_n <- nrow(dk)
    s_i <- s_i + 1L
    sum_rows[[s_i]] <- row
  }
  list(summary = do.call(rbind, sum_rows), draws = d)
}

.mgcfa_boot_recovery_ci <- function(boot_b, boot_a, conf_level = 0.95) {
  if (is.null(boot_b) || is.null(boot_a) || nrow(boot_b) == 0L || nrow(boot_a) == 0L) {
    return(data.frame())
  }
  key <- c("boot", "step", "comparison", "outcome")
  m <- merge(
    boot_b[, c(key, "dmacs", "dmacs_signed", "dmacs_true", "sdi2", "udi2", "sudi2"), drop = FALSE],
    boot_a[, c(key, "dmacs", "dmacs_signed", "dmacs_true", "sdi2", "udi2", "sudi2"), drop = FALSE],
    by = key,
    suffixes = c("_biased", "_adjusted")
  )
  if (nrow(m) == 0L) {
    return(data.frame())
  }
  m$dmacs_change <- m$dmacs_adjusted - m$dmacs_biased
  m$dmacs_signed_change <- m$dmacs_signed_adjusted - m$dmacs_signed_biased
  m$dmacs_true_change <- m$dmacs_true_adjusted - m$dmacs_true_biased
  m$sdi2_change <- m$sdi2_adjusted - m$sdi2_biased
  m$udi2_change <- m$udi2_adjusted - m$udi2_biased
  m$sudi2_change <- m$sudi2_adjusted - m$sudi2_biased

  alpha <- (1 - conf_level) / 2
  qlow <- alpha
  qhigh <- 1 - alpha
  metrics <- c("dmacs_change", "dmacs_signed_change", "dmacs_true_change", "sdi2_change", "udi2_change", "sudi2_change")
  keys <- unique(m[, c("step", "comparison", "outcome"), drop = FALSE])
  out <- list()
  out_i <- 0L
  for (k in seq_len(nrow(keys))) {
    key_k <- keys[k, , drop = FALSE]
    mk <- m[
      m$step == key_k$step &
        m$comparison == key_k$comparison &
        m$outcome == key_k$outcome,
      ,
      drop = FALSE
    ]
    row <- key_k
    for (mm in metrics) {
      vals <- as.numeric(mk[[mm]])
      vals <- vals[is.finite(vals)]
      row[[paste0(mm, "_ci_low")]] <- if (length(vals) > 0L) stats::quantile(vals, probs = qlow, na.rm = TRUE, names = FALSE) else NA_real_
      row[[paste0(mm, "_ci_high")]] <- if (length(vals) > 0L) stats::quantile(vals, probs = qhigh, na.rm = TRUE, names = FALSE) else NA_real_
    }
    row$boot_n <- nrow(mk)
    out_i <- out_i + 1L
    out[[out_i]] <- row
  }
  do.call(rbind, out)
}

.mgcfa_group_d_from_stats <- function(stats_df, ref_group, focal_group) {
  ref <- stats_df[stats_df$group == ref_group, , drop = FALSE]
  focal <- stats_df[stats_df$group == focal_group, , drop = FALSE]
  if (nrow(ref) < 1L || nrow(focal) < 1L) {
    return(data.frame())
  }
  merge_df <- merge(
    ref[, c("composite", "mean", "sd"), drop = FALSE],
    focal[, c("composite", "mean", "sd"), drop = FALSE],
    by = "composite",
    suffixes = c("_ref", "_focal")
  )
  pooled_sd <- sqrt((merge_df$sd_ref^2 + merge_df$sd_focal^2) / 2)
  d <- (merge_df$mean_focal - merge_df$mean_ref) / pooled_sd
  data.frame(
    composite = merge_df$composite,
    mean_diff = merge_df$mean_focal - merge_df$mean_ref,
    d = d,
    stringsAsFactors = FALSE
  )
}

.mgcfa_observed_bias_delta <- function(
  biased_stats,
  adjusted_stats,
  ref_group,
  focal_group,
  step
) {
  d_biased <- .mgcfa_group_d_from_stats(biased_stats, ref_group = ref_group, focal_group = focal_group)
  d_adjusted <- .mgcfa_group_d_from_stats(adjusted_stats, ref_group = ref_group, focal_group = focal_group)
  if (nrow(d_biased) == 0L || nrow(d_adjusted) == 0L) {
    return(data.frame())
  }
  m <- merge(
    d_biased,
    d_adjusted,
    by = "composite",
    suffixes = c("_biased", "_adjusted")
  )
  data.frame(
    step = step,
    comparison = paste(focal_group, "vs", ref_group),
    composite = m$composite,
    d_biased = as.numeric(m$d_biased),
    d_adjusted = as.numeric(m$d_adjusted),
    dmacs_signed = as.numeric(m$d_adjusted - m$d_biased),
    mean_diff_biased = as.numeric(m$mean_diff_biased),
    mean_diff_adjusted = as.numeric(m$mean_diff_adjusted),
    stringsAsFactors = FALSE
  )
}

.mgcfa_observed_delta_against_biased <- function(
  biased_stats,
  compare_stats,
  ref_group,
  focal_group,
  step,
  class_label
) {
  d_biased <- .mgcfa_group_d_from_stats(biased_stats, ref_group = ref_group, focal_group = focal_group)
  d_cmp <- .mgcfa_group_d_from_stats(compare_stats, ref_group = ref_group, focal_group = focal_group)
  if (nrow(d_biased) == 0L || nrow(d_cmp) == 0L) {
    return(data.frame())
  }
  m <- merge(d_biased, d_cmp, by = "composite", suffixes = c("_biased", "_class"))
  data.frame(
    step = step,
    class = class_label,
    comparison = paste(focal_group, "vs", ref_group),
    composite = m$composite,
    d_biased = as.numeric(m$d_biased),
    d_class_partial = as.numeric(m$d_class),
    dmacs_signed = as.numeric(m$d_class - m$d_biased),
    stringsAsFactors = FALSE
  )
}

.mgcfa_model_latent_stats <- function(fit, model_label) {
  mean_lv <- tryCatch(lavaan::lavInspect(fit, "mean.lv"), error = function(e) NULL)
  cov_lv <- tryCatch(lavaan::lavInspect(fit, "cov.lv"), error = function(e) NULL)
  gl <- tryCatch(as.character(lavaan::lavInspect(fit, "group.label")), error = function(e) NULL)
  nobs <- tryCatch(as.numeric(lavaan::lavInspect(fit, "nobs")), error = function(e) NULL)
  if (is.null(mean_lv) || is.null(cov_lv) || is.null(gl)) {
    return(data.frame())
  }
  if (!is.list(mean_lv)) {
    mean_lv <- list(mean_lv)
  }
  if (!is.list(cov_lv)) {
    cov_lv <- list(cov_lv)
  }
  if (length(mean_lv) != length(gl) || length(cov_lv) != length(gl)) {
    return(data.frame())
  }
  if (is.null(nobs) || length(nobs) != length(gl) || any(!is.finite(nobs))) {
    nobs <- rep(1, length(gl))
  }

  out <- list()
  out_i <- 0L
  lv_names <- colnames(mean_lv[[1L]]) %||% names(mean_lv[[1L]])
  for (g in seq_along(gl)) {
    mu <- as.numeric(mean_lv[[g]])
    names(mu) <- lv_names
    sigma <- cov_lv[[g]]
    sdv <- sqrt(pmax(diag(sigma), 0))
    names(sdv) <- colnames(sigma)
    common <- intersect(names(mu), names(sdv))
    for (lv in common) {
      out_i <- out_i + 1L
      out[[out_i]] <- data.frame(
        model = model_label,
        latent = lv,
        group = gl[[g]],
        mean = mu[[lv]],
        sd = sdv[[lv]],
        n = nobs[[g]],
        stringsAsFactors = FALSE
      )
    }
  }
  if (out_i == 0L) {
    return(data.frame())
  }
  df <- do.call(rbind, out)

  pooled <- do.call(rbind, lapply(split(df, df$latent), function(dfi) {
    n_tot <- sum(dfi$n)
    w <- dfi$n / n_tot
    m <- sum(w * dfi$mean)
    v <- if (n_tot > 1) {
      sum((dfi$n - 1) * (dfi$sd^2) + dfi$n * (dfi$mean - m)^2) / (n_tot - 1)
    } else {
      sum(w * dfi$sd^2)
    }
    data.frame(
      model = unique(dfi$model)[[1L]],
      latent = unique(dfi$latent)[[1L]],
      group = "Pooled",
      mean = m,
      sd = sqrt(max(v, 0)),
      n = n_tot,
      stringsAsFactors = FALSE
    )
  }))
  rbind(df, pooled)
}

.mgcfa_latent_group_d <- function(stats_df, ref_group, focal_group) {
  ref <- stats_df[stats_df$group == ref_group, , drop = FALSE]
  focal <- stats_df[stats_df$group == focal_group, , drop = FALSE]
  if (nrow(ref) == 0L || nrow(focal) == 0L) {
    return(data.frame())
  }
  m <- merge(
    ref[, c("latent", "mean", "sd"), drop = FALSE],
    focal[, c("latent", "mean", "sd"), drop = FALSE],
    by = "latent",
    suffixes = c("_ref", "_focal")
  )
  pooled_sd <- sqrt((m$sd_ref^2 + m$sd_focal^2) / 2)
  data.frame(
    latent = m$latent,
    mean_diff = m$mean_focal - m$mean_ref,
    d = (m$mean_focal - m$mean_ref) / pooled_sd,
    stringsAsFactors = FALSE
  )
}

.mgcfa_latent_bias_delta <- function(biased_stats, adjusted_stats, ref_group, focal_group, step) {
  if (nrow(biased_stats) == 0L || nrow(adjusted_stats) == 0L) {
    return(data.frame())
  }
  d_b <- .mgcfa_latent_group_d(biased_stats, ref_group = ref_group, focal_group = focal_group)
  d_a <- .mgcfa_latent_group_d(adjusted_stats, ref_group = ref_group, focal_group = focal_group)
  if (nrow(d_b) == 0L || nrow(d_a) == 0L) {
    return(data.frame())
  }
  m <- merge(d_b, d_a, by = "latent", suffixes = c("_biased", "_adjusted"))
  data.frame(
    step = step,
    comparison = paste(focal_group, "vs", ref_group),
    latent = m$latent,
    d_biased = as.numeric(m$d_biased),
    d_adjusted = as.numeric(m$d_adjusted),
    dmacs_signed = as.numeric(m$d_adjusted - m$d_biased),
    mean_diff_biased = as.numeric(m$mean_diff_biased),
    mean_diff_adjusted = as.numeric(m$mean_diff_adjusted),
    stringsAsFactors = FALSE
  )
}

.mgcfa_term_classes <- function(fit, terms) {
  terms <- .mgcfa_normalize_terms(terms)
  if (length(terms) == 0L) {
    return(stats::setNames(character(), character()))
  }
  pt <- tryCatch(lavaan::parTable(fit), error = function(e) NULL)
  if (is.null(pt) || nrow(pt) == 0L) {
    return(stats::setNames(rep("unknown", length(terms)), terms))
  }
  ov <- tryCatch(lavaan::lavNames(fit, type = "ov"), error = function(e) character())
  lv <- tryCatch(lavaan::lavNames(fit, type = "lv"), error = function(e) character())
  rows <- vector("list", nrow(pt))
  out_i <- 0L
  for (i in seq_len(nrow(pt))) {
    term <- .mgcfa_term_from_parts(
      lhs = as.character(pt$lhs[[i]]),
      op = as.character(pt$op[[i]]),
      rhs = as.character(pt$rhs[[i]])
    )
    cls <- .mgcfa_constraint_class(
      op = as.character(pt$op[[i]]),
      lhs = as.character(pt$lhs[[i]]),
      rhs = as.character(pt$rhs[[i]]),
      ov_names = ov,
      lv_names = lv
    )
    if (is.na(cls)) {
      cls <- "unknown"
    }
    out_i <- out_i + 1L
    rows[[out_i]] <- data.frame(term = term, class = cls, stringsAsFactors = FALSE)
  }
  map <- do.call(rbind, rows[seq_len(out_i)])
  map <- map[!duplicated(map$term), , drop = FALSE]
  cls <- map$class[match(terms, map$term)]
  cls[is.na(cls) | !nzchar(cls)] <- "unknown"
  stats::setNames(cls, terms)
}

.mgcfa_refit_with_group_partial <- function(base_fit, group_partial_terms) {
  gp <- .mgcfa_normalize_terms(group_partial_terms)
  tryCatch(
    stats::update(base_fit, group.partial = gp),
    error = function(e) NULL
  )
}
