#' Launch the Interactive AutomaticMGCFA App
#'
#' Opens a Shiny-based interface for running `mgcfa_auto()` with guided inputs,
#' diagnostics, stage summaries, and fit plots.
#'
#' The app supports two workflows:
#' - Raw data upload (`.csv` or `.rds` data.frame), with group/variable selection.
#' - Summary input upload (`.rds` list), either the direct argument list for
#'   `mgcfa_auto()` matrix mode or an object returned by `mgcfa_make_summary()`.
#'
#' @return Invisibly returns the app object after launching.
#' @export
mgcfa_launch_app <- function() {
  app <- mgcfa_app()
  shiny::runApp(app)
  invisible(app)
}

#' Build the Interactive AutomaticMGCFA App
#'
#' Creates and returns the Shiny app object used by
#' `mgcfa_launch_app()`.
#'
#' @return A `shiny.appobj`.
#' @export
mgcfa_app <- function() {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package `shiny` is required for the interactive app. Install it with install.packages('shiny').", call. = FALSE)
  }

  step_choices <- c(
    "configural", "metric", "scalar", "strict",
    "lv.variances", "lv.covariances", "residual.covariances", "regressions", "means"
  )

  ui <- shiny::fluidPage(
    shiny::titlePanel("AutomaticMGCFA Interactive Runner"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::h4("Input Type"),
        shiny::radioButtons(
          "data_mode",
          label = NULL,
          choices = c("Raw data" = "raw", "Summary stats (RDS list)" = "summary"),
          selected = "raw"
        ),

        shiny::conditionalPanel(
          "input.data_mode == 'raw'",
          shiny::fileInput("raw_file", "Upload Raw Data (.csv, .rds)", accept = c(".csv", ".rds")),
          shiny::checkboxInput("raw_header", "CSV has header", value = TRUE),
          shiny::selectInput("raw_sep", "CSV separator", choices = c("," = ",", "Semicolon" = ";", "Tab" = "\t"
          )),
          shiny::selectInput("group_var", "Group variable", choices = character()),
          shiny::uiOutput("raw_vars_ui")
        ),

        shiny::conditionalPanel(
          "input.data_mode == 'summary'",
          shiny::fileInput("summary_file", "Upload Summary Object (.rds)", accept = c(".rds")),
          shiny::helpText("Expected fields: sample_cov, sample_nobs, optional sample_mean, sample_sd, group_labels, matrices_are_cor.")
        ),

        shiny::hr(),
        shiny::h4("Model"),
        shiny::selectInput(
          "model_type",
          "Model type",
          choices = c("Custom" = "custom", "Single factor" = "single_factor", "EFA-based" = "efa", "PCA-based" = "pca"),
          selected = "custom"
        ),
        shiny::conditionalPanel(
          "input.model_type == 'custom'",
          shiny::textAreaInput(
            "model_syntax",
            "Lavaan model syntax",
            value = "F =~ x1 + x2 + x3 + x4",
            rows = 5
          )
        ),
        shiny::conditionalPanel(
          "input.model_type != 'custom'",
          shiny::numericInput("n_factors", "Number of factors/components", value = 1, min = 1, step = 1),
          shiny::numericInput("loading_threshold", "Loading threshold", value = 0.30, min = 0, max = 1, step = 0.01),
          shiny::textInput("rotation", "Rotation", value = "oblimin"),
          shiny::checkboxInput("allow_cross_loadings", "Allow cross-loadings above threshold", value = FALSE)
        ),

        shiny::hr(),
        shiny::h4("Invariance + Search"),
        shiny::checkboxGroupInput(
          "include_steps",
          "Stages to run",
          choices = step_choices,
          selected = step_choices
        ),
        shiny::checkboxInput("stop_early", "Stop at first unacceptable stage", value = TRUE),
        shiny::selectInput(
          "partial_auto_search",
          "Automatic partial search",
          choices = c("always", "never", "prompt"),
          selected = "always"
        ),
        shiny::selectInput(
          "partial_failure_criterion",
          "Failure criterion",
          choices = c(
            "chisq_pvalue",
            "delta_cfi",
            "aic_weight",
            "bic_weight",
            "aic_bic_weight",
            "measure_change",
            "none"
          ),
          selected = "chisq_pvalue"
        ),
        shiny::checkboxInput("use_failure_rules", "Use multi-rule failure decision", value = FALSE),
        shiny::conditionalPanel(
          "input.use_failure_rules == false",
          shiny::numericInput("partial_failure_threshold", "Failure threshold (optional)", value = NA, step = 0.01),
          shiny::conditionalPanel(
            "input.partial_failure_criterion == 'measure_change'",
            shiny::textInput(
              "partial_failure_measure",
              "Failure fit measure (lavaan::fitMeasures name)",
              value = "aic",
              placeholder = "e.g., aic, bic, cfi, tli, rmsea, srmr, chisq"
            ),
            shiny::selectInput(
              "partial_failure_direction",
              "Failure direction",
              choices = c("decrease", "increase"),
              selected = "decrease"
            )
          )
        ),
        shiny::conditionalPanel(
          "input.use_failure_rules == true",
          shiny::numericInput("failure_n_rules", "Number of failure rules", value = 3, min = 1, max = 8, step = 1),
          shiny::uiOutput("failure_rules_ui"),
          shiny::selectInput(
            "failure_rule_policy",
            "Failure rule policy",
            choices = c("all", "majority", "any", "at_least"),
            selected = "majority"
          ),
          shiny::conditionalPanel(
            "input.failure_rule_policy == 'at_least'",
            shiny::numericInput("failure_rule_min", "Failure minimum passing rules", value = 2, min = 1, step = 1)
          )
        ),
        shiny::selectInput(
          "partial_search_criterion",
          "Partial-search criterion",
          choices = c(
            "chisq_pvalue",
            "delta_cfi",
            "aic_weight",
            "bic_weight",
            "aic_bic_weight",
            "measure_change"
          ),
          selected = "chisq_pvalue"
        ),
        shiny::checkboxInput("use_search_rules", "Use multi-rule partial-search decision", value = FALSE),
        shiny::conditionalPanel(
          "input.use_search_rules == false",
          shiny::numericInput("partial_search_threshold", "Partial-search threshold (optional)", value = NA, step = 0.01),
          shiny::conditionalPanel(
            "input.partial_search_criterion == 'measure_change'",
            shiny::textInput(
              "partial_search_measure",
              "Search fit measure (lavaan::fitMeasures name)",
              value = "aic",
              placeholder = "e.g., aic, bic, cfi, tli, rmsea, srmr, chisq"
            ),
            shiny::selectInput(
              "partial_search_direction",
              "Search direction",
              choices = c("decrease", "increase"),
              selected = "decrease"
            )
          )
        ),
        shiny::conditionalPanel(
          "input.use_search_rules == true",
          shiny::numericInput("search_n_rules", "Number of partial-search rules", value = 3, min = 1, max = 8, step = 1),
          shiny::uiOutput("search_rules_ui"),
          shiny::selectInput(
            "search_rule_policy",
            "Partial-search rule policy",
            choices = c("all", "majority", "any", "at_least"),
            selected = "majority"
          ),
          shiny::conditionalPanel(
            "input.search_rule_policy == 'at_least'",
            shiny::numericInput("search_rule_min", "Partial-search minimum passing rules", value = 2, min = 1, step = 1)
          )
        ),
        shiny::numericInput(
          "partial_ic_bic_weight",
          "BIC weight for AIC/BIC-weight criteria (AIC = 1 - BIC)",
          value = 0.50,
          min = 0,
          max = 1,
          step = 0.05
        ),
        shiny::helpText(
          "For AIC/BIC-weight criteria, the candidate model is compared to the previous model using pairwise model weights.",
          "AIC-only = aic_weight, BIC-only = bic_weight, or use aic_bic_weight with custom BIC weight.",
          "For measure_change, any valid lavaan fit measure name can be used (including those provided by installed extensions)."
        ),

        shiny::actionButton("run_mgcfa", "Run MGCFA", class = "btn-primary")
      ),
      shiny::mainPanel(
        shiny::tabsetPanel(
          shiny::tabPanel(
            "Input Preview",
            shiny::h4("Current Input Summary"),
            shiny::verbatimTextOutput("input_summary_text"),
            shiny::tableOutput("input_preview_table")
          ),
          shiny::tabPanel(
            "Results",
            shiny::h4("Run Status"),
            shiny::verbatimTextOutput("run_status_text"),
            shiny::h4("Stage Overview"),
            shiny::tableOutput("overview_table"),
            shiny::h4("Failures / Not Applicable"),
            shiny::tableOutput("failures_table"),
            shiny::h4("Decision Trace"),
            shiny::tableOutput("decision_trace_table")
          ),
          shiny::tabPanel(
            "Plots",
            shiny::fluidRow(
              shiny::column(
                width = 4,
                shiny::selectInput("plot_measure", "Fit measure", choices = c("aic", "bic", "cfi", "rmsea", "srmr", "chisq"), selected = "aic")
              ),
              shiny::column(
                width = 4,
                shiny::radioButtons("plot_type", "Plot type", choices = c("Delta" = "delta", "Raw" = "raw"), inline = TRUE)
              )
            ),
            shiny::plotOutput("fit_plot", height = "520px")
          )
        )
      )
    )
  )

  server <- function(input, output, session) {
    rv <- shiny::reactiveValues(result = NULL, report = NULL, error = NULL, warnings = character())

    raw_data <- shiny::reactive({
      shiny::req(input$data_mode == "raw")
      shiny::req(!is.null(input$raw_file))
      .mgcfa_app_read_raw_data(
        path = input$raw_file$datapath,
        header = isTRUE(input$raw_header),
        sep = input$raw_sep %||% ","
      )
    })

    summary_args <- shiny::reactive({
      shiny::req(input$data_mode == "summary")
      shiny::req(!is.null(input$summary_file))
      obj <- readRDS(input$summary_file$datapath)
      .mgcfa_app_unpack_summary_input(obj)
    })

    output$failure_rules_ui <- shiny::renderUI({
      shiny::req(isTRUE(input$use_failure_rules))
      n <- as.integer(input$failure_n_rules %||% 1L)
      n <- max(1L, min(8L, if (is.na(n)) 1L else n))
      rows <- lapply(seq_len(n), function(i) .mgcfa_app_rule_row_ui("failure", i))
      do.call(shiny::tagList, rows)
    })

    output$search_rules_ui <- shiny::renderUI({
      shiny::req(isTRUE(input$use_search_rules))
      n <- as.integer(input$search_n_rules %||% 1L)
      n <- max(1L, min(8L, if (is.na(n)) 1L else n))
      rows <- lapply(seq_len(n), function(i) .mgcfa_app_rule_row_ui("search", i))
      do.call(shiny::tagList, rows)
    })

    shiny::observeEvent(raw_data(), {
      dat <- raw_data()
      shiny::updateSelectInput(session, "group_var", choices = names(dat), selected = names(dat)[[1L]])
    }, ignoreInit = TRUE)

    output$raw_vars_ui <- shiny::renderUI({
      shiny::req(input$data_mode == "raw")
      dat <- raw_data()
      grp <- input$group_var %||% ""
      cols <- setdiff(names(dat), grp)
      shiny::selectInput("raw_variables", "Indicators (for non-custom model types)", choices = cols, selected = cols, multiple = TRUE)
    })

    output$input_summary_text <- shiny::renderPrint({
      if (identical(input$data_mode, "raw")) {
        dat <- raw_data()
        cat("Raw data loaded\n")
        cat("Rows:", nrow(dat), " Columns:", ncol(dat), "\n")
        if (!is.null(input$group_var) && nzchar(input$group_var)) {
          cat("Group variable:", input$group_var, "\n")
        }
      } else {
        s <- summary_args()
        cat("Summary input loaded\n")
        cat("Groups:", length(s$sample_cov), " Variables:", ncol(s$sample_cov[[1L]]), "\n")
        cat("Matrices are correlations:", isTRUE(s$matrices_are_cor), "\n")
        cat("Means supplied:", !is.null(s$sample_mean), "\n")
      }
    })

    output$input_preview_table <- shiny::renderTable({
      if (identical(input$data_mode, "raw")) {
        utils::head(raw_data(), 10)
      } else {
        s <- summary_args()
        as.data.frame(round(s$sample_cov[[1L]], 3))
      }
    }, rownames = TRUE)

    shiny::observeEvent(input$run_mgcfa, {
      rv$error <- NULL
      rv$warnings <- character()
      rv$result <- NULL
      rv$report <- NULL

      fit_or_err <- tryCatch({
        args <- list(
          model_type = input$model_type,
          include_steps = input$include_steps,
          partial_auto_search = input$partial_auto_search,
          partial_failure_criterion = "chisq_pvalue",
          partial_search_criterion = "chisq_pvalue",
          stop_at_first_unacceptable = isTRUE(input$stop_early)
        )
        if (isTRUE(input$use_failure_rules)) {
          n_fail <- as.integer(input$failure_n_rules %||% 1L)
          n_fail <- max(1L, min(8L, if (is.na(n_fail)) 1L else n_fail))
          fail_criteria <- vapply(seq_len(n_fail), function(i) as.character(input[[paste0("failure_rule_criterion_", i)]] %||% "chisq_pvalue"), character(1L))
          fail_thresholds <- vapply(seq_len(n_fail), function(i) as.numeric(input[[paste0("failure_rule_threshold_", i)]] %||% NA_real_), numeric(1L))
          fail_measures <- vapply(seq_len(n_fail), function(i) as.character(input[[paste0("failure_rule_measure_", i)]] %||% "aic"), character(1L))
          fail_dirs <- vapply(seq_len(n_fail), function(i) as.character(input[[paste0("failure_rule_direction_", i)]] %||% "decrease"), character(1L))
          fail_w <- vapply(seq_len(n_fail), function(i) as.numeric(input[[paste0("failure_rule_ic_bic_weight_", i)]] %||% input$partial_ic_bic_weight %||% 0.5), numeric(1L))
          args$partial_failure_rules <- .mgcfa_app_build_rules(
            criteria = fail_criteria,
            thresholds = fail_thresholds,
            measures = fail_measures,
            directions = fail_dirs,
            ic_bic_weights = fail_w
          )
          args$partial_failure_rule_policy <- as.character(input$failure_rule_policy %||% "majority")
          if (identical(args$partial_failure_rule_policy, "at_least")) {
            min_fail <- as.integer(input$failure_rule_min %||% NA_integer_)
            if (!is.na(min_fail) && min_fail > 0L) {
              args$partial_failure_rule_min <- min_fail
            }
          }
        } else {
          fail_cfg <- .mgcfa_app_parse_criterion(
            criterion = input$partial_failure_criterion,
            threshold = input$partial_failure_threshold,
            measure = input$partial_failure_measure %||% "aic",
            direction = input$partial_failure_direction %||% "decrease",
            ic_bic_weight = input$partial_ic_bic_weight,
            allow_none = TRUE
          )
          args$partial_failure_criterion <- fail_cfg$criterion
          if (!is.null(fail_cfg$threshold)) {
            args$partial_failure_threshold <- fail_cfg$threshold
          }
          if (!is.null(fail_cfg$measure)) {
            args$partial_failure_measure <- fail_cfg$measure
            args$partial_failure_direction <- fail_cfg$direction
          }
          if (!is.null(fail_cfg$ic_bic_weight)) {
            args$partial_ic_bic_weight <- fail_cfg$ic_bic_weight
          }
        }

        if (isTRUE(input$use_search_rules)) {
          n_search <- as.integer(input$search_n_rules %||% 1L)
          n_search <- max(1L, min(8L, if (is.na(n_search)) 1L else n_search))
          search_criteria <- vapply(seq_len(n_search), function(i) as.character(input[[paste0("search_rule_criterion_", i)]] %||% "chisq_pvalue"), character(1L))
          search_thresholds <- vapply(seq_len(n_search), function(i) as.numeric(input[[paste0("search_rule_threshold_", i)]] %||% NA_real_), numeric(1L))
          search_measures <- vapply(seq_len(n_search), function(i) as.character(input[[paste0("search_rule_measure_", i)]] %||% "aic"), character(1L))
          search_dirs <- vapply(seq_len(n_search), function(i) as.character(input[[paste0("search_rule_direction_", i)]] %||% "decrease"), character(1L))
          search_w <- vapply(seq_len(n_search), function(i) as.numeric(input[[paste0("search_rule_ic_bic_weight_", i)]] %||% input$partial_ic_bic_weight %||% 0.5), numeric(1L))
          args$partial_search_rules <- .mgcfa_app_build_rules(
            criteria = search_criteria,
            thresholds = search_thresholds,
            measures = search_measures,
            directions = search_dirs,
            ic_bic_weights = search_w
          )
          args$partial_search_rule_policy <- as.character(input$search_rule_policy %||% "majority")
          if (identical(args$partial_search_rule_policy, "at_least")) {
            min_search <- as.integer(input$search_rule_min %||% NA_integer_)
            if (!is.na(min_search) && min_search > 0L) {
              args$partial_search_rule_min <- min_search
            }
          }
        } else {
          search_cfg <- .mgcfa_app_parse_criterion(
            criterion = input$partial_search_criterion,
            threshold = input$partial_search_threshold,
            measure = input$partial_search_measure %||% input$partial_failure_measure %||% "aic",
            direction = input$partial_search_direction %||% input$partial_failure_direction %||% "decrease",
            ic_bic_weight = input$partial_ic_bic_weight,
            allow_none = FALSE
          )
          args$partial_search_criterion <- search_cfg$criterion
          if (!is.null(search_cfg$threshold)) {
            args$partial_search_threshold <- search_cfg$threshold
          }
          if (!is.null(search_cfg$measure)) {
            args$partial_search_measure <- search_cfg$measure
            args$partial_search_direction <- search_cfg$direction
          }
          if (!is.null(search_cfg$ic_bic_weight)) {
            args$partial_ic_bic_weight <- search_cfg$ic_bic_weight
          }
        }

        if (identical(input$model_type, "custom")) {
          model_txt <- trimws(input$model_syntax %||% "")
          if (!nzchar(model_txt)) {
            stop("Custom model syntax is required.", call. = FALSE)
          }
          args$model <- model_txt
        } else {
          args$model <- NULL
          args$n_factors <- as.integer(input$n_factors)
          args$loading_threshold <- as.numeric(input$loading_threshold)
          args$rotation <- as.character(input$rotation)
          args$allow_cross_loadings <- isTRUE(input$allow_cross_loadings)
          if (identical(input$data_mode, "raw")) {
            vars <- as.character(input$raw_variables %||% character())
            if (length(vars) < 2L) {
              stop("Select at least two indicator variables for non-custom model types.", call. = FALSE)
            }
            args$variables <- vars
          }
        }

        if (identical(input$data_mode, "raw")) {
          dat <- raw_data()
          grp <- as.character(input$group_var %||% "")
          if (!nzchar(grp) || !(grp %in% names(dat))) {
            stop("Choose a valid group variable from the raw data.", call. = FALSE)
          }
          args$data <- dat
          args$group <- grp
        } else {
          s <- summary_args()
          args$sample_cov <- s$sample_cov
          args$sample_nobs <- s$sample_nobs
          args$sample_mean <- s$sample_mean
          args$sample_sd <- s$sample_sd
          args$group_labels <- s$group_labels
          args$matrices_are_cor <- isTRUE(s$matrices_are_cor)
        }

        warn_buf <- character()
        fit <- withCallingHandlers(
          do.call(mgcfa_auto, args),
          warning = function(w) {
            warn_buf <<- c(warn_buf, conditionMessage(w))
            tryInvokeRestart("muffleWarning")
          }
        )
        rv$warnings <- unique(as.character(warn_buf))
        fit
      }, error = function(e) e)

      if (inherits(fit_or_err, "error")) {
        rv$error <- conditionMessage(fit_or_err)
        shiny::showNotification(rv$error, type = "error", duration = NULL)
      } else {
        rv$result <- fit_or_err
        rv$report <- mgcfa_report(fit_or_err, include_plots = FALSE)
        shiny::showNotification("MGCFA run completed.", type = "message", duration = 3)
      }
    })

    output$run_status_text <- shiny::renderPrint({
      if (!is.null(rv$error)) {
        cat("Run failed:\n")
        cat(rv$error, "\n")
        return(invisible(NULL))
      }
      shiny::req(rv$result)
      out <- rv$result
      cat("Run successful\n")
      cat("Fitted stages:", paste(names(out$fits), collapse = ", "), "\n")
      if (length(out$recovered_steps %||% character()) > 0L) {
        cat("Recovered with partial invariance:", paste(out$recovered_steps, collapse = ", "), "\n")
      }
      failed <- names(Filter(function(z) isTRUE(z$failed), out$step_failures %||% list()))
      if (length(failed) > 0L) {
        cat("Failed non-partial stages:", paste(failed, collapse = ", "), "\n")
      }
      if (length(rv$warnings) > 0L) {
        cat("Warnings:\n")
        for (w in rv$warnings) cat("-", w, "\n")
      }
    })

    output$overview_table <- shiny::renderTable({
      shiny::req(rv$report)
      rv$report$overview
    }, rownames = FALSE)

    output$failures_table <- shiny::renderTable({
      shiny::req(rv$report)
      rv$report$failures
    }, rownames = FALSE)

    output$decision_trace_table <- shiny::renderTable({
      shiny::req(rv$report)
      rv$report$decision_trace
    }, rownames = FALSE)

    output$fit_plot <- shiny::renderPlot({
      shiny::req(rv$result)
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Install ggplot2 to show plots.")
        return(invisible(NULL))
      }
      p <- mgcfa_plot_fit(
        rv$result,
        measures = input$plot_measure,
        include_non_partial = TRUE,
        plot_type = input$plot_type
      )
      print(p)
    })
  }

  shiny::shinyApp(ui = ui, server = server)
}

.mgcfa_app_read_raw_data <- function(path, header = TRUE, sep = ",") {
  ext <- tolower(tools::file_ext(path %||% ""))
  if (identical(ext, "csv")) {
    dat <- utils::read.csv(path, header = isTRUE(header), sep = sep, stringsAsFactors = FALSE)
  } else if (identical(ext, "rds")) {
    dat <- readRDS(path)
  } else {
    stop("Raw input must be a .csv or .rds file.", call. = FALSE)
  }
  if (!is.data.frame(dat)) {
    stop("Raw upload must evaluate to a data.frame.", call. = FALSE)
  }
  dat
}

.mgcfa_app_unpack_summary_input <- function(obj) {
  if (!is.list(obj)) {
    stop("Summary input .rds must contain a list.", call. = FALSE)
  }

  args <- obj$mgcfa_args %||% obj

  if (!is.null(args$sample_cor) && is.null(args$sample_cov)) {
    args$sample_cov <- args$sample_cor
    args$matrices_are_cor <- TRUE
  }

  need <- c("sample_cov", "sample_nobs")
  miss <- need[!vapply(need, function(nm) !is.null(args[[nm]]), logical(1L))]
  if (length(miss) > 0L) {
    stop("Summary input is missing required fields: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  list(
    sample_cov = args$sample_cov,
    sample_mean = args$sample_mean %||% NULL,
    sample_nobs = args$sample_nobs,
    sample_sd = args$sample_sd %||% NULL,
    group_labels = args$group_labels %||% NULL,
    matrices_are_cor = isTRUE(args$matrices_are_cor)
  )
}

.mgcfa_app_rule_row_ui <- function(prefix, i) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package `shiny` is required.", call. = FALSE)
  }
  id <- function(stem) paste0(prefix, "_rule_", stem, "_", i)
  shiny::tagList(
    shiny::h5(sprintf("Rule %d", i)),
    shiny::selectInput(
      id("criterion"),
      "Criterion",
      choices = c(
        "chisq_pvalue",
        "delta_cfi",
        "aic_weight",
        "bic_weight",
        "aic_bic_weight",
        "measure_change"
      ),
      selected = "chisq_pvalue"
    ),
    shiny::numericInput(id("threshold"), "Threshold (optional)", value = NA, step = 0.01),
    shiny::textInput(
      id("measure"),
      "Fit measure (used by measure_change)",
      value = "aic",
      placeholder = "e.g., aic, bic, cfi, tli, rmsea, srmr, chisq"
    ),
    shiny::selectInput(
      id("direction"),
      "Direction (used by measure_change)",
      choices = c("decrease", "increase"),
      selected = "decrease"
    ),
    shiny::numericInput(
      id("ic_bic_weight"),
      "BIC weight (used by aic_bic_weight)",
      value = 0.50,
      min = 0,
      max = 1,
      step = 0.05
    ),
    shiny::hr()
  )
}

.mgcfa_app_build_rules <- function(
  criteria,
  thresholds = NULL,
  measures = NULL,
  directions = NULL,
  ic_bic_weights = NULL
) {
  criteria <- as.character(criteria %||% character())
  n <- length(criteria)
  if (n < 1L) {
    stop("At least one rule is required.", call. = FALSE)
  }

  expand <- function(x, default) {
    if (is.null(x)) {
      return(rep(default, n))
    }
    if (length(x) == 1L) {
      return(rep(x, n))
    }
    if (length(x) != n) {
      stop("Rule vectors must be length 1 or match the number of criteria.", call. = FALSE)
    }
    x
  }

  thresholds <- expand(thresholds, NA_real_)
  measures <- expand(measures, "aic")
  directions <- expand(directions, "decrease")
  ic_bic_weights <- expand(ic_bic_weights, 0.5)

  out <- vector("list", n)
  for (i in seq_len(n)) {
    cfg <- .mgcfa_app_parse_criterion(
      criterion = criteria[[i]],
      threshold = thresholds[[i]],
      measure = measures[[i]],
      direction = directions[[i]],
      ic_bic_weight = ic_bic_weights[[i]],
      allow_none = FALSE
    )
    rule <- list(criterion = cfg$criterion)
    if (!is.null(cfg$threshold)) {
      rule$threshold <- cfg$threshold
    }
    if (!is.null(cfg$measure)) {
      rule$measure <- cfg$measure
      rule$direction <- cfg$direction
    }
    if (!is.null(cfg$ic_bic_weight)) {
      rule$ic_bic_weight <- cfg$ic_bic_weight
    }
    out[[i]] <- rule
  }
  out
}

.mgcfa_app_parse_criterion <- function(
  criterion,
  threshold = NULL,
  measure = "aic",
  direction = c("decrease", "increase"),
  ic_bic_weight = 0.5,
  allow_none = FALSE
) {
  criterion <- as.character(criterion %||% "")
  direction <- match.arg(as.character(direction)[[1L]], choices = c("decrease", "increase"))

  if (identical(criterion, "none")) {
    if (!isTRUE(allow_none)) {
      stop("`none` is not allowed for this criterion slot.", call. = FALSE)
    }
    return(list(
      criterion = "none",
      threshold = NULL,
      measure = NULL,
      direction = NULL,
      ic_bic_weight = NULL
    ))
  }

  if (!criterion %in% c("chisq_pvalue", "delta_cfi", "aic_weight", "bic_weight", "aic_bic_weight", "measure_change")) {
    stop("Unsupported criterion: ", criterion, call. = FALSE)
  }

  if (!is.null(threshold) && is.finite(threshold)) {
    threshold_out <- as.numeric(threshold)
  } else {
    threshold_out <- NULL
  }

  if (identical(criterion, "aic_weight")) {
    return(list(
      criterion = "aic_bic_weight",
      threshold = threshold_out,
      measure = NULL,
      direction = NULL,
      ic_bic_weight = 0
    ))
  }
  if (identical(criterion, "bic_weight")) {
    return(list(
      criterion = "aic_bic_weight",
      threshold = threshold_out,
      measure = NULL,
      direction = NULL,
      ic_bic_weight = 1
    ))
  }
  if (identical(criterion, "aic_bic_weight")) {
    w <- as.numeric(ic_bic_weight)[[1L]]
    if (!is.finite(w) || w < 0 || w > 1) {
      stop("BIC weight must be in [0, 1].", call. = FALSE)
    }
    return(list(
      criterion = "aic_bic_weight",
      threshold = threshold_out,
      measure = NULL,
      direction = NULL,
      ic_bic_weight = w
    ))
  }
  if (identical(criterion, "measure_change")) {
    m <- trimws(as.character(measure %||% ""))
    if (!nzchar(m)) {
      stop("A fit measure name is required for `measure_change`.", call. = FALSE)
    }
    return(list(
      criterion = "measure_change",
      threshold = threshold_out,
      measure = m,
      direction = direction,
      ic_bic_weight = NULL
    ))
  }

  list(
    criterion = criterion,
    threshold = threshold_out,
    measure = NULL,
    direction = NULL,
    ic_bic_weight = NULL
  )
}
