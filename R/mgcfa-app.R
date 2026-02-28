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
          choices = c("chisq_pvalue", "delta_cfi", "aic_bic_weight", "measure_change", "none"),
          selected = "chisq_pvalue"
        ),
        shiny::numericInput("partial_failure_threshold", "Failure threshold (optional)", value = NA, step = 0.01),
        shiny::selectInput(
          "partial_search_criterion",
          "Partial-search criterion",
          choices = c("chisq_pvalue", "delta_cfi", "aic_bic_weight", "measure_change"),
          selected = "chisq_pvalue"
        ),
        shiny::numericInput("partial_search_threshold", "Partial-search threshold (optional)", value = NA, step = 0.01),

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
          partial_failure_criterion = input$partial_failure_criterion,
          partial_search_criterion = input$partial_search_criterion,
          stop_at_first_unacceptable = isTRUE(input$stop_early)
        )

        if (is.finite(input$partial_failure_threshold)) {
          args$partial_failure_threshold <- input$partial_failure_threshold
        }
        if (is.finite(input$partial_search_threshold)) {
          args$partial_search_threshold <- input$partial_search_threshold
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
