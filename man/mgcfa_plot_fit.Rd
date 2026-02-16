\name{mgcfa_plot_fit}
\alias{mgcfa_plot_fit}
\title{Plot MGCFA Fit Evolution}
\usage{
mgcfa_plot_fit(
  x,
  measures = c("aic"),
  include_non_partial = TRUE,
  baseline_step = "configural",
  plot_type = c("delta", "raw"),
  color_values = c(selected = "#1768AC", non_partial = "#B9473E"),
  linetype_values = c(selected = "solid", non_partial = "dashed"),
  point_size = 2.2,
  line_size = 0.9,
  facet_scales = "free_y",
  facet_ncol = NULL,
  title = NULL
)
}
\arguments{
\item{x}{An object returned by \code{mgcfa_auto()}.}

\item{measures}{Fit measures to plot.}

\item{include_non_partial}{Include failed non-partial models as dashed lines.}

\item{baseline_step}{Step used for delta plotting.}

\item{plot_type}{One of \code{"delta"} or \code{"raw"}.}

\item{color_values}{Named colors for \code{"selected"} and
\code{"non_partial"}. You may also provide \code{"selected_mean"}.}

\item{linetype_values}{Named linetypes for \code{"selected"} and
\code{"non_partial"}. You may also provide \code{"selected_mean"}.}

\item{point_size}{Point size.}

\item{line_size}{Line width.}

\item{facet_scales}{Facet scale behavior passed to
\code{ggplot2::facet_wrap()}.}

\item{facet_ncol}{Optional facet column count.}

\item{title}{Optional custom plot title. If \code{NULL}, a default title is
generated from the selected fit measure.}
}
\value{
A \code{ggplot} object.
}
\description{
Creates ggplot2-based visuals for fit-index evolution across invariance steps.
Selected models are shown with solid lines and failed non-partial models with
dashed lines. Default labels use title-cased invariance stages and
measure-aware axis text (for example Change in AIC From Configural).
}
