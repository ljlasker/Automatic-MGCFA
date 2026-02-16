\name{mgcfa_tidy_fit}
\alias{mgcfa_tidy_fit}
\title{Tidy MGCFA Fit Results}
\usage{
mgcfa_tidy_fit(
  x,
  measures = c("chisq", "df", "cfi", "rmsea", "srmr", "aic", "bic"),
  include_non_partial = TRUE,
  baseline_step = "configural",
  add_delta = TRUE,
  digits = 3L,
  rounding = c("signif", "round", "none")
)
}
\arguments{
\item{x}{An object returned by \code{mgcfa_auto()}.}

\item{measures}{Fit measures to extract.}

\item{include_non_partial}{If \code{TRUE}, include failed non-partial models
for failed steps when available.}

\item{baseline_step}{Step used for delta calculations.}

\item{add_delta}{If \code{TRUE}, include changes from \code{baseline_step}.}

\item{digits}{Number of digits used for numeric formatting in the returned
table.}

\item{rounding}{Numeric rounding mode. One of \code{"signif"},
\code{"round"}, or \code{"none"}.}
}
\value{
A long-format \code{data.frame} with one row per
step/measure/model-variant.
}
\description{
Converts \code{mgcfa_result} objects into a tidy table for reporting and
plotting, including selected models and optional failed non-partial models.
}
