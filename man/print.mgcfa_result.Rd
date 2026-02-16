\name{print.mgcfa_result}
\alias{print.mgcfa_result}
\title{Print Method for MGCFA Results}
\usage{
\method{print}{mgcfa_result}(x, digits = 3L, rounding = c("signif", "round", "none"), show_freed_parameters = FALSE, verbose = FALSE, ...)
}
\arguments{
\item{x}{An object returned by \code{mgcfa_auto()}.}

\item{digits}{Number of digits used for numeric formatting.}

\item{rounding}{Numeric rounding mode. One of \code{"signif"},
\code{"round"}, or \code{"none"}.}

\item{show_freed_parameters}{If \code{TRUE}, print the constraints freed by
automatic partial-model discovery for each recovered step.}

\item{verbose}{If \code{TRUE}, print full diagnostic detail in addition to the
simplified overview table.}

\item{...}{Unused.}
}
\value{
Returns \code{x} invisibly.
}
\description{
Prints a compact summary of mode, model type, fitted invariance steps,
fit statistics, failed non-partial outputs for failed constrained steps,
and optional automatic partial-search diagnostics.
}
