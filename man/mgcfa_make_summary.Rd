\name{mgcfa_make_summary}
\alias{mgcfa_make_summary}
\title{Build Summary Inputs from Raw Data}
\usage{
mgcfa_make_summary(
  data,
  group,
  variables,
  matrix_type = c("cor", "cov"),
  cor_digits = NULL,
  sd_digits = NULL,
  mean_digits = NULL,
  cov_digits = NULL,
  format = c("none", "truncate", "round"),
  drop_incomplete = TRUE
)
}
\arguments{
\item{data}{Data frame containing indicators and a grouping variable.}

\item{group}{Grouping variable name in \code{data}.}

\item{variables}{Character vector of indicator names.}

\item{matrix_type}{One of \code{"cor"} or \code{"cov"}.}

\item{cor_digits}{Decimal places for correlation matrices when
\code{matrix_type = "cor"}. Ignored when \code{NULL}.}

\item{sd_digits}{Decimal places for group SD vectors when
\code{matrix_type = "cor"}. Ignored when \code{NULL}.}

\item{mean_digits}{Decimal places for group means. Ignored when \code{NULL}.}

\item{cov_digits}{Decimal places for covariance matrices when
\code{matrix_type = "cov"}. Ignored when \code{NULL}.}

\item{format}{Decimal formatting method: \code{"none"},
\code{"truncate"}, or \code{"round"}.}

\item{drop_incomplete}{Logical; if \code{TRUE}, uses complete cases for
\code{variables} within each group.}
}
\value{
A list with \code{sample_cov}, \code{sample_mean}, \code{sample_nobs}, and
\code{group_labels}, plus \code{matrices_are_cor} and \code{sample_sd}
(\code{NULL} in covariance mode).
}
\description{
Creates group-wise summary inputs compatible with \code{mgcfa_auto()} from
raw data. Optionally applies paper-style decimal formatting (truncate or
round) to correlations/covariances, SDs, and means.
}
