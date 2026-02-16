\name{mgcfa_build_model}
\alias{mgcfa_build_model}
\title{Build Lavaan Syntax for MGCFA}
\usage{
mgcfa_build_model(
  model_type = c("custom", "single_factor", "efa", "pca"),
  model = NULL,
  variables = NULL,
  data = NULL,
  input_matrix = NULL,
  matrix_nobs = NULL,
  n_factors = 1L,
  loading_threshold = 0.30,
  rotation = "oblimin",
  allow_cross_loadings = FALSE
)
}
\arguments{
\item{model_type}{One of \code{"custom"}, \code{"single_factor"},
\code{"efa"}, or \code{"pca"}.}

\item{model}{Lavaan syntax when \code{model_type = "custom"}.}

\item{variables}{Observed indicator names.}

\item{data}{Optional raw data used to build EFA/PCA syntax.}

\item{input_matrix}{Optional covariance/correlation matrix used to build
EFA/PCA syntax.}

\item{matrix_nobs}{Optional sample size for matrix-based EFA/PCA.}

\item{n_factors}{Number of factors/components for EFA/PCA.}

\item{loading_threshold}{Absolute loading threshold for retaining paths.}

\item{rotation}{Rotation passed to \code{psych::fa()} or
\code{psych::principal()}.}

\item{allow_cross_loadings}{If \code{TRUE}, keep all loadings above threshold.
If \code{FALSE}, keep the strongest loading per indicator.}
}
\value{
A character string containing lavaan syntax.
}
\description{
Creates SEM/CFA syntax for \code{mgcfa_auto()} from either user syntax,
a single-factor template, or thresholded exploratory loadings from EFA/PCA.
}
