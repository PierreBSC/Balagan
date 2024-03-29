% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Get_MVS_threshold_beta_distribution.R
\name{Get_MVS_threshold_beta_distribution}
\alias{Get_MVS_threshold_beta_distribution}
\title{Computing the Minimum Variance Stratification (MVS) thresholding for a beta distribution of known parameters.}
\usage{
Get_MVS_threshold_beta_distribution(
  L = 6,
  alpha_parameter = 1,
  beta_parameter = 2
)
}
\arguments{
\item{L}{number of strata. Should be between 4 and 6}

\item{x}{numeric vector to be stratified}

\item{distribution}{Probability distribution model used for the thresholding. So far only the beta distribution is available.}
}
\value{
Returns a plot as well as tau and No estimate for each of the parameter value.
}
\description{
Subfunction used to compute the MVS thresholding for a theoretical beta distribution of known parameters. For more details please have a look the paper "Minimum Variance Stratification" by Tore Dalenius and Joseph L. Hodge in 1959
}
