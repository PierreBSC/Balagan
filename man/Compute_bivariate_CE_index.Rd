% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Compute_bivariate_CE_index.R
\name{Compute_bivariate_CE_index}
\alias{Compute_bivariate_CE_index}
\title{Compute the Clark-Evans index for a pair of point patterns}
\usage{
Compute_bivariate_CE_index(X_pattern, Y_pattern, k_parameter = 1)
}
\arguments{
\item{X_pattern}{a ppp object.}

\item{Y_pattern}{a ppp object.}

\item{k_parameter}{a ppp object.}
}
\value{
Returns the bivariate Clark-Evans Index for this pair of point pattern.
}
\description{
Compute the Clark-Evans index for a pair of point patterns. Internal function called by the CE_interaction_tensor() function
}
\examples{
CE_interaction = Compute_bivariate_CE_index(X_pattern,Y_pattern,k = k_parameter)
}
