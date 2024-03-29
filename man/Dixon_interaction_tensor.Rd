% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Dixon_interaction_tensor.R
\name{Dixon_interaction_tensor}
\alias{Dixon_interaction_tensor}
\title{Compute cell type segregation using Dixon's approach}
\usage{
Dixon_interaction_tensor(sce, type_output = "Z-score")
}
\arguments{
\item{sce}{a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.}
}
\value{
Returns three-order tensor with the log-odds ratio of each pair of cell type
}
\description{
Compute cell type segregation using Dixon's approach described in papers by Dixon (1994 and 2002)
}
\examples{
Dixon_tensor = Dixon_interaction_tensor(sce)
}
