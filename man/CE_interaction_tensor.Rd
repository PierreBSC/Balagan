% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CE_interaction_tensor.R
\name{CE_interaction_tensor}
\alias{CE_interaction_tensor}
\title{Compute the Clark-Evans index for each pair of cell type in each image and store them in a three-order tensor}
\usage{
CE_interaction_tensor(sce, type_output = "Index", Perform_symmetrization = T)
}
\arguments{
\item{sce}{a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.}
}
\value{
Returns three-order tensor with the log2 CE index of each pair of cell type in each image
}
\description{
Compute the Clark-Evans index for each pair of cell type in each image and store them in a three-order tensor
}
\examples{
CE_tensor = CE_interaction_tensor(sce)
}
