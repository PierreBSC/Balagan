% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Random_spatial_sampling.R
\name{Random_spatial_sampling}
\alias{Random_spatial_sampling}
\title{Random spatial sampling of square areas}
\usage{
Random_spatial_sampling(
  sce,
  width_FOV = 400,
  height_FOV = 400,
  N_samplings = 10,
  Selected_image = 1,
  plot_result = TRUE
)
}
\arguments{
\item{sce}{a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.}

\item{width_FOV}{height of the individual rectangles}

\item{height_FOV}{width of the individual rectangles}

\item{N_samplings}{number of individual rectangles samples}

\item{Selected_image}{of which image/FOV should the sampling be performed ?}
}
\value{
Returns a list with the fields "List_sampled_cells" describing which cells have been selected and the field "List_sampled_cluster" where their associated cluster/label have been sampled
}
\description{
Perform a random spatial sampling with a given number of rectangles of with a given width and height
}
\examples{
Simple_sampling = Random_spatial_sampling(sce,width_FOV=400,height_FOV=400,N_samplings=10,Selected_image = 1,plot_result=TRUE)
}
