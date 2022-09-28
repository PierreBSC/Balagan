#' @rdname Gamma_transformation
#' @title Normalisation and scaling of marker expression based on Wilson-Hilferty transform
#'
#' @description This function normalize and transform the intensity using the Wilson-Hilferty transform (cube root transform) and then scale the intensity
#'
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param clip_quantile the percentile used to clip intensity values.
#' @return Returns a \code{\link[SingleCellExperiment]{SingleCellExperiment}} object with a new assay slot called "Count_normalised_intensity"
#'
#' @examples
#' sce = Gamma_transformation(sce,clip_quantile=0.001)
#' @import SingleCellExperiment
#' @export


Gamma_transformation = function(sce,clip_quantile=NULL) {
  cat("Transforming values... \n")
  
  Transformed_data = t(assays(sce)[["Raw_intensity"]])
  Transformed_data = Transformed_data^(0.33)
  
  if (!is.null(clip_quantile)) {
    cat("Scaling/clipping values.. \n")
    Transformed_data = apply(Transformed_data,MARGIN = 2,FUN = function(x) {
      x_min = quantile(x,clip_quantile)
      x_max = quantile(x,1-clip_quantile)
      x[x<x_min] = x_min
      x[x>x_max] = x_max
      return(x)
    })
  }
  assay(sce, "Gamma_transformed_intensity",withDimnames = FALSE) <- t(Transformed_data)
  return(sce)
}
  