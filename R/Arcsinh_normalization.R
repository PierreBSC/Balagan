#' @rdname Arcsinh_normalization
#' @title Transformation of raw intensity using arcsinh transform 
#'
#' @description This function transform the raw intensity 
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param compute_cofactor a boolean describing if a global lambda cofactor should be computed 
#' @param cofactor_value numeric value describing the arcsinh cofactor. If defined, overide the "compute_cofactor" argument.
#' @param trimming_value numeric value describing the quantile paramater used for trimming. Set to NULL by default.


#' @return Returns a \code{\link[SingleCellExperiment]{SingleCellExperiment}} object with a new assay slot called "Arcsinh_transformed_intensity"
#'
#' @examples
#' sce = Arcsinh_normalization(sce,cofactor_value = 50)
#' @import SingleCellExperiment
#' @export

Arcsinh_normalization = function(sce,compute_cofactor=TRUE, cofactor_value = NULL) {

  if (!compute_cofactor & is.null(cofactor_value)) {
    stop("Please provide an adapted strategy : either compute lambda or provide a value !")
  }
  
  if (compute_cofactor) {
    cat("Computing the optimal lambda parameter...")
    #Computing the relation between mean and variance to optimize the transformation
    Mean_channel_intensity = apply(assays(sce)[["Raw_intensity"]],MARGIN = 1,FUN = mean)
    Mean_channel_intensity_squared = Mean_channel_intensity^2
    Variance_channel_intensity = apply(assays(sce)[["Raw_intensity"]],MARGIN = 1,FUN = var)
    
    nonlin_mod = nls(Variance_channel_intensity ~ Mean_channel_intensity + lambda*Mean_channel_intensity_squared ,start=list(lambda=0))
    Lambda_parameter = coef(nonlin_mod)
    par(las=1,bty="l")
    plot(log10(Mean_channel_intensity),log10(Variance_channel_intensity),xlab="Mean intensity (log10)",
         ylab="Intensity variance (log10)",pch=21,bg='red3',cex=1.3,cex.lab=1.2)
    m = lm(log10(Variance_channel_intensity)~log10(Mean_channel_intensity))
    abline(coef(m),lwd=2,lty=2)
    cat(" done ! \n")
  }
  
  if (!is.null(cofactor_value)) {
    Lambda_parameter = cofactor_value
  }
  cat(paste("Lambda parameter value :",as.character(Lambda_parameter)))
  #For justification : see the excellent paper from Bartlett "The Use of Transformations"
  Transformed_data = 1/Lambda_parameter * asinh(Lambda_parameter*sqrt(assays(sce)[["Raw_intensity"]]))
  
  #Adding the matrix as a new assay slot 
  assay(sce, "Arcsinh_transformed_intensity") <- (Transformed_data)
  
  return(sce)
}

