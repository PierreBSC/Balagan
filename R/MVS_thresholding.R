#' @rdname MVS_thresholding
#' @title Minimum Variance Stratification thresholding
#'
#' @description Function that computes the optimal thresholding of a numeric vector into L strata using the Minimum Variance Stratification (MVS) thresholding. For a complete description of the MVS thresholding strategy please have a look the paper "Minimum Variance Stratification" by Tore Dalenius and Joseph L. Hodge in 1959
#'
#' @param x numeric vector to be stratified
#' @param L number of strata. Should be between 4 and 6 
#' @param distribution Probability distribution model used for the thresholding. So far only the beta distribution is available.
#' @return Returns a vector containing the L-1 thresholding values and produces a plot showing the threshold values.
#' @examples
#' x = rbeta(n = 1000,shape1 = 2,shape2 = 5)
#' MVS_thresholding(x,L=6)
#' @import fitdistrplus 
#' @import zipfR
#' @export


MVS_thresholding = function(x,L=6,distribution="beta") {
  if (!distribution %in% c("beta")) {
    stop("Please select an appropriate distribution. Currently only the beta distribution is available. Other distributions such as the Gamma distribution will be implemented")
  }
  
  if (distribution=="beta") {
    x = x[x!=0]
    x = x[x!=1]
    
    fit_marker_distribution = fitdist(x,distr = "beta",method = "mme",start = list(shape1=5,shape2=5))
    alpha_parameter = coef(fit_marker_distribution)[1]
    beta_parameter = coef(fit_marker_distribution)[2]
    
    Vector_threshold = Get_MVS_threshold_beta_distribution(L = L,alpha_parameter,beta_parameter)
    hist(as.numeric(x),100,xaxs='i',yaxs='i',cex.lab=1.3,freq=F,xlab="Marker intensity",main="MVS thresholding")
    curve(dbeta(x = x,shape1 = coef(fit_marker_distribution)[1],shape2 = coef(fit_marker_distribution)[2]),add = T,lwd=2,from=0,to=1)
    abline(v=Vector_threshold,lwd=2,lty=2,col="red")
    
  }
  return(Vector_threshold)
}
