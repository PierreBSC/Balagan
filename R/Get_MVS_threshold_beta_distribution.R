#' @rdname Get_MVS_threshold_beta_distribution
#' @title Computing the Minimum Variance Stratification (MVS) thresholding for a beta distribution of known parameters.
#'
#' @description Subfunction used to compute the MVS thresholding for a theoretical beta distribution of known parameters. For more details please have a look the paper "Minimum Variance Stratification" by Tore Dalenius and Joseph L. Hodge in 1959
#'
#' @param x numeric vector to be stratified
#' @param L number of strata. Should be between 4 and 6 
#' @param distribution Probability distribution model used for the thresholding. So far only the beta distribution is available.
#' @return Returns a plot as well as tau and No estimate for each of the parameter value.
#' @export


Get_MVS_threshold_beta_distribution = function(L=6,alpha_parameter = 1, beta_parameter = 2) {
  #L : number of strata
  if (L<2){
    stop("There must be at least 2 strata !")
  }
  Vector_threshold = c()
  for (h in 1:(L-1)) {
    m = optimize(f = function(x) {(Rbeta(x = x,a = (alpha_parameter+1)/2,b=(beta_parameter+1)/2)-h/L)^2},lower = 0,upper = 1)
    Vector_threshold = c(Vector_threshold,m$minimum)
  }
  Vector_threshold = c(0,Vector_threshold,1)
  return(Vector_threshold)
}
