#' @rdname Strata_thresholding
#' @title Computing thresholding values for 
#'
#' @description Function that computes the optimal thresholding of a numeric vector into L strata using dufferent approaches including the 'Cumsqrt' or Minimum Variance Stratification method 
#' (for a detailed review see the paper "Stratification of Skewed Populations: A Review")
#'
#' @param x numeric vector to be stratified
#' @param L number of strata. Should be between 4 and 6 
#' @param type_thresholding how should stratum be computed based on the intensity of the covariate ? Can be "Cumsqrt" or "Geometric"
#' @return Returns a vector containing the L-1 thresholding values and produces a plot showing the threshold values.
#' @examples
#' x = rbeta(n = 1000,shape1 = 2,shape2 = 5)
#' MVS_thresholding(x,L=6)
#' @export


Strata_thresholding = function(x,L=6,type_thresholding="Cumsqrt") {
  
  density_temp = density(x)
  
  if (type_thresholding=="Cumsqrt") {

    cumsum_temp = cumsum(sqrt(density_temp$y))
    max_value_temp = max(cumsum_temp)
    Reached_values = seq(0,max_value_temp,length.out=L+1)
    
    Thresholding_values = c()
    
    for (k in 1:L){
      u = Reached_values[k+1]
      Thresholding_values = c(Thresholding_values,density_temp$x[min(which(cumsum_temp>=u))])
    }
    
    
  }
  
  if (type_thresholding=="Geometric") {
    
    Max_value = max(x)
    Min_value = 1
    Ratio_parameter = (Max_value/Min_value)^(1/L)
    Thresholding_values = Min_value * (Ratio_parameter^(1:L))
  }
  
  plot(density_temp,xlim=c(0,quantile(x,0.999)),main="",xlab="x",xaxs='i',yaxs='i')
  abline(v=Thresholding_values,lwd=0.8,lty=2,col="red")
  Thresholding_values[L] =max(x)
  Thresholding_values = c(0,Thresholding_values)

  return(Thresholding_values)
}
