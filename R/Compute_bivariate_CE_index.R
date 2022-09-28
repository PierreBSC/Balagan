#' @rdname Compute_bivariate_CE_index
#' @title Compute the Clark-Evans index for a pair of point patterns
#'
#' @description Compute the Clark-Evans index for a pair of point patterns. Internal function called by the CE_interaction_tensor() function
#'
#' @param X_pattern a ppp object.
#' @param Y_pattern a ppp object.
#' @param k_parameter a ppp object.
#' @return Returns the bivariate Clark-Evans Index for this pair of point pattern.
#' @examples
#' CE_interaction = Compute_bivariate_CE_index(X_pattern,Y_pattern,k = k_parameter)
#' @import spatstat
#' @export



Compute_bivariate_CE_index = function(X_pattern,Y_pattern,k_parameter = 1) {
  
  NN_distance = nncross(X_pattern,Y_pattern,k = k_parameter)
  
  #Observed values 
  Mean_distance = mean(NN_distance$dist)
  
  #Theoretical disrtibution
  Lambda_2 = Y_pattern$n / spatstat::area(Y_pattern)
  Mean_distance_theo = 1/(2*sqrt(Lambda_2))
  Variance_theo = (4-pi)/(4*Lambda_2*pi*X_pattern$n)
  
  #Statistical test
  Z_score = (Mean_distance-Mean_distance_theo)/Variance_theo
  P_value = pnorm(Z_score,mean = 0,sd = 1,lower.tail = F,log.p = T)
  
  #Index
  Bivariate_CE_index = Mean_distance * 2 *sqrt(Lambda_2)
  return(c(Bivariate_CE_index,Z_score,P_value))
}

