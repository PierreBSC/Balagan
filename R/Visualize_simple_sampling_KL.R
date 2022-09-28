#' @rdname Visualize_sampling_KL
#' @title Visualization of simple sampling analysis with a focus on the cell composition. 
#'
#' @description Function that plot and analyze the results from the Perform_sampling_analysis() and estimate the theta, KL_0 and KL_b parameters (Kullback-Leibler divergence analysis).
#'
#' @param Sampling_data_frame a data.frame object resulting from the Perform_sampling_analysis() function
#' @return Returns a plot as well as estimate values for theta, KL_0 and KL_b parameters. 
#' @examples
#' #Simple case where only the number of sample regions changes :
#' Visualize_sampling_KL(Simple_sampling_analysis)
#' @export


Visualize_simple_sampling_KL = function (Sampling_data_frame) 
{
  Max_KL = round(1.2 * max(Sampling_data_frame$Mean_KL_divergence))
  Max_sampled_FoV = round(1.2 * max(Sampling_data_frame$N_sampling))
  par(las = 1, bty = "l")
  plot(Sampling_data_frame$N_sampling, Sampling_data_frame$Mean_KL_divergence, 
       xlim = c(0, Max_sampled_FoV), ylim = c(0, Max_KL), 
       xaxs = "i", yaxs = "i", xlab = "Number of sampled regions", 
       ylab = "Mean KL divergence", cex.lab = 1.2, 
       pch = 21, bg = "red3", cex = 2)
  
  
  y = Sampling_data_frame$Mean_KL_divergence
  x = Sampling_data_frame$N_sampling
  x = x-1
  
  expo_baseline_KL = nls(y ~ KL_0 *  exp(-x/theta)+KL_b, 
                         start = list(KL_0 = 1, theta = 2,KL_b=0.01))
  
  
  curve(expr = coef(expo_baseline_KL)[1]*exp(-(x-1)/coef(expo_baseline_KL)[2])+coef(expo_baseline_KL)[3], add = T, col = "red", lwd = 2, from = 0, to = 100, lty = 2)  
  Result_vector = coef(expo_baseline_KL)
  R_squared = cor(predict(expo_baseline_KL, newdata = x), 
                  y)^2
  Result_vector = c(Result_vector, R_squared)
  names(Result_vector)[4] = "R_squared"
  abline(h = Result_vector[1], lwd = 1.5, lty = 2, col = "grey")
  legend("right", legend = c(paste("KL_0 =", round(Result_vector[1], 2)), paste("Theta =", round(Result_vector[2], 2)),paste("KL_b =", round(Result_vector[3], 2)), paste("R2 =", round(Result_vector[4], 3))), bty = "n", cex = 1.2)
  return(Result_vector)
}





