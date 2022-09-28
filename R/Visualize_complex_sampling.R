#' @rdname Visualize_complex_sampling_KL
#' @title Visualization of complex sampling analysis focusing on cell composition
#'
#' @description Function that plot and analyze the results from the Perform_sampling_analysis() and estimate the and estimate the theta, KL_0 and KL_b parameters (Kullback-Leibler divergence analysis) when various parameters value have been used .
#'
#' @param Sampling_data_frame a data.frame object resulting from the Perform_sampling_analysis() function
#' @param Parameter_table a table that describe the parameter 
#' @return Returns a plot as well as tau and No estimate for each of the parameter value.
#' @examples
#' #Simple case where only the number of sample regions changes :
#' Visualize_complex_sampling_KL(Simple_sampling_analysis)
#' @import grDevices 
#' @import RColorBrewer 
#' @export



Visualize_complex_sampling_KL   = function (Sampling_data_frame, Parameter_table) 
{
  Max_KL = round(1.2 * max(Sampling_data_frame$Mean_KL_divergence))
  Max_sampled_FoV = round(1.2 * max(Sampling_data_frame$N_sampling))
  X = c()
  for (k in 1:ncol(Parameter_table)) {
    X = paste(X, Parameter_table[, k], sep = "_")
  }
  X = substr(X, start = 2, stop = 500)
  Parameter_values = unique(X)
  par(las = 1, bty = "l")
  plot(Sampling_data_frame$N_sampling, Sampling_data_frame$Mean_KL_divergence, 
       xlim = c(0, Max_sampled_FoV), ylim = c(0, Max_KL), 
       xaxs = "i", yaxs = "i", xlab = "Number of sampled regions", 
       ylab = "Mean KL divergence", cex.lab = 1.2, 
       pch = 21, bg = cluster_to_color(as.numeric(as.factor(X))), 
       cex = 2)
  segments(x0 = Sampling_data_frame$N_sampling, x1 = Sampling_data_frame$N_sampling, 
           y0 = Sampling_data_frame$Mean_KL_divergence - Sampling_data_frame$Sd_KL_divergence, 
           y1 = Sampling_data_frame$Mean_KL_divergence + Sampling_data_frame$Sd_KL_divergence, 
           lwd = 0.2)
  Result_table = c()
  for (k in 1:length(Parameter_values)) {
    y = Sampling_data_frame$Mean_KL_divergence[X == Parameter_values[k]]
    x = Sampling_data_frame$N_sampling[X == Parameter_values[k]]
    x = x-1
    
    expo_baseline_KL = nls(y ~ KL_0 *  exp(-x/theta)+KL_b, 
                           start = list(KL_0 = 1, theta = 2,KL_b=0.01))
    
    
    
    curve(expr = coef(expo_baseline_KL)[1]*exp(-(x-1)/coef(expo_baseline_KL)[2])+coef(expo_baseline_KL)[3], add = T, 
          col = "red", lwd = 2, from = 0, to = 100, lty = 2)
    Result_vector = coef(expo_baseline_KL)
    R_squared = cor(predict(expo_baseline_KL, newdata = x), 
                    y)^2
    Result_vector = c(Result_vector, R_squared)
    names(Result_vector)[4] = "R_squared"
    Result_table = rbind(Result_table, Result_vector)
  }
  return(Result_table)
}
