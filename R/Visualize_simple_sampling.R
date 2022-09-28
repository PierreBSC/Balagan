#' @rdname Visualize_sampling
#' @title Visualization of simple sampling analysis
#'
#' @description Function that plot and analyze the results from the Perform_sampling_analysis() and estimate the tau and No parameters.
#'
#' @param Sampling_data_frame a data.frame object resulting from the Perform_sampling_analysis() function
#' @return Returns a plot as well as tau and No estimate for each of the parameter value.
#' @examples
#' #Simple case where only the number of sample regions changes :
#' Visualize_sampling(Simple_sampling_analysis)
#' @export


Visualize_simple_sampling = function(Sampling_data_frame) {
  Max_cluster_recovered = round(1.2*max(Sampling_data_frame$Mean_number_cluster))
  Max_sampled_FoV = round(1.2*max(Sampling_data_frame$N_sampling))
  
  par(las=1,bty="l")
  plot(Sampling_data_frame$N_sampling,Sampling_data_frame$Mean_number_cluster,
       xlim=c(0,Max_sampled_FoV),ylim=c(0,Max_cluster_recovered),xaxs='i',yaxs='i',
       xlab="Number of sampled regions",ylab="Mean number of recovered clusters",cex.lab=1.2,
       pch=21,bg="red3",cex=2)
  
  y = Sampling_data_frame$Mean_number_cluster
  x = Sampling_data_frame$N_sampling
  expo_model_number_cluster = nls(y~N*(1-exp(-x/tau)),start=list(N=20,tau=5)) 
  curve(expr = coef(expo_model_number_cluster)[1]*(1-exp(-x/coef(expo_model_number_cluster)[2])) ,add = T,col="red",lwd=2,from = 0,to=100,lty=2)
  Result_vector = coef(expo_model_number_cluster)
  R_squared = cor(predict(expo_model_number_cluster,newdata = x),y)^2
  Result_vector = c(Result_vector,R_squared)
  names(Result_vector)[3] = "R_squared"
  
  abline(h=Result_vector[1],lwd=1.5,lty=2,col="grey")
  
  legend("right",legend = c(paste("N =",round(Result_vector[1],2)),
                            paste("tau =",round(Result_vector[2],2)),
                            paste("R2 =",round(Result_vector[3],3))),bty="n",cex=1.2)
  return(Result_vector)
}
