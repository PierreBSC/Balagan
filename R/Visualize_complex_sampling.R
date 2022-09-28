#' @rdname Visualize_complex_sampling
#' @title Visualization of complex sampling analysis
#'
#' @description Function that plot and analyze the results from the Perform_sampling_analysis() and estimate the tau and No parameters when various parameters value have been used .
#'
#' @param Sampling_data_frame a data.frame object resulting from the Perform_sampling_analysis() function
#' @param Parameter_table a table that describe the parameter 
#' @return Returns a plot as well as tau and No estimate for each of the parameter value.
#' @examples
#' #Simple case where only the number of sample regions changes :
#' Visualize_sampling(Simple_sampling_analysis)
#' @import grDevices 
#' @import RColorBrewer 
#' @export


Visualize_complex_sampling = function(Sampling_data_frame,Parameter_table) {
  
  Max_cluster_recovered = round(1.2*max(Sampling_data_frame$Mean_number_cluster))
  Max_sampled_FoV = round(1.2*max(Sampling_data_frame$N_sampling))
  
  X = c()
  for (k in 1:ncol(Parameter_table)) {
    X = paste(X,Parameter_table[,k],sep = "_")
  }
  X = substr(X,start = 2,stop = 500)
  Parameter_values = unique(X)
  
  
  #Plot
  par(las=1,bty="l")
  plot(Sampling_data_frame$N_sampling,Sampling_data_frame$Mean_number_cluster,
       xlim=c(0,Max_sampled_FoV),ylim=c(0,Max_cluster_recovered),xaxs='i',yaxs='i',
       xlab="Number of sampled regions",ylab="Mean number of recovered clusters",cex.lab=1.2,
       pch=21,bg=.cluster_to_color(as.numeric(as.factor(X))),cex=2)
  
  segments(x0 = Sampling_data_frame$N_sampling,x1 = Sampling_data_frame$N_sampling,
           y0 = Sampling_data_frame$Mean_number_cluster-Sampling_data_frame$Sd_number_cluster,
           y1 = Sampling_data_frame$Mean_number_cluster+Sampling_data_frame$Sd_number_cluster,lwd=0.2)

  Result_table = c()
  
  for (k in 1:length(Parameter_values)) {
    y = Sampling_data_frame$Mean_number_cluster[X==Parameter_values[k]]
    x = Sampling_data_frame$N_sampling[X==Parameter_values[k]]
    expo_model_number_cluster = nls(y~N*(1-exp(-x/tau)),start=list(N=20,tau=5)) 
    curve(expr = coef(expo_model_number_cluster)[1]*(1-exp(-x/coef(expo_model_number_cluster)[2])) ,add = T,col="red",lwd=2,from = 0,to=100,lty=2)
    Result_vector = coef(expo_model_number_cluster)
    R_squared = cor(predict(expo_model_number_cluster,newdata = x),y)^2
    Result_vector = c(Result_vector,R_squared)
    names(Result_vector)[3] = "R_squared"
    Result_table = rbind(Result_table,Result_vector)
    
  }
  abline(h=max(Result_table[,1]),lwd=2,lty=2,col="grey")


  return(Result_table)
}
