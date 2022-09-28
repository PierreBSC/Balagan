#' @rdname Compute_CP_decomposition
#' @title Compute the Canonical Polyadic decomposition of a tensor object and reorder the data in a more understandable way
#'
#' @description Compute the Canonical Polyadic decomposition of a tensor object and reorder the data in a more understandable way
#'
#' @param Interaction_tensor a three-order tensor object,
#' @param num_components number of components to compute
#' @param max_iter number of iterations
#' @param tol relative error tolerance 
#' @param Center_data center the data through third dimension before CP 
#' @param Scale_data scale the data through third dimension before CP 
#' @return Returns a list containing the score matrix (score for each sample in each CP dimension), the list of typical matrices (product of the left and right CP vectors) and a vector with the lambda values
#' @examples
#' CP_decomposition = Compute_CP_decomposition(sce)
#' @import rTensor
#' @export

Compute_CP_decomposition = function(Interaction_tensor,num_components=5,max_iter = 1000,tol = 1e-6,Center_data = T, Scale_data=F) {
  
  if (Center_data) {
    cat("Centering the data \n !")
    N_col = dim(Interaction_tensor)[1]
    N_row = dim(Interaction_tensor)[2]
    
    for (i in 1:N_col) {
      for (j in 1:N_row) {
        Interaction_tensor[i,j,]@data = Interaction_tensor[i,j,]@data - mean(Interaction_tensor[i,j,]@data)
      }
    }
  }
  
  if (Scale_data) {
    cat("Scaline the data \n !")
    N_col = dim(Interaction_tensor)[1]
    N_row = dim(Interaction_tensor)[2]
    
    for (i in 1:N_col) {
      for (j in 1:N_row) {
        Interaction_tensor[i,j,]@data = Interaction_tensor[i,j,]@data / sd(Interaction_tensor[i,j,]@data)
      }
    }
    
  }
  
  
  cat("Computing the CP decomposition : \n")
  Resulting_decomposition = cp(Interaction_tensor,num_components = num_components,max_iter = max_iter,tol = tol)
  
  par(las=1,bty="l",mfcol=c(1,2))
  barplot(Resulting_decomposition$lambdas,ylab=c("Lambda value"),xlab="CP dimension",cex.lab=1.3)
  plot(as.numeric(Resulting_decomposition$est@data),as.numeric(Interaction_tensor@data),pch=21,bg="orange",
       ylim=range(c(as.numeric(Resulting_decomposition$est@data),as.numeric(Interaction_tensor@data))),
       xlim= range(c(as.numeric(Resulting_decomposition$est@data),as.numeric(Interaction_tensor@data))),
       xaxs="i",yaxs='i',xlab="Estimated values",ylab="Observed values",cex.lab=1.3)
  abline(0,1,lwd=2,lty=2)
  R_coef = cor(as.numeric(Resulting_decomposition$est@data),as.numeric(Interaction_tensor@data))
  legend("topleft",paste("R=",as.character(round(R_coef,2))),bty="n",cex=1.4)
  
  
  cat("Reorganizing the result of the CP decomposition")
  
  List_typical_matrix = c()
  for (k in 1:num_components) {
    List_typical_matrix[[k]] = Resulting_decomposition$U[[1]][,k]%*%t(Resulting_decomposition$U[[2]][,k])
    colnames(List_typical_matrix[[k]]) = 1:ncol(List_typical_matrix[[k]])
    rownames(List_typical_matrix[[k]]) = 1:nrow(List_typical_matrix[[k]])
   }
  
  Score_matrix = Resulting_decomposition$U[[3]]
  colnames(Score_matrix) = paste("CP_dimension_",as.character(1:ncol(Score_matrix)),sep = "")
  return(list(Score_matrix = Score_matrix,List_typical_matrix= List_typical_matrix,Lambda_vector = Resulting_decomposition$lambdas,R_coef = R_coef))
  
  
  
}
