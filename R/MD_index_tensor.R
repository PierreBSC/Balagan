#' @rdname MD_index_tensor
#' @title Compute the mean direction index
#'
#' @description Compute the mean direction (MD) index for each pair of cell type in each image and store them in a three-order tensor
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param k_parameter number of NN to use to compute the mean direction index
#' @return Returns three-order tensor with the log2 CE index of each pair of cell type in each image
#' @examples
#' MD_tensor = MD_index_tensor(sce)
#' @import SingleCellExperiment
#' @import spatstat.geom
#' @import rTensor
#' @import doParallel
#' @import foreach
#' @export

Compute_bivariate_MD_index = function(X_pattern,Y_pattern,k_parameter = 4,normalize_MD = TRUE) {
  
  if (!k%in%c(2,3,4,5,6)) {
    stop("Please select a k parameter among 2, 3, 4, 5 or 6 !")
  }
  
  NN_ID = nncross(X_pattern,Y_pattern,k = 1:(k_parameter),what = "which")
  
  
  ###
  Delta_X = c()
  for (k in 1:k_parameter) {
    Delta_X = cbind(Delta_X,Y_pattern$x[NN_ID[,k]]-X_pattern$x)
  }
  
  Delta_Y = c()
  for (k in 1:k_parameter) {
    Delta_Y = cbind(Delta_Y,Y_pattern$y[NN_ID[,k]]-X_pattern$y)
  }
  
  #Normalising to have vectors of similar size and only keeping the direction
  Norm_vector = sqrt(Delta_X^2 + Delta_Y^2)
  Delta_X = Delta_X/Norm_vector
  Delta_Y = Delta_Y/Norm_vector
  
  #Computing the length of the summed vectors
  MD_index =rowSums(Delta_X)^2+rowSums(Delta_Y)^2
  Mean_MD_index = mean(MD_index)
  

  
  #Observed values 
  return(Mean_MD_index)
}


Compute_MD_index = function(X_pattern,k_parameter = 4) {
  
  if (!k%in%c(2,3,4,5,6)) {
    stop("Please select a k parameter among 2, 3, 4, 5 or 6 !")
  }
  
  NN_ID = nncross(X_pattern,X_pattern,k = 2:(k_parameter+1),what = "which")
  
  
  ###
  Delta_X = c()
  for (k in 1:k_parameter) {
    Delta_X = cbind(Delta_X,Y_pattern$x[NN_ID[,k]]-X_pattern$x)
  }
  
  Delta_Y = c()
  for (k in 1:k_parameter) {
    Delta_Y = cbind(Delta_Y,Y_pattern$y[NN_ID[,k]]-X_pattern$y)
  }
  
  #Normalising to have vectors of similar size and only keeping the direction
  Norm_vector = sqrt(Delta_X^2 + Delta_Y^2)
  Delta_X = Delta_X/Norm_vector
  Delta_Y = Delta_Y/Norm_vector
  
  #Computing the length of the summed vectors
  MD_index =rowSums(Delta_X)^2+rowSums(Delta_Y)^2
  Mean_MD_index = mean(MD_index)
  
  
  
  #Observed values 
  return(Mean_MD_index)
}


MD_interaction_tensor = function(sce,type_output="Index",Perform_symmetrization=T) {
  

  N_cluster = max(colLabels(sce))
  N_image = max(sce$ImageNumber)
  List_MD_matrix = c()
  
  cat(paste("Creating parallel backend using"), as.character(metadata(sce)$N_core), "cores \n")
  registerDoParallel(metadata(sce)$N_core)
  
  List_MD_matrix=foreach(k = 1:N_image) %dopar% {
    
    cat(paste("Computing interactions for image",as.character(k),"\n"))
    MD_index_matrix = matrix(data = 0,nrow = max(colLabels(sce)),ncol= N_cluster)
    X_range = range(sce$Location_Center_X[sce$ImageNumber==k])
    Y_range = range(sce$Location_Center_Y[sce$ImageNumber==k])
    
    
    for (i in 1:N_cluster) {
      
      for (j in 1:N_cluster) {
        
        X_pattern = data.frame(X = sce$Location_Center_X[sce$ImageNumber==k & colLabels(sce)==i ],
                               Y = sce$Location_Center_Y[sce$ImageNumber==k & colLabels(sce)==i ])
        X_pattern = ppp(X_pattern$X,X_pattern$Y,window = owin(X_range,yrange = Y_range))
        
        Y_pattern = data.frame(X = sce$Location_Center_X[sce$ImageNumber==k & colLabels(sce)==j ],
                               Y = sce$Location_Center_Y[sce$ImageNumber==k & colLabels(sce)==j ])
        Y_pattern = ppp(Y_pattern$X,Y_pattern$Y,window = owin(X_range,yrange = Y_range))
        
        if (i!=j & npoints(X_pattern)>10 & npoints(Y_pattern)>10) {
          Temp_MD_anlaysis = Compute_bivariate_MD_index(X_pattern,Y_pattern,k_parameter = 4)
          MD_index_matrix[i,j] = Temp_MD_anlaysis
        }
        
        if (i==j & npoints(X_pattern)>10) {
          Temp_MD_anlaysis = Compute_MD_index(X_pattern,k_parameter = 4)
          MD_index_matrix[i,j] = Temp_MD_anlaysis
          
        }
        
        
      }
    }
    #Cleaning and symmetrizing the CE index matrix
    

      
    rownames(MD_index_matrix) = 1:nrow(MD_index_matrix)
    colnames(MD_index_matrix) = 1:ncol(MD_index_matrix)
      
    MD_index_matrix
  }
  
  cat("Storing the data into a three-order tensor \n")
  
  Tensor_object = rand_tensor(modes = c(N_cluster, N_cluster, N_image), drop = FALSE)
  for (k in 1:N_image) {
    Tensor_object[,,k] = List_MD_matrix[[k]]
  }
  
  return(Tensor_object)
  
}
