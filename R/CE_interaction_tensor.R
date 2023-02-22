#' @rdname CE_interaction_tensor
#' @title Compute the Clark-Evans index for each pair of cell type in each image and store them in a three-order tensor
#'
#' @description Compute the Clark-Evans index for each pair of cell type in each image and store them in a three-order tensor
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @return Returns three-order tensor with the log2 CE index of each pair of cell type in each image
#' @examples
#' CE_tensor = CE_interaction_tensor(sce)
#' @import SingleCellExperiment
#' @import spatstat.geom
#' @import rTensor
#' @import doParallel
#' @import foreach
#' @export


CE_interaction_tensor = function(sce,type_output="Index",Perform_symmetrization=T) {
  
  if (!type_output%in%c("Index","Z-score")) {
    stop("Please select an appropriate type of output (Index or Z-score)")
  }
  
  N_cluster = max(colLabels(sce))
  N_image = max(sce$ImageNumber)
  List_CE_matrix = c()
  
  cat(paste("Creating parallel backend using"), as.character(metadata(sce)$N_core), "cores \n")
  registerDoParallel(metadata(sce)$N_core)
  
  List_CE_matrix=foreach(k = 1:N_image) %dopar% {
    
    cat(paste("Computing interactions for image",as.character(k),"\n"))
    CE_Z_score_matrix = matrix(data = 0,nrow = max(colLabels(sce)),ncol= N_cluster)
    CE_index_matrix = matrix(data = 1,nrow = max(colLabels(sce)),ncol= N_cluster)
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
          Temp_CE_anlaysis = Compute_bivariate_CE_index(X_pattern,Y_pattern,k_parameter = 1)
          CE_Z_score_matrix[i,j] = Temp_CE_anlaysis[2]
          CE_index_matrix[i,j] = Temp_CE_anlaysis[1]
        }
        
        if (i==j & npoints(X_pattern)>10) {
          Temp_CE_anlaysis = Compute_bivariate_CE_index(X_pattern,Y_pattern,k_parameter = 2)
          CE_Z_score_matrix[i,j] = Temp_CE_anlaysis[2]
          CE_index_matrix[i,j] = Temp_CE_anlaysis[1]
          
        }
        

      }
    }
    #Cleaning and symmetrizing the CE index matrix
    
    if (type_output=="Index") {
      CE_index_matrix[is.na(CE_index_matrix)]=1
      CE_index_matrix[CE_index_matrix==0]=1
      CE_index_matrix = log2(CE_index_matrix)
      
      if (Perform_symmetrization) {
        CE_index_matrix = (CE_index_matrix + t(CE_index_matrix))/2
      }
      
      
      rownames(CE_index_matrix) = 1:nrow(CE_index_matrix)
      colnames(CE_index_matrix) = 1:ncol(CE_index_matrix)
      
      x = CE_index_matrix
    }
    
    
    if (type_output=="Z-score") {
      CE_Z_score_matrix[is.na(CE_Z_score_matrix)]=0
      CE_Z_score_matrix = (CE_Z_score_matrix + t(CE_Z_score_matrix))/2
      rownames(CE_Z_score_matrix) = 1:nrow(CE_Z_score_matrix)
      colnames(CE_Z_score_matrix) = 1:ncol(CE_Z_score_matrix)
      
      x = CE_Z_score_matrix
    }
    
    x
  }
  
  cat("Storing the data into a three-order tensor \n")
  
  Tensor_object = rand_tensor(modes = c(N_cluster, N_cluster, N_image), drop = FALSE)
  for (k in 1:N_image) {
    Tensor_object[,,k] = List_CE_matrix[[k]]
  }
  
  return(Tensor_object)
  
}
