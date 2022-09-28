#' @rdname Dixon_interaction_tensor
#' @title Compute cell type segregation using Dixon's approach
#'
#' @description Compute cell type segregation using Dixon's approach described in papers by Dixon (1994 and 2002)
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @return Returns three-order tensor with the log-odds ratio of each pair of cell type
#' @examples
#' Dixon_tensor = Dixon_interaction_tensor(sce)
#' @import SingleCellExperiment
#' @import spatstat
#' @import rTensor
#' @import doParallel
#' @import foreach
#' @import Matrix

#' @export


Dixon_interaction_tensor = function(sce,type_output="Z-score") {
  
  if (!type_output%in%c("logg-odds","Z-score")) {
    stop("Please select an appropriate type of output (logg-odds or Z-score)")
  }
  
  N_cluster = max(colLabels(sce))
  N_image = max(sce$ImageNumber)
  List_CE_matrix = c()
  
  cat(paste("Creating parallel backend using"), as.character(metadata(sce)$N_core), "cores \n")
  registerDoParallel(metadata(sce)$N_core)
  
  List_Dixon_matrix=foreach(k = 1:N_image) %dopar% {
    
    cat(paste("Computing interactions for image",as.character(k),"\n"))
    X_range = range(sce$Location_Center_X[sce$ImageNumber==k])
    Y_range = range(sce$Location_Center_Y[sce$ImageNumber==k])
    
    Levels_clustering = levels(factor(colLabels(sce)))
        
    Temp_pattern = data.frame(X = sce$Location_Center_X[sce$ImageNumber==k ],
                               Y = sce$Location_Center_Y[sce$ImageNumber==k ])
    Temp_pattern = ppp(Temp_pattern$X,Temp_pattern$Y,window = owin(X_range,yrange = Y_range),
                        marks = factor(colLabels(sce)[sce$ImageNumber==k],levels = Levels_clustering))
    
    NN_distance = nncross(Temp_pattern,Temp_pattern,k = 2)
    NN_mark_table = data.frame(Point_mark = marks(Temp_pattern),NN_mark = marks(Temp_pattern)[NN_distance$which.2] )
    
    Contingency_table = table(NN_mark_table$Point_mark,NN_mark_table$NN_mark)
    Contingency_vector = table(marks(Temp_pattern))
    N = npoints(Temp_pattern)
    
    ##W_table : table where wij tells if j is the NN of i (Sparse matrix)
    
    W_table = matrix(0,nrow = N,ncol=N)
    W_table = as(W_table,"dgCMatrix")
    
    for (i in 1:N) {
      j = NN_distance[i,2]
      W_table[i,j] = 1
    }
    
    ##Computing R (number of reflexive NN ) and Q (number of points sharing a neighbor) factors
    
    R = sum(as.numeric((W_table)*t(W_table)))
    Q = table(table(NN_distance$which.2))
    Q = sum(Q[-1]*(2:length(Q)))
    
    ##Computing the logg-odd ratio
    
    if (type_output%in%"logg-odds") {
      
      #Computing the diagonal
      S_diagonal = log2(diag(Contingency_table)*(N-Contingency_vector)/((Contingency_vector-diag(Contingency_table))*(Contingency_vector-1)))
      
      #Computing the non diagonal terms
      Association_matrix = Contingency_table
      for (i in 1:N_cluster) {
        for (j in 1:N_cluster) {
          Association_matrix[i,j] = Association_matrix[i,j]*(N-Contingency_vector[j]-1)/(Contingency_vector[i]-Association_matrix[i,j])/Contingency_vector[i]
        }
      }
      Association_matrix = log2(Association_matrix)
      diag(Association_matrix) = S_diagonal
      
      colnames(Association_matrix) = 1:N_cluster
      rownames(Association_matrix) = 1:N_cluster
      Association_matrix[is.infinite(Association_matrix)] = 0
      x = Association_matrix
    }
    
    ##Computing the Z-score
    
    if (type_output%in%"Z-score") {
      
      #
      Z_score_matrix = Contingency_table
      Association_matrix = matrix(0,nrow = N_cluster,ncol=N_cluster)
      for (i in 1:N_cluster) {
        for (j in 1:N_cluster) {
          
          if (i==j) {
            
            N_i = Contingency_vector[i]
            Pii =  N_i/N * (N_i-1)/(N-1) 
            Piii =  Pii * (N_i-2)/(N-2)
            Piiii =  Piii * (N_i-3)/(N-3)
            
            Expected_NN = N * Pii
            Variance_NN = (N+R) * Pii + (2*N-2*R+Q)* Piii + (N^2-3*N-Q+R)*Piiii - (N*Pii)^2
            
            Z_score_NN = (Contingency_table[i,i]-Expected_NN)/sqrt(Variance_NN)
            Association_matrix[i,j] = Z_score_NN
          }
          
          if (i!=j) {
            
            N_i = Contingency_vector[i]
            N_j = Contingency_vector[j]
            
            Pij =  N_i * N_j / N / (N-1)
            Piij = N_i/N * (N_i-1)/(N-1) * N_j/(N-2)
            Piijj = Piij * (N_j-1)/(N-3)
            
            Expected_NN = Pij * N 
            Variance_NN = N*Pij + Q*Piij +(N^2-3*N-Q+R)*Piijj -(N*Pij)^2
            
            Z_score_NN = (Contingency_table[i,j]-Expected_NN)/sqrt(Variance_NN)
            Association_matrix[i,j] = Z_score_NN
          }
          
        }
      }
      

      colnames(Association_matrix) = 1:N_cluster
      rownames(Association_matrix) = 1:N_cluster
      Association_matrix[is.infinite(Association_matrix)] = 0
      Association_matrix[is.na(Association_matrix)]=0
      x = Association_matrix
    }
    


    x
  }
  
  cat("Storing the data into a three-order tensor \n")
  
  Tensor_object = rand_tensor(modes = c(N_cluster, N_cluster, N_image), drop = FALSE)
  for (k in 1:N_image) {
    Tensor_object[,,k] = List_Dixon_matrix[[k]]
  }
  
  return(Tensor_object)
  
}
