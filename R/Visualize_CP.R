#' @rdname Visualize_CP
#' @title Visualize the result of a Canonical Polyadic (CP) decomposition
#'
#' @description Basic function allowing to visualize the result of a Canonical Polyadic (CP) decomposition
#'
#' @param CP_object the output of the cp (rTensor package) function
#' @param Plot_type the type of plot to produce. Has to be chosen among "all", "loading matrix" or
#' @return 
#' @examples
#' Visualize_CP(CP_object)
#' @import pheatmap
#' @export


Visualize_CP = function(CP_decomposition,k=1) {
  
  
  N_rank = ncol(CP_decomposition$U[[1]])
  
  List_prototypical_matrix = c()
  
    List_prototypical_matrix[[k]] = matrix(CP_decomposition$U[[1]][,k],ncol = 1)%*%matrix(CP_decomposition$U[[2]][,k],nrow = 1)
    x = CP_decomposition$U[[3]][,k]
    List_prototypical_matrix[[k]] = List_prototypical_matrix[[k]] / mean(x)
    rownames(List_prototypical_matrix[[k]]) = paste("Cluster",1:nrow(List_prototypical_matrix[[k]]))
    colnames(List_prototypical_matrix[[k]]) = paste("Cluster",1:ncol(List_prototypical_matrix[[k]]))
    
    pheatmap(List_prototypical_matrix[[k]],main = paste("Prototypical matrix of rank",as.character(k)),cluster_rows = F,cluster_cols = F)
  
  
  
}
