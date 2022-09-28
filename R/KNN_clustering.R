#' @rdname KNN_clustering
#' @title Perform cell clustering on a computed KNN graph 
#'
#' @description Perform cell clustering on a computed KNN graph 
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param assay_type name of the data slot on which to perform the clustering (Raw_intensity, Arcsinh_transformed_intensity....). By default the regression-normalized data are used
#' @param K number of neighbors for the KNN graph computation
#' @param clustering_method method used for graph clustering. Has to be chose among "Louvain","Greedy" and "Infomap"
#' @return Returns an updated \code{\link[SingleCellExperiment]{SingleCellExperiment}} object with an updated colLabels slot
#'
#' @examples 
#' sce = KNN_clustering(sce,K = 30,clustering_method = "Louvain")
#'#'
#' @import SingleCellExperiment 
#' @import igraph 
#' @import N2R 

#' @export

KNN_clustering = function(sce,K=30,clustering_method = "Louvain",assay_type="Count_normalised_intensity",metric="angular") {
  
  if (!assay_type%in%names(assays(sce))) {
    stop("The slot required does not exist. Please select an existing slot !")
  }
  
  if (!clustering_method %in%c("Louvain","Greedy","Infomap")) {
    stop("The clustering method required does not exist. Please choose among Louvain,Greedy and Infomap !")
  }
    
  
  if (!metric %in%c("L2","angular")) {
    stop("The distance metric required does not exist. Please choose among angular or L2!")
  }
  
  
  data_to_cluster =assay(sce,assay_type)
  data_to_cluster = t(data_to_cluster)
  
  Channel_for_clustering = rowData(sce)$Used_for_clustering
  
  #If no selection of the channels : using all channels
  if (sum(Channel_for_clustering)==0) {
    Channel_for_clustering = rep(TRUE,ncol(data_to_cluster))
  }
  
  cat(paste("Clustering of the data using",as.character(sum(Channel_for_clustering)),"channels from the",assay_type,"slot ! \n"))
  
  data_to_cluster = data_to_cluster[,Channel_for_clustering]
  

  KNN_graph_matrix =  Knn(as.matrix(data_to_cluster), K, nThreads=metadata(sce)$N_core, verbose=F, indexType=metric)
  KNN_graph_matrix = KNN_graph_matrix + t(KNN_graph_matrix)
  Final_graph <- graph_from_adjacency_matrix(KNN_graph_matrix,mode='undirected',weighted=TRUE)
  cat("KNN computed \n")
  
  if (clustering_method == "Louvain") {
    Clustering <- cluster_louvain(Final_graph)
  }
  
  if (clustering_method == "Greedy") {
    Clustering <- cluster_fast_greedy(Final_graph,modularity = TRUE)
  }
  
  if (clustering_method == "Infomap") {
    Clustering <- cluster_infomap(Final_graph,modularity = TRUE)
  }
  
  Clustering_group = membership(Clustering)
  colLabels(sce) = Clustering_group
  cat(paste(as.character(length(unique(Clustering_group))),"clusters have been identified \n"))
  return(sce)
}
