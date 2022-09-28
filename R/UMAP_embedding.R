#' @rdname UMAP_embedding
#' @title Computation of a UMAP embedding 
#'
#' @description Computes a UMAP embedding 
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param assay_type name of the data slot on which to perform the clustering (Raw_intensity, Arcsinh_transformed_intensity....). By default the regression-normalized data are used
#'
#' @return Return an updated \code{\link[SingleCellExperiment]{SingleCellExperiment}} object with an updated reducedDim slot.
#' @examples 
#' sce = UMAP_embedding(sce,assay_type = "Count_normalised_intensity", metric="correlation",n_neighbors = 30)
#'
#'
#' @import uwot 
#' @export

UMAP_embedding = function(sce,assay_type = "Count_normalised_intensity",metric="correlation",n_neighbors = 30,use) {
  
  if (!assay_type%in%names(assays(sce))) {
    stop("The slot required does not exist. Please select an existing slot !")
  }

  
  Channel_for_clustering = rowData(sce)$Used_for_clustering
  
  #If no selection of the channels : using all channels
  if (sum(Channel_for_clustering)==0) {
    Channel_for_clustering = rep(TRUE,ncol(data_to_cluster))
  }
  

  
  data_to_project =assay(sce,assay_type)
  data_to_project = t(data_to_project)
  data_to_project = data_to_project[,Channel_for_clustering]
  
  
  cat("Computing UMAP embedding...")
  umap_embedding = umap(data_to_project,n_neighbors = n_neighbors,pca = NULL,metric = metric,verbose=F,fast_sgd = TRUE,n_threads = metadata(sce)$N_core)
  cat(" done ! \n")
  
  reducedDim(sce, "UMAP") <- umap_embedding
  return(sce)
}

