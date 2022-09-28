#' @rdname Plot_cluster_spatial
#' @title Visualization of cluster spatial distribution
#'
#' @description Plot the spatial distribution of the identified clusters or of any other provided grouping of the cells
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}. Clustering of the cell should have been previously performed.
#' @param Image_number Number of name of the image/ROI to be plotted
#' @param Cex_parameter Scaling factor for the size of the cells 
#' @param Specific_cluster Colors points based based on their belonging to a specific cluster 
#' @param Provided_group User-provided clustering
#' @return Return a plot 
#' @examples 
#' Plot_cluster_spatial(sce,Image_number = 1)
#'
#'
#' @import SingleCellExperiment 
#' @export

Plot_cluster_spatial = function(sce,Image_number = 1,Cex_parameter=10,Specific_cluster=NULL,Provided_cluster=NULL) {
  if (is.null(colLabels(sce))) {
    stop("Please compute clustering first ! \n")
  }
  
  
  Temp_location_data = data.frame(X=sce$Location_Center_X[sce$ImageNumber==Image_number],
                                  Y=sce$Location_Center_Y[sce$ImageNumber==Image_number])
  Temp_cluster_data = colLabels(sce)[sce$ImageNumber==Image_number]
  
  if (!is.null(Provided_cluster)) {
    Temp_cluster_data = Provided_cluster
  }
  
  Dimension = metadata(sce)$dimension
  
  if ("Cell_size"%in%colnames(colData(sce))) {
    Temp_size_data = sce$Cell_size[sce$ImageNumber==Image_number]
  }
  
  if (!"Cell_size"%in%colnames(colData(sce))) {
    Temp_size_data = rep(10,sum(sce$ImageNumber==Image_number))
  }
  
  
  if (Dimension == "2D") {
    Temp_size_data = sqrt(Temp_size_data)
  }
  
  if (Dimension == "3D") {
    Temp_size_data = (Temp_size_data)^1/3
  }
  
  
  
  
  par(las=1,bty="l")
  
  if (is.null(Specific_cluster)) {
    color_temp_vector = .cluster_to_color(Temp_cluster_data)
  }
  
  if (!is.null(Specific_cluster)) { 
    if (Specific_cluster%in%unique(Temp_cluster_data)) {
      color_temp_vector = (as.numeric(Temp_cluster_data==Specific_cluster))
    }
  }
  
  

  plot(Temp_location_data,cex=Temp_size_data/Cex_parameter,pch=21,bg=color_temp_vector)
  
}
