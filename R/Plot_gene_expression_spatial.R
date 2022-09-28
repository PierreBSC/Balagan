#' @rdname Plot_gene_expression_spatial
#' @title Visualization of spatial gene expression
#'
#' @description Plot the spatial distribution of a given gene normalised expression
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param assay_type name of the data slot on which to perform the clustering (Raw_intensity, Arcsinh_transformed_intensity....). By default the regression-normalized data are used
#' @param Image_number Number of name of the image/ROI to be plotted
#' @param Gene Name of the gene to be plotted
#' @param Cex_parameter Scaling factor for the size of the cells 
#'
#' @return Return a plot 
#' @examples 
#' Plot_gene_expression_spatial(sce,Gene = "CD45")
#'
#'
#' @import SingleCellExperiment 
#' @export

Plot_gene_expression_spatial = function(sce,Image_number = 1,Cex_parameter=5,assay_type = "Count_normalised_intensity",Gene=NULL) {
  
  if (is.null(Gene) | !Gene%in%rownames(sce@assays@data@listData$Raw_intensity) ) {
    stop("Please select a correct gene to plot ! \n")
  }
  
  if (!assay_type%in%names(assays(sce))) {
    stop("The slot required does not exist. Please select an existing slot !")
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
  

  Temp_location_data = data.frame(X=sce$Location_Center_X[sce$ImageNumber==Image_number],
                                  Y=sce$Location_Center_Y[sce$ImageNumber==Image_number])
  Temp_expression_data = as.numeric(assay(sce,assay_type)[Gene,])
  Temp_expression_data = Temp_expression_data - min(Temp_expression_data)
  par(las=1,bty="l")
  plot(Temp_location_data,cex=Temp_size_data/Cex_parameter,pch=21,bg=.color_convertion(Temp_expression_data)[sce$ImageNumber==Image_number])
  
}
