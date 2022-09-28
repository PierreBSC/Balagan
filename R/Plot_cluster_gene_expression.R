#' @rdname Plot_cluster_gene_expression
#' @title Visualization of a specific gene expression across clusters
#'
#' @description Plot a boxplot of gene expression across clusters
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param assay_type name of the data slot on which to perform the clustering (Raw_intensity, Arcsinh_transformed_intensity....). By default the regression-normalized data are used
#' @param Gene Name of the gene to be plotted
#'
#' @return Return a plot 
#' @examples 
#' Plot_cluster_gene_expression(sce,Gene = "CD45")
#'
#'
#' @import SingleCellExperiment 
#' @export

Plot_cluster_gene_expression = function(sce,assay_type = "Count_normalised_intensity",Gene=NULL) {
  
  if (is.null(Gene) | !Gene%in%rownames(sce@assays@data@listData$Raw_intensity) ) {
    stop("Please select a correct gene to plot ! \n")
  }
  
  if (!assay_type%in%names(assays(sce))) {
    stop("The slot required does not exist. Please select an existing slot !")
  }
  
  Temp_expression_data = as.numeric(assay(sce,assay_type)[Gene,])
  
  par(las=1,bty="l")
  Color_vector_temp = .cluster_to_color(unique(colLabels(sce)))
  boxplot(Temp_expression_data~colLabels(sce),outline=F,xlab="Cluster",ylab="Gene expression value",main=Gene,col=Color_vector_temp)
}
