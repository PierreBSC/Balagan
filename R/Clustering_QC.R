#' @rdname Clustering_QC
#' @title Comparison and QC of clustering result 
#'
#' @description This function controls the quality of the clustering by comparing cluster 
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param assay_type name of the data slot on which to perform the clustering (Raw_intensity, Arcsinh_transformed_intensity....). By default the regression-normalized data are used
#' @return Returns a \code{\link[SingleCellExperiment]{SingleCellExperiment}} object with a new assay slot called "Count_normalised_intensity"
#'
#' @examples
#' sce = Count_normalization(sce,residual_normalisation = "Anscombe")
#' @import SingleCellExperiment
#' @import pheatmap

#' @export

Clustering_QC = function(sce,assay_type="Raw_intensity") {
  
  if (is.null(colLabels(sce))) {
    stop("Please perform the clustering step before !")
    
    Clustering_temp = as.character(colLabels(sce))
    Mean_cluster_profile = aggregate(as.data.frame(t(assay(sce,assay_type))),FUN = median,by=list(Clustering_temp))
    rownames(Mean_cluster_profile) = Mean_cluster_profile$Group.1
    Mean_cluster_profile = Mean_cluster_profile[,-1]
    
    Meta_clustering = pheatmap(cor(t(Mean_cluster_profile),method = "spearman"),clustering_method = "ward.D")
    Meta_clustering = Meta_clustering$tree_col
    
    pheatmap(t(Mean_cluster_profile),cluster_cols = Meta_clustering,scale = "row",clustering_method = "ward")
    
  }
}