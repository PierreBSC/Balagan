#' @rdname Global_alpha_estimation
#' @title Estimation of alpha using many small FoVs
#'
#' @description Function that estimates the alpha parameter of a given tissue using many small FoVs instead of a unique large FoVs. 
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param List_image a vector with the number of the images to use for the estimate. If not provided, all images will be used.
#' @param Threshold_detection minimal number of cells required to consider a cell type/cluster as detected. By default set to 20.
#' @param Value_split vector containing the fraction of the FoV that are taken, i.e numbers bigger than 1. If not provided the following vector will be used : 1, 1.1, 1.2, 1.5, 1.8, 2, 2.5 and 3.
#' @return Returns a data.frame with the estimate of alpha for each FoV as well as the quality of each fitting (R-squared).
#' @examples
#' #Simple case where only the number of sample regions changes :
#' Alpha_estimate = Global_alpha_estimation(sce)
#' @import SingleCellExperiment
#' @export


Global_alpha_estimation = function(sce,List_image=NULL,Threshold_detection=20) {
  

  if (is.null(Value_split)) {
    Value_split = c(1,1.1,1.2,1.5,1.8,2,2.5,3) 
  }
  
  if (is.null(List_image)) {
    List_image = unique(sce$ImageNb)
  }
  
  ####
  List_coef = c()
  R_squared = c()
  
  
  for (k in List_image) {
    
    sce_temp = sce[,sce$ImageNb==k]
    x_range = max(sce_temp$Location_Center_X)
    y_range = max(sce_temp$Location_Center_Y)
    
    center_x = x_range/2
    center_y = y_range/2
    N_clusters_detected = c()
    
    for (i in Value_split) {
      x_range_temp = c(center_x-x_range/(i*2),center_x+x_range/(i*2))
      y_range_temp = c(center_y-y_range/(i*2),center_y+y_range/(i*2))
      sce_temp_filtered = sce_temp[,sce_temp$Location_Center_X<x_range_temp[2] & sce_temp$Location_Center_X>x_range_temp[1] & 
                                     sce_temp$Location_Center_Y>y_range_temp[1] & sce_temp$Location_Center_Y<y_range_temp[2] ]
      x = table(factor(colLabels(sce_temp_filtered),levels = unique(colLabels(sce))))
      N_clusters_detected = c(N_clusters_detected,sum(x>Threshold_detection))
    }
    
    m = lm(log10(N_clusters_detected)~log10(Value_split),subset = N_clusters_detected>0)
    List_coef = c(List_coef,coef(m)[2])
    R_squared = c(R_squared,summary(m)$r.squared)
  }
  Table_estimation = data.frame(Alpha = -List_coef,R_squared = R_squared,row.names = List_image)
  return(Table_estimation)
}


