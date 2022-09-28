#' @rdname Random_spatial_sampling
#' @title Random spatial sampling of square areas  
#'
#' @description Perform a random spatial sampling with a given number of rectangles of with a given width and height
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param width_FOV height of the individual rectangles 
#' @param height_FOV width of the individual rectangles
#' @param N_samplings number of individual rectangles samples
#' @param N_samplings number of individual rectangles samples
#' @param Selected_image of which image/FOV should the sampling be performed ?
#' @return Returns a list with the fields "List_sampled_cells" describing which cells have been selected and the field "List_sampled_cluster" where their associated cluster/label have been sampled
#'
#' @examples
#' Simple_sampling = Random_spatial_sampling(sce,width_FOV=400,height_FOV=400,N_samplings=10,Selected_image = 1,plot_result=TRUE)
#' @import SingleCellExperiment
#' @export
#' 
Random_spatial_sampling = function(sce,width_FOV=400,height_FOV=400,N_samplings=10,Selected_image = 1,plot_result=TRUE) {
  
  sce = sce[,sce$ImageNumber==Selected_image]
  
  Equivalent_radius = sqrt((width_FOV^2)/4+(height_FOV^2)/4) *2
  
  x_range_sampling = range(sce$Location_Center_X)
  y_range_sampling = range(sce$Location_Center_Y)
  
  List_center = c()
  
  for (k in 1:N_samplings) {
    Is_in_empty_space = FALSE
    #print(k)
    while(!Is_in_empty_space) {
      position_temp = c(runif(n = 1,min =x_range_sampling[1],max = x_range_sampling[2] ),
                        runif(n = 1,min =y_range_sampling[1],max = y_range_sampling[2] ))
      Dist_matrix = dist(rbind(position_temp,List_center),method = "manhattan")
      Dist_matrix =as.matrix(Dist_matrix)
      if (nrow(Dist_matrix)==1) {
        Is_in_empty_space = T
        List_center = rbind(List_center,position_temp)
      }
      if (nrow(Dist_matrix)!=1) {
        List_distance_temp = Dist_matrix[1,]
        if (min(List_distance_temp[-1])>Equivalent_radius) {
          Is_in_empty_space = T
          List_center = rbind(List_center,position_temp)
          
        }
      }
      
      
    }
  }
  
  List_sampled_cells = c()
  List_sampled_cluster = c()
  List_sample = c()
  
  for (k in 1:nrow(List_center)) {
    center_temp = List_center[k,]
    Selected_cells = which(sce$Location_Center_X> (center_temp[1]-width_FOV/2) &  sce$Location_Center_X < center_temp[1]+width_FOV/2 & sce$Location_Center_Y > center_temp[2]-height_FOV/2 & sce$Location_Center_Y < center_temp[2]+height_FOV/2 )
    List_sampled_cells[[k]] = Selected_cells
    List_sampled_cluster[[k]] = colLabels(sce)[Selected_cells]
    List_sample = c(List_sample,rep(paste("Sample",k,sep = "_"),length(Selected_cells)))
  }
  
  if (plot_result) {
    
    par(bty="n",las=1)
    plot(sce$Location_Center_X,sce$Location_Center_Y,pch=21,bg=.cluster_to_color(colLabels(sce)))

    for (k in 1:nrow(List_center)) {
      center_temp = List_center[k,]
      rect(xleft = center_temp[1]-width_FOV/2,ybottom = center_temp[2]-height_FOV/2,
           xright = center_temp[1]+width_FOV/2,ytop = center_temp[2]+height_FOV/2,col = "black",density = 40)
    }
  }
  
  List_result = list(List_sampling = List_sample,
                     List_sampled_cells = unlist(List_sampled_cells),
                     List_sampled_cluster = unlist(List_sampled_cluster))
  
  return(List_result)
}
