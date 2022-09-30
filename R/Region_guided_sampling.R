#' @rdname Region_guided_sampling
#' @title Random sampling within a specific part of the tissue.
#'
#' @description Random sampling within a specific part of the tissue using a thresholded image, typically done by MVS. 
#'
#' @param Segmentation_image a cimg object 
#' @param Selected_regions number of the stratum to be sampled.
#' @param N_FoV number of Field of Views (FoVs) to be sampled within the region of interest
#' @param FoV_size size of Field of Views (FoVs) to be sampled within the region of interest
#' @param show_plot boolean object. Should the result of the sampling be shown ? 
#' @return Returns the vector of all FoV center.
#' @examples
#' Region_guided_sampling(MVS_thresholded_image,Selected_region = 1,N_FoV = 10,FoV_size = 200,show_plot = T)
#' @import imager 
#' @export



Region_guided_sampling = function(Segmentation_image,Selected_regions=1,N_FoV=1,FoV_size=100,show_plot=FALSE) {
  
  
  List_coordinate = c()
  
  Im_width = imager::width(raw_Image)
  Im_height = imager::height(raw_Image)
  Segmentation_image_temp = Segmentation_image==Selected_regions
  Segmentation_image_temp = as.cimg(Segmentation_image_temp)
  
  for (i in 1:N_FoV) {
    
    #Computing the Distance Transform image
    DT_temp = distance_transform(Segmentation_image_temp,value = min(Segmentation_image_temp),metric = 0 )
    
    #Where can we take a FOV of a sufficient size ?
    DT_prob_map = DT_temp
    DT_prob_map[DT_prob_map<FoV_size/2]=0
    DT_prob_map = as.matrix(DT_prob_map)
    DT_prob_map = DT_prob_map/sum(DT_prob_map)
    #Randomly sampling a location for the FOV
    
    Pixel_in_border = TRUE
    n = 1
    while (Pixel_in_border) {
      
      Selected_pixel = sample(x = 1:(nrow(DT_prob_map)*ncol(DT_prob_map)),size = 1,prob = as.numeric(DT_prob_map))
      m = coord.index(DT_temp,Selected_pixel)
      
      if ( m$x > (FoV_size/2) & m$x<(Im_width-FoV_size/2) & m$y> (FoV_size/2) & m$y < (Im_height-FoV_size/2)) {
        Pixel_in_border=FALSE
      }
      n = n +1
    }
    #print(n)
    List_coordinate = rbind(List_coordinate,c(m$x,m$y))  
    
    #Updating the original map
    x_min = m$x-FoV_size/2
    x_max = m$x+FoV_size/2
    y_min = m$y-FoV_size/2
    y_max = m$y+FoV_size/2
    
    List_pixels_to_update = expand.grid(x = x_min:x_max ,y = y_min:y_max)
    
    at(Segmentation_image_temp,x=List_pixels_to_update[,1],y = List_pixels_to_update[,2]) = 0
  }
  
  colnames(List_coordinate) = c("x","y")
  
  if (show_plot) {
    plot(Segmentation_image_temp)
    points(List_coordinate,pch=21,bg="red3")
    
  }
  return(List_coordinate)
  
}
