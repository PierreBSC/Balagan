#' @rdname Perform_stratified_sampling_simulation
#' @title Perform stratified sampling simulation for a list of markers
#'
#' @description Allows to perform stratified sampling simulation for a list of markers.
#'
#' @param sce  a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param tiff_file path toward the raw tiff file
#' @param panel_file path toward a csv file associating each channel of the tiff file to gene. The file should have column with name "Target"
#' @param L number of strata. Should ideally be between 4 and 6 
#' @param N_FoV Number of fields of views
#' @param Fov_size Size (i.e width/height) of the square FoVs
#' @param sigma smoothing parameter  
#' @param specific_channel which channels should be tested for stratified sampling ? If not provided all channels will be used
#' @param perform_rotation performing rotation between the sce object and tiff image coordinate 
#' @param type_stratification how should the allocation of FoVs across strata be computed ? Can be "Neyman" or "Proportional" for Neyman/optimal and proportional allocation respectively
#' @param type_thresholding how should stratum be computed based on the intensity of the covariate ? Can be "Cumsqrt" or "Geometric"
#' @param Parallel_computing should the different sampling be performed in parallel ? Boolean.
#' @param show_plot should the intermediate plots be shown, i.e results of tissue stratification. Boolean.
#' @return A list of table, each table corresponding to a different stratified sampling 
#' @import SingleCellExperiment 
#' @import imager 
#' @import tiff 
#' @export


Perform_stratified_sampling_simulation = function(sce,tiff_file,panel_file,N_simulations=30,N_FoV=10,FoV_size=100,type_stratification="Neyman",type_thresholding="Cumsqrt",
                                                  L=6,sigma=100,specific_channels=NULL,perform_rotation=TRUE,show_plot=TRUE,Parallel_computing=TRUE) {
  
  #Checking the different parameter values
  if (L<2){
    stop("There must be at least 2 strata !")
  }
  
  if (!file.exists(tiff_file)) {
    stop("Please provide a correct path to the tiff file")
  }
  
  if (!file.exists(panel_file)) {
    stop("Please provide a correct path to the panel file")
  }
  
  if (!type_stratification %in% c("Proportional","Neyman")) {
    stop("Please provide a correct type of stratification (Either Neyman or Proportional)")
  }
  
  if (!type_thresholding %in% c("Cumsqrt","Geometric")) {
    stop("Please provide a correct type of thresholding (Either Cumsqrt or Geometric)")
  }
  
  
  List_sampling = c()
  
  #Loading the files
  cat("Loading the raw tiff file... ")
  Panel_data = read.delim(panel_file,sep=",")
  raw_Image = readTIFF(tiff_file,all = TRUE)
  cat("done ! \n")
  N_channels = length(raw_Image)
  
  if (is.null(specific_channels)) {
    specific_channels = 1:N_channels
  }
  
  for (k in specific_channels) {
    cat(paste("Performing stratified sampling for"),Panel_data$Target[k],"\n")
    
    Image_temp = raw_Image[[k]]
    Image_temp = as.cimg(Image_temp)
    if (perform_rotation) {
      Image_temp= imager::imrotate(Image_temp,angle = -90)
    }
    
    #Smoothing the data
    Image_temp_smoothed = vanvliet(vanvliet(Image_temp,sigma = sigma,order = 0,axis = "x"),sigma = sigma,axis = "y",order = 0)
    Image_temp_smoothed = Image_temp_smoothed/max(Image_temp_smoothed)
    
    
    #Performing the tissue stratification
    cat("Computing the optimal thresholds...")
    temp_thresholding = Strata_thresholding(as.numeric(Image_temp_smoothed),L = L,type_thresholding=type_thresholding)
    title(paste("Stratification thresholds for gene",Panel_data$Target[k]))
    MVS_thresholded_image = Image_temp_smoothed
    for (h in 1:L) {
      MVS_thresholded_image[Image_temp_smoothed>=temp_thresholding[h] & Image_temp_smoothed<=temp_thresholding[h+1]]=h
    }
    
    if (show_plot) {
      plot(as.cimg(MVS_thresholded_image),main=paste("Stratification for channel",Panel_data$Target[k]))
    }
    cat(" done ! \n")
    
    #Computing the area of each stratum
    Pixel_number_per_strata = table(factor(MVS_thresholded_image,levels = 1:L))
    Pixel_fraction_per_strata = as.numeric(Pixel_number_per_strata)/sum(Pixel_number_per_strata)
    
    #Estimating variance
    if (type_stratification=="Neyman") {
      Estimated_variance = aggregate(as.numeric(Image_temp_smoothed),by=list(as.numeric(MVS_thresholded_image)),FUN = var)
      Estimated_variance = Estimated_variance$x
    }
    
    #Assigning the number of sampled FoV for each stratum
    
    if (type_stratification=="Neyman") {
      Proportion_FoV_assignments = sqrt(Estimated_variance)*Pixel_fraction_per_strata
      Proportion_FoV_assignments = Proportion_FoV_assignments/sum(Proportion_FoV_assignments)
    }
    
    if (type_stratification=="Proportional") {
      Proportion_FoV_assignments = Pixel_fraction_per_strata
    }
    
    N_FoV_per_region = smart_round(x = Proportion_FoV_assignments*N_FoV)
    #Checking that each stratum is sampled at least once
    if (sum(N_FoV_per_region==0)>0) {
      Empty_slot = which(N_FoV_per_region==0)
      N_FoV_per_region_unrounded = Proportion_FoV_assignments*N_FoV
      for (i in 1:length(Empty_slot)) {
        N_FoV_per_region[Empty_slot[i]] = 1
        N_FoV_per_region[which.max(N_FoV_per_region_unrounded)] = N_FoV_per_region[which.max(N_FoV_per_region_unrounded)]-1
        N_FoV_per_region_unrounded[which.max(N_FoV_per_region_unrounded)] = N_FoV_per_region_unrounded[which.max(N_FoV_per_region_unrounded)]-1
        
      }
    }
    
    #Finally perfoming the stratification
    Sampling_temp = Stratified_sampling(Thresholded_image = MVS_thresholded_image,sce = sce,N_FoV_per_region =N_FoV_per_region ,FoV_size = FoV_size,
                                        N_sampling = N_simulations ,Weight_strata =Pixel_fraction_per_strata,Parallel_computing=Parallel_computing,name_stratifying_marker=Panel_data$Target[k])
    colnames(Sampling_temp) = paste("Cluster_",colnames(Sampling_temp),sep = "")
    rownames(Sampling_temp) = paste("Simulation",1:nrow(Sampling_temp),sep = "_")
    
    List_sampling[[Panel_data$Target[k]]] = Sampling_temp
  }
  
  return(List_sampling)
}
