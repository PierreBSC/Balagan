#' @rdname Stratified_sampling
#' @title Simulation of a Stratified Sampling 
#'
#' @description Perform several simulation of a stratified sampling using a stratified image 
#'
#' @param Thresholded_image a cimg object 
#' @param sce SingleCellExperiment object processed by balagan
#' @param N_FoV_per_region numeric vector. Number of FoV to sample in each stratum.
#' @param FoV_size size of Field of Views (FoVs) to be sampled.
#' @param N_sampling number of simulated sampling to perform.
#' @param Weight_strata numeric vector that must sum to 1, corresponding to the contribution of each stratum to the general tissue composition.
#' @param List_strata  (optional) character or numeric vector. Ordered list of the stratum names.
#' @return Returns a table with the estimated cell density for each cell type and each simulation
#' @examples
#'Stratified_sampling(Thresholded_image = MVS_thresholded_image,sce = sce,N_FoV_per_region = c(1,2,3,1,2),FoV_size = 200,N_sampling = 50 ,Weight_strata =c(0.1,0.2,0.4,0.2,0.1),List_strata = 1:6)
#' @import imager 
#' @export


Stratified_sampling = function(Thresholded_image,sce,N_FoV_per_region,FoV_size = 100,N_sampling = 50,
                               Weight_strata=NULL,List_strata = NULL) {
  
  
  #Checking inputs and defining variables
  
  if (is.null(List_strata))  {
    List_strata = unique(Thresholded_image)
  }
  
  N_strata = length(List_strata)
  

  #If no weight vector provided : 
  if (is.null(Weight_strata)) {
    cat("No weight vector was provided. Performing proportional allocation using total surface ! \nThis is not recommended and a proper weighting vector should be provided by the user ! \n")
    Weight_strata = table(Thresholded_image)
    Weight_strata = Weight_strata/sum(Weight_strata)
    Weight_strata = as.numeric(Weight_strata)
  }
  
  ##
  List_cell_types = unique(colLabels(sce))
  
  
  Stratified_sampling_estimation = c()
  
  for (i in 1:N_sampling) {
    
    Table_estimation_mu = c()
    
    #If one region is not sampled : problem...
    
    #We go region by region
    List_sampled_cell_per_region = c()
    print(i)
    
    for (k in 1:length(List_strata)) {
      Position_sampling = Region_guided_sampling(Thresholded_image,
                                                 Selected_region = List_strata[k],
                                                 N_FoV = N_FoV_per_region[k],FoV_size = FoV_size,show_plot = F)
      
      Count_cell_temp = c()
      
      for (j in 1:N_FoV_per_region[k]) {
        Selected_cells = sce$Location_Center_X>(Position_sampling[j,1]-FoV_size/2) & sce$Location_Center_X<(Position_sampling[j,1]+FoV_size/2) & sce$Location_Center_Y>(Position_sampling[j,2]-FoV_size/2) & sce$Location_Center_Y<(Position_sampling[j,2]+FoV_size/2)
        Sampled_cell_cluster_cluster = colLabels(sce)[Selected_cells]
        Count_cell_temp = rbind(Count_cell_temp,table(factor(Sampled_cell_cluster_cluster,levels = List_cell_types)))
      }
      Mu_estimate_temp = colMeans(Count_cell_temp)
      Table_estimation_mu  = rbind(Table_estimation_mu,Mu_estimate_temp)
      
    }
    
    Stratified_sampling_estimation_temp = Weight_strata%*%Table_estimation_mu #Computing the final values/estimates by taking the weighted mean of each estimator
    
    Stratified_sampling_estimation = rbind(Stratified_sampling_estimation,Stratified_sampling_estimation_temp)
  }
  
  
  #Generating the QC plot
  Real_cell_density = table(factor(colLabels(sce),unique(colLabels(sce))))
  Real_cell_density = Real_cell_density/(ncol(Thresholded_image)*nrow(Thresholded_image)) * FoV_size^2
  
  par(las=1,bty="l")
  boxplot(Stratified_sampling_estimation,outline=F,xlab="Cell phenotype",ylab="Estimated cell density",
          main="Results of stratified sampling")
  points(1:length(List_cell_types),Real_cell_density,pch=21,bg="red3",cex=2)
  
  return(Stratified_sampling_estimation)
}
