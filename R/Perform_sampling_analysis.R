#' @rdname Perform_sampling_analysis
#' @title Sampling analysis
#'
#' @description Perform a sampling analysis where one or various sampling parameters are varying   
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param Selected_image of which image/FOV should the sampling analysis be performed ?
#' @param N_times number of times each type of sampling is performed (typically between 20 and 100)
#' @param N_sampling_region_vector vector (or integer) describing the different values taken by the number of sampled regions
#' @param width_FOV_vector vector (or real number) describing the width of the FOV
#' @param height_FOV_vector vector (or real number) describing the height of the FOV
#' @param Threshold_detection_cluster real number corresponding to the minimal number of a given cell type to be considered as detected. Typically around 50 or 100.
#' @return Returns a data.frame containing various summary statistics such as the mean number of detected clusters, the KL divergence with the real cell type composition etc...
#' @examples
#' #Simple case where only the number of sample regions changes :
#' Perform_sampling_analysis(sce,Selected_image=1,N_times=50,N_sampling_region_vector=1:20,width_FOV_vector=400,height_FOV_vector=400,Threshold_detection=50)
#' @import pheatmap
#' @export


Perform_sampling_analysis = function(sce,Selected_image=1,N_times=50,
                                     N_sampling_region_vector=1:20,
                                     width_FOV_vector=400,height_FOV_vector=400,
                                     Threshold_detection_cluster=50) {
  
  #Checking/Modifiying the shape of the parameter vector
  Length_vector = c(length(N_sampling_region_vector),length(width_FOV_vector),length(height_FOV_vector))
  Length_vector = unique(Length_vector)
  Length_vector = Length_vector[Length_vector!=1]
  
  if (length(Length_vector)>1) {
    stop("multiple parameter vectors of with a size bigger than 1 have been provided. Please use appropriate parameters !")
  }
  
  if (length(N_sampling_region_vector)==1) {
    N_sampling_region_vector = rep(N_sampling_region_vector,Length_vector)
  }
  
  if (length(width_FOV_vector)==1) {
    width_FOV_vector = rep(width_FOV_vector,Length_vector)
  }
  
  if (length(height_FOV_vector)==1) {
    height_FOV_vector = rep(height_FOV_vector,Length_vector)
  }
  
  ##Computing the global cell composition
  
  Global_composition = table(factor(colLabels(sce)))
  Global_composition_normalised =Global_composition/sum(Global_composition)
  Global_composition = log(Global_composition/prod(Global_composition)^(1/length(Global_composition)))

  #Doing the sampling
  
  Mean_number_cluster_identified = c()
  Sd_number_cluster_identified = c()
  
  Mean_divergence_global_composition = c()
  Sd_divergence_global_composition = c()
  
  Mean_correlation_composition = c()
  Sd_correlation_composition = c()
  
  Mean_KL_score = c()
  Sd_KL_score = c()
  
  
  ##Performing the sampling in itself
  
  for (i in 1:Length_vector) {
    
    print(i)
    Number_cluster_identified_temp = c()
    Divergence_temp = c()
    Correlation_temp = c()
    KL_temp = c()
    
    for (j in 1:N_times) {
      
      x = Random_spatial_sampling(sce,Selected_image = Selected_image,width_FOV = width_FOV_vector[i],height_FOV = height_FOV_vector[i],
                                  N_samplings = N_sampling_region_vector[i],plot_result = F)
      
      table_sampled_clusters = table(factor(x$List_sampled_cluster,levels = levels(factor(colLabels(sce)))))
      table_sampled_clusters_normalized = table_sampled_clusters/sum(table_sampled_clusters)
      
      Number_cluster_identified_temp = c(Number_cluster_identified_temp,sum(table_sampled_clusters>Threshold_detection_cluster))
      KL_temp = c(KL_temp,sum(table_sampled_clusters_normalized*log(table_sampled_clusters_normalized/Global_composition_normalised),na.rm = T))
      
      sampled_composition = table_sampled_clusters+1 #Adding a pseudo count of 1 for the Atchison distance
      sampled_composition = log(sampled_composition/(prod(sampled_composition)^(1/length(sampled_composition))))
      Aitchison_distance = sqrt(sum((sampled_composition-Global_composition)^2))
      Divergence_temp = c(Divergence_temp,Aitchison_distance)
      
      table_sampled_cluster_normalised = table_sampled_clusters/sum(table_sampled_clusters)
      R_temp = cor(table_sampled_cluster_normalised,Global_composition_normalised)
      Correlation_temp = c(Correlation_temp,R_temp)
    }
    
    Mean_divergence_global_composition = c(Mean_divergence_global_composition,mean(Divergence_temp))
    Sd_divergence_global_composition = c(Sd_divergence_global_composition,sd(Divergence_temp))
    
    Mean_number_cluster_identified = c(Mean_number_cluster_identified,mean(Number_cluster_identified_temp))
    Sd_number_cluster_identified = c(Sd_number_cluster_identified,sd(Number_cluster_identified_temp))
    
    Mean_correlation_composition = c(Mean_correlation_composition,mean(Correlation_temp))
    Sd_correlation_composition = c(Sd_correlation_composition,sd(Correlation_temp))
    
    Mean_KL_score = c(Mean_KL_score,mean(KL_temp))
    Sd_KL_score = c(Sd_KL_score,sd(KL_temp))
    
  }
  
  Statistic_data_frame = data.frame(N_sampling = N_sampling_region_vector,
                                    Mean_number_cluster = Mean_number_cluster_identified,
                                    Sd_number_cluster = Sd_number_cluster_identified,
                                    Mean_KL_divergence = Mean_KL_score,
                                    Sd_KL_divergence = Sd_KL_score,
                                    Mean_correlation_composition = Mean_correlation_composition,
                                    Sd_correlation_composition = Sd_correlation_composition,
                                    Mean_Aitchison_distance = Mean_divergence_global_composition,
                                    Sd_Aitchison_distance = Sd_divergence_global_composition)
  
  return(Statistic_data_frame)
}
