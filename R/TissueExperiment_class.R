library(spatstat)
library(FactoMineR)
library(foreach)
library(doParallel)

setClass("TissueExperiment", representation(r_range = "numeric",
                                            Contingency_table = "table",
                                            CA_analysis = "list",
                                            List_global_K_function = "data.frame",
                                            List_global_pcf_function = "data.frame",
                                            List_cluster_K_function = "list",
                                            List_cluster_pcf_function = "list"))


Create_TissueExperiment = function(sce,r_range = seq(0,300,length.out =100)) {
  
  N_Images = length(unique(sce$ImageNumber))
  #Computing zero/first order statistics and CA object
  
  cat("Computing Correspondance Analysis object...")
  Contingency_table = table(sce$ImageNumber,colLabels(sce))
  colnames(Contingency_table) = paste("Cluster_",colnames(Contingency_table),sep = "")
  rownames(Contingency_table) = paste("Image_",rownames(Contingency_table),sep = "")
  
  CA_analysis = CA(Contingency_table,graph = F)
  
  cat("done ! \n")
  
  
  #Computing global Ripley K functions 
  cat("Computing global Ripley's K functions..")
  registerDoParallel(cores = metadata(sce)$N_core)
  List_Kest = foreach(k=1:N_Images) %dopar% {
    
    Temp_location_data = data.frame(X = sce$Location_Center_X[sce$ImageNumber==k],
                                    Y = sce$Location_Center_Y[sce$ImageNumber==k])

    
    Temp_ppp_object = ppp(x = Temp_location_data$X,y=Temp_location_data$Y,
                          window = owin(xrange = range(Temp_location_data$X),yrange = range(Temp_location_data$Y)))
    Kest_temp = Kest(Temp_ppp_object,correction = "isotropic",r = r_range)
    Kest_temp
  }
  
  Merged_Kest = List_Kest[[1]]
  
  for (k in 2:length(List_Kest)) {
    Merged_Kest = suppressWarnings(bind.fv(Merged_Kest,List_Kest[[k]]$iso))
    
  }
  names(Merged_Kest) = c("r","theo",paste("Kest_Image_",1:length(List_Kest),sep = ""))
  fvlabels(Merged_Kest) = c("r","Kest[pois](r)",paste("Kest[Image_",1:length(List_Kest),"](r)",sep = ""))
  cat(" done ! \n")
  
  #Computing global pcf function 
  cat("Computing global pair Correlation Functions (pcf)..")
  
  List_pcf = foreach(i=unique(sce$ImageNumber)) %dopar% {
    #Estimating the pcf from the K function using the following strategy (Spatstat doc)
    # apply smoothing to Y(r) = K(r)/(2 * pi * r) constraining Y(0) = 0, estimate the derivative of Y, and solve
    
    pcf_temp = pcf(List_Kest[[i]],method = "b")
    
  }
  
  Merged_pcf = List_pcf[[1]]
  
  for (k in 2:length(List_pcf)) {
    Merged_pcf = suppressWarnings(bind.fv(Merged_pcf,List_pcf[[k]]$pcf))
    
  }
  names(Merged_pcf) = c("r","theo",paste("Pcf_Image_",1:length(List_pcf),sep = ""))
  fvlabels(Merged_pcf) = c("r","pcf[pois](r)",paste("pcf[Image_",1:length(List_pcf),"](r)",sep = ""))
  cat(" done ! \n")
  
  #Computing the cluster specific K functions
  
  cat("Computing the cluster-specific K and pair correlation functions in each image... \n")
  List_cluster_K_functions = c()
  List_cluster_pcf = c()
  Ordered_list_cluster = as.numeric(levels(factor(colLabels(sce))))
  
  for (k in 1:N_Images) {
    
    #Computing the K function for each cluster in a parallelised way..
    
    cat(paste("Processing of Image"),as.character(k),"...")
    List_Cluster_K_function_temp = foreach(i=Ordered_list_cluster) %dopar% {
      
      Selected_cells = sce$ImageNumber==k & colLabels(sce)==i

      if (sum(Selected_cells)>0) {
        Temp_ppp_object = ppp(x = sce$Location_Center_X[Selected_cells],
                              y = sce$Location_Center_Y[Selected_cells],
                              window = owin(xrange = range(sce$Location_Center_X[Selected_cells]),
                                            yrange = range(sce$Location_Center_Y[Selected_cells])))
        
        Kest_temp = Kest(Temp_ppp_object,correction = "isotropic",r = r_range)
      }
      
      if (sum(Selected_cells)<=1) {
        Kest_temp = fv(data.frame(r = r_range,iso = NA),argu = "r",valu = "iso")
      }
      Kest_temp
    }
    cat(" done !\n")
    
    
    #Aggregating the K functions from the same image 
    
    Cluster_K_function_temp_merged = List_Cluster_K_function_temp[[1]]
    Cluster_K_function_temp_merged$theo = NULL
    
    for (i in 2:length(List_Cluster_K_function_temp)) {
      Cluster_K_function_temp_merged = suppressWarnings(bind.fv(Cluster_K_function_temp_merged,List_Cluster_K_function_temp[[i]]$iso))
      
    }
    names(Cluster_K_function_temp_merged) = c("r",paste("Kest_cluster",1:length(List_Cluster_K_function_temp),sep = ""))
    fvlabels(Cluster_K_function_temp_merged) = c("r",paste("pcf[Cluster_",1:length(List_Cluster_K_function_temp),"](r)",sep = ""))
    List_cluster_K_functions[[k]] = Cluster_K_function_temp_merged
    
    
    #Generating and aggregating the pcf from the same image 
    
    Cluster_pcf_temp_merged = data.frame(r = r_range)
    
    for (i in 1:length(List_Cluster_K_function_temp)) {
      
      X = List_Cluster_K_function_temp[[i]]
      
      if (sum(is.na(X$iso))==length(X$iso)) {
        Cluster_pcf_temp_merged[,paste("Cluster",i,sep = "_")]= NA
      }
      
      if (sum(is.na(X$iso))!=length(X$iso)) {
        Temp_pcf = pcf(List_Cluster_K_function_temp[[i]],method = "b")
        
        if (length(Temp_pcf$pcf)!=length(r_range)) {
          Cluster_pcf_temp_merged[,paste("Cluster",i,sep = "_")]= NA
        }
        
        
        if (length(Temp_pcf$pcf)==length(r_range)) {
          Cluster_pcf_temp_merged[,paste("Cluster",i,sep = "_")]= Temp_pcf$pcf
        }
      }
      
    }
    
    
    List_cluster_pcf[[k]] = Cluster_pcf_temp_merged
  }
  
  #Computing the spatial simpson index alpha 
  
  Spatial_alpha_index = foreach(k=1:N_Images,.combine = cbind) %dopar% {
    
    
    Temp_location_data = data.frame(X = sce$Location_Center_X[sce$ImageNumber==k],
                                    Y = sce$Location_Center_Y[sce$ImageNumber==k])
    
    #Computing the lambda parameters
    Area_size = spatstat::area(owin(xrange = range(Temp_location_data$X),yrange = range(Temp_location_data$Y)))
    Lambda_parameter_global = nrow(Temp_location_data)/Area_size
    Lambda_parameter_cluster = table(factor(colLabels(sce))[sce$ImageNumber==k])/Area_size
    
    #Extracting the K and pair correlation function
    Temp_global_K_function = as.data.frame(Merged_Kest)[,k+2]
    
    Temp_cluster_K_function = as.data.frame(List_cluster_K_functions[[k]])
    Temp_cluster_K_function = Temp_cluster_K_function[,c(-1,-2)]
    
    #na_columns = colSums(is.na(Temp_cluster_K_function))==nrow(Temp_cluster_K_function)
    #Temp_cluster_K_function[,na_columns] = 0
        
    Temp_cluster_K_function = apply(X = Temp_cluster_K_function,MARGIN = 1,FUN = function(x) {x*Lambda_parameter_cluster^2})
    Temp_cluster_K_function = t(Temp_cluster_K_function)
    Temp_cluster_K_function = apply(X = Temp_cluster_K_function, MARGIN = 2,FUN = function(x) {x/(Lambda_parameter_global^2*Temp_global_K_function)})
    
    Alpha_index = 1- rowSums(Temp_cluster_K_function,na.rm = T)
    Alpha_index      
  }
  
  
  #Computing the spatial simpson index beta 
  
  Spatial_beta_index = foreach(k=1:N_Images,.combine = cbind) %dopar% {
    
    
    Temp_location_data = data.frame(X = sce$Location_Center_X[sce$ImageNumber==k],
                                    Y = sce$Location_Center_Y[sce$ImageNumber==k])
    
    #Computing the lambda parameters
    Area_size = spatstat::area(owin(xrange = range(Temp_location_data$X),yrange = range(Temp_location_data$Y)))
    Lambda_parameter_global = nrow(Temp_location_data)/Area_size
    Lambda_parameter_cluster = table(factor(colLabels(sce))[sce$ImageNumber==k])/Area_size
    
    #Extracting the K and pair correlation function
    Temp_global_pcf = as.data.frame(Merged_pcf)[,k+2]
    
    Temp_cluster_pcf = as.data.frame(List_cluster_pcf[[k]])
    Temp_cluster_pcf = Temp_cluster_pcf[,c(-1,-2)]
    

    Temp_cluster_pcf = apply(X = Temp_cluster_pcf,MARGIN = 1,FUN = function(x) {x*Lambda_parameter_cluster^2})
    Temp_cluster_pcf = t(Temp_cluster_pcf)
    Temp_cluster_pcf = apply(X = Temp_cluster_pcf, MARGIN = 2,FUN = function(x) {x/(Lambda_parameter_global^2*Temp_global_pcf)})
    
    Beta_index = 1- rowSums(Temp_cluster_pcf,na.rm = T)
    Beta_index      
  }
  
  
  #Computing the R index 
  
  CE_index_matrix = matrix(data = 1,nrow = length(unique(sce$ImageNumber)),ncol = length(unique(colLabels(sce))))
  colnames(CE_index_matrix) = paste("Cluster",levels(factor(unique(colLabels(sce)))))
  rownames(CE_index_matrix) = paste("Image",unique(unique(sce$ImageNumber)))
  
  CE_p_value_matrix = matrix(data = 1,nrow = length(unique(sce$ImageNumber)),ncol = length(unique(colLabels(sce))))
  colnames(CE_p_value_matrix) = paste("Cluster",levels(factor(unique(colLabels(sce)))))
  rownames(CE_p_value_matrix) = paste("Image",unique(unique(sce$ImageNumber)))
  
  
  N_minimal_point = 15
  
  for (k in unique(sce$ImageNumber)) {
    for (i in as.numeric(levels(factor(unique(colLabels(sce)))))) {
      
      Selected_images = colLabels(sce)==i & sce$ImageNumber==k
      
      if (sum(Selected_images)>N_minimal_point) {
        Temp_ppp_object = ppp(sce$Location_Center_X[Selected_images],sce$Location_Center_Y[Selected_images],
                              window = owin(xrange = range(sce$Location_Center_X[Selected_images]),yrange = range(sce$Location_Center_Y[Selected_images])))
  
        CE_test_temp = clarkevans.test(Temp_ppp_object,alternative = "clustered")
        CE_index_matrix[k,i] = CE_test_temp$statistic
        CE_p_value_matrix[k,i] = CE_test_temp$p.value
      }
    }
  }
  
  CE_p_value_matrix_corrected = matrix(data = p.adjust(CE_p_value_matrix),
                                       nrow = length(unique(sce$ImageNumber)),ncol = length(unique(colLabels(sce))))
  colnames(CE_p_value_matrix_corrected) = paste("Cluster",levels(factor(unique(colLabels(sce)))))
  rownames(CE_p_value_matrix_corrected) = paste("Image",unique(unique(sce$ImageNumber)))

  Image_number = 4
  Cluster_number = 13
  plot(r_range,(List_cluster_pcf[[Image_number]])[[Cluster_number+2]],xaxs="i",type="o",bg="orange",pch=21)
  Plot_cluster_spatial(sce,Image_number = Image_number,Cex_parameter = 4,Specific_cluster =Cluster_number )
  
  
  #Computing the mingling index 
  library(RANN)
  K_parameter = 30
  Mean_mingling_index = c()
  for (k in unique(sce$ImageNumber)) {
    
    Selected_cells = sce$ImageNumber==k
    Selected_position = data.frame(Location_X = sce$Location_Center_X[Selected_cells],
                                   Location_Y = sce$Location_Center_Y[Selected_cells])
    Selected_labels = colLabels(sce)[Selected_cells]
    
    KNN_matrix_temp = nn2(Selected_position,k = K_parameter+1)
    KNN_matrix_temp = KNN_matrix_temp$nn.idx
    
    KNN_cluster_temp = apply(KNN_matrix_temp,MARGIN = 1,FUN = function(x) {Selected_labels[x]})
    KNN_cluster_temp = t(KNN_cluster_temp)
    Mingling_index_list = apply(KNN_cluster_temp,MARGIN = 1,FUN = function(x) {
      x = x==x[1]
      return(sum(x[-1])/K_parameter)
    })
    Mean_mingling_index = c(Mean_mingling_index,median(Mingling_index_list))
  }
  
  
}
