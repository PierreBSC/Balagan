#' @rdname Patch_dectection
#' @title Cellular patch detection 
#'
#' @description This function identifies patches of cells, i.e physically neighbor cells that share a similar cellular environement  
#'
#' @param sce a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' @param graph_type the type of graph to use to compute cellular neighborhood. Can be "KNN" for a k-nearest neighbor , "Radius" for epsilon-radius graph or  "Gabriel for a Gabriel graph.
#' @param graph_parameter the parameter value for the graph construction. Corresponds to the number of neighbor for the KNN and the radius value for the epsilon-radius graph.
#' @param neighbor_degree the degree of the neigbor used to compute the similarity of cellular environement between two cells. Typically 1 or 2. 
#' @param distance_composition type of distance to use to compare neighbor composition. Can be "cosine" for cosine similarity, "JS" for Jensen-Shanon divergence
#' @param clustering_method method used to clusters the cells once similarity/distances have been computed 
#' @param min_points minimal number of points to consider a community of points a patch
#' @param image_list images/FOV/location on which to perform the analysis. If not provided the analysis is done on all images.
#' @return 
#'
#' @examples
#' sce = Patch_detection(sce,graph_type="KNN",graph_parameter=30,neighbor_degree=2)
#' @import SingleCellExperiment
#' @import igraph
#' @import dbscan  
#' @import spatgraphs
#' @import spatstat
#' @import doParallel
#' @import foreach
#' @import N2R
#' @export

Patch_detection = function(sce,graph_type = "KNN",graph_parameter = NULL,neighbor_degree=2,clustering_method="saturation",distance_composition="Bray-Curtis",min_points = 15,image_list=NULL) {
  

  if (!graph_type%in%c("KNN","Radius","Gabriel")) {
    stop("Please provided a correct type of graph to compute !")
  }
  
  if (graph_type=="KNN" & is.null(graph_parameter)) {
    cat("No paramater provided for KNN graph computation. K parameter set to 15 \n !")
    graph_parameter = 15
  }
  
  if (graph_type=="Radius" & is.null(graph_parameter)) {
    cat("No paramater provided for Epsilon-radius graph computation. Epsilon parameter set to 50 \n !")
    graph_parameter = 50
  }
  
  if (is.null(image_list)) {
    image_list = unique(sce$ImageNumber)
  }
  
  cat(paste("Creating parallel backend using"),as.character(metadata(sce)$N_core),"cores \n")
  registerDoParallel(metadata(sce)$N_core)
  
  List_patch_analysis =foreach(selected_image=image_list) %dopar%  {
    
    #Data loading
    Selected_cells = sce$ImageNumber==selected_image
    Location_points = data.frame(Location_X = sce$Location_Center_X[Selected_cells],
                                 Location_Y = sce$Location_Center_Y[Selected_cells])
    Point_label = colLabels(sce)[Selected_cells]
    List_values = levels(factor(Point_label))
    
    #Graph computation 
    
    if (graph_type=="KNN") {
      KNN_graph_matrix = N2R::Knn(as.matrix(Location_points), k = graph_parameter,quiet = T, nThreads = 1, 
                                  verbose = F, indexType = "L2")
      KNN_graph_matrix = KNN_graph_matrix + t(KNN_graph_matrix)
      Final_graph <- graph_from_adjacency_matrix(KNN_graph_matrix, 
                                                 mode = "undirected", weighted = TRUE)
    }

    if (graph_type=="Radius") {
      RNN_list = dbscan::frNN(as.matrix(Location_points), sort=F,eps = graph_parameter)
      Final_graph <- graph_from_adj_list(RNN_list$id,mode = "all")
    }
    
    if (graph_type=="Gabriel") {
      Gabriel_graph_adj_list = spatgraph(ppp(x = Location_points$Location_X,y = Location_points$Location_Y,
                                             window = owin(range(Location_points$Location_X),yrange = range(Location_points$Location_Y))),type = "gabriel")
      Final_graph = graph_from_adj_list(Gabriel_graph_adj_list$edges,"all")
    }
    
    
    #Computing interactions scores between each pair of neighbor cells
    
    List_interactions = as_edgelist(Final_graph)
    List_weights = c()
    
    for (k in 1:nrow(List_interactions)) {
      
      node_i = List_interactions[k,1]
      node_j = List_interactions[k,2]
      
      Neighbor_i = as.numeric(neighborhood(graph = Final_graph,nodes = node_i,order = neighbor_degree)[[1]])
      Neighbor_j = as.numeric(neighborhood(graph = Final_graph,nodes = node_j,order = neighbor_degree)[[1]])
      
      Neighbor_i_specific = Neighbor_i[!Neighbor_i%in%Neighbor_j]
      Neighbor_j_specific = Neighbor_j[!Neighbor_j%in%Neighbor_i]
      
      composition_i = Point_label[Neighbor_i_specific]
      composition_j = Point_label[Neighbor_j_specific]
      
      
      composition_i_normalized = table(factor(composition_i,levels = List_values))/length(List_values)
      composition_j_normalized = table(factor(composition_j,levels = List_values))/length(List_values)
      
      composition_i = table(factor(composition_i,levels = List_values))
      composition_j = table(factor(composition_j,levels = List_values))
      
      if (distance_composition=="cosine") {
        Distance_score = 1-sum((composition_i_normalized*composition_j_normalized))/(sqrt(sum(composition_i_normalized^2))*sqrt(sum(composition_j_normalized^2)))
      }

      if (distance_composition=="JS") {
        composition_mixed = (composition_i_normalized+composition_j_normalized)/2
        Distance_score = composition_i_normalized*log(composition_i_normalized/composition_mixed) +composition_j_normalized*log(composition_j_normalized/composition_mixed)
        Distance_score = sum(Distance_score,na.rm = T)
      }
      
      if (distance_composition=="Bray-Curtis") {
        
        U = rbind(composition_i,composition_j)
        Min_comp = apply(U,MARGIN = 2,FUN = min)
        Distance_score = 1 - sum(Min_comp)/sum(U)
      }
      
      
      List_weights = c(List_weights,Distance_score)
    }
    
    if (distance_composition%in%c("cosine","Bray-Curtis")) {
      List_weights[is.na(List_weights)] = 1
    }
    
    ###Going from distance to similarity 
    
    List_weights_transformed = exp(-List_weights/sd(List_weights))
    E(Final_graph)$weight = List_weights_transformed*100
    
    ## Clustering on the graph
    
    if (clustering_method=="greedy") {
      Clustering_final = cluster_fast_greedy(Final_graph,modularity = T)
      Clustering_final = membership(Clustering_final)
    }

    if (clustering_method=="saturation") {
      
      Similarity_matrix = as_adjacency_matrix(Final_graph,sparse = T,attr = "weight")
      List_modularity_score = c()
      List_number_components = c()
      
      #Finding the best cut
      for (i in seq(min(E(Final_graph)$weight),quantile(E(Final_graph)$weight,0.99),length.out=30)) {
        print(i)
        Similarity_matrix_sat = Similarity_matrix
        Similarity_matrix_sat[Similarity_matrix_sat<i] = 0
        graph_temp = graph_from_adjacency_matrix(Similarity_matrix_sat,mode = "undirected",weighted = T)
        
        connected_component = components(graph_temp)
        connected_component = membership(connected_component)
        
        modularity_score = modularity(x = Final_graph,membership = connected_component)
        List_modularity_score = c(List_modularity_score,modularity_score)
        List_number_components = c(List_number_components,length(unique(connected_component)))
      }
      
      #Applying the cut in itself
  
      Threshold_similarity = seq(0,1,length.out=50)[which.max(List_modularity_score)]
      Similarity_matrix_sat = Similarity_matrix
      Similarity_matrix_sat[Similarity_matrix_sat<Threshold_similarity] = 0
      Graph_cuted = graph_from_adjacency_matrix(Similarity_matrix_sat,mode = "undirected",weighted = T)
      Louvain_clustering = fastgreedy.community(Graph_cuted)
      Clustering_final = Louvain_clustering$membership
      
      #Cleaning the resulting clusters
      N_point_cluster = table(factor(Clustering_final))
      clusters_to_filter = which(N_point_cluster<min_points)
      Clustering_final[Clustering_final%in%clusters_to_filter] = 0
      
      plot.igraph(Graph_cuted,layout = as.matrix(Location_points),vertex.label=NA,vertex.size=2,
                  edge.width = E(Graph_cuted)$weight,vertex.color=.cluster_to_color(Clustering_final))
      

    }

    Clustering_final
  }
  

  return(List_patch_analysis)
}

