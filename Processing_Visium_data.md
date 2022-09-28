# Processing and analysis of Visium datasets 

While Balagan was initially designed to work on Multiplexed Imaging (MI) data, sampling analyses can also be performed on spatial transcriptomic datasets. Here, we describe how to process raw Visium data (i.e output of CellRanger) so that they can be analyzed by Balagan.

Example datasets can be found on the [10X website](https://www.10xgenomics.com/resources/datasets?menu%5Bproducts.name%5D=Spatial%20Gene%20Expression&query=&page=1&configure%5Bfacets%5D%5B0%5D=chemistryVersionAndThroughput&configure%5Bfacets%5D%5B1%5D=pipeline.version&configure%5BhitsPerPage%5D=500).

## Loading and pre-processing of the raw data

Here we will use the [pagoda2 package](https://github.com/kharchenkolab/pagoda2) from the Kharchenko Lab to load the data and later process them :

We start by loading the different packages and the expression matrix :

```r
library(pagoda2)
library(CountClust)
library(doParallel)
library(balagan)
data_raw = read.10x.matrices("Visium_Data/raw_feature_bc_matrix//")
l = colnames(data_raw)
l = substr(l,start = 5,stop = 600)
colnames(data_raw) = l
```

We then load the location data (i.e correspondance between the beads barcode and their spatial location) and normalize data size so that we get the exact location of each bead (µm unit):

```r
data_location = read.delim("Desktop/TLO_analysis/Breast_cancer/Data_1/spatial/tissue_positions_list.csv",sep=",",header = F,row.names = 1)
data_location = data_location[colnames(data_raw),]
data_location = data_location[,c(2,3)]
colnames(data_location) = c("Location_X","Location_Y")
data_location$Location_X = data_location$Location_X*6400/78
data_location$Location_Y = data_location$Location_Y*6400/128
```

We can then remove low quality beads/location and genes :

```r
Lib_size = colSums(data_raw)
hist(log10(Lib_size),100)
Gene_size = rowSums(data_raw)
hist(log10(Gene_size+1),100)

#We remove beads with less than 1000 UMIs and genes with less than 10 genes
data_count = data_raw[Gene_size>10,Lib_size>1000]
Multiple_gene = names(which(table(rownames(data_count))>1))
data_count = data_count[!rownames(data_count)%in%Multiple_gene,]
Location_count = data_location[Lib_size>1000,]
```


## Analysis of the data

We can now perform the real analysis : first we look for the most variable genes.

```r
r <- Pagoda2$new(data_count,log.scale=T)
r$adjustVariance(plot=T,gam.k=5,)
```

We now take the 1500 most variables genes (can be adjusted based on the dataset) and perform Latent Dirichlet Allocation reduction dimension. To improve quality of the analysis we remove Immunoglobulin and mitochondrial genes :

```r
Selected_genes = r$getOdGenes(1500)
Selected_genes = Selected_genes[!grepl(Selected_genes,pattern = "IGH")]
Selected_genes = Selected_genes[!grepl(Selected_genes,pattern = "IGK")]
Selected_genes = Selected_genes[!grepl(Selected_genes,pattern = "IGL")]
Selected_genes = Selected_genes[!grepl(Selected_genes,pattern = "MT-")]
```

Then we perform LDA with different numbers of topics : this step is computationnaly heavy and is performed here in a parallel manner :

```r
registerDoParallel(cores=4)

Model_LDA_merged = foreach(k=c(10,15,20,25,30)) %dopar% {
  Model_LDA_temp = FitGoM(as.matrix(t(data_count[Selected_genes,])),K = k,tol=100,options="BIC")
  Model_LDA_temp
}
```

The BIC value of each model is then extracted and plotted. Similarly to a Scree plot for a PCA, we select the number of topic that display a clear elbow. Here we select the model with 20 topics :

```r
Selected_models = Model_LDA_merged[[3]]
Mixing = Selected_models$fit$omega
Contribution = Selected_models$fit$theta
Marginal_distribution = rowSums(data_count[Selected_genes,])/sum(data_count[Selected_genes,])
```

We now inject the LDA reduction to the pagoda2 object and perform the single-cell clustering :

```r
r$reductions$LDA = Mixing
r$makeKnnGraph(k=15,type='LDA',distance = "cosine")
r$getKnnClusters(method=multilevel.community,type='LDA',name = "LDA_cluster")
r$getDifferentialGenes(type = "LDA",clusterType = "LDA_cluster",verbose = T,z.threshold = 3)
```

Finally, we create the SingleCellExperiment experiment object :

```r
sce = SingleCellExperiment(assays = list(Raw_intensity = as.matrix(t(Mixing))), 
                           metadata = list(dimension = "2D", 
                                           N_core = 8, Is_nuc_cyt = F))
colLabels(sce) = as.numeric(r$clusters$LDA$LDA_cluster)
sce$Location_Center_X = Location_count$Location_X
sce$Location_Center_Y = Location_count$Location_Y
sce$ImageNumber =  1
```

The sce object can now be used for sampling analysis with balagan !