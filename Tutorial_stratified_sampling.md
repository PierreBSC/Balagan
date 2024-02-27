# Stratified spatial sampling analysis

This tutorial describes how to perform a simulated stratified spatial sampling on a piece of tissue. I strongly recommend the users to familiarize themselves with the stratified sampling approach by reading the excellent book by Cochran "Sampling Techniques, Third edition". 

To perform a stratified sampling, two elements will be needed :
- A SingleCellExperiment object, obtained by using our package balagan to process MI data.
- A .tiff file, corresponding to the raw IMC data.
- A .csv file that describes the antibody panel used to generate the IMC data

The order of the .csv file rows must be coherent with the order of the channels for the tiff file ! We also assume that the panel file contains a column called "Target" that corresponds to the gene name of each channel.


We first load the required packages and the SingleCellExperiment object 

```r
library(tiff)
library(balagan)
library(imager)
sce = readRDS("path/to/LN_sce.rds)
```


To be sure that the data are correct we first load the tiff file, plot channel number 10 (CD20, a B-cell marker that should clearly delimitate the different regions of the lymph node) and project the individual cell location colored by cell type from the SCE object.


```r
x =readTIFF("LN_quadrat_analysis/LN_panorama_raw_file.tiff",all = TRUE,as.is = TRUE)
plot(log(1+as.cimg(x[[10]])))
points(sce$Location_Center_X,sce$Location_Center_Y,col=string.to.colors(colLabels(sce)),cex=0.1)
```

<img src="Screenshot/LN_structure.png" alt="Stratified_sampling_estimation" width='500'> 


We observe a really nice overlap: the data are good and we can launch the stratified sampling.

```r
Systemic_stratified_sampling = Perform_stratified_sampling_simulation(sce,tiff_file = "LN_quadrat_analysis/LN_panorama_raw_file.tiff",
                                                                      panel_file="LN_quadrat_analysis/Panel_panorama_1.csv",N_simulations=30,
                                                                      N_FoV=10,FoV_size=100,type_thresholding = "Cumsqrt",type_stratification="Neyman",specific_channels = NULL,
                                                                      L=6,sigma=50,perform_rotation=FALSE,show_plot=TRUE,Parallel_computing=TRUE) 
```
Here many options are available :

- The sce file
- The path to the .tiff file (**tiff_file**)
- The path to the panel .csv file (**panel_file**)
- The number of simulations (N_simulations>=30)
- The number (**N_FoV**) and size (**FoV_size**) of the Fields of View
- The type of thresholding (**type_thresholding**), can be "cumsqrt" or "geometric"
- The allocation strategy of FoVs across strata (**type_stratification**), can be "Neyman" or "Proportional"
- Which channel to use (**specific_channels**). If not specified or set to NULL all channels will be tested
- The number of strata (**L**), should not be higher to 6 due to the limited gain and increased computational cost
- The smoothing parameter to be applied to the raw image (**sigma**)
- The activation of the parallisation (Parallel_computing). If set to TRUE the number of parallel process corresponds to the N_core parameter of the sce metadata slot. For instance to set it to 10 simply type :
```r
metadata(sce)$N_core = 10
```
- The plotting option (**show_plot**) which when set to TRUE allows to see all the plots produced during the sampling simulation



Once the sampling is finished (this can take some time...) we can compute the mean and variance of the estimated cell density :

```r
Mean_population_estimation = c()
Var_population_estimation = c()

for (k in 1:length(Systemic_stratified_sampling_2)) {
  u = Systemic_stratified_sampling_2[[k]]
  Mean_population_estimation = rbind(Mean_population_estimation,colMeans(u))
  Var_population_estimation = rbind(Var_population_estimation,apply(u,MARGIN = 2,FUN = var))
}
rownames(Mean_population_estimation) = Panel_data$Target
rownames(Var_population_estimation) = Panel_data$Target
```

Now we can compare a regular random sampling 

```r
Random_sampling = c()
for (k in 1:50) {
  m = Random_spatial_sampling(sce,width_FOV = 100,height_FOV = 100,
                              N_samplings = 10,Selected_image = 1,plot_result = FALSE)
  count_cells = table(factor(m$List_sampled_cluster,levels = unique(colLabels(sce))))
  count_cells = count_cells/sum(count_cells)
  Random_sampling = rbind(Random_sampling,count_cells)
}
Random_sampling = Random_sampling * 100

Mean_random = colMeans(Random_sampling)
Var_random = apply(Random_sampling,MARGIN = 2,FUN = var)

```

Finally we can compute the log2 variance ratio between random and stratified sampling :

```r
Var_ratio = apply(Var_population_estimation,MARGIN = 1,FUN = function(x) {log2(x/Var_random)})
Var_ratio_mean  = apply(Var_ratio,MARGIN = 2,FUN = mean)
```
The obtained matrix can be interpretated as following: negative values corresponds to a lower variance using stratified sampling compared to random sampling and thus stratified sampling is more efficient than random sampling.

Finally we can visualize the efficiency of the stratified sampling for each cell type and for each marker used :

```r
library(pheatmap)
pheatmap(-Var_ratio_mean,clustering_method = "ward",scale = "none")
```
