<img src="Screenshot/Balagan_name.jpeg" alt="Balagan_name.jpeg" width='500'>             <img src="Screenshot/Logo_v1.jpeg" alt="Logo_v1.jpeg" width='200' align="right">

Balagan is an R-package dedicated to the study of large-scale Multiplexed Imaging (MI) datasets. It contains several tools allowing to :

- Normalize and process the original data, i.e transforming and clustering the data.
- Infer the best parameters for an optimal spatial sampling strategy.
- Vizualise the data with various plotting functions.
- Study and quantify cell-cell interactions using various interaction scores derived from the spatial point pattern theory and tensor decomposition approaches.

Balagan is based on the **SingleCellExperiment** object structure and is therefore compatible with a variety of other single-cell analysis tools.
This package is aimed to provide advanced and robust statistical tool for the analysis of MI and therefore I **strongly recommend** the user to read all the mathematical papers mentionned in the documentation.

Scripts initially written for our [Nature Method paper](https://www.nature.com/articles/s41592-022-01692-z) have been integrated to this package.

# Installation

Balagan can be installed from the source file :

```r
install.packages(c("devtools","pheatmap","rTensor","spatgraphs","statmod","imager"))
devtools::install_local("Path/to/balagan_0.1.tar.gz",dependencies = T)
```

It can also be installed using devtools :

```r
devtools::install_github("PierreBSC/Balagan")
```

# List of tutorials

As mentionned above, Balagan can be used to perform various tasks and comprehensive tutorials have been written for each of them:
- [Basic MI data processing](https://github.com/PierreBSC/Balagan/blob/main/Tutorial_data_processing.md)
- [Basic spatial sampling analysis](https://github.com/PierreBSC/Balagan/blob/main/Tutorial_sampling.md) and the corresponding guide for [Visium data processing](https://github.com/PierreBSC/Balagan/blob/main/Processing_Visium_data.md)
- [Stratified spatial sampling tutorial](https://github.com/PierreBSC/Balagan/blob/main/Tutorial_stratified_sampling.md)


