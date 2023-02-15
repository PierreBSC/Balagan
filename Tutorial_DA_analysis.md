# Differential Abundance (DA) analysis of cell count data 

The final stage of multiplexed imaging data analysis usually consists in the comparison of specific cell type abundance between groups of samples. To do so, one can use Negative-Binomial (NB) or Beta-Binomial (BB) based regressions to detect difference of abundances. Balagan implements basic differential tests but for more complicated cases, the original functions contained in the packages [aod](https://cran.r-project.org/web/packages/aod/) and [MASS](https://cran.r-project.org/web/packages/MASS/) should be used.

## Basic differential abundance test using NB test

We start by a very simple usecase where two groups of observations have been sampled. In the case of multiplexed imaging, this will corresponds to Fields of View (FoVs) coming from two biological samples. We assume that all FoVs have the same size.

We start by loading the balagan library and by creating simulated data :

```r
library(balagan)
Sample_1 = rnbinom(n = 100,mu = 10,size = 2)
Sample_2 = rnbinom(n = 100,mu = 20,size = 2)
Count_vector = c(Sample_1,Sample_2)
Condition_vector = c(rep("Sample 1",100),rep("Sample 2",100))
```
We can now perform the DA analysis based on a NB distribution with **Likelihood Ratio Test** (LRT) or a **Wald's test**:

```r
NB_DA_test(Count_vector,Condition_vector,type_test = "LRT")
NB_DA_test(Count_vector,Condition_vector,type_test = "Wald")
```

When only two groups of samples are compared, LRT and Wald's test will usually display similar results. However when more than two groups of samples are studied, and that the statistical significance of more than one term can be tested, the two tests differ :

- The LRT compares the full model to the null model, similarly to an ANOVA, and the significance of each term (i.e the effect of each group) is not tested. Thus only one p-value is provided.
- The Wald's test identifies terms that are significantly different from zero, thus if N groups of samples are compared, N-1 p-values will be computed.


## Basic differential abundance test using BB test

DA analysis based on the BB distribution can be performed similarly to NB test, the only difference being that we do not need to assume similar FoV size and that we have to provide the total number of cells present in the FoV. Again we simulate two groups of sample :

```r
P_sample_1 = rbeta(n = 100,shape1 = 2,shape2 = 20)
P_sample_2 = rbeta(n = 100,shape1 = 6,shape2 = 10)

Total_cell_number = 400
Sample_1 = rbinom(n = 100,size = Total_cell_number,prob =P_sample_1 )
Sample_2 = rbinom(n = 100,size = Total_cell_number,prob =P_sample_2 )

Count_vector = c(Sample_1,Sample_2)
Condition_vector = c(rep("Sample 1",100),rep("Sample 2",100))
```

Then, a BB based test can be performed, either using a LRT or Wald's test :

```r
BB_DA_test(Count_vector,Total_cell_number,Condition_vector,type_test = "LRT")
BB_DA_test(Count_vector,Total_cell_number,Condition_vector,type_test = "Wald")
```

## Advanced differential abundance test

While the build-in balagan functions can be applied easily, they only work for very simple experimental design : one categorical factor with multiple groups. If one wants to include continuous covariate or multiple categorical factors, then native functions from the packages [aod](https://cran.r-project.org/web/packages/aod/) and [MASS](https://cran.r-project.org/web/packages/MASS/) have to be used. We strongly advise the user to familiarize themselve with the **Generalized Linear Model** (GLM) mathematical and R framework. A very nice introduction can be found in the excellent book ["An Introduction to Statistical Learning"](https://www.statlearning.com) (Chapter 4).

### NB based approach
Complex DA analysis based on NB distribution are performed through the **glm.nb()** function from the **MASS** package.
Here we will simulate cell count data coming from FoVs of variable size 

```r
FoV_size_1 = (rnorm(n = 100,mean = 0.8,sd = 0.1))^2
FoV_size_2 = (rnorm(n = 100,mean = 0.8,sd = 0.1))^2

Sample_1 = rnbinom(n = 100,mu = FoV_size_1*20,size = 2)
Sample_2 = rnbinom(n = 100,mu = FoV_size_2*100,size = 2)

Count_vector = c(Sample_1,Sample_2)
FoV_size_vector = c(FoV_size_1,FoV_size_2)
Condition_vector = c(rep("Sample 1",100),rep("Sample 2",100))
```

We can now 

## Mixed model for clinical cohort analysis

