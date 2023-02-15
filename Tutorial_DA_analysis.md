#Differential Abundance (DA) analysis of cell count data 

The final stage of multiplexed imaging data analysis usually consists in the comparison of specific cell type abundance between groups of samples. To do so, one can use Negative-Binomial (NB) or Beta-Binomial (BB) based regressions to detect difference of abundances. Balagan implements basic differential tests but for more complicated cases, the original functions contained in the packages [aod](https://cran.r-project.org/web/packages/aod/) and [MASS](https://cran.r-project.org/web/packages/MASS/) should be used.

##Basic differential abundance test using NB test

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

 


##Advanced differential abundance test

##Mixed model for clinical cohort analysis

