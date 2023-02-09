#' @rdname BB_DA_test
#' @title Differential abundance analysis using a Beta Binomial (BB) based model
#'
#' @description Performs a differential abundance analysis between groups of FoVs. Unlike the NB based test, the FoVs do not need to have the same size.
#'
#' @param Count_vector a numerical vector containing the list of count values for a given cell of interest.
#' @param Total_cell Total number of cells in each FoV 
#' @param Condition_vector a character vector of factor of the same size as Count_vector that describes to which group each sample belongs.
#' @param type_test the type of statistical test to perform in order to compute the statistical significance. Has to be either LRT (Likelihood Ratio Test) or Wald (Wald's test)
#' @return  P-value for differential abundance between the groups of sample. If LRT has been selected only one p-value will be provided. If Wald test was selected one p-value per condition - 1 will be provided
#' @examples
#' Count_vector = c(rnbinom(n=100,mu=2,size=1),rnbinom(n=50,mu=6,size=1))
#' Condition_vector = c(rep("Control",100),rep("Treatment",50))
#' Total_cell = rep(200,150)
#'BB_DA_test(Count_vector,Total_cell,Condition_vector)
#' @import aod 
#' @export


BB_DA_test = function(Count_vector,Total_cell,Condition_vector,type_test = "LRT") {
  
  if (!type_test%in%c("LRT","Wald")) {
    warning("No proper name for the type of test has been provided. By default a Likelihood Ratio Test (LRT) will be perfomred !")
    type_test = "LRT"
  }
  
  data_BB =data.frame(n=Total_cell,y=Count_vector,Condition_vector=Condition_vector)
  
  model_1 =betabin(cbind(y,n-y)~Condition_vector,data=data_BB,~Condition_vector)
  model_0 =betabin(cbind(y,n-y)~1,data=data_BB,~Condition_vector)
  
  
  if (type_test=="LRT"){
    LRT_test_BB = anova(model_0,model_1)
    Pvalue_test = LRT_test_BB@anova.table$`P(> Chi2)`[2]
  }
  
  if (type_test=="Wald"){
    Pvalue_test = c()
    for (k in 2:length(unique(Condition_vector))) {
      Wald_test_BB = wald.test(b = coef(model_1),Sigma = vcov(model_1,comlete=F),Terms = k)
      Pvalue_test = c(Pvalue_test,Wald_test_BB$result$chi2[3])
    }
    names(Pvalue_test) = unique(Condition_vector)[-1]
  }
  return(Pvalue_test)
}
