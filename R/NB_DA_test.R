#' @rdname NB_DA_test
#' @title Differential abundance analysis using a Negative Binomial (NB) based model
#'
#' @description Performs a differential abundance analysis between groups of FoVs. It is important to note that the size of the FoV is assumed to be equal here ! 
#'
#' @param Count_vector a numerical vector containing the list of count values
#' @param Condition_vector a character vector of factor of the same size as Count_vector that describes to which group each sample belongs.
#' @param type_test the type of statistical test to perform in order to compute the statistical significance. Has to be either LRT (Likelihood Ratio Test) or Wald (Wald's test)
#' @return  P-value for differential abundance between the groups of sample. If LRT has been selected only one p-value will be provided. If Wald test was selected one p-value per condition - 1 will be provided
#' @examples
#' Count_vector = c(rnbinom(n=100,mu=2,size=1),rnbinom(n=50,mu=6,size=1))
#' Condition_vector = c(rep("Control",100),rep("Treatment",50))
#'NB_DA_test(Count_vector,Condition_vector)
#' @import MASS 
#' @export


NB_DA_test = function(Count_vector,Condition_vector,type_test = "LRT") {
  
  m =glm.nb(Count_vector~Condition_vector)
  m_0 =glm.nb(Count_vector~1)
  
  if (!type_test%in%c("LRT","Wald")) {
    warning("No proper name for the type of test has been provided. By default a Likelihood Ratio Test (LRT) will be perfomred !")
    type_test = "LRT"
  }
  
  if (type_test=="LRT"){
    Anova_NB=anova(m,m_0,test="Chisq")
    Pvalue_test = Anova_NB$`Pr(Chi)`[2]
  }
  
  if (type_test=="Wald"){
    summary_m = coef(summary(m))
    Wald_test_NB = summary_m[-1,1]^2/summary_m[-1,2]^2
    Pvalue_test = pchisq(df = 1,q = Wald_test_NB,lower.tail = F)
  }
  return(Pvalue_test)
}
