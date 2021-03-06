
###
### Sample size calculations for detecting a specific alternative 
### 
###



#' Sample size calculations for detecting a specific alternative 
#' 
#' This function calculates the sample size needed to detect a specific alternative hypothesis with a given power at a given significance level.
#' For the details of the method implemented by this function, see the
#' references.
#' 
#' @param data A data frame containing the relevant variables. The names for the variables should be ``Z'' for the treatment assignment, ``Y'' for the treatment outcome, ``A'' for the treatment assignment mechanism, and ``id'' for the cluster ID. The variable for the cluster ID should be a factor.

#' @param mu The effect size (i.e. the largest direct effect across treatment assignment mechanisms).
#' @param qa The proportions of different treatment assignment mechanisms.
#' @param alpha The given significance level (default 0.05).
#' @param beta The given power level (default 0.2).
#' 
#' @return A list of class \code{sampleSRE} which contains the following item:
#' \item{samplesize}{ A list of the calculated necessary nubmer of clusters for each assignment mechanism in order to detect a specific alternative with a given power at a given significance level. }
#'
#' 
#' 
#' @author Kosuke Imai, Department of Statistics, Harvard University
#' \email{imai@harvard.edu}, \url{https://imai.fas.harvard.edu/};
#' Zhichao Jiang, School of Public Health and Health Sciences, University of Massachusetts Amherst
#' \email{zhichaojiang@umass.edu};
#' Karissa Huang, Department of Statistics, Harvard College
#' \email{krhuang@college.harvard.edu}
#' @references Zhichao Jiang, Kosuke Imai (2020).
#' \dQuote{Statistical Inference and Power Analysis for Direct and Spillover Effects in Two-Stage Randomized Experiments}, \emph{Technical Report}.
#' @keywords two-stage randomized experiments
#' @name Calsamplesize
#' @import stats
#' @import quadprog
#' 
#' 
#' @export Calsamplesize
#' 







Calsamplesize <- function (data, mu, qa, alpha=0.05, beta=0.2){  
  
  var_data <- calpara(data)
  r <- var_data$r
  sigma <- var_data$sigma.tot
  pa <- sort(unique(data$A))
  n <- round(var_data$n.avg)
  
  
  m <-  length(qa)  
  D0 <-  array(0,dim=c(2*m,2*m))
  
  for (i in 1:m){
    D0[2*i-1,2*i-1] <-  (r+(1-pa[i])*(1-r)/n/pa[i])/qa[i]
    D0[2*i,2*i] <-  (r+pa[i]*(1-r)/n/(1-pa[i]))/qa[i]
    D0[2*i-1,2*i] <-  0
    D0[2*i,2*i-1] <-  0
  }
  
  ## contrast matrices
  C1 <-  array(0,dim=c(m,2*m))
  C2 <-  rep(0,2*m)
  
  for (i in 1:m){
    C1[i,2*i-1] <-  1
    C1[i,2*i] <-  -1
    C2[2*i-1] <-  qa[i]
    C2[2*i] <-  -qa[i]
  }
  
  C3 <-  array(0,dim=c(2*m-2,2*m))
  for ( i in 1:(m-1)){
    C3[i,2*i-1] <- 1
    C3[i,2*i+1] <- -1
    C3[m-1+i,2*i] <- 1
    C3[m-1+i,2*i+2] <- -1
  }
  
  #### sample size formulas
  
  ## direct effect
  s1 <-   Cals(m,alpha,beta)
  J.DE <-      s1*max(diag(C1 %*% D0 %*% t(C1)))*sigma/mu^2
  
  ## marginal direct effect
  
  s2 <-   Cals(1,alpha,beta)
  J.MDE <-      s2*t(C2) %*% D0 %*% (C2)*sigma/mu^2
  
  
  ##### spillover effect
  s3  <-   Cals(2*m-2,alpha,beta)
  
  J.SE <-  s3*sigma/quadprogSE(C3 %*% D0 %*% t(C3))/mu^2
  
  samplesize <-  c(J.DE,J.MDE,J.SE)
  class(samplesize) <- "sample"
  return (samplesize)   
}
