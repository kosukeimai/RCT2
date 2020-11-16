
###
### Hypothesis testing for three null hypotheses
### 
###



#' Sample size calculations for 
#'
#' 
#' This function tests the null hypotheses of no direct effect, no marginal direct effect, and no spillover effect.
#' 
#' For the details of the method implemented by this function, see the
#' references.
#' 
#' @param data  A data frame containing the relevant variables. The names for the variables should be: ``Z'' for the treatment assignment,  ``D''  for the actual received treatment, ``Y'' for the outcome, ``A'' for the treatment assignment mechanism and ``id'' for the cluster ID. The variable for the cluster id should be a factor.
#' @param effect Specify which null hypothesis to be tested. ``DE'' for direct effect, ``ME'' for marginal effect, and ``SE'' for spillover effect.
#' @param alpha The level of significance at which the test is to be run (default is 0.05).
#' @return A list of class \code{Test2SRE} which contains the following item:
#' \item{rej}{ Rejection region for test conducted. }
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
#' 
#' 
#' @export Test2SRE







#######  Sample size formula
#### Calculate s^2 in the paper
library(stats)

Funs = function(x,m,alpha=0.05, beta=0.2 ){
  chi <-  qchisq(1-alpha, m)
  return (  pchisq(chi,m,ncp=x)  - beta    )
}

Cals = function(m, alpha=0.05, beta=0.2){
  uniroot(Funs,m=m, alpha=alpha,beta=beta,lower=0,upper=100)$root
}


####  minimize  s'(C3 D.hat C3)^{-1}s under max |ASE(z,a,a')|=1
## The function enumerates all the possible situations when max |ASE(z;a,a')|=1
## For example, one situtation with m = 3  is  -1 \leq ASE(z;a,a') \leq 1 for all z a a' and
##   ASE(0;1,2) =1

### this may not be efficient 
quadprogSE = function(D.se){
  m <-  (dim(D.se)[1]+2)/2
  
  # quadratic function
  D <-  2*solve(D.se)
  
  # construct the constraint  A s \geq s0
  A.sub <-  NULL
  for ( i in 1:(m-1)){
    for (j in i:(m-1)){
      A.sub <- rbind(A.sub, c(rep(0,i-1), rep(1,j+1-i),rep(0,m-1-j)))
    }
  }
  zeros <-  array(0,dim=c(m*(m-1)/2,m-1))
  A1 <-  cbind(A.sub,zeros)
  A1 <-  rbind(A1,-A1)
  A0 <-  cbind(zeros,A.sub)
  A0 <-  rbind(A0,-A0)
  A <-  rbind(A1,A0)
  min <-  100
  for ( i in 1:(m*(m-1)/2)){
    s0 <-  rep(-1,2*m*(m-1))
    s0[i] <-  s0[i]+2
    res <- solve.QP(Dmat=D, dvec = rep(0,2*m-2),Amat = t(A),bvec= s0)
    if (res$value<min){min <-  res$value}
    s0 <-  rep(-1,2*m*(m-1))
    s0[m*(m-1)+i] <-  s0[m*(m-1)+i]+2
    res <- solve.QP(Dmat=D, dvec = rep(0,2*m-2),Amat = t(A),bvec= s0)
    if (res$value<min){min <-  res$value}
  }
  return(min)
}



###  sigma: the total variance
### r is the intra-class correlationn coefficient
### mu is the effect size
### pa is the treated propotion under different treatment assignment mechnisms
### qa is the proportion of different treatment assignment mechnisms

Calsamplesize = function (mu,n,qa, pa, r, sigma=1, alpha=0.05, beta=0.2){  
  m <-  length(qa)  
  D0 <-  array(0,dim=c(2*m,2*m))
  
  for (i in 1:m){
    D0[2*i-1,2*i-1] <-  (r+(1-pa[i])*(1-r)/n/pa[i])/qa[i]
    D0[2*i,2*i] <-  (r+pa[i]*(1-r)/n/(1-pa[i]))/qa[i]
    D0[2*i-1,2*i] <-  0
    D0[2*i,2*i-1] <-  0
  }
  
  ## constrast matrices
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
  return (samplesize)	    
}
