
###
### Hypothesis testing for three null hypotheses
### 
###



#' Hypothesis testing for three null hypotheses
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


### Testing the hypotheses of DE=0,MDE=0, SE=0 
Test2SRE <- function(data,effect = "DE", alpha = 0.05){
  ### read in and format the data
  cluster.id=unique(data$id)	
  n.cluster=length(cluster.id)	
  Z=vector("list", n.cluster) 	
  Y=vector("list", n.cluster)
  A=rep(0,n.cluster)
  for (i in 1:n.cluster){
    Z[[i]]=as.numeric(data$Z[data$id==cluster.id[i]])
    Y[[i]]=data$Y[data$id==cluster.id[i]]
    if (length(unique(data$A[data$id==cluster.id[i]]))!=1){
      stop( paste0('The assignment mechanism in cluster ',i,' should be the same.'))
    }
    A[i]=data$A[data$id==cluster.id[i]][1]
  }
  
  Ja <- table(A)
  J <- sum(Ja)
  qa <- Ja/J
  est <-  CalAPO(Z,A,Y)
  m <- length(Ja)
  C1 <-  array(0,dim=c(m,2*m))
  C2 <-  rep(0,2*m)
  
  for (a in 1:m){
    C1[a,2*a-1] <-  1
    C1[a,2*a] <-  -1
    C2[2*a-1] <-  qa[a]
    C2[2*a] <-  -qa[a]
  }
  
  C3 = array(0,dim=c(2*m-2,2*m))
  for ( a in 1:(m-1)){
    C3[a,2*a-1] <- 1
    C3[a,2*a+1] <- -1
    C3[m-1+a,2*a] <- 1
    C3[m-1+a,2*a+2] <- -1
  }
  
  
  if ( effect == "DE"){ 
    rej <-  (t(C1%*%est$Y.hat)%*%solve(C1%*%est$cov.hat%*%t(C1))%*%(C1%*%est$Y.hat)) > qchisq(1-alpha,m)
  }
  if (effect == "MDE"){
    rej <-  (sum(C2*est$Y.hat))^2/(t(C2)%*%est$cov.hat%*%(C2)) > qchisq(1-alpha,1)
  }
  if (effect == "SE"){
    rej <-  (t(C3%*%est$Y.hat)%*%solve(C3%*%est$cov.hat%*%t(C3))%*%(C3%*%est$Y.hat)) >  qchisq(1-alpha,2*m-2)
  }
  
  return(drop(rej))
}