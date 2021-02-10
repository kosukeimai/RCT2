
###
### Sample size parameter calculations for detecting a specific alternative 
### 
###



#' Sample size parameter calculations for detecting a specific alternative 
#' 
#' This function calculates the parameters needed for the method to calculate sample size
#' references.
#' 
#' 
#' @param Z A vector of the treatment assignments.
#' @param A A vector of treatment assignment mechanisms. 
#' @param Y A vector of potential outcomes.
#' 
#' @return A list of class \code{calpara} which contains the following item:
#' \item{sigmaw}{ The within-cluster variance of the potential outcomes, with the assumption that the all of the variances the same. }
#' \item{sigmab}{ The between-cluster variance of the potential outcomes, with the assumption that all of the variances are the same. }
#' \item{r}{ The intraclass correlation coefficient with respect to the potential outcomes. }
#' \item{sigma.tot}{ The total variance of the potential outcomes. }
#' \item{n.avg}{ The mean of the number of treated observations by cluster. }
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
#' @name calpara
#' @import stats
#' @import quadprog
#' 
#' 
#' @export calpara
#' 


calpara <- function(Z,A,Y){
  # number of clusters
  n.lea <- length(A)

  est.Yj <-  array(dim=c(n.lea,2))
  est.sigmaj <- array(dim=c(n.lea,2))

  for (j in 1:n.lea){
    Z.sub <-  Z[[j]]
    Y.sub <-  Y[[j]]
    n1.sub <- sum(Z.sub)
    n0.sub <- sum(1-Z.sub)
    est.Yj[j,2] <- sum(Z.sub*Y.sub)/n1.sub
    est.Yj[j,1] <- sum((1-Z.sub)*Y.sub)/n0.sub
    est.sigmaj [j,2] <-   sum( (Y.sub-est.Yj[j,2])^2*Z.sub)/(n1.sub-1)
    est.sigmaj [j,1] <-   sum( (Y.sub-est.Yj[j,1])^2*(1-Z.sub))/(n0.sub-1)
  }

  Ja <- table(A)

  sigmab1 <- rep(-1,3)
  sigmab0 <- rep(-1,3)
  est.Y1 <-  rep(-1,3)
  est.Y0 <-  rep(-1,3)

  for ( a in 1:3){
    est.Y1[a] <- sum(est.Yj[,2]*(A==a))/Ja[a]
    sigmab1[a] <-     sum((est.Yj[,2]-est.Y1[a])^2*(A==a))/(Ja[a]-1)
    est.Y0[a] <- sum(est.Yj[,1]*(A==a))/Ja[a]
    sigmab0[a] <-     sum((est.Yj[,1]-est.Y0[a])^2*(A==a))/(Ja[a]-1)
  }

  n1 <- sapply(Z,sum)
  n <- sapply(Z,length)


  sigmaw <- mean(est.sigmaj)
  sigmab <- mean(c(sigmab1,sigmab0))- mean(1/n1-1/n)* sigmaw

  return(list(sigmaw = sigmaw,sigmab = sigmab,r=sigmab/(sigmab+sigmaw), sigma.tot = sigmab+sigmaw, n.avg = mean(n)))

}

