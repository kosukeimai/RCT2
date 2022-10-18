
###
### Sample size parameter calculations for detecting a specific alternative 
### 
###



#' Sample size parameter calculations for detecting a specific alternative 
#' 
#' This function calculates the parameters needed for the method to calculate sample size
#' references.
#' 
#' @param data A data frame containing the relevant variables. The names for the variables should be ``Z'' for the treatment assignment, ``Y'' for the treatment outcome, ``A'' for the treatment assignment mechanism, and ``id'' for the cluster ID. The variable for the cluster ID should be a factor.
#' 
#' @return A list of class \code{calpara} which contains the following item:
#' \item{sigmaw}{ The within-cluster variance of the potential outcomes, with the assumption that the all of the variances the same. }
#' \item{sigmab}{ The between-cluster variance of the potential outcomes, with the assumption that all of the variances are the same. }
#' \item{r}{ The intraclass correlation coefficient with respect to the potential outcomes. }
#' \item{sigma.tot}{ The total variance of the potential outcomes. }
#' \item{n.avg}{ The mean of the number of treated observations by cluster. }
#' 
#' @examples 
#' data(jd)
#' data_LTFC <- data.frame(jd$assigned, jd$pct0, jd$cdd6m, jd$anonale)
#' colnames(data_LTFC) <- c("Z", "A", "Y", "id")
#' var.LTFC <- calpara(data_LTFC)
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


calpara <- function(data){
  data <- data
  clusters <- unique(data$id)
  # number of clusters
  n.clusters <- length(clusters)
  

  A <- tapply(data$A, data$id, mean)
  
  A <- as.numeric(factor(A))
  uniqueA <- length(unique(A))
  
  est.Yj <-  array(dim=c(n.clusters,2))
  est.sigmaj <- array(dim=c(n.clusters,2))
  
  n1 <- tapply((1-data$Z), data$id, sum)
  n0 <- tapply(data$Z, data$id, sum)
  
  est.Yj[,1] <- tapply((1-data$Z)*data$Y, data$id, sum)/n1
  est.Yj[,2] <- tapply(data$Z*data$Y, data$id, sum)/n0
  
  cluster.length <- tapply(data$Y, data$id, length)
  est.Yj1.rep <- rep(est.Yj[,1], times = cluster.length)
  est.Yj2.rep <- rep(est.Yj[,2], times = cluster.length)
  est.sigmaj[,2] <- tapply( (data$Y-est.Yj2.rep)^2*data$Z, data$id, sum)/(n0-1)
  est.sigmaj[,1] <- tapply( (data$Y-est.Yj1.rep)^2*(1-data$Z), data$id, sum)/(n1-1)
  
  Ja <- table(A)
  
  sigmab1 <- rep(-1,uniqueA)
  sigmab0 <- rep(-1,uniqueA)
  est.Y1 <-  rep(-1,uniqueA)
  est.Y0 <-  rep(-1,uniqueA)
  
  Yj2 <- data.frame(cbind(est.Yj[,2], A))
  Yj2 <- Yj2[order(Yj2$A), ]
  Yj1 <- data.frame(cbind(est.Yj[,1], A))
  Yj1 <- Yj1[order(Yj1$A), ]
  est.Y1 <- tapply(Yj2$V1, Yj2$A, sum)/Ja
  est.Y0 <- tapply(Yj1$V1, Yj1$A, sum)/Ja
  
  A.count <- table(A)
  est.Y1.rep <- rep(est.Y1, times = A.count)
  est.Y0.rep <- rep(est.Y0, times = A.count)
  sigmab1 <- tapply( (Yj2$V1-est.Y1.rep)^2, Yj2$A, sum)/(Ja-1)
  sigmab0 <- tapply( (Yj1$V1-est.Y0.rep)^2, Yj1$A, sum)/(Ja-1)
  n <- tapply(data$Z, data$id, length)
  
  
  sigmaw <- mean(est.sigmaj)
  sigmab <- mean(c(sigmab1,sigmab0))- mean(1/n0-1/n)* sigmaw

  return(list(sigmaw = sigmaw,sigmab = sigmab,r=sigmab/(sigmab+sigmaw), sigma.tot = sigmab+sigmaw, n.avg = mean(n)))

}

