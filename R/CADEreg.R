###
### Regression-based method for the complier average direct effect
### 
###



#' Regression-based method for the complier average direct effect
#'
#' 
#' This function computes the point estimates of the complier average direct effect (CADE) and four
#'  different variance estimates: the HC2 variance, the cluster-robust variance, the cluster-robust HC2
#'  variance and the variance proposed in the reference. The estimators calculated using this function
#'  are cluster-weighted, i.e., the weights are equal for each cluster. To obtain the indivudal-weighted
#'  estimators, please multiply the recieved treatment and the outcome by \code{n_jJ/N}, where
#'  \code{n_j} is the number of individuals in cluster \code{j}, \code{J} is the number of clusters and 
#'  \code{N} is the total number of individuals. 
#'
#' 
#' For the details of the method implemented by this function, see the
#' references.
#' 
#' @param data  A data frame containing the relevant variables. The names for the variables should be: ``Z'' for the treatment assignment,  ``D''  for the actual received treatment, ``Y'' for the outcome, ``A'' for the treatment assignment mechanism and ``id'' for the cluster ID. The variable for the cluster id should be a factor.
#' @param ci.level A double between 0 and 1 specifying the confidence interval level to be output.
#' @return A list of class \code{CADEreg} which contains the following items:
#' \item{CADE1}{ The point estimate of CADE(1).  } \item{CADE0}{ The point estimate of CADE(0).  } 
#' \item{var1.clu}{ The cluster-robust variance of CADE(1).   } \item{var0.clu}{ The cluster-robust variance of CADE(0).  }
#' \item{var1.clu.hc2}{ The cluster-robust HC2 variance of CADE(1).   } 
#' \item{var0.clu.hc2}{ The cluster-robust HC2 variance of CADE(0).    } 
#' \item{var1.hc2}{ The  HC2 variance of CADE(1).    } 
#' \item{var0.hc2}{ The  HC2 variance of CADE(0).    } 
#' \item{var1.ind}{ The  individual-robust variance of CADE(1).    } 
#' \item{var0.ind}{ The  individual-robust variance of CADE(0).    } 
#' \item{var1.reg}{ The  proposed variance of CADE(1).    } 
#' \item{var0.reg}{ The  proposed variance of CADE(0).    } 
#' @author Kosuke Imai, Department of Statistics, Harvard University
#' \email{imai@harvard.edu}, \url{https://imai.fas.harvard.edu/};
#' Zhichao Jiang, School of Public Health and Health Sciences, University of Massachusetts Amherst
#' \email{zhichaojiang@umass.edu};
#' Karissa Huang, Department of Statistics, Harvard College
#' \email{krhuang@college.harvard.edu}
#' @references Kosuke Imai, Zhichao Jiang and Anup Malani (2018).
#' \dQuote{Causal Inference with Interference and Noncompliance in the Two-Stage Randomized Experiments}, \emph{Technical Report}. Department of Politics, Princeton
#' University.
#' @keywords two-stage randomized experiments
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @export CADEreg


CADEreg=function(data, ci.level=0.95){
  data <- data
  ## validate ci.level
  if(ci.level>1){stop('Please specify the confidence interval as a decimal between 0 and 1.')}
  if(!is.factor(data$id)){stop('The cluster_id should be a factor variable.')}
  
  ## pre-processing & calculating the weights
  data <- data[order(data$id), ]
  cluster.id <- unique(data$id)
  n.cluster <- length(cluster.id)
  n <- tapply(data$Z, data$id, length)
  
  
  data_sub <- data %>% select(.data$id, .data$A)
  A <- data_sub[!duplicated(data_sub$id), ]
  A <- A[order(A$id),]
  A <- A$A
  
  W <- rep(0, length(data$id))
  J1 <- sum(A)
  J0 <- n.cluster-J1
  n1 <- tapply(data$Z, data$id, sum)
  n0 <- n-n1
  
  index.l <- as.numeric( c(1, cumsum(n)+1) )
  index.l <- head(index.l, -1)
  index.r <- as.numeric( cumsum(n) )
  
  for(j in 1:n.cluster){
    index <- index.l[j]:index.r[j]
    Zj <- data$Z[data$id == cluster.id[j]]
    W[index] <- ifelse(A[j] == 1, 1/J1, 1/(n.cluster-J1)) * ifelse(Zj == 1, 1/n1[j], 1/n0[j])
  }

  
  # ## Design matrix in the first stage
  A.reg <- data$A
  Z.reg <- data$Z
  D.reg <- data$D
  Y.reg <- data$Y
  
  
  ### analysis

  X <- cbind(A.reg, 1-A.reg,  Z.reg*A.reg, Z.reg*(1-A.reg)  )

  reg1s <- lm(D.reg~0+X,weights=W)
  D.hat <- X%*%reg1s$coefficients

  M <- cbind(A.reg, 1-A.reg,  D.hat*A.reg, D.hat*(1-A.reg)  )
  reg2s <- lm(Y.reg~0+M,weights=as.vector(W))
  res <- Y.reg-cbind(A.reg, 1-A.reg,  D.reg*A.reg, D.reg*(1-A.reg)  )%*%reg2s$coefficients
  ###  variance

  ## cluster robust variance
  MM <- t(M)%*%diag(W)%*%M

  var.cluster.med <- array(0,dim=c(4,4))
  for( j in 1:n.cluster){
    index <- index.l[j]:index.r[j]
    Mj <- M[index,]
    if (A[j]==1){
      Sj <- cbind(W[index],   0, W[index]*D.hat[index],0)

      var.cluster.med <- var.cluster.med+t(Sj)%*%res[index] %*%t(t(Sj)%*%res[index])
    }else{
      Sj <- cbind(0,W[index],   0, W[index]*D.hat[index])

      var.cluster.med <- var.cluster.med+t(Sj)%*%res[index] %*%t(t(Sj)%*%res[index])

    }
  }
  var.cluster <- solve(MM)%*%var.cluster.med%*%solve(MM)


  ## cluster robust hc2 variance
  MM <- t(M)%*%diag(W)%*%M

  var.cluster.med <- array(0,dim=c(4,4))
  for( j in 1:n.cluster){
    index <- index.l[j]:index.r[j]
    Mj <- M[index,]
    if (A[j]==1){
      Sj <- cbind(W[index],   0, W[index]*D.hat[index],0)*sqrt(J1/(J1-1))

      var.cluster.med <- var.cluster.med+t(Sj)%*%res[index] %*%t(t(Sj)%*%res[index])
    }else{
      Sj <- cbind(0,W[index],   0, W[index]*D.hat[index])*sqrt((J0)/(J0-1))

      var.cluster.med <- var.cluster.med+t(Sj)%*%res[index] %*%t(t(Sj)%*%res[index])

    }
  }

  var.cluster.hc2 <- solve(MM)%*%var.cluster.med%*%solve(MM)


  ### individual robust hc2
  res.ind <- rep(0,sum(n))
  var.ind.med <- array(0,dim=c(4,4))
  for (j in 1:n.cluster){
    Zj <- data$Z[data$id == cluster.id[j]]
    index <- index.l[j]:index.r[j]
    adj1 <- sum(res[index]*Zj/sum(Zj))
    adj0 <- sum(res[index]*(1-Zj)/sum(1-Zj))

    res.ind[index] <- res[index] - ifelse(Zj==1,adj1,adj0)
  }

  for (j in 1:n.cluster){
    for (i in 1:n[j]){
      index <- index.l[j]-1+i
      var.ind.med <- var.ind.med+(M[index,])%*% t( M[index,])  *W[index]^2 * ifelse(Z.reg[index]==1, n1[j]/(n1[j]-1),(n0[j])/(n0[j]-1))*res.ind[index]^2
    }
  }
  var.ind <- solve(MM)%*%var.ind.med%*%solve(MM)


  ### traditional hc2 variance
  var.hc2.med <- array(0,dim=c(4,4))
  for (j in 1:n.cluster){
    for (i in 1:n[j]){
      index <- index.l[j]-1+i
      if (A[j]==1){
        constant <- ifelse(Z.reg[index]==1, J1*n1[j]/(J1*n1[j]-1),J1*n0[j]/(J1*n0[j]-1))
      }else{
        constant <- ifelse(Z.reg[index]==1, J0*n1[j]/(J1*n1[j]-1),J0*n0[j]/(J1*n0[j]-1))
      }

      var.hc2.med <- var.hc2.med+(M[index,])%*% t( M[index,])  *W[index]^2 * constant*res[index]^2
    }
  }
  var.hc2 <- solve(MM)%*%var.hc2.med%*%solve(MM)

  ## results
  est.CADE1 <- reg2s$coefficients[3]
  est.CADE0 <- reg2s$coefficients[4]
  var1.cluster <- var.cluster[3,3]
  var0.cluster <- var.cluster[4,4]
  var1.cluster.hc2 <- var.cluster.hc2[3,3]
  var0.cluster.hc2 <- var.cluster.hc2[4,4]
  var1.ind <- var.ind[3,3]
  var0.ind <- var.ind[4,4]
  var1.reg <- (1-J1/n.cluster)*var.cluster.hc2[3,3]+(J1/n.cluster)*var.ind[3,3]
  var0.reg <- (J1/n.cluster)*var.cluster.hc2[4,4]+(1-J1/n.cluster)*var.ind[4,4]
  var1.hc2 <- var.hc2[3,3]
  var0.hc2 <- var.hc2[4,4]

  ## confidence intervals
  ci.tail <- (1-ci.level)/2
  qnorm <- qnorm(ci.level+ci.tail)

  est.CADE1.CI95 <- c(est.CADE1-qnorm*var1.reg, est.CADE1+qnorm*var1.reg)
  est.CADE0.CI95 <- c(est.CADE1-qnorm*var0.reg, est.CADE1+qnorm*var0.reg)
  
  colnames(est.CADE1) <- NULL
  colnames(est.CADE0) <- NULL
  colnames(est.CADE1.CI95) <- NULL
  colnames(est.CADE0.CI95) <- NULL

  output <- list(CADE1=est.CADE1,CADE0=est.CADE0, var1.clu=var1.cluster,var0.clu=var0.cluster,var1.clu.hc2=var1.cluster.hc2,var0.clu.hc2=var0.cluster.hc2,
                var1.ind=var1.ind,var0.ind=var0.ind,var1.reg=var1.reg,var0.reg=var0.reg,var1.hc2=var1.hc2,var0.hc2=var0.hc2,
                CADE1.CI=est.CADE1.CI95, CADE0.CI=est.CADE0.CI95)
  class(output) <- "regression"

  return(output)
}










