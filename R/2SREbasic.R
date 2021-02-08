##### basic functions for use in 2SRE


get_unit_vec <- function(i, len){
  ei <- rep(0, len)
  ei[i] <- 1
  return(ei)
}

#######  Sample size formula
#### Calculate s^2 in the paper
library(stats)
library(quadprog)

Funs = function(x,m,alpha=0.05, beta=0.2 ){
  chi <-  qchisq(1-alpha, m)
  return (  pchisq(chi,m,ncp=x)  - beta    )
}

Cals = function(m, alpha=0.05, beta=0.2){
  uniroot(Funs,m=m, alpha=alpha,beta=beta,lower=0,upper=100)$root
}


####  minimize  s'(C3 D.hat C3)^{-1}s under max |ASE(z,a,a')|=1
## The function enumerates all the possible situations when max |ASE(z;a,a')|=1
## For example, one situation with m = 3  is  -1 \leq ASE(z;a,a') \leq 1 for all z a a' and
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


Calpara <- function(Z,A,Y){
  # number of clusters
  n.lea <- length(A)
  m <- length(unique(A))
  
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
  
  sigmab1 <- rep(-1,m)
  sigmab0 <- rep(-1,m)
  est.Y1 <-  rep(-1,m)
  est.Y0 <-  rep(-1,m)
  
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
  
  return(list(sigmaw = sigmaw,sigmab = sigmab,r=sigmab/(sigmab+sigmaw), n.avg = mean(n)))
  
}






### sample size calculations for the parameters
# Calpara <- function(Z,A,Y){
#   # number of clusters
#   n.lea <- length(A)
#   
#   est.Yj <-  array(dim=c(n.lea,2))
#   est.sigmaj <- array(dim=c(n.lea,2))
#   
#   for (j in 1:n.lea){
#     Z.sub <-  Z[[j]]
#     Y.sub <-  Y[[j]]
#     n1.sub <- sum(Z.sub)
#     n0.sub <- sum(1-Z.sub)
#     est.Yj[j,2] <- sum(Z.sub*Y.sub)/n1.sub
#     est.Yj[j,1] <- sum((1-Z.sub)*Y.sub)/n0.sub
#     est.sigmaj [j,2] <-   sum( (Y.sub-est.Yj[j,2])^2*Z.sub)/(n1.sub-1)
#     est.sigmaj [j,1] <-   sum( (Y.sub-est.Yj[j,1])^2*(1-Z.sub))/(n0.sub-1)
#   }
#   
#   Ja <- table(A)
#   
#   sigmab1 <- rep(-1,3)
#   sigmab0 <- rep(-1,3)
#   est.Y1 <-  rep(-1,3)
#   est.Y0 <-  rep(-1,3)
#   
#   for ( a in 1:3){
#     est.Y1[a] <- sum(est.Yj[,2]*(A==a))/Ja[a]
#     sigmab1[a] <-     sum((est.Yj[,2]-est.Y1[a])^2*(A==a))/(Ja[a]-1)
#     est.Y0[a] <- sum(est.Yj[,1]*(A==a))/Ja[a]
#     sigmab0[a] <-     sum((est.Yj[,1]-est.Y0[a])^2*(A==a))/(Ja[a]-1)
#   }
#   
#   n1 <- sapply(Z,sum)
#   n <- sapply(Z,length)
#   
#   
#   sigmaw <- mean(est.sigmaj)
#   sigmab <- mean(c(sigmab1,sigmab0))- mean(1/n1-1/n)* sigmaw
#   
#   return(list(sigmaw = sigmaw,sigmab = sigmab,r=sigmab/(sigmab+sigmaw), n.avg = mean(n)))
#   
# }
