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

