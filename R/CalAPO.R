
###
### Point Estimation and Variance for the unit-level direct effect (ADE), marginal direct effect (MDE), and unit level spillover effect (ASE)
### 
###



#' Point Estimation and Variance for the unit-level direct effect (ADE), marginal direct effect (MDE), and unit level spillover effect (ASE)
#'
#' 
#' This function calculates the estimated average potential outcomes Y(z,a), point estimates for the ADE, MDE, and ASE, and conservative covariance matrix estimates.
#'
#' 
#' For the details of the method implemented by this function, see the
#' references.
#' 
#' @param Z A vector of the treatment assignments.
#' @param A A vector of treatment assignment mechanisms. 
#' @param Y A vector of potential outcomes.

#' @return A list of class \code{CalAPO} which contains the following items:
#' \item{Y.hat}{ Estimate of the average potential outcomes. }
#' \item{ADE.est}{ Estimate of the unit level direct effect. }
#' \item{MDE.est}{ Estimate of the marginal direct effect. }
#' \item{ASE.est}{ Estimate of the unti level spillover effect. }
#' \item{cov.hat}{ Conservative covariance matrix for the estimated potential outcomes. }
#' \item{var.hat.ADE}{ Estimated variance of the ADE. }
#' \item{var.hat.MDE}{ Estimated variance of the MDE. }
#' \item{var.hat.ASE}{ Estimated variance of the ASE. }
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
#' @export CalAPO


CalAPO <- function (Z, A, Y){
  ### format the data
  Ja <- table(A)
  J <- sum(Ja)
  m <- length(Ja)
  
  #######  point estimators
  n1 <- sapply(Z,sum)
  n <- sapply(Z,length)
  n0 <- n-n1
  
  Y1j.hat <- rep(0,J)
  Y0j.hat <- rep(0,J)
  
  
  for (j in 1:J){
    Y1j.hat[j] <- sum(Y[[j]]*Z[[j]])/n1[j]
    Y0j.hat[j] <- sum(Y[[j]]*(1-Z[[j]]))/n0[j]
    
  }
  
  
  ###  estiamted APO   (Y(1,1),Y(0,1)...,Y(1,m),Y(0,m))
  Y.hat <-  rep(0,2*m)
  for (a in 1:m){
    Y.hat[2*a-1] <-  sum(Y1j.hat* (A==a))/Ja[a]
    Y.hat[2*a] <-  sum(Y0j.hat* (A==a))/Ja[a]
  }
  
  ### estimated ADE
  C1 <- matrix(data = NA, nrow = m, ncol = 2*m)
  for (i in 1:m){
    ind <- 2*i-1
    eind <- get_unit_vec(ind, 2*m)
    eindplus1 <- get_unit_vec(ind + 1, 2*m)
    C1[i, ] <- eind-eindplus1
  }
  C1 <- t(C1)
  
  ADE <- C1%*%Y.hat
  
  ### estimated MDE
  C2 <- rep(0, 2*m)
  j <- 1
  for (i in 1:(2*m)){
    
    if(i %% 2 == 1){
      C2[i] <- Ja[j]
      j <- j+1
    }else if(i%%2 == 0){
      C2[i] <- -C2[i-1]
    }

  }
  
  MDE <- C2/J*Y.hat
  
  ### estimated ASE
  C3 <- matrix(0, nrow = 2*m-2, ncol = 2*m)
  for(i in 1:(2*m-2)){
    if(i%%2 == 1){
      C3[i, ] <- get_unit_vec(i, 2*m)-get_unit_vec(i+2, 2*m)
    }else if(i%%2 == 0){
      C3[i, ] <- get_unit_vec(i, 2*m)-get_unit_vec(i+2, 2*m)
    }
  }
  C30 <- matrix(0, nrow = m-1, ncol = 2*m)
  C31 <- matrix(0, nrow = m-1, ncol = 2*m)
  j <- 1 
  for(i in 1:(2*m-2)){
    
    if(i%%2 == 1){
      C30[j, ] <- C3[i,]
      j <- j+1
    }
  }
  k <- 1
  for(i in 1:(2*m-2)){
    if(i%%2 == 0){
      C31[k, ] <- C3[i,]
      k <- k+1
    }
  }
  
  C3 <- rbind(C30, C31)
  ASE <- C3 %*% Y.hat
  
  
  
  ## covariance matrix est
  cov.hat <-  array(0,dim=c(2*m,2*m))
  
  
  for (a in 1:m){
    cov.hat[2*a-1,2*a-1] <-   1/Ja[a]*sum((Y1j.hat-Y.hat[2*a-1] )^2*(A==a))/(Ja[a]-1)
    cov.hat[2*a,2*a] <-   1/Ja[a]*sum((Y0j.hat-Y.hat[2*a] )^2*(A==a))/(Ja[a]-1)
    cov.hat[2*a,2*a-1] <-   1/Ja[a]*sum((Y1j.hat-Y.hat[2*a-1] )*(Y0j.hat-Y.hat[2*a])*(A==a))/(Ja[a]-1)
    cov.hat[2*a-1,2*a] <- cov.hat[2*a,2*a-1]
  }
  
  ### covariance matrix of ADE, MDE, ASE
  hat.D <- cov.hat*J
  var.hat.ADE <- C1%*%hat.D%*%t(C1)/J
  var.hat.MDE <- C2%*%hat.D%*%t(C2)/J
  var.hat.ASE <- C3%*%hat.D%*%t(C3)/J
  
  return(list(Y.hat=Y.hat, ADE.est = ADE, MDE.est = MDE, ASE.est = ASE, cov.hat = cov.hat, var.hat.ADE = var.hat.ADE, var.hat.MDE = var.hat.MDE, var.hat.ASE = var.hat.ASE))
}
