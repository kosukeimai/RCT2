###
### Randomization-based method for the complier average direct effect and the complier average spillover effect
### 
###



#' Randomization-based method for the complier average direct effect and the complier average spillover effect
#'
#' 
#' This function computes the point estimates and variance estimates of the complier average direct effect (CADE)  and the complier average spillover effect (CASE).
#' The estimators calculated using this function are either individual weighted or cluster-weighted. The point estimates and variances of ITT effects are also included. 
#'
#' 
#' For the details of the method implemented by this function, see the
#' references.
#' 
#' @param data  A data frame containing the relevant variables. The names for the variables should be: ``Z'' for the treatment assignment,  ``D''  for the actual received treatment, ``Y'' for the outcome, ``A'' for the treatment assignment mechanism and ``id'' for the cluster ID. The variable for the cluster id should be a factor.
#' @param individual  A binary variable with TRUE for  individual-weighted estimators and FALSE for cluster-weighted estimators.
#' @param ci A numeric variable between 0 and 1 for the level of the confidence interval to be returned.
#' @return A list of class \code{CADErand} which contains the following items:
#' \item{CADE1}{ The point estimate of CADE(1).  } \item{CADE0}{ The point estimate of CADE(0).  } 
#'\item{CASE1}{ The point estimate of CASE(1).  } \item{CASE0}{ The point estimate of CASE(0).  } 
#'\item{var.CADE1}{ The  variance estimate of CADE(1).   } 
#'\item{var.CADE0}{ The  variance estimate of CADE(0).   } 
#'\item{var.CASE1}{ The  variance estimate of CASE(1).   } 
#'\item{var.CASE0}{ The  variance estimate of CASE(0).   } 
#'\item{DEY1}{ The point estimate of DEY(1).  } \item{DEY0}{ The point estimate of DEY(0).  } 
#'\item{DED1}{ The point estimate of DED(1).  } \item{DED0}{ The point estimate of DED(0).  } 
#'\item{var.DEY1}{ The  variance estimate of DEY(1).   } 
#'\item{var.DEY0}{ The  variance estimate of DEY(0).   } 
#'\item{var.DED1}{ The  variance estimate of DED(1).   } 
#'\item{var.DED0}{ The  variance estimate of DED(0).   } 
#'\item{SEY1}{ The point estimate of SEY(1).  } \item{SEY0}{ The point estimate of SEY(0).  } 
#'\item{SED1}{ The point estimate of SED(1).  } \item{SED0}{ The point estimate of SED(0).  } 
#'\item{var.SEY1}{ The  variance estimate of SEY(1).   } 
#'\item{var.SEY0}{ The  variance estimate of SEY(0).   } 
#'\item{var.SED1}{ The  variance estimate of SED(1).   } 
#'\item{var.SED0}{ The  variance estimate of SED(0).   } 
#'\item{NADE1}{ The point estimate of NADE(1).  } \item{NADE0}{ The point estimate of NADE(0).  }
#'\item{NASE1}{ The poinst estimate of NASE(1).  } \item{NASE0}{ The point estimate of NADE(0).  }
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
#' @export CADErand


CADErand<-function(data,individual=1, ci = 0.95){
	## transform the data into list 	
  if(!is.factor(data$id)){stop('The cluster_id should be a factor variable.')}
  cluster.id<-unique(data$id)	
  n.cluster<-length(cluster.id)	
  Z<-vector("list", n.cluster) 	
  D<-vector("list", n.cluster) 
  Y<-vector("list", n.cluster) 
  A<-rep(0,n.cluster)
  for (i in 1:n.cluster){
  	Z[[i]]<-as.numeric(data$Z[data$id==cluster.id[i]])
  	D[[i]]<-as.numeric(data$D[data$id==cluster.id[i]])
  	Y[[i]]<-data$Y[data$id==cluster.id[i]]
  	if (length(unique(data$A[data$id==cluster.id[i]]))!=1){
  		stop( paste0('The assignment mechanism in cluster ',i,' should be the same.'))
  	}
  	A[i]<-data$A[data$id==cluster.id[i]][1]
  }

	
	
  n<-sapply(Z,length)
  N<-sum(n)
  J<-length(n)
  if (individual==1){
    for ( i in 1:n.cluster){
      Y[[i]]<-Y[[i]]*n[i]*J/N
      D[[i]]<-D[[i]]*n[i]*J/N
    }
  }

  A2 <- factor(A)
  A_mat <- model.matrix(~A2-1)
  init <- sapply(Difflist(D,Productlist(D,Z)),sum)/(n-sapply(Z,sum))
  
  
  Dj0 <- init*A_mat # this is est.Dj00 and est.Dj01
  Dj0 <- as.data.frame(Dj0)
  
  # we do the same thing for the case z=1
  Dj1init <- sapply(Productlist(D,Z),sum)/(sapply(Z,sum))
  
  Dj1 <- Dj1init*A_mat # this is est.Dj10 and est.Dj11
  Dj1 <- as.data.frame(Dj1)
  
  # next step up averaging analogous
  # z = 0
  D0 <- t(as.data.frame(colSums(Dj0)/colSums(A_mat))) # this is est.D00 and est.D01
  D1 <- t(as.data.frame(colSums(Dj1)/colSums(A_mat))) # this is est.D10 and est.D11
  D0 <- rbind(D0, D0[rep(1, nrow(A_mat)-1), ]) # change nrow(DEDj) to something more general
  D1 <- rbind(D1, D1[rep(1, nrow(A_mat)-1), ])
  
  # analogous DED calculations 
  DEDj <- as.data.frame(Dj1-Dj0 )# this is est.DEDj0 and est.DEDj1
  DED <-as.data.frame(D1-D0) # this is est.DED0 and est.DED1
  
  
  # make sure that this is correct
  lenSED <- ncol(A_mat)*(ncol(A_mat)-1)/2
  SED0 <- rep(0, lenSED)
  SED1 <- rep(0, lenSED)
  for(i in 2:lenSED){
    SED0[i-1] <- D0[1, i]-D0[1, i-1]
    SED1[i-1] <- D1[1, i]-D1[1, i-1]
  }
  
  est.SED <- c(SED0, SED1)
  
  # analogous variance calculations
  # DED <- rbind(DED, DED[rep(1, nrow(DEDj)-1), ])
  est.xiDE <- colSums((DEDj - DED)^2*(A_mat))/(colSums(A_mat)-1) # est.xiDE0 and est.xiDE1
  est.xib0 <- colSums((Dj0-D0)^2*A_mat)/(colSums(A_mat)-1) # est.xib00 and est.xib01
  est.xib1 <- colSums((Dj1-D1)^2*A_mat)/(colSums(A_mat)-1) # est.xib10 and est.xib11
  
  
  est.xij0 <- matrix(0, nrow=ncol(A_mat), ncol=J) # est.xij00 and est.xij01
  est.xij1 <- matrix(0, nrow=ncol(A_mat), ncol=J) # est.xj10 and est.xij11
  
  for (j in 1:J){
    tmp1 <- Dj0[j, ]
    tmp2 <- rbind(tmp1, tmp1[rep(1,length(D[[j]])-1), ])
    est.xij0[, j] <- unlist(colSums((D[[j]]-tmp2)^2*(1-Z[[j]]))/(sum(1-Z[[j]])-1)) # est.xij00 and est.xij01
    
    tmp3 <- Dj1[j, ]
    tmp4 <- rbind(tmp3, tmp3[rep(1, length(D[[j]])-1), ])
    est.xij1[, j] <- unlist(colSums((D[[j]]-tmp4)^2*(Z[[j]]))/(sum(Z[[j]])-1)) # est.xj10 and est.xij11
  }
  
  est.xij0 <- as.data.frame(est.xij0)
  est.xij1 <- as.data.frame(est.xij1)
  
  
  Z_sum <- sapply(Z, sum)
  n_Z_sum <- n-sapply(Z, sum)
  denom1 <- t(matrix( Z_sum , length(Z_sum) , nrow(est.xij0) ))
  denom2 <- t(matrix( n_Z_sum , length(n_Z_sum) , nrow(est.xij0) ))
  
  var.DED <- est.xiDE*(1/colSums(A_mat)-1/J)+rowSums((est.xij0/denom2+est.xij1/denom1)*t(A_mat))/J/colSums(A_mat) # var.DED0 and var.DED1
  var.SED1 <- sum(est.xib1/colSums(A_mat))
  var.SED0 <- sum(est.xib0/colSums(A_mat))
  var.SED <- c(var.SED0, var.SED1)
  est.Yj0 <- sapply(Difflist(Y,Productlist(Y,Z)),sum)/(n-sapply(Z,sum))*A_mat # est.Yj00 and est.Yj01
  est.Yj1 <- sapply(Productlist(Y,Z),sum)/(sapply(Z,sum))*A_mat # est.Yj10 and est.Yj11
  est.Y0 <- colSums(est.Yj0*A_mat)/colSums(A_mat) # est.Y00 and est.Y01
  est.Y1 <- colSums(est.Yj1*A_mat)/colSums(A_mat) # est.Y10 and est. Y11
  est.DEYj <- est.Yj1-est.Yj0 # est.DEYj0 and est.DEYj1
  est.DEY <- est.Y1-est.Y0 # est.DEY0 and est.DEY1
  lenSEY <- ncol(A_mat)*(ncol(A_mat)-1)/2
  SEY0 <- rep(0, lenSEY) # this is est.SEY0
  SEY1 <- rep(0, lenSEY) # this is est.SEY1
  for(i in 2:lenSEY){
    SEY0[i-1] <- est.Y0[i]-est.Y0[i-1]
    SEY1[i-1] <- est.Y1[i]-est.Y1[i-1]
  }
  est.SEY <- c(SEY0, SEY1)
  
  
  
  
  # variance for spillover effects
  est.DEY <- t(matrix( est.DEY , length(est.DEY) , nrow(est.DEYj) ))
  est.Y0 <- t(matrix( est.Y0 , length(est.Y0) , nrow(est.Yj0) ))
  est.Y1 <- t(matrix( est.Y1 , length(est.Y1) , nrow(est.Yj1) ))
  
  est.sigmaDE <- colSums((est.DEYj-est.DEY)^2*A_mat)/(colSums(A_mat)-1) # est.sigmaDEO and est.sigmaDE1
  est.sigmab0 <- colSums((est.Yj0-est.Y0)^2*A_mat)/(colSums(A_mat)-1) # est.sigmab00 and sigmab01
  est.sigmab1 <- colSums((est.Yj1-est.Y1)^2*A_mat)/(colSums(A_mat)-1) # est.sigmab10 and est.sigmab11
  est.sigmaj0 <- matrix(0, 2, J) # est.sigmaj00 and est.sigmaj01
  est.sigmaj1 <- matrix(0, 2, J) # est.sigmaj10 and est.sigma11
  
  for(j in 1:J){
    tmp1 <- est.Yj0[j, ]
    tmp2 <- t( matrix( tmp1 , length(tmp1) , length(Y[[j]]) ) )
    est.sigmaj0[, j] <- colSums((Y[[j]]-tmp2)^2*(1-Z[[j]]))/(sum(1-Z[[j]])-1)
    
    tmp3 <- est.Yj1[j, ]
    tmp4 <- t( matrix( tmp3 , length(tmp3) , length(Y[[j]]) ) )
    est.sigmaj1[, j] <- colSums((Y[[j]]-tmp4)^2*(Z[[j]]))/(sum(Z[[j]])-1)
  }
  
  est.sigmaj0 <- as.data.frame(est.sigmaj0)
  est.sigmaj1 <- as.data.frame(est.sigmaj1)
  
  # var.DEY0 and var.DEY1
  var.DEY <- est.sigmaDE*(1/colSums(A_mat)-1/J)+rowSums( ( t(apply( est.sigmaj0, 1, "/", n-sapply(Z, sum)))+t(apply( est.sigmaj1, 1, "/", sapply(Z, sum))))*t(A_mat))/J/colSums(A_mat)
  # var.SEY0 and var.SEY1
  var.SEY <- c(sum(est.sigmab0/colSums(A_mat)), sum(est.sigmab1/colSums(A_mat)))
  
  # analogous covariance calculations
  est.zetaDE <- colSums( (est.DEYj-est.DEY)*(DEDj-DED)*(A_mat) )/(colSums(A_mat)-1) # est.zetaDE0 and est.zetaDE1
  est.zetab0 <- colSums( (est.Yj0-est.Y0)*(Dj0-D0)*(A_mat) )/(colSums(A_mat)-1) # est.zetab00 and est.zetab01
  est.zetab1 <- colSums( (est.Yj1-est.Y1)*(Dj1-D1)*(A_mat) )/(colSums(A_mat)-1) # est.zetab10 and est.zetab11
  
  est.zetaj0 <- matrix(0, nrow=ncol(A_mat), ncol=J) # est.zetaj01 and est.zetaj01
  est.zetaj1 <- matrix(0, nrow=ncol(A_mat), ncol=J) # est.zetaj10 and est.zetaj11
  
  for(j in 1:J){
    tmp1 <- Dj0[j, ]
    tmp2 <- est.Yj0[j, ]
    tmp3 <- rbind(tmp1, tmp1[rep(1,length(D[[j]])-1), ])
    tmp4 <- t( matrix( tmp2 , length(tmp2) , length(Y[[j]]) ) )
    est.zetaj0[, j] <- colSums((Y[[j]]-tmp4)*(D[[j]]-tmp3)*(1-Z[[j]]))/(sum(1-Z[[j]])-1)
    
    tmp5 <- Dj1[j, ]
    tmp6 <- est.Yj1[j, ]
    tmp7 <- rbind(tmp5, tmp5[rep(1,length(D[[j]])-1), ])
    tmp8 <- t( matrix( tmp6 , length(tmp6) , length(Y[[j]]) ) )
    est.zetaj1[, j] <- colSums((Y[[j]]-tmp8)*(D[[j]]-tmp7)*(Z[[j]]))/(sum(Z[[j]])-1)
    
  }
  
  est.zeta <- est.zetaDE*(1/colSums(A_mat)-1/J) + rowSums( (t(apply( est.zetaj0, 1, "/", n-sapply(Z, sum))) + t(apply( est.zetaj1, 1, "/",  sapply(Z, sum))))*t(A_mat) )/J/colSums(A_mat) 
  est.zetab <- c(sum(est.zetab0/colSums(A_mat)) ,sum( est.zetab1/colSums(A_mat))) # est.zetab0 and est.zetab1
  
  
  
  #### CADE and CASE
  est.CADE <- (est.DEY/DED)[1,] # est.CADE0 est.CADE1
  est.CASE <- c(SEY0/SED0, SEY1/SED1) # est. CASE0 and est.CASE1 , may look different when there are >= 3 assignment mechanisms
  est.varCADE <- (var.DEY-2*est.CADE*est.zeta + est.CADE^2*var.DED)/DED[1,]^2 # est.varCADE0 and est.varCADE1
  est.varCASE <- (var.SEY-2*est.CASE*est.zetab + est.CASE^2*var.SED)/est.SED^2 # est.varCASE0 and est.varCASE1
  
  
  
  # standard deviations
  est.stdCADE <- sqrt(est.varCADE)
  est.stdCASE <- sqrt(est.varCASE)
  
  std.DEY <- sqrt(var.DEY)
  std.SEY <- sqrt(var.SEY)
  
  std.DED <- sqrt(var.DED)
  std.SED <- sqrt(var.SED)
  
  #### original code #####
  
  # right confidence intervals 
  level <- qnorm((1-ci)/2, 0, 1)
  rci.CADE <- round(est.CADE-level*est.stdCADE, 3)
  rci.CASE <- round(est.CASE-level*est.stdCASE, 3)
  
  rci.DEY <- round(est.DEY[1,]-level*std.DEY, 3)
  rci.SEY <- round(est.SEY-level*std.SEY, 3)
  
  rci.DED <- round(DED[1,]-level*std.DED, 3)
  rci.SED <- round(est.SED-level*std.SED, 3)
  
  # left confidence intervals 
  lci.CADE <- round(est.CADE+level*est.stdCADE, 3)
  lci.CASE <- round(est.CASE+level*est.stdCASE, 3)
  
  lci.DEY <- round(est.DEY[1, ]+level*std.DEY, 3)
  lci.SEY <- round(est.SEY+level*std.SEY, 3)
  
  lci.DED <- round(DED[1,]+level*std.DED, 3)
  lci.SED <- round(est.SED+level*std.SED, 3)
  
# 
#   class(output) = "random"
  
  return(est.CADE)
}

