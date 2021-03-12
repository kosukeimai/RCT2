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
#' \item{CADE}{ The point estimates of the CADE for each assignment mechanism.  }
#'\item{CASE}{ The point estimate of CASE for each assignment mechanism.  } 
#'\item{var.CADE1}{ The  variance estimate of CADE for each assignment mechanism.   } 
#'\item{var.CASE1}{ The  variance estimate of CASE for each assignment mechanism.   } 
#'\item{DEY1}{ The point estimate of DEY for each assignment mechanism.  }
#'\item{DED1}{ The point estimate of DED for each assignment mechanism.  }
#'\item{var.DEY1}{ The  variance estimate of DEY for each assignment mechanism.   } 
#'\item{var.DED1}{ The  variance estimate of DED for each assignment mechanism.   } 
#'\item{SEY1}{ The point estimate of SEY for each pairwise groups of assignment mechanisms.  } 
#'\item{SED1}{ The point estimate of SED for each pairwise groups of assignment mechanisms.  }
#'\item{var.SEY1}{ The  variance estimate of SEY for each pairwise groups of assignment mechanisms.   } 
#'\item{var.SED1}{ The  variance estimate of SED for each pairwise groups of assignment mechanisms.   } 
#'\item{lci.CADE}{ The left endpoint for the confidence intervals for the CADE from each assignment mechanism. }
#'\item{rci.CADE}{ The right endpoint for the confidence intervals for the CADE from each assignment mechanism. }
#'\item{lci.CASE}{ The left endpoint for the confidence intervals for the CASE from each assignment mechanism. }
#'\item{rci.CASE}{ The left endpoint for the confidence intervals for the CASE from each assignment mechanism. }
#'\item{lci.DEY}{ The left endpoint for the confidence intervals for the DEY from each assignment mechanism. }
#'\item{rci.DEY}{ The left endpoint for the confidence intervals for the DEY from each assignment mechanism. }
#'\item{lci.SEY}{ The left endpoint for the confidence intervals for the SEY from each pairwise groups of assignment mechanisms. }
#'\item{rci.SEY}{ The left endpoint for the confidence intervals for the SEY from each pairwise groups of assignment mechanism. }
#'\item{lci.DED}{ The left endpoint for the confidence intervals for the DED from each assignment mechanism. }
#'\item{rci.DED}{ The left endpoint for the confidence intervals for the DED from each assignment mechanism. }
#'\item{lci.SED}{ The left endpoint for the confidence intervals for the SED from each pairwise groups of assignment mechanism. }
#'\item{rci.SED}{ The left endpoint for the confidence intervals for the SED from each pairwise groups of assignment mechanism. }
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
#' @importFrom utils head
#' @export CADErand


CADErand<-function(data,individual=1, ci = 0.95){
  id <- NULL
  data$id <- factor(data$id)
  data <- data[order(data$id),]
  
  N <- length(data$id)
  
  cluster.id <- unique(data$id)	
  n.cluster <- length(cluster.id)	
  J <- n.cluster
  
  
  individual <- 1
  
  # data frame of id and frequency of the id
  freq <- data.frame( table(data$id) )
  n <- freq$Freq
  
  merged <- merge(data, freq, by.x = "id", by.y = "Var1")
  if(individual == 1){
    merged$Y <- merged$Y * merged$Freq * n.cluster / N
    merged$D <- merged$D * merged$Freq * n.cluster / N
  }
  
  n.mech <- unique(merged$A)
  n.amech <- length(n.mech)
  
  # full dataframe A
  testA <- data.frame(factor(merged$A) )
  A_mat <- model.matrix(~ . + 0, data=testA, contrasts.arg = lapply(testA, contrasts, contrasts=FALSE))
  merged <- cbind(merged, A_mat)
  
  # A by cluster
  A_index <- match(unique(merged$id), merged$id)
  A_cluster <- factor(merged$A[A_index])
  A_cluster_mat <- model.matrix(~A_cluster-1)
  A_sums <- colSums(A_cluster_mat)
  A3 <- tapply(merged$A, merged$id, mean )

  A4 <- matrix(0, nrow = n.cluster, ncol = n.amech)
  for(i in 1:n.cluster){
    mech <- A3[i]
    A4[i, mech+1] <- 1
  }
  
  
  
  # common variable
  Z_sum <- tapply(merged$Z, merged$id, sum)
  
  
  
  Dj0_init <- A_mat*(merged$D - merged$D*merged$Z)
  Dj1_init <- A_mat*(merged$D*merged$Z)
  Dj0 <- matrix(data = NA, nrow = n.cluster, ncol = n.amech)
  Dj1 <- matrix(data = NA, nrow = n.cluster, ncol = n.amech)
  for(i in 1:length(n.mech)){
    Dj0[, i] <- tapply(Dj0_init[, i], merged$id, sum)/(n-Z_sum)
    Dj1[, i] <- tapply(Dj1_init[, i], merged$id, sum)/Z_sum
  }
  
  
  # DED calculations
  D0 <- t(as.data.frame(colSums(Dj0)/A_sums)) # this is est.D00 and est.D01
  D1 <- t(as.data.frame(colSums(Dj1)/A_sums))
  
  
  
  DEDj <- Dj1-Dj0 # this is est.DEDj0 and est.DEDj1
  DED <- D1-D0
  
  
  # spillover calculations
  lenSED <- ncol(A_mat)*(ncol(A_mat)-1)/2
  SED0 <- rep(0, lenSED)
  SED1 <- rep(0, lenSED)
  for(i in 2:lenSED){
    SED0[i-1] <- D0[1, i]-D0[1, i-1]
    SED1[i-1] <- D1[1, i]-D1[1, i-1]
  }
  
  est.SED <- c(SED0, SED1)
  
  
  # analogous variance calculations
  est.xiDE <- colSums( (sweep(DEDj, 2, DED))^2*A4)/(A_sums-1) # est.xiDE0 and est.xiDE1
  est.xib0 <- colSums( (sweep(Dj0, 2, D0))^2*A4)/(A_sums-1) # est.xib00 and est.xib01
  est.xib1 <- colSums( (sweep(Dj1, 2, D1))^2*A4)/(A_sums-1) # est.xib10 and est.xib11
  
  est.xij0 <- matrix(0, nrow=n.amech, ncol=J) # est.xij00 and est.xij01
  est.xij1 <- matrix(0, nrow=n.amech, ncol=J) # est.xj10 and est.xij11
  
  for (j in 1:J){
    tmp1 <- Dj0[j, ]
    tmp2 <- t( matrix( tmp1 , length(tmp1) , n[j] ) )
    tmp_subset <- subset(merged, id %in% freq$Var1[j])
    D_temp <- tmp_subset$D
    Z_temp <- tmp_subset$Z
    est.xij0[, j] <- colSums((D_temp - tmp2)^2*(1-Z_temp))/(sum(1-Z_temp)-1) # est.xij00 and est.xij01
    
    tmp3 <- Dj1[j, ]
    tmp4 <- t( matrix( tmp3, length(tmp3), n[j] ) )
    est.xij1[, j] <- colSums((D_temp - tmp4)^2*(Z_temp))/(sum(Z_temp) - 1)  
  }
  
  denom1 <- t(matrix( Z_sum , length(Z_sum) , nrow(est.xij0) ))
  denom2 <- t(matrix( freq$Freq - Z_sum , length(Z_sum) , nrow(est.xij0) ))
  
  var.DED <- est.xiDE*(1/A_sums-1/J)+rowSums((est.xij0/denom2+est.xij1/denom1)*t(A4))/J/A_sums # var.DED0 and var.DED1
  var.SED1 <- sum(est.xib1/A_sums)
  var.SED0 <- sum(est.xib0/A_sums)
  var.SED <- c(var.SED0, var.SED1)
  
  
  est.Yj0 <- as.vector(tapply( (merged$Y - merged$Y*merged$Z), merged$id, sum)/(freq$Freq-Z_sum) )*A4
  est.Yj1 <- as.vector( tapply( (merged$Y*merged$Z ), merged$id, sum) / Z_sum )*A4
  
  est.Y0 <- colSums(est.Yj0*A4)/A_sums # est.Y00 and est.Y01
  est.Y1 <- colSums(est.Yj1*A4)/A_sums # est.Y10 and est. Y11
  est.DEYj <- est.Yj1-est.Yj0 # est.DEYj0 and est.DEYj1
  est.DEY <- est.Y1-est.Y0 # est.DEY0 and est.DEY1
  lenSEY <- ncol(A4)*(ncol(A4)-1)/2
  SEY0 <- rep(0, lenSEY) # this is est.SEY0
  SEY1 <- rep(0, lenSEY) # this is est.SEY1
  for(i in 2:lenSEY){
    SEY0[i-1] <- est.Y0[i]-est.Y0[i-1]
    SEY1[i-1] <- est.Y1[i]-est.Y1[i-1]
  }
  est.SEY <- c(SEY0, SEY1)
  
  # variance for spillover effects
  est.sigmaDE <- colSums((sweep(est.DEYj, 2, est.DEY ))^2*A4)/(A_sums-1) # est.sigmaDEO and est.sigmaDE1
  est.sigmab0 <- colSums((sweep( est.Yj0, 2, est.Y0 )) ^2*A4)/(A_sums-1) # est.sigmab00 and sigmab01
  est.sigmab1 <- colSums((sweep( est.Yj1, 2, est.Y1))^2*A4)/(A_sums-1) # est.sigmab10 and est.sigmab11
  est.sigmaj0 <- matrix(0, 2, J) # est.sigmaj00 and est.sigmaj01
  est.sigmaj1 <- matrix(0, 2, J) # est.sigmaj10 and est.sigma11
  
  for(j in 1:J){
    tmp1 <- est.Yj0[j, ]
    tmp2 <- t( matrix( tmp1 , length(tmp1) , n[[j]] ) )
    tmp_subset <- subset(merged, id %in% freq$Var1[j])
    Y_temp <- tmp_subset$Y
    Z_temp <- tmp_subset$Z
    est.sigmaj0[, j] <- colSums((Y_temp-tmp2)^2*(1-Z_temp))/(sum(1-Z_temp)-1)
    
    tmp3 <- est.Yj1[j, ]
    tmp4 <- t( matrix( tmp3 , length(tmp3) , n[[j]]) )
    est.sigmaj1[, j] <- colSums((Y_temp-tmp4)^2*(Z_temp))/(sum(Z_temp)-1)
  }
  
  # var.DEY0 and var.DEY1
  var.DEY <- est.sigmaDE*(1/A_sums-1/J)+rowSums((est.sigmaj0/denom2+est.sigmaj1/denom1)*t(A4))/J/A_sums # var.DED0 and var.DED1
  
  # var.SEY0 and var.SEY1
  var.SEY <- c(sum(est.sigmab0/A_sums), sum(est.sigmab1/A_sums))
  
  # analogous covariance calculations
  est.zetaDE <- colSums( (est.DEYj-est.DEY)*( sweep(DEDj, 2, DED) )*A4 )/(A_sums-1) # est.zetaDE0 and est.zetaDE1
  est.zetab0 <- colSums( (est.Yj0-est.Y0)*( sweep(Dj0, 2, D0))*A4 )/(A_sums-1) # est.zetab00 and est.zetab01
  est.zetab1 <- colSums( (est.Yj1-est.Y1)*(sweep(Dj1, 2, D1 ))*A4 )/(A_sums-1) # est.zetab10 and est.zetab11
  
  est.zetaj0 <- matrix(0, nrow=ncol(A_mat), ncol=J) # est.zetaj01 and est.zetaj01
  est.zetaj1 <- matrix(0, nrow=ncol(A_mat), ncol=J) # est.zetaj10 and est.zetaj11
  
  for(j in 1:J){
    tmp_subset <- subset(merged, id %in% freq$Var1[j])
    Y_temp <- tmp_subset$Y
    Z_temp <- tmp_subset$Z
    D_temp <- tmp_subset$D
    tmp1 <- Dj0[j, ]
    tmp2 <- est.Yj0[j, ]
    tmp3 <- t( matrix(tmp1, length(tmp1), n[j] ) )
    tmp4 <- t( matrix( tmp2 , length(tmp2) , n[j] ) )
    est.zetaj0[, j] <- colSums((Y_temp-tmp4)*(D_temp-tmp3)*(1-Z_temp))/(sum(1-Z_temp)-1)
    
    tmp5 <- Dj1[j, ]
    tmp6 <- est.Yj1[j, ]
    tmp7 <- t( matrix(tmp5, length(tmp5), n[j]) )
    tmp8 <- t( matrix( tmp6 , length(tmp6) , n[j] ) )
    est.zetaj1[, j] <- colSums((Y_temp-tmp8)*(D_temp-tmp7)*(Z_temp))/(sum(Z_temp)-1)
    
  }
  
  est.zeta <- est.zetaDE*(1/A_sums-1/J) + rowSums((est.zetaj0/denom2+est.zetaj1/denom1)*t(A4))/J/A_sums
  est.zetab <- c(sum(est.zetab0/A_sums) ,sum( est.zetab1/A_sums)) # est.zetab0 and est.zetab1
  
  #### CADE and CASE
  est.CADE <- (est.DEY/DED) # est.CADE0 est.CADE1
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
  ci <- 0.05
  
  
  level <- qnorm((1-ci)/2, 0, 1)
  rci.CADE <- round(est.CADE-level*est.stdCADE, 3)
  rci.CASE <- round(est.CASE-level*est.stdCASE, 3)
  
  rci.DEY <- round(est.DEY-level*std.DEY, 3)
  rci.SEY <- round(est.SEY-level*std.SEY, 3)
  
  rci.DED <- round(DED-level*std.DED, 3)
  rci.SED <- round(est.SED-level*std.SED, 3)
  
  # left confidence intervals 
  lci.CADE <- round(est.CADE+level*est.stdCADE, 3)
  lci.CASE <- round(est.CASE+level*est.stdCASE, 3)
  
  lci.DEY <- round(est.DEY+level*std.DEY, 3)
  lci.SEY <- round(est.SEY+level*std.SEY, 3)
  
  lci.DED <- round(DED+level*std.DED, 3)
  lci.SED <- round(est.SED+level*std.SED, 3)
  
  rownames(est.CADE) <- NULL
  rownames(est.varCADE) <- NULL
  rownames(est.DEY) <- NULL
  rownames(DED) <- NULL
  rownames(var.DED) <- NULL
  rownames(var.DEY) <- NULL


  
  output <- list(CADE = est.CADE, CASE = est.CASE, var.CADE = est.varCADE, var.CASE = est.varCASE,
                 DEY = est.DEY, DED = DED, var.DED = var.DED,
                 var.DEY = var.DEY, var.DED = var.DED,
                 SEY = est.SEY, SED = est.SED,
                 var.SEY = var.SEY, var.SED = var.SED,
                 lci.CADE = lci.CADE, lci.CASE = lci.CASE,
                 lci.DEY = lci.DEY, lci.SEY = lci.SEY,
                 lci.DED = lci.DED, lci.SED = lci.SED,
                 rci.CADE = rci.CADE, rci.CASE = rci.CASE, 
                 rci.DEY = rci.DEY, rci.SEY = rci.SEY,
                 rci.DED = rci.DED, rci.SED = rci.SED)
  
  class(output) = "random"

  
  return(output)
}

