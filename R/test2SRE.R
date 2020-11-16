

### Testing the hypotheses of DE=0,MDE=0, SE=0 
Test2SRE <- function(Z,A,Y,effect = "DE", alpha = 0.05){
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
    rej <-  (t(C3%*%est$Y.hat)%*%solve(C3%*%est$cov.hat%*%t(C3))%*%(C3%*%est$Y.hat)) >  qchisq(10alpha,2*m-2)
  }
  
  return(drop(rej))
}