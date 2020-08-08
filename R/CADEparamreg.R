

CADEparamreg=function(data, assign.prob, ci.level=0.95){
  ## validate ci.level
  if(ci.level>1){stop('Please specify the confidence interval as a decimal between 0 and 1.')}
  
  ## transform the data into list 	
  if(!is.factor(data$id)){stop('The cluster_id should be a factor variable.')}
  cluster.id=unique(data$id)	
  n.cluster=length(cluster.id)	
  Z=vector("list", n.cluster) 	
  D=vector("list", n.cluster) 
  Y=vector("list", n.cluster) 
  A=rep(0,n.cluster)
  for (i in 1:n.cluster){
    Z[[i]]=as.numeric(data$Z[data$id==cluster.id[i]])
    D[[i]]=as.numeric(data$D[data$id==cluster.id[i]])
    Y[[i]]=data.hosp$Y[data$id==cluster.id[i]]
    if (length(unique(data$A[data$id==cluster.id[i]]))!=1){
      stop( paste0('The assignment mechanism in cluster ',i,' should be the same.'))
    }
    A[i]=data$A[data$id==cluster.id[i]][1]
  }
  
  
  data=data[!is.na(data$Y),]
  n=sapply(Z,length)
  id.remove= cluster.id[n-sapply(Z,sum)<=1 | sapply(Z,sum)<=1]
  data=data[!(data$id %in% id.remove),]
  
  
  
  ### model-based analysis
  ## ITT effects
  
  data.ITT=data
  data.ITT$A[data.ITT$A==1]=assign.prob
  data.ITT$A[data.ITT$A==0]=1-assign.prob
  lmITT=lm(Y~Z+A+Z:A, data=data.ITT)
  var.lmITT=vcovHC(lmITT, type = "HC2", cluster = "id", adjust = T)
  matrix1=array(0,dim=c(4,4))
  matrix1[1,]=c(0,1,0,1)
  matrix1[2,2]=1
  matrix1[3,c(3,4)]=1
  matrix1[4,3]=1
  effect=as.numeric(matrix1%*%lmITT$coefficients)
  var.effect=matrix1%*%var.lmITT%*%t(matrix1)
  sd.effect=sqrt(diag(var.effect))
  
  

  
  ## iv estimate
  
  data.iv=data
  data.iv$propD=data.iv$D
  ### create a new variable for the proportion of D
  cluster.id=unique(data.iv$id)	
  n.cluster=length(cluster.id)	
  for (i in 1:n.cluster){
    data.iv$propD[data.iv$id==cluster.id[i]]= mean(data.iv$D[data.iv$id==cluster.id[i]])
  }
  
  lmiv=  ivreg(Y~D+propD+D:propD | Z+A+Z:A,data=data.iv)
  var.lmiv=vcovHC(lmiv, type = "HC2", cluster = "id", adjust = T)
  
  effect.iv=as.numeric(matrix1%*%lmiv$coefficients)
  var.effect.iv=matrix1%*%var.lmiv%*%t(matrix1)
  sd.effect.iv=sqrt(diag(var.effect.iv))  
  
  

  
  
  ## results
  ci.tail=(1-ci.level)/2
  qnorm=qnorm(ci.level+ci.tail)
  
  ITT.effect=round(effect)
  ITT.LCI=round(effect-qnorm*sd.effect)
  ITT.RCI=round(effect+qnorm*sd.effect)
  
  
  IV.effect=round(effect.iv)
  IV.LCI=round(effect.iv-qnorm*sd.effect.iv)
  IV.RCI=ITT.RCI=round(effect+qnorm*sd.effect)
  
  t.stat.ITT=summary(lmITT)$coefficients[,3]
  p.vals.ITT=2*pt(t.stat.ITT, sum(n)-2, lower.tail = F)
  
  t.stat.iv=summary(lmiv)$coefficients[,3]
  p.vals.iv=2*pt(t.stat.iv, sum(n)-2, lower.tail = F)
  
  return(list(ITT.DE=ITT.effect[c(1,2)], ITT.SE=ITT.effect[c(3,4)], ITT.DE.CI=list(c(ITT.LCI[1], ITT.RCI[1]), c(ITT.LCI[2], ITT.RCI[2])),
              ITT.SE.CI=list(c(ITT.LCI[3], ITT.RCI[3]), c(ITT.LCI[4], ITT.RCI[4])),
              IV.DE=IV.effect[c(1, 2)], IV.SE=IV.effect[c(3, 4)], IV.DE.CI=list(c(IV.LCI[1], IV.RCI[1]), c(IV.LCI[2], IV.RCI[2])),
              IV.SE.CI=list(c(IV.LCI[3], IV.RCI[3]), c(IV.LCI[4], IV.RCI[4])),
              ITT.tstat=t.stat.ITT, ITT.pvals=p.vals.ITT, iv.tstat=t.stat.iv, iv.pvals=p.vals.iv))
}









