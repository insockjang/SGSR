require(doParallel)

parallel_stepwiseDecision<-function(featureData,responseData,groups,coreNum = 1, iterations = 10){
  penalty<-rep(1,ncol(featureData))
  path<-c()                               
  mse<-list()
  k1 <- 0
  while(k1 <= iterations){
    k1<-k1+1
    MSE<-c()
    set.seed(2)
    fit<-cv.glmnet(featureData,responseData,alpha=1,nfolds=5,penalty.factor = penalty)
    a<-min(fit$cvm)
    
    
    grouping<-function(kkk){
      group <- groups[[kkk]]
      b<-match(group,colnames(featureData))
      if(length(which(is.na(b)==0))>0){
        penalty2<-penalty
        penalty2[b[which(is.na(b)==0)]]<-0
        set.seed(2)
        fit<-cv.glmnet(featureData,responseData,alpha=1,nfolds=5,penalty.factor = penalty2)
        return(min(fit$cvm))
      }else{
        return(a)
      }
    }
    
    
    M<-mclapply(1:length(groups),function(x)grouping(x),mc.cores= coreNum)
           
    MSE<-do.call("c",M)                                 
    
    if(min(MSE)<a){
      group1<-groups[[which.min(MSE)]]
      path<-rbind(path,c(which.min(MSE),min(MSE)))
      b<-match(group1,colnames(featureData))
      penalty[b[which(is.na(b)==0)]]<-0
      mse[[k1]]<-MSE
    }else{
      break
    }
    a<-min(MSE)
  }
  return(list(penalty = penalty,path = path,MSE = mse))
}