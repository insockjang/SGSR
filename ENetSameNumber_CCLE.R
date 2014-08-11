ENetSameNumber_CCLE<-function(pathwayName,dataCombine,ALPHA = 0.5){
  ### Cross training and testing
  library(predictiveModeling)
  library(synapseClient)
  synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
  load("~/SGSR_01/graphite_pathways_structure.Rdata")
  
  source("~/PredictiveModel_pipeline/myData_CCLE_new.R")
  
  dataSets<-myData_CCLE_new(dataCombine,"ActArea")
    
  source("~/PredictiveModel_pipeline/R5/myEnetModel1.R")
  source("~/PredictiveModel_pipeline/R5/crossValidatePredictiveModel1.R")
  
  # testfunction<-function(kk){
  allProcess<-function(kk){  
    filename1 = paste("~/Result_priorIncorporateLassoNew_filterVar02/",dataCombine,"/CCLE/",pathwayName,"/PriorENet_alpha_",ALPHA,"_cvDrug_",kk,".Rdata",sep = "")
    if(!file.exists(filename1)){
      #########################################################################################################
      ######## Training and Testing data are scaled(normalized) vs. raw(unnormalized) #######################
      #########################################################################################################
      
      # data preprocessing for preselecting features      
      filteredData<-filterPredictiveModelData(dataSets$featureData,dataSets$responseData[,kk,drop=FALSE],featureVarianceThreshold = 0.2)
      
      # filtered feature and response data
      filteredFeatureData <- filteredData$featureData
      filteredResponseData <- filteredData$responseData
      
      ## scale these data
      filteredFeatureDataScaled <- scale(filteredFeatureData)
      filteredResponseDataScaled <- scale(filteredResponseData)
      
      load(paste("~/Result_priorIncorporateLassoNew_filterVar02/",dataCombine,"/CCLE/",pathwayName,"/PriorIncorporated_cvDrug_",kk,".Rdata",sep = ""))
      
      set.seed(2)
      foldIndices <- createFolds(filteredFeatureDataScaled[,1], k = 5, list = TRUE)
      
      processer<-function(k){                
        fit.sgsr<-cv.glmnet(filteredFeatureDataScaled[-foldIndices[[k]],], filteredResponseDataScaled[-foldIndices[[k]]], alpha = 1,nfolds = 5, penalty.factor = resultSTEP[[k]]$penalty)
        num.of.features<-fit.sgsr$nzero[which.min(fit.sgsr$cvm)]    
        fit.enet.1<-cv.glmnet(filteredFeatureDataScaled[-foldIndices[[k]],], filteredResponseDataScaled[-foldIndices[[k]]], alpha = ALPHA,nfolds = 5)    
        lambda.1<-fit.enet.1$lambda[which.min(abs(fit.enet.1$nzero - num.of.features))]    
        P.resp.enet.1<-predict(fit.enet.1,newx = filteredFeatureDataScaled[foldIndices[[k]],], s= lambda.1)
        return(list(testPredictions = P.resp.enet.1,testObservations = filteredResponseDataScaled[foldIndices[[k]]]))
      }  
      require(multicore)
      resultsScale<-mclapply(1:5,function(x)processer(x),mc.cores=5)
      
      save(resultsScale,file = filename1)
    }
  }
  require(multicore)
  resultsScale<-mclapply(1:24,function(x)allProcess(x),mc.cores=2)
  
}