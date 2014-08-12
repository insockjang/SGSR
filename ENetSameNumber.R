ENetSameNumber_CCLE<-function(pathwayName,dataCombine,ALPHA = 0.5,KK){
  ### Cross training and testing
  library(predictiveModeling)
  library(synapseClient)
  require(devtools)
  
  a<-synGet("syn2604222")
  load(a@filePath)
  
  source_url("https://raw.githubusercontent.com/Sage-Bionetworks/PredictiveModel_pipeline/master/R5/myEnetModel1.R")
  source_url("https://raw.githubusercontent.com/Sage-Bionetworks/PredictiveModel_pipeline/master/R5/crossValidatePredictiveModel1.R")
  source_url("https://raw.githubusercontent.com/Sage-Bionetworks/PredictiveModel_pipeline/master/myData_CCLE_new.R")
    
  dataSets<-myData_CCLE_new(dataCombine,"ActArea")
    
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
  resultsScale<-mclapply(KK,function(x)allProcess(x),mc.cores=2)
  
}


ENetSameNumber_Sanger<-function(pathwayName,dataCombine,ALPHA = 0.5,KK,mcCoreNum = 32){
  ### Cross training and testing
  library(predictiveModeling)
  library(synapseClient)
  require(devtools)
  
  a<-synGet("syn2604222")
  load(a@filePath)
  
  source_url("https://raw.githubusercontent.com/Sage-Bionetworks/PredictiveModel_pipeline/master/R5/myEnetModel1.R")
  source_url("https://raw.githubusercontent.com/Sage-Bionetworks/PredictiveModel_pipeline/master/R5/crossValidatePredictiveModel1.R")
  source_url("https://raw.githubusercontent.com/Sage-Bionetworks/PredictiveModel_pipeline/master/myData_Sanger.R")
    
  dataSets<-myData_Sanger(dataCombine,"IC50")
    
  # testfunction<-function(kk){
  allProcess<-function(kk){  
    filename1 = paste("~/Result_priorIncorporateLassoNew_filterVar02/",dataCombine,"/Sanger/",pathwayName,"/PriorENet_alpha_",ALPHA,"_cvDrug_",kk,".Rdata",sep = "")
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
      
      load(paste("~/Result_priorIncorporateLassoNew_filterVar02/",dataCombine,"/Sanger/",pathwayName,"/PriorIncorporated_cvDrug_",kk,".Rdata",sep = ""))
      
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
  
  #KK=c(1,2,4,6,9,11,12,13,14,17,39,41,75,82,84,87,90,94,95,103,108,110,122,123,124,126,127,129)
  require(multicore)
  resultsScale<-mclapply(KK,function(x)allProcess(x),mc.cores= mcCoreNum)
  
}