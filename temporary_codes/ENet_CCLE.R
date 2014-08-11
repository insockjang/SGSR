ENet_CCLE<-function(dataCombine){
  ### DEMO Stepwise grouping Lasso
  require(predictiveModeling)
  require(synapseClient)
  # synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
  source("~/PredictiveModel_pipeline/R5/crossValidatePredictiveModel1.R")
  source("~/PredictiveModel_pipeline/R5/myEnetModel1.R")
  source("~/PredictiveModel_pipeline/myData_CCLE_new.R")
  source("~/SGSR/parallel_stepwiseDecision.R")
  
  
  ###################################################
  #### Load CCLE Molecular Feature Data from Synapse ####
  ###################################################
  dataSets<-myData_CCLE_new(dataCombine,"ActArea")
  
  
  alphas =unique(createENetTuneGrid()[,1])  
  processFunction<-function(kk){  
    filename = paste("~/Result_priorIncorporateLasso/",dataCombine,"/CCLE/ENet/cvDrug_",kk,".Rdata",sep = "")
    if(!file.exists(filename)){
      #########################################################################################################
      ######## Training and Testing data are scaled(normalized) vs. raw(unnormalized) #######################
      #########################################################################################################
      
      # data preprocessing for preselecting features
      filteredData<-filterPredictiveModelData(dataSets$featureData,dataSets$responseData[,kk,drop=FALSE])
      
      # filtered feature and response data
      filteredFeatureData <- filteredData$featureData
      filteredResponseData <- filteredData$responseData
      
      ## scale these data
      filteredFeatureDataScaled <- scale(filteredFeatureData)
      filteredResponseDataScaled <- scale(filteredResponseData)
      
      resultsScale<-crossValidatePredictiveModel1(filteredFeatureDataScaled, filteredResponseDataScaled, model = myEnetModel1$new(), alpha=alphas, numFolds=5, nfolds = 5)
      save(resultsScale,file = filename)
    }
  }
  require(multicore)
  ppp<-mclapply(1:24,function(x)processFunction(x),mc.cores=10)
}
