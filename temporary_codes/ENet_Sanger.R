ENet_Sanger<-function(dataCombine){
  ### DEMO Stepwise grouping ENet
  require(predictiveModeling)
  require(synapseClient)
  # synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
  source("~/PredictiveModel_pipeline/R5/crossValidatePredictiveModel1.R")
  source("~/PredictiveModel_pipeline/R5/myEnetModel1.R")
  source("~/PredictiveModel_pipeline/myData_Sanger.R")
  source("~/SGSR/parallel_stepwiseDecision.R")
  
  
  ###################################################
  #### Load CCLE Molecular Feature Data from Synapse ####
  ###################################################
  dataSets<-myData_Sanger(dataCombine,"IC50")
  
  alphas =unique(createENetTuneGrid()[,1])  
  processFunction<-function(kk){  
    filename = paste("~/Result_priorIncorporateLasso/",dataCombine,"/Sanger/ENet/cvDrug_",kk,".Rdata",sep = "")
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
  ppp<-mclapply(1:138,function(x)processFunction(x),mc.cores=10)
}

