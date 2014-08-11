PriorRandomLasso_CCLE<-function(pathwayName,dataCombine,KK = c(1:24),mcCoreNum = 32){
  ### DEMO Stepwise grouping Lasso
  require(predictiveModeling)
  require(synapseClient)
  # synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
  source("~/PredictiveModel_pipeline/R5/crossValidatePredictiveModel1.R")
  source("~/PredictiveModel_pipeline/R5/myEnetModel1.R")
  source("~/PredictiveModel_pipeline/myData_CCLE_new.R")
#   source("~/SGSR/parallel_stepwiseDecision.R")
  
  ###################################################
  #### Load Pathways                             ####
  ###################################################
  GRAPHITE<-synGet("syn2135029")
  
  load(GRAPHITE@filePath)
  pathwayName<-toupper(pathwayName)
  if(is.element(pathwayName,"BIOCARTA")){
    allPathways <- BIOCARTA
  }
  if(is.element(pathwayName,"KEGG")){
    allPathways <- KEGG
  }
  if(is.element(pathwayName,"REACTOME")){
    allPathways <- REACTOME
  }
  if(is.element(pathwayName,"NCI")){
    allPathways <- NCI
  }    
  
  ###################################################
  #### Load CCLE Molecular Feature Data from Synapse ####
  ###################################################
  dataSets<-myData_CCLE_new(dataCombine,"ActArea")
  
  for(kk in KK){
    filename = paste("~/Result_priorIncorporateLassoNew_filterVar02/",dataCombine,"/CCLE/",pathwayName,"/restoredRandom_cvDrug_",kk,".Rdata",sep = "")
    if(!file.exists(filename)){
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
      
      resultsScale <- foreach(k = 1:length(foldIndices)) %dopar% {
        if(length(which(resultSTEP[[k]]$penalty==0))==0){
          load(paste("~/Result_priorIncorporateLassoNew_filterVar02/",dataCombine,"/CCLE/",pathwayName,"/restoredPriorIncorporated_cvDrug_",kk,".Rdata",sep = ""))
          results<-list()
          for(pp in 1:50){
            results[[pp]]<-resultsScale[[k]]
          }
          return(results)  
        }else{
          myRand<-function(pp){
            set.seed(pp)
            foldModel <- myEnetModel1$new()      
            foldModel$customTrain(filteredFeatureDataScaled[-foldIndices[[k]],], filteredResponseDataScaled[-foldIndices[[k]]], alpha = 1, nfolds = 5,penalty.factor = sample(resultSTEP[[k]]$penalty))
            
            res <- list(trainPredictions = foldModel$customPredict(filteredFeatureDataScaled[-foldIndices[[k]],]), 
                        trainObservations = filteredResponseDataScaled[-foldIndices[[k]]],
                        testPredictions = foldModel$customPredict(filteredFeatureDataScaled[foldIndices[[k]],]),
                        testObservations = filteredResponseDataScaled[foldIndices[[k]]])
            
            return(res)           
          }
          
          results<-mclapply(1:50,function(x)myRand(x),mc.cores= mcCoreNum)
          return(results)
        }
      }   
      save(resultsScale,file = filename)
    }
  }
}