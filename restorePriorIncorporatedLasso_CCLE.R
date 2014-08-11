restorePriorIncorporatedLasso_CCLE<-function(pathwayName,dataCombine){
  ### DEMO Stepwise grouping Lasso
  require(predictiveModeling)
  require(synapseClient)
  # synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
  source("~/PredictiveModel_pipeline/R5/myEnetModel1.R")
  source("~/PredictiveModel_pipeline/myData_CCLE_new.R")
  
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
  MSigDB<-synGet("syn2227979")
  load(MSigDB@filePath)
  if(is.element(pathwayName,"GO_BP")){
    allPathways <- MSigDB$C5.GO_BP
  }    
  if(is.element(pathwayName,"GO_MF")){
    allPathways <- MSigDB$C5.GO_MF
  }    
  
  ###################################################
  #### Load CCLE Molecular Feature Data from Synapse ####
  ###################################################
  dataSets<-myData_CCLE_new(dataCombine,"ActArea")
  
  #   require(graphite)
  groups=list()
  for(k in 1:length(allPathways)){
    a=allPathways[[k]]
    a1<-paste(a,"_expr",sep="")
    a2<-paste(a,"_copy",sep="")
    a3<-paste(a,"_mut",sep="")
    aa<-union(a1,union(a2,a3))
    groups[[k]]<-aa
  }
  
  #   KK<-c(1,3,4,5,9,10,12,13,14,15,16,17,20,21,6,23)
  KK<-c(1:24)
  for(kk in KK){
    filename1 = paste("~/Result_priorIncorporateLassoNew_filterVar02/",dataCombine,"/CCLE/",pathwayName,"/PriorIncorporated_cvDrug_",kk,".Rdata",sep = "")
    load(filename1)  
    filename = paste("~/Result_priorIncorporateLassoNew_filterVar02/",dataCombine,"/CCLE/",pathwayName,"/restoredPriorIncorporated_cvDrug_",kk,".Rdata",sep = "")
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
      
      set.seed(2)
      foldIndices <- createFolds(filteredFeatureDataScaled[,1], k = 5, list = TRUE)
      
      resultsScale <- foreach(k = 1:length(foldIndices)) %dopar% {      
        foldModel <- myEnetModel1$new()      
        foldModel$customTrain(filteredFeatureDataScaled[-foldIndices[[k]],], filteredResponseDataScaled[-foldIndices[[k]]], alpha = 1, nfolds = 5,penalty.factor = resultSTEP[[k]]$penalty)
        
        res <- list(trainPredictions = foldModel$customPredict(filteredFeatureDataScaled[-foldIndices[[k]],]), 
                    trainObservations = filteredResponseDataScaled[-foldIndices[[k]]],
                    testPredictions = foldModel$customPredict(filteredFeatureDataScaled[foldIndices[[k]],]),
                    testObservations = filteredResponseDataScaled[foldIndices[[k]]])
        
        return(res)           
      }    
      
      
      save(resultsScale,file = filename)
      
    }
  }  
}