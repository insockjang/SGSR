PriorIncorporatedLasso_Sanger<-function(pathwayName,dataCombine){
  ### DEMO Stepwise grouping Lasso
  require(predictiveModeling)
  require(synapseClient)
  # synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
  source("~/PredictiveModel_pipeline/R5/crossValidatePredictiveModel1.R")
  source("~/PredictiveModel_pipeline/R5/myEnetModel1.R")
  source("~/PredictiveModel_pipeline/myData_Sanger.R")
  source("~/SGSR/parallel_stepwiseDecision.R")
  
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
  dataSets<-myData_Sanger(dataCombine,"IC50")
  
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
  
  KK<-c(129,127,95,94,108,39,11,17,90,2,122,124,126,1,87,84,6,75)
  for(kk in KK){
    filename = paste("~/Result_priorIncorporateLasso_filterVar02//",dataCombine,"/Sanger/",pathwayName,"/PriorIncorporated_cvDrug_",kk,".Rdata",sep = "")
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
      STEP<-parallel_stepwiseDecision(filteredFeatureDataScaled,filteredResponseDataScaled,groups,coreNum = 32,100)
      
      
      
      set.seed(2)
      resultsScale<-crossValidatePredictiveModel1(filteredFeatureDataScaled, filteredResponseDataScaled, model = myEnetModel1$new(), alpha=1, numFolds=5, nfolds = 5,penalty.factor = STEP$penalty)
      save(resultsScale,STEP,file = filename)
    }
  }
}

