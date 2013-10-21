PriorRandomLasso_CCLE<-function(pathwayName,dataCombine){
  ### DEMO Stepwise grouping Lasso
  require(predictiveModeling)
  require(synapseClient)
  # synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
  source("~/PredictiveModel_pipeline/R5/crossValidatePredictiveModel1.R")
  source("~/PredictiveModel_pipeline/R5/myEnetModel1.R")
  source("~/PredictiveModel_pipeline/myData_CCLE_new.R")
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
  
  ###################################################
  #### Load CCLE Molecular Feature Data from Synapse ####
  ###################################################
  dataSets<-myData_CCLE_new(dataCombine,"ActArea")
  
  
  
  for(kk in 1:24){
    filename = paste("~/Result_priorIncorporateLasso/",dataCombine,"/CCLE/",pathwayName,"/PriorRandom_cvDrug_",kk,".Rdata",sep = "")
    if(!file.exists(filename)){
      #########################################################################################################
      ######## Training and Testing data are scaled(normalized) vs. raw(unnormalized) #######################
      #########################################################################################################
      
      # data preprocessing for preselecting features
      filteredData<-filterPredictiveModelData(dataSets$featureData,dataSets$responseData[,kk,drop=FALSE], featureVarianceThreshold = 0.01, corPValThresh = 0.1)
      
      # filtered feature and response data
      filteredFeatureData <- filteredData$featureData
      filteredResponseData <- filteredData$responseData
      
      ## scale these data
      filteredFeatureDataScaled <- scale(filteredFeatureData)
      filteredResponseDataScaled <- scale(filteredResponseData)
      
      load(paste("~/Result_priorIncorporateLasso/",dataCombine,"/CCLE/",pathwayName,"/PriorIncorporated_cvDrug_",kk,".Rdata",sep = ""))
      aaa<-unique(STEP$path[,1])
      
      bbb<-c()    
      for(k.1 in 1:length(aaa)){
        bbb = union(bbb,allPathways[[aaa[k.1]]])
      }
      
      bbb1<-paste(bbb,"_expr",sep="")
      bbb2<-paste(bbb,"_copy",sep="")
      bbb3<-paste(bbb,"_mut",sep="")
      bbbb<-union(bbb1,union(bbb2,bbb3))
      
      
      b.1<-intersect(bbbb,colnames(filteredFeatureDataScaled))
      
      
      
      testfunction<-function(k.2){  
        set.seed(k.2)
        penalty_vector<-rep(1,ncol(filteredFeatureDataScaled))    
        penalty_vector[sample(ncol(filteredFeatureDataScaled),length(b.1))]<-0
        
        set.seed(2)  
        resultsScale<-crossValidatePredictiveModel1(filteredFeatureDataScaled, filteredResponseDataScaled, model = myEnetModel1$new(), alpha=1, numFolds=5, nfolds = 5,penalty.factor = penalty_vector)      
        return(resultsScale)
      }
      require(multicore)
      resultsScale<-mclapply(1:100,function(x)testfunction(x),mc.cores=5)
      save(resultsScale,file = filename)
      
    }
  }
}