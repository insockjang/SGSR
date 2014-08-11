randomPathway_Sanger<-function(pathwayName,dataCombine,KK = NA,mcCoreNum = 32){
  ### Cross training and testing
  library(predictiveModeling)
  library(synapseClient)
  synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
  
  source("~/PredictiveModel_pipeline/myData_Sanger.R")
  dataSets<-myData_Sanger(dataCombine,"IC50")
  
  a<-synGet("syn2604222")
  load(a@filePath)
  
  pathwayName<-toupper(pathwayName)
  if(is.element(pathwayName,"BIOCARTA")){
    allPathways <- structure.BIOCARTA
  }
  if(is.element(pathwayName,"KEGG")){
    allPathways <- structure.KEGG
  }  
  if(is.element(pathwayName,"NCI")){
    allPathways <- structure.NCI
  }    
  if(is.element(pathwayName,"GO_BP")){
    allPathways <- structure.GO_BP
  }    
  if(is.element(pathwayName,"GO_MF")){
    allPathways <- structure.GO_MF
  }    
  
  
  
  
  source("~/PredictiveModel_pipeline/R5/myEnetModel1.R")
  source("~/PredictiveModel_pipeline/R5/crossValidatePredictiveModel1.R")
  
  # testfunction<-function(kk){
  for(kk in KK){  
    filename = paste("~/Result_priorIncorporateLassoNew_filterVar02/",dataCombine,"/Sanger/",pathwayName,"/PriorRandomPathwayStructure_cvDrug_",kk,".Rdata",sep = "")
    if(!file.exists(filename)){
      #########################################################################################################
      ######## Training and Testing data are scaled(normalized) vs. raw(unnormalized) #######################
      #########################################################################################################
      
      # data preprocessing for preselecting features
      filteredData<-filterPredictiveModelData(dataSets$featureData,dataSets$responseData[,kk,drop=FALSE],featureVarianceThreshold=0.2)
      
      # filtered feature and response data
      filteredFeatureData <- filteredData$featureData
      filteredResponseData <- filteredData$responseData
      
      ## scale these data
      filteredFeatureDataScaled <- scale(filteredFeatureData)
      filteredResponseDataScaled <- scale(filteredResponseData)
      
      # Unique genes in input data
      NN<-colnames(filteredFeatureDataScaled)
      NN1<-strsplit(NN,"_")
      MM<-c()
      for(pp in 1:length(NN1)){
        MM<-c(MM,NN1[[pp]][1])
      }
      MM1<-unique(MM)
      OO1<-rownames(allPathways)
      qq<-match(OO1,MM1)
      qq1<-which(!is.na(qq))
      modPathway<-allPathways[qq1,]
      
      load(paste("~/Result_priorIncorporateLassoNew_filterVar02/",dataCombine,"/Sanger/",pathwayName,"/PriorIncorporated_cvDrug_",kk,".Rdata",sep = ""))
      
      set.seed(2)
      foldIndices <- createFolds(filteredFeatureDataScaled[,1], k = 5, list = TRUE)
      
      resultsScale <- foreach(k = 1:length(foldIndices)) %dopar% {
        aa<-as.numeric(unique(resultSTEP[[k]]$path[,1]))
        if(length(which(resultSTEP[[k]]$penalty==0))==0){
          load(paste("~/Result_priorIncorporateLassoNew_filterVar02/",dataCombine,"/Sanger/",pathwayName,"/restoredPriorIncorporated_cvDrug_",kk,".Rdata",sep = ""))
          results<-list()
          for(pp in 1:50){
            results[[pp]]<-resultsScale[[k]]
          }
          return(results)  
        }else{
          myRand<-function(pp){
            set.seed(pp)            
            
            random.pathway<-modPathway
            rownames(random.pathway)<-sample(rownames(modPathway))
            if(length(aa)>1){
              aaa<-apply(random.pathway[,aa],1,sum)
            }
            if(length(aa)==1){
              aaa<-random.pathway[,aa]
            }
            N1<-N2<-N3<-c()
            if(length(grep("E",dataCombine))>0){
              N1<-paste(names(which(aaa>0)),"_expr",sep="")
            }
            if(length(grep("C",dataCombine))>0){
              N2<-paste(names(which(aaa>0)),"_copy",sep="")
            }
            if(length(grep("Mh",dataCombine))>0 | length(grep("Mo",dataCombine))>0){
              N3<-paste(names(which(aaa>0)),"_mut",sep="")
            }        
            N<-union(N1,union(N2,N3))
            
            penalty_vector<-rep(1,ncol(filteredFeatureDataScaled))    
            a4<-match(N,colnames(filteredFeatureDataScaled))    
            penalty_vector[a4[which(!is.na(a4))]]<-0
            
            foldModel <- myEnetModel1$new()      
            foldModel$customTrain(filteredFeatureDataScaled[-foldIndices[[k]],], filteredResponseDataScaled[-foldIndices[[k]]], alpha = 1, nfolds = 5,penalty.factor = penalty_vector)
            
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