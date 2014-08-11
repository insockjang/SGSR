randomPathway_CCLE<-function(pathwayName,dataCombine){
  ### Cross training and testing
  library(predictiveModeling)
  library(synapseClient)
  synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
  
  source("~/PredictiveModel_pipeline/myData_CCLE_new.R")
  dataSets<-myData_CCLE_new(dataCombine,"ActArea")
  
  load("~/SGSR/graphite_pathways_structure.Rdata")
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
  for(kk in 1:24){  
    filename = paste("~/Result_priorIncorporateLasso_filterVar02/",dataCombine,"/CCLE/",pathwayName,"/PriorRandomPathwayStructure_cvDrug_",kk,".Rdata",sep = "")
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
      
      load(paste("~/Result_priorIncorporateLasso_filterVar02/",dataCombine,"/CCLE/",pathwayName,"/PriorIncorporated_cvDrug_",kk,".Rdata",sep = ""))
      aa<-as.numeric(unique(STEP$path[,1]))
      if(length(aa)==0){
        testfunction<-function(k){
          return(resultsScale)
        }
        resultsScale<-mclapply(1:100,function(x)testfunction(x),mc.cores=32)
      }else{      
        
        testfunction<-function(kk.1){
          set.seed(kk.1)
          random.pathway<-modPathways
          rownames(random.pathway)<-sample(rownames(modPathways))
          if(length(aa)>1){
            aaa<-apply(random.pathway[,aa],1,sum)
          }
          if(length(aa)==1){
            aaa<-random.pathway[,aa]
          }
          if(is.element(dataCombine,"E")){
            N<-paste(names(which(aaa>0)),"_expr",sep="")
          }
          if(is.element(dataCombine,"C")){
            N<-paste(names(which(aaa>0)),"_copy",sep="")
          }
          if(is.element(dataCombine,"Mh") | is.element(dataCombine,"Mo")){
            N<-paste(names(which(aaa>0)),"_mut",sep="")
          }        
          penalty_vector<-rep(1,ncol(filteredFeatureDataScaled))    
          a4<-match(N,colnames(filteredFeatureDataScaled))    
          penalty_vector[a4[which(!is.na(a4))]]<-0
          
          resultsScale<-crossValidatePredictiveModel1(filteredFeatureDataScaled, filteredResponseDataScaled, model = myEnetModel1$new(), alpha=1, numFolds=5, nfolds = 5,penalty.factor = penalty_vector)      
          return(resultsScale)
        }
        
        require(multicore)
        resultsScale<-mclapply(1:100,function(x)testfunction(x),mc.cores=32)
      }
      save(resultsScale,file = filename)
    }
  }
}