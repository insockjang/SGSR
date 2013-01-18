### DEMO Stepwise grouping Lasso
library(predictiveModeling)
library(synapseClient)
library(graphite)
#synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
library(multicore)
library(doMC)
registerDoMC()

source("~/DrugResponse/R5/crossValidatePredictiveModel1.R")
source("~/DrugResponse/R5/myEnetModel.R")
###################################################
#### Load CCLE Molecular Feature Data from Synapse ####
###################################################
id_copyLayer <- "269019"     
layer_copy <- loadEntity(id_copyLayer)
eSet_copy <- layer_copy$objects$eSet_copy

id_oncomapLayer <- "1528027"  
layer_oncomap <- loadEntity(id_oncomapLayer)
eSet_oncomap <- layer_oncomap$objects$eSet_hybrid

id_exprLayer <- "269056" 
layer_expr <- loadEntity(id_exprLayer)
eSet_expr <- layer_expr$objects$eSet_expr

id_drugLayer <- "269024" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$adf_drug


#featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr, copy = eSet_copy, mut = eSet_oncomap))
featureData <- exprs(eSet_expr)

# NA filter for training set
featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "rows")
dataSets_ccle <- createFeatureAndResponseDataList(t(featureData_filtered),adf_drug)


load("~/DrugResponse/pathway_analysis/graphite_pathways.Rdata")
groups=list()
# for(k in 1:length(KEGG)){
#   groups[[k]]=nodes(KEGG[[k]])
# }
for(k in 1:length(BIOCARTA)){
  groups[[k]]=nodes(BIOCARTA[[k]])
}
# for(k in 1:length(REACTOME)){
#   groups[[k]]=nodes(REACTOME[[k]])
# }

for(kk in 1:ncol(dataSets_ccle$responseData)){  
  
  
  #########################################################################################################
  ########  Training and Testing data are scaled(normalized) vs. raw(unnormalized)  #######################
  #########################################################################################################
  
  # data preprocessing for preselecting features
  filteredData<-filterPredictiveModelData(dataSets_ccle$featureData,dataSets_ccle$responseData[,kk,drop=FALSE], featureVarianceThreshold = 0.01, corPValThresh = 0.1)
  
  # filtered feature and response data
  filteredFeatureData  <- filteredData$featureData
  filteredResponseData <- filteredData$responseData
  
  ## scale these data    
  filteredFeatureDataScaled <- scale(filteredFeatureData)
  filteredResponseDataScaled <- scale(filteredResponseData)  
  
  
  # 5 fold cross validation 
  
  alphas  = unique(createENetTuneGrid()[,1])
  lambdas = createENetTuneGrid(alphas = 1)[,2]
  
  penalty<-rep(1,ncol(filteredFeatureDataScaled))
  mse<-list()
  k1 <- 0
  while(k1 <= 10){
    k1<-k1+1
    MSE<-c()
    set.seed(2)
    fit<-cv.glmnet(filteredFeatureDataScaled,filteredResponseDataScaled,alpha=1,lambda = lambdas,nfolds=5,penalty.factor = penalty)
    a<-min(fit$cvm)
    M<-foreach(kkk = 1:length(groups)) %dopar% {
      group <- groups[[kkk]]
      b<-match(group,colnames(filteredFeatureDataScaled))
      if(length(which(is.na(b)==0))>0){
        penalty2<-penalty
        penalty2[b[which(is.na(b)==0)]]<-0
        set.seed(2)
        fit<-cv.glmnet(filteredFeatureDataScaled,filteredResponseDataScaled,alpha=1,lambda = lambdas,nfolds=5,penalty.factor = penalty2)
        return(min(fit$cvm))
      }else{
        return(a)
      }     
    }
    MSE<-do.call("c",M)
    
    print(k1)
    
    if(min(MSE)<a){
      group1<-groups[[which.min(MSE)]]
      b<-match(group1,colnames(filteredFeatureDataScaled))
      penalty[b[which(is.na(b)==0)]]<-0
      mse[[k1]]<-MSE
    }else{
      break
    }
  }
  
}
