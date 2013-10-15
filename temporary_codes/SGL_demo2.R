### DEMO Stepwise grouping Lasso
library(predictiveModeling)
library(synapseClient)

source("~/PredictiveModel_pipeline/R5/crossValidatePredictiveModel1.R")
source("~/SGSR/stepwiseDecision.R")
source("~/PredictiveModel_pipeline/R5/myEnetModel1.R")

###################################################
#### Load CCLE Molecular Feature Data from Synapse ####
###################################################
id_exprLayer <- "syn1757082" 
layer_expr <- loadEntity(id_exprLayer)
eSet_expr <- layer_expr$objects$eSet_expr

id_copyLayer <- "syn1757086"     
layer_copy <- loadEntity(id_copyLayer)
eSet_copy <- layer_copy$objects$eSet_copy

id_hybridLayer <-  "syn1757084" 
layer_hybrid <- loadEntity(id_hybridLayer)
eSet_hybrid <- layer_hybrid$objects$eSet_hybrid

id_drugLayer <- "syn1757078" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$drugCCLE_ActArea


# featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr, copy = eSet_copy, mut = eSet_oncomap))
featureData <- exprs(eSet_hybrid)

# NA filter for training set
featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "rows")
dataSets <- createFeatureAndResponseDataList(t(featureData_filtered),adf_drug)


load("~/PathwayCuration/InteractDB_newly_curated_subnetwork_per_gene.Rdata")
group<-controled

for(k in 1:length(group)){
  group[[k]]=union(group[[k]],names(group)[k])
}
kk=3

for(kk in 1:24){
  
  #########################################################################################################
  ########  Training and Testing data are scaled(normalized) vs. raw(unnormalized)  #######################
  #########################################################################################################
  
  # data preprocessing for preselecting features
  filteredData<-filterPredictiveModelData(dataSets$featureData,dataSets$responseData[,kk,drop=FALSE], featureVarianceThreshold = 0.01, corPValThresh = 0.1)
  
  # filtered feature and response data
  filteredFeatureData  <- filteredData$featureData
  filteredResponseData <- filteredData$responseData
  
  ## scale these data    
  filteredFeatureDataScaled <- scale(filteredFeatureData)
  filteredResponseDataScaled <- scale(filteredResponseData)  
  
  set.seed(2)
  STEP<-stepwiseDecision(filteredFeatureDataScaled,filteredResponseDataScaled,groups,iterations= 50)
  
  set.seed(2)
  resultsScale<-crossValidatePredictiveModel1(filteredFeatureDataScaled, filteredResponseDataScaled, model = myEnetModel$new(), alpha=1, lambda = lambdas, numFolds=3, nfolds = 3,penalty.factor = STEP$penalty
  save(resultsScale,file = paste("~/SGSR/numFolds_3_nfolds_3/SGL_all_",kk,".Rdata",sep = ""))
  
  
  set.seed(2)
  resultsScale<-crossValidatePredictiveModel1(filteredFeatureDataScaled, filteredResponseDataScaled, model = myEnetModel$new(), alpha=1, lambda = lambdas, numFolds=3,nfolds =3)
  save(resultsScale,file = paste("~/SGSR/numFolds_3_nfolds_3/Lasso_",kk,".Rdata",sep = ""))
  
}