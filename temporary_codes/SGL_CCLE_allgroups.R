### DEMO Stepwise grouping Lasso
library(predictiveModeling)
library(synapseClient)
library(graphite)
source("~/DrugResponse/R5/crossValidatePredictiveModel1.R")
source("~/SGSR/mySGLModel.R")
source("~/DrugResponse/R5/myEnetModel.R")

library(multicore)
library(doMC)
registerDoMC()

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


# featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr, copy = eSet_copy, mut = eSet_oncomap))
featureData <- exprs(eSet_expr)

# NA filter for training set
featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "rows")
dataSets_ccle <- createFeatureAndResponseDataList(t(featureData_filtered),adf_drug)


load("~/DrugResponse/pathway_analysis/graphite_pathways.Rdata")
alphas  = unique(createENetTuneGrid()[,1])
lambdas = createENetTuneGrid(alphas = 1)[,2]


# 5 fold cross validation 
a1<-length(NCI)
a2<-length(KEGG)
a3<-length(BIOCARTA)
a4<-length(REACTOME)

groups=list()
for(k in 1:length(NCI)){
  groups[[k]]=nodes(NCI[[k]])
}

for(k in (a1+1):(a1+length(KEGG))){
  groups[[k]]=nodes(KEGG[[k-a1]])
}

for(k in (a1+a2+1):(a1+a2+length(BIOCARTA))){
  groups[[k]]=nodes(BIOCARTA[[k-(a1+a2)]])
}

for(k in (a1+a2+a3+1):(a1+a2+a3+length(REACTOME))){
  groups[[k]]=nodes(REACTOME[[k-(a1+a2+a3)]])
}

for(kk in 1:24){
  
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
  filteredResponseDataScaled <- (filteredResponseData)  
  
  
  set.seed(2)
  resultsScale<-crossValidatePredictiveModel1(filteredFeatureDataScaled, filteredResponseDataScaled, model = mySGLModel$new(), alpha=1, lambda = lambdas, numFolds=3, nfolds = 3,iterations=length(groups),groups = groups)
  save(resultsScale,file = paste("~/SGSR/numFolds_3_nfolds_3/SGL_all_",kk,".Rdata",sep = ""))
  
  
  set.seed(2)
  resultsScale<-crossValidatePredictiveModel1(filteredFeatureDataScaled, filteredResponseDataScaled, model = myEnetModel$new(), alpha=1, lambda = lambdas, numFolds=3,nfolds =3)
  save(resultsScale,file = paste("~/SGSR/numFolds_3_nfolds_3/Lasso_",kk,".Rdata",sep = ""))
  
}