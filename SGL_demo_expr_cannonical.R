### DEMO Stepwise grouping Lasso
library(predictiveModeling)
library(synapseClient)
source("~/DrugResponse/R5/crossValidatePredictiveModel1.R")
source("~/SGSR/mySGLModel.R")

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


kk=2

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


# 5 fold cross validation 

groups=list()
for(k in 1:length(NCI)){
  groups[[k]]=nodes(NCI[[k]])
}

set.seed(2)
resultsScale<-crossValidatePredictiveModel1(filteredFeatureDataScaled, filteredResponseDataScaled, model = mySGLModel$new(), alpha=1, lambda = lambdas, numFolds=10, nfolds = 5,iterations=30,groups = groups)
save(resultsScale,file = paste("SGL_NCI_",kk,".Rdata",sep = ""))



groups=list()
for(k in 1:length(BIOCARTA)){
  groups[[k]]=nodes(BIOCARTA[[k]])
}
set.seed(2)
resultsScale<-crossValidatePredictiveModel1(filteredFeatureDataScaled, filteredResponseDataScaled, model = mySGLModel$new(), alpha=1, lambda = lambdas, numFolds=10, nfolds = 5,iterations=30,groups = groups)
save(resultsScale,file = paste("SGL_BIOCARTA_",kk,".Rdata",sep = ""))


groups=list()
for(k in 1:length(KEGG)){
  groups[[k]]=nodes(KEGG[[k]])
}
set.seed(2)
resultsScale<-crossValidatePredictiveModel1(filteredFeatureDataScaled, filteredResponseDataScaled, model = mySGLModel$new(), alpha=1, lambda = lambdas, numFolds=10, nfolds = 5,iterations=30,groups = groups)
save(resultsScale,file = paste("SGL_KEGG_",kk,".Rdata",sep = ""))

groups=list()
for(k in 1:length(REACTOME)){
  groups[[k]]=nodes(REACTOME[[k]])
}
set.seed(2)
resultsScale<-crossValidatePredictiveModel1(filteredFeatureDataScaled, filteredResponseDataScaled, model = mySGLModel$new(), alpha=1, lambda = lambdas, numFolds=10, nfolds = 5,iterations=30,groups = groups)
save(resultsScale,file = paste("SGL_REACTOME_",kk,".Rdata",sep = ""))



corPearsonScale<-c()
corSpearmanScale<-c()

trPred <- foreach(k = 1:5) %do%{resultsScale[[k]]$trainPredictions}
tePred <- foreach(k = 1:5) %do%{resultsScale[[k]]$testPredictions}
trObsr <- foreach(k = 1:5) %do%{resultsScale[[k]]$trainObservations}
teObsr <- foreach(k = 1:5) %do%{resultsScale[[k]]$testObservations}

allTrPred<-do.call("c",trPred)
allTePred<-do.call("c",tePred)
allTrObsr<-do.call("c",trObsr)
allTeObsr<-do.call("c",teObsr)

c(cor(allTrPred,allTrObsr),cor(allTePred,allTeObsr))
c(cor(allTrPred,allTrObsr,method = "spearman"),cor(allTePred,allTeObsr,method = "spearman"))
# 
