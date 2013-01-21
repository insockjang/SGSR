### DEMO Stepwise grouping Lasso
library(predictiveModeling)
library(synapseClient)
library(graphite)

library(multicore)
library(doMC)
registerDoMC()

source("~/DrugResponse/R5/crossValidatePredictiveModel1.R")
source("~/SGR/SGL/mySGLModel.R")
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


featureData <- createAggregateFeatureDataSet(list(expr = eSet_expr, copy = eSet_copy, mut = eSet_oncomap))

# NA filter for training set
featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "rows")
dataSets_ccle <- createFeatureAndResponseDataList(t(featureData_filtered),adf_drug)


load("~/DrugResponse/pathway_analysis/graphite_pathways.Rdata")
groups=list()
for(k in 1:length(BIOCARTA)){
  a=nodes(BIOCARTA[[k]])
  a1<-paste(a,"_expr",sep="")
  a2<-paste(a,"_copy",sep="")
  a3<-paste(a,"_mut",sep="")
  aa<-union(a1,union(a2,a3))
  groups[[k]]<-aa
}

kk=1  

#########################################################################################################
########  Training and Testing data are scaled(normalized) vs. raw(unnormalized)  #######################
#########################################################################################################

# data preprocessing for preselecting features
filteredData<-filterPredictiveModelData(dataSets_ccle$featureData,dataSets_ccle$responseData[,kk,drop=FALSE])

# filtered feature and response data
filteredFeatureData  <- filteredData$featureData
filteredResponseData <- filteredData$responseData

## scale these data    
filteredFeatureDataScaled <- scale(filteredFeatureData)
#filteredResponseDataScaled <- scale(filteredResponseData)  
filteredResponseDataScaled <- (filteredResponseData)  


# 5 fold cross validation 

alphas  = unique(createENetTuneGrid()[,1])
lambdas = createENetTuneGrid(alphas = 1)[,2]

set.seed(2)
resultsScale<-crossValidatePredictiveModel1(filteredFeatureDataScaled, filteredResponseDataScaled, model = mySGLModel$new(), alpha=1, lambda = lambdas, nfolds = 5,iterations=5,groups = groups)

set.seed(2)
resultsScale1<-crossValidatePredictiveModel1(filteredFeatureDataScaled, filteredResponseDataScaled, model = myEnetModel$new(), alpha=1, lambda = lambdas, nfolds = 5)



corPearsonScale<-c()
corSpearmanScale<-c()

trPred <- foreach(k = 1:5) %do%{resultsScale1[[k]]$trainPredictions}
tePred <- foreach(k = 1:5) %do%{resultsScale1[[k]]$testPredictions}
trObsr <- foreach(k = 1:5) %do%{resultsScale1[[k]]$trainObservations}
teObsr <- foreach(k = 1:5) %do%{resultsScale1[[k]]$testObservations}

allTrPred<-do.call("c",trPred)
allTePred<-do.call("c",tePred)
allTrObsr<-do.call("c",trObsr)
allTeObsr<-do.call("c",teObsr)

c(cor(allTrPred,allTrObsr),cor(allTePred,allTeObsr))
c(cor(allTrPred,allTrObsr,method = "spearman"),cor(allTePred,allTeObsr,method = "spearman"))

