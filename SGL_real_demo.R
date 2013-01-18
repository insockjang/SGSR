### DEMO Stepwise grouping Lasso with spliting training and testing data 
library(predictiveModeling)
library(synapseClient)
library(graphite)
library(multicore)
library(doMC)
registerDoMC()

source("~/DrugResponse/R5/crossValidatePredictiveModel1.R")
source("~/DrugResponse/R5/myEnetModel.R")
###################################################
#### Load CCLE Expression Molecular Feature Data from Synapse ####
###################################################
id_exprLayer <- "269056" 
layer_expr <- loadEntity(id_exprLayer)
eSet_expr <- layer_expr$objects$eSet_expr

id_drugLayer <- "269024" 
layer_drug <- loadEntity(id_drugLayer)
adf_drug <- layer_drug$objects$adf_drug


featureData <- exprs(eSet_expr)

# NA filter for training set
featureData_filtered <- filterNasFromMatrix(featureData, filterBy = "rows")
dataSets_ccle <- createFeatureAndResponseDataList(t(featureData_filtered),adf_drug)


# Biocarta is used for testing
load("~/DrugResponse/pathway_analysis/graphite_pathways.Rdata")
groups=list()
for(k in 1:length(BIOCARTA)){
  groups[[k]]=nodes(BIOCARTA[[k]])
}

kk=2


# data preprocessing for preselecting features
filteredData<-filterPredictiveModelData(dataSets_ccle$featureData,dataSets_ccle$responseData[,kk,drop=FALSE], featureVarianceThreshold = 0.01, corPValThresh = 0.1)


set.seed(2)
trainNum<-sample(length(filteredData$responseData),round(length(filteredData$responseData)*2/3))
testNum<-setdiff(c(1:length(filteredData$responseData)),trainNum)
# filtered feature and response data
filteredFeatureData  <- filteredData$featureData[trainNum,]
filteredResponseData <- filteredData$responseData[trainNum]

## scale these data    
filteredFeatureDataScaled <- scale(filteredFeatureData)
filteredResponseDataScaled <- scale(filteredResponseData)  

testFeatureData  <- filteredData$featureData[testNum,]
testResponseData <- filteredData$responseData[testNum]

testFeatureDataScaled <- scale(testFeatureData)
testResponseDataScaled <- scale(testResponseData)  

alphas  = unique(createENetTuneGrid()[,1])
lambdas = createENetTuneGrid(alphas = 1)[,2]

# main algorithm to find the best combination of pathway 
penalty<-rep(1,ncol(filteredFeatureDataScaled))
mse<-list()
k1 <- 0
while(k1 <= 20){ # you can increase the recursive round but I set 20 rounds by default
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

# now check the improvement by selected priors with testing dataset
SGLmodel<-myEnetModel$new()
set.seed(2)
SGLmodel$customTrain(filteredFeatureDataScaled,filteredResponseDataScaled,alpha= 1, lambda = lambdas,penalty.factor = finalPenalty)
SGL<-SGLmodel$customPredict(testFeatureDataScaled)
cor(SGL,testResponseScaled)

# benchmarking Lasso model without priors with testing dataset
Lassomodel<-myEnetModel$new()
set.seed(2)
Lassomodel$customTrain(filteredFeatureDataScaled,filteredResponseDataScaled,alpha= 1, lambda = lambdas)
Lasso<-Lassomodel$customPredict(testFeatureDataScaled)
cor(Lasso,testResponseScaled)


# Done
