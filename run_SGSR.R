# run SGSR
library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
K<-setdiff(c(1:138),c(103,14,129,1,39,90,127,123,13,4,110,6,9,12,95,94,108,11,17,2,122,124,126,87,84,75,41,82))
source("~/SGSR_01/PriorIncorporatedLasso_Sanger.R")
PathwayName<-c("NCI")#,"BIOCARTA","NCI","GO_BP","GO_MF")
DataCombine<-c("E")
for(k1 in PathwayName){
  for(k2 in DataCombine){
    PriorIncorporatedLasso_Sanger(k1,k2,KK=K)
  }
}

# run SGSR
library(synapseClient)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
K<-setdiff(c(1:138),c(103,14,129,1,39,90,127,123,13,4,110,6,9,12,95,94,108,11,17,2,122,124,126,87,84,75,41,82))
source("~/SGSR_01/restorePriorIncorporatedLasso_Sanger.R")
source("~/SGSR_01/PriorRandomLasso_Sanger.R")
source("~/SGSR_01/randomPathways_Sanger.R")
source("~/SGSR_01/ENetSameNumber_Sanger.R")

PathwayName<-c("NCI")#,"BIOCARTA","NCI","GO_BP","GO_MF")
DataCombine<-c("E")
for(k1 in PathwayName){
  for(k2 in DataCombine){
    restorePriorIncorporatedLasso_Sanger(k1,k2,KK=K)
    PriorRandomLasso_Sanger(k1,k2,KK=K)
    randomPathways_Sanger(k1,k2,KK=K)    
    ENetSameNumber_Sanger(k1,k2,KK=K)
  }
}


