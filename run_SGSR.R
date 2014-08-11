
# run SGSR
source("/home/ubuntu/ijang/SGSR/PriorIncorporatedLasso_Sanger.R")
source("/home/ubuntu/ijang/SGSR/ENetSameNumber_Sanger.R")
source("/home/ubuntu/ijang/SGSR/randomPathways_Sanger.R")
source("/home/ubuntu/ijang/SGSR/PriorRandomLasso_Sanger.R")
PathwayName<-c("KEGG")
DataCombine<-c("C")
ALPHAS=c(0.1)
for(k1 in PathwayName){
  for(k2 in DataCombine){
    PriorRandomLasso_Sanger(k1,k2)
    randomPathway_Sanger(k1,k2)
    for(k3 in ALPHAS){
      ENetSameNumber_Sanger(k1,k2,ALPHA = k3)
    }
  }
}




# run SGSR
source("/home/ubuntu/ijang/SGSR/PriorIncorporatedLasso_Sanger.R")
PathwayName<-c("NCI")
DataCombine<-c("C")
for(k1 in PathwayName){
  for(k2 in DataCombine){
    PriorIncorporatedLasso_Sanger(k1,k2)
  }
}



source("~/SGSR/temporary_codes/Lasso_CCLE.R")
source("~/SGSR/temporary_codes/ENet_CCLE.R")
DataCombine<-c("Mh","E","C")
for(k2 in DataCombine){
  Lasso_CCLE(k2)
  Lasso_ENet(k2)
}

source("~/SGSR/temporary_codes/Lasso_Sanger.R")
source("~/SGSR/temporary_codes/ENet_Sanger.R")
DataCombine<-c("E","C","Mh")
for(k2 in DataCombine){  
  Lasso_Sanger(k2)  
  ENet_Sanger(k2)
}




# run SGSR
source("/home/ubuntu/ijang/SGSR/PriorIncorporatedLasso_Sanger.R")
PathwayName<-c("GO_BP","GO_MF")
DataCombine<-c("Mh")
for(k1 in PathwayName){
  for(k2 in DataCombine){
    PriorIncorporatedLasso_Sanger(k1,k2)
  }
}


# run SGSR
source("/home/ubuntu/ijang/SGSR/PriorIncorporatedLasso_CCLE.R")
PathwayName<-c("GO_BP","GO_MF")
DataCombine<-c("Mh")
for(k1 in PathwayName){
  for(k2 in DataCombine){
    PriorIncorporatedLasso_CCLE(k1,k2)
  }
}


library(predictiveModeling)
synapseLogin("in.sock.jang@sagebase.org","tjsDUD@")
