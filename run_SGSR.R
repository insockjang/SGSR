# run SGSR
# example: NCI database with Mh dataset in CCLE

# Step 01 : stepwise prior selection
source("~/SGSR/SGSR_prior.R")
PathwayName<-c("NCI")
DataCombine<-c("Mh")
for(k1 in PathwayName){
  for(k2 in DataCombine){
    SGSR_CCLE(k1,k2)
  }
}


# Step 02: prediction is separately run with restoreSGSR.R

source("~/SGSR/restoreSGSR.R")
PathwayName<-c("NCI")
DataCombine<-c("Mh")
for(k1 in PathwayName){
  for(k2 in DataCombine){
    restoreSGSR_CCLE(k1,k2)
  }
}

# Step 03 : Elastic Net with fixed sparsity with SGSR

source("~/SGSR/ENetSameNumber.R")
PathwayName<-c("NCI")
DataCombine<-c("Mh")
for(k1 in PathwayName){
  for(k2 in DataCombine){
    ENetSameNumber_CCLE(k1,k2)
  }
}

# Step 04 : Random Genes approach : null distribution for SGSR

source("~/SGSR/randomGene.R")
PathwayName<-c("NCI")
DataCombine<-c("Mh")
for(k1 in PathwayName){
  for(k2 in DataCombine){
    randomGene_CCLE(k1,k2)
  }
}

# Step 05 : Random pathway structure approach : null distribution for SGSR

source("~/SGSR/randomPathways.R")
PathwayName<-c("NCI")
DataCombine<-c("Mh")
for(k1 in PathwayName){
  for(k2 in DataCombine){
    randomPathway_CCLE(k1,k2)
  }
}
