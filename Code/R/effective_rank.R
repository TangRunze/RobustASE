rm(list = ls())

setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")

dataNameVec <- c("m2g", "ndmg2")
dataNameDisplayVec <- c("m2g", "ndmg2")

source("mylibrary.R")
require(Matrix)
eRank <- list()
for (iData in 1:length(dataNameVec)) {
  dataName <- dataNameVec[iData]
  tmpList <- ReadData(dataName)
  AList <- tmpList[[1]]
  n <- tmpList[[2]]
  M <- tmpList[[3]]
  rm(tmpList)
  Abar <- add(AList)/M
  D0 <- Diagonal(n, x = rowSums(Abar)/(n-1))
  singularValueVec <- abs(eigen(Abar + D0)$values)
  pVec <- singularValueVec/sum(singularValueVec)
  eRank[[iData]] <- exp(-sum(pVec*log(pVec)))
}

