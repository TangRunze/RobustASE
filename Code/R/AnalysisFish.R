rm(list = ls())
# setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
setwd("E:/GitHub/RobustASE/Code/R")

q <- 0.9
dataName <- "fish"
isSVD <- 0

train1Vec <- c("D2")
train2Vec <- c("G4")
testVec <- c("D3", "G5")
# train1Vec <- c("D1")
# train2Vec <- c("G2")
# testVec <- c("D2", "D3", "G4", "G5")
labelTrain1 <- substr(train1Vec, 1, 1)
labelTrain2 <- substr(train2Vec, 1, 1)

if (isSVD) {
  fileName = paste("../../Result/result_", dataName, "_brute_",
                   paste(train1Vec, collapse=""), "_",
                   paste(train2Vec, collapse=""), "_",
                   paste(testVec, collapse=""), "_q_", q, "_svd.RData", sep="")
} else {
  fileName = paste("../../Result/result_", dataName, "_brute_",
                   paste(train1Vec, collapse=""), "_",
                   paste(train2Vec, collapse=""), "_",
                   paste(testVec, collapse=""), "_q_", q, "_eig.RData", sep="")
}

load(fileName)

source("function_collection.R")

mean(errorMLEVec)
mean(errorMLqEVec)
mcnemar.test(BuildMcNemarMatrix(errorMLEVec, errorMLqEVec))

mean(errorMLEVec)
rowMeans(errorMLEASEVec)

mean(errorMLEVec)
mean(errorMLEASE_ZGVec)
mcnemar.test(BuildMcNemarMatrix(errorMLEVec, errorMLEASE_ZGVec))

mean(errorMLqEVec)
rowMeans(errorMLqEASEVec)

mean(errorMLqEVec)
mean(errorMLqEASE_ZGVec)
mcnemar.test(BuildMcNemarMatrix(errorMLqEVec, errorMLqEASE_ZGVec))

