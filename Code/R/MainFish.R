rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")

# # load("../../Data/gg-v2-14_02_01-deltat50-th0.Rbin")
# load("../../Data/cov-v2-14_02_01-deltat50-th0.RBin")
# load("../../Data/tlab.Rbin")
# fileName <- paste("../../Data/fish.RData")
# save(dd, tlab, file=fileName)

fileName <- paste("../../Data/fish.RData")
load(fileName)

# i <- 1
# j <- 4
# tmp <- sapply(1:length(dd), function(ind) {dd[[ind]][i, j]})
# hist(tmp)

# Consider the odd id of graphs to make them independent
dd <- dd[seq(1, length(dd), 2)]
# Make sure every time frame of 50 time belong to the same category
labelVec <- rep("0", 1, length(dd))
for (i in 1:(length(dd))) {
  if (all(tlab[((i - 1)*length(tlab)/length(dd) + 1):(i*length(tlab)/length(dd))]
          == tlab[(i*length(tlab)/length(dd))])) {
    labelVec[i] <- tlab[(i*length(tlab)/length(dd))]
  }
}

n <- dim(dd[[1]])[1]
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

nvTrain1 <- sapply(1:length(labelVec), function(i) {labelVec[i] %in% train1Vec})
nvTrain2 <- sapply(1:length(labelVec), function(i) {labelVec[i] %in% train2Vec})
nvTest <- sapply(1:length(labelVec), function(i) {labelVec[i] %in% testVec})
AListTrain1 <- dd[nvTrain1]
AListTrain2 <- dd[nvTrain2]
AListTest <- dd[nvTest]
labelTest <- labelVec[nvTest]

source("function_collection.R")
require(parallel)

dVec <- seq(1, n, 5)
if (dVec[length(dVec)] != n) {
  dVec <- c(dVec, n)
}
nD <- length(dVec)

out <- ExpAllDimClassify(AListTrain1, labelTrain1, AListTrain2, labelTrain2,
                         AListTest, labelTest, dVec, q, isSVD)

errorMLE <- out[[1]]/length(AListTest)
errorMLqE <- out[[2]]/length(AListTest)
errorMLEASE <- out[[3]]/length(AListTest)
errorMLEASE_ZG <- out[[4]]/length(AListTest)
errorMLqEASE <- out[[5]]/length(AListTest)
errorMLqEASE_ZG <- out[[6]]/length(AListTest)

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

save(errorMLE, errorMLqE, errorMLEASE, errorMLEASE_ZG,
     errorMLqEASE, errorMLqEASE_ZG, n, dVec, file=fileName)

