rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")

# load("../../Data/gg-v2-14_02_01-deltat50-th0.Rbin")
# load("../../Data/tlab.Rbin")
# fileName <- paste("../../Data/fish.RData")
# save(dd, tlab, file=fileName)

fileName <- paste("../../Data/fish.RData")
load(fileName)

# i <- 2
# j <- 5
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
labelTrain1 <- substr(train1Vec, 1, 1)
train2Vec <- c("G4")
labelTrain2 <- substr(train2Vec, 1, 1)
testVec <- c("D3", "G5")

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

errorABar <- matrix(0, nD, nIter)
errorPHat <- matrix(0, nD, nIter)
errorABarASE <- matrix(0, nD, nIter)
errorPHatASE <- matrix(0, nD, nIter)

out <- mclapply(1:nIter, function(x) ExpAllDim(M, m, dVec, AList, ASum, q, isSVD), 
                mc.cores=nCores)
out = array(unlist(out), dim = c(2*nD+6, nIter))

errorABar = out[1,]
errorABarASE = out[1+(1:nD),]
errorPHat = out[nD+2,]
errorPHatASE = out[nD+2+(1:nD),]
dimZGABar = out[2*nD+3,]
dimUSVTABar = out[2*nD+4,]
dimZGPHat = out[2*nD+5,]
dimUSVTPHat = out[2*nD+6,]

errorABarZG = rep(0, length(dimZGABar))
errorABarUSVT = rep(0, length(dimUSVTABar))
errorPHatZG = rep(0, length(dimZGPHat))
errorPHatUSVT = rep(0, length(dimUSVTPHat))
for (i in 1:length(dimZGABar)) {
  errorABarZG[i] = errorABarASE[dimZGABar[i], i]
  errorABarUSVT[i] = errorABarASE[dimUSVTABar[i], i]
  errorPHatZG[i] = errorPHatASE[dimZGPHat[i], i]
  errorPHatUSVT[i] = errorPHatASE[dimUSVTPHat[i], i]
}

if (isSVD) {
  fileName = paste("../../Result/result_", dataName, "_new_brute_", "m_", m, "_q_", q, "_svd.RData", sep="")
} else {
  fileName = paste("../../Result/result_", dataName, "_new_brute_", "m_", m, "_q_", q, "_eig.RData", sep="")
}

save(errorABar, errorABarASE, errorPHat, errorPHatASE,
     errorABarZG, errorABarUSVT, errorPHatZG, errorPHatUSVT, 
     dimZGABar, dimUSVTABar, dimZGPHat, dimUSVTPHat,
     n, M, m, dVec, nIter, file=fileName)
