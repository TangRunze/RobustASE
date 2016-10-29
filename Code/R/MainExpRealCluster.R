rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")

# mVec <- c(2, 5)
mVec <- 2
q <- 0.9
nIter <- 200
nCores <- 2
# dataName <- "desikan"
dataName <- "CPAC200"

isSVD <- 0

source("function_collection.R")
require(parallel)

tmpList <- ReadDataWeighted(dataName, DA = F)
AList <- tmpList[[1]]
n <- tmpList[[2]]
M <- tmpList[[3]]
rm(tmpList)


# source("getElbows.R")
# # require(rARPACK)
# require(irlba)
# dMax = n-10
# elbMat = matrix(0, M, 3)
# evalMat = matrix(0, M, dMax)
# eval3Mat = matrix(0, M, 3)
# plot(1:dMax, type="n", ylim=c(0,50), xlab="embedding dimension", ylab="eigen value")
# for (i in 1:M) {
#   A = AList[[i]]
#   vecs = irlba(A, dMax, dMax)$d
#   # A = as.matrix(A_all[[i]])
#   # vecs = eigs_sym(A, dMax, which = "LM")$values
#   evalMat[i,] = vecs
#   elb = getElbows(vecs, plot=F)
#   elbMat[i,] = elb
#   eval3Mat[i,] = vecs[elb]
#   points(vecs, type="l", col="grey")
#   points(elb, vecs[elb], pch=19, col=2:4, cex=0.5)
# }
# save(evalMat, file="eval.RData")

load("eval.RData")

library(mclust)
GMM = Mclust(evalMat[,1:2], 2)

plot(1:(n-10), type="n", ylim=c(0,50), xlab="embedding dimension", ylab="eigen value")
for (i in 1:M) {
  if (GMM$classification[i] == 1) {
    colStr = "grey72"
  } else {
    colStr = "grey52"
  }
  points(evalMat[i,], type="l", col=colStr)
#   points(elbMat[i,], eval3Mat[i,], pch=19, col=2:4, cex=0.5)
  GMM$classification == 1
}
title("2 Clusters using the first 2 eigenvalues")

cl <- GMM$classification
AList1 <- list()
AList2 <- list()
M1 <- 0
M2 <- 0

for (i in 1:M) {
  if (cl[[i]] == 1) {
    M1 <- M1 + 1
    AList1[[M1]] <- AList[[i]]
  } else {
    M2 <- M2 + 1
    AList2[[M2]] <- AList[[i]]
  }
}


###### Choose Class 1 ######
AList <- AList1
M <- M1

###### Choose Class 2 ######
# AList <- AList2
# M <- M2


dVec <- 1:n
nD <- length(dVec)
ASum <- add(AList)

for (m in mVec) {
  print(c(m, isSVD))
  
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
    fileName = paste("../../Result/result_", dataName, "_brute_", "m_", m, "_q_", q, "_svd_cluster1.RData", sep="")
  } else {
    fileName = paste("../../Result/result_", dataName, "_brute_", "m_", m, "_q_", q, "_eig_cluster1.RData", sep="")
  }
  
  save(errorABar, errorABarASE, errorPHat, errorPHatASE,
       errorABarZG, errorABarUSVT, errorPHatZG, errorPHatUSVT, 
       dimZGABar, dimUSVTABar, dimZGPHat, dimUSVTPHat,
       n, M, m, dVec, nIter, file=fileName)
}