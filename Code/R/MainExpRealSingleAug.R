rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")

mVec <- c(2, 5)
# mVec <- 5
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
AList0 <- AList
n <- tmpList[[2]]
M <- tmpList[[3]]
rm(tmpList)

ASum0 <- add(AList0)

for (i in 1:M) {
  print(i)
  d <- 55
  # d <- 20
  AASE <- ase(AList[[i]], d, isSVD)
  if (d == 1) {
    AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
  } else {
    AHat <- AASE[[3]][ ,1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][ ,1:d])
  }
  AList[[i]] <- regularize(AHat)
}

dVec <- 1:n
nD <- length(dVec)
ASum <- add(AList)

for (m in mVec) {
  print(c(m, isSVD))
  
  errorABar <- matrix(0, nD, nIter)
  errorPHat <- matrix(0, nD, nIter)
  errorABarASE <- matrix(0, nD, nIter)
  errorPHatASE <- matrix(0, nD, nIter)
  errorABar0 <- matrix(0, nD, nIter)
  errorPHat0 <- matrix(0, nD, nIter)
  errorABarASE0 <- matrix(0, nD, nIter)
  errorPHatASE0 <- matrix(0, nD, nIter)
  
  out <- mclapply(1:nIter, function(x) ExpAllDimSingleAug(M, m, dVec, AList, ASum,
                                                          AList0, ASum0, q, isSVD),
                  mc.cores=nCores)
  out = array(unlist(out), dim = c(4*nD+8, nIter))
  
  errorABar = out[1,]
  errorABarASE = out[1+(1:nD),]
  errorPHat = out[nD+2,]
  errorPHatASE = out[nD+2+(1:nD),]
  dimZGABar = out[2*nD+3,]
  dimUSVTABar = out[2*nD+4,]
  dimZGPHat = out[2*nD+5,]
  dimUSVTPHat = out[2*nD+6,]
  errorABar0 = out[2*nD+7,]
  errorABarASE0 = out[2*nD+7+(1:nD),]
  errorPHat0 = out[3*nD+8,]
  errorPHatASE0 = out[3*nD+8+(1:nD),]
  
  errorABarZG = rep(0, length(dimZGABar))
  errorABarUSVT = rep(0, length(dimUSVTABar))
  errorPHatZG = rep(0, length(dimZGPHat))
  errorPHatUSVT = rep(0, length(dimUSVTPHat))
  errorABarZG0 = rep(0, length(dimZGABar))
  errorABarUSVT0 = rep(0, length(dimUSVTABar))
  errorPHatZG0 = rep(0, length(dimZGPHat))
  errorPHatUSVT0 = rep(0, length(dimUSVTPHat))
  for (i in 1:length(dimZGABar)) {
    errorABarZG[i] = errorABarASE[dimZGABar[i], i]
    errorABarUSVT[i] = errorABarASE[dimUSVTABar[i], i]
    errorPHatZG[i] = errorPHatASE[dimZGPHat[i], i]
    errorPHatUSVT[i] = errorPHatASE[dimUSVTPHat[i], i]
    
    errorABarZG0[i] = errorABarASE0[dimZGABar[i], i]
    errorABarUSVT0[i] = errorABarASE0[dimUSVTABar[i], i]
    errorPHatZG0[i] = errorPHatASE0[dimZGPHat[i], i]
    errorPHatUSVT0[i] = errorPHatASE0[dimUSVTPHat[i], i]
  }
  
  if (isSVD) {
    fileName = paste("../../Result/result_", dataName, "_brute_", "m_", m,
                     "_q_", q, "_svd_single_aug_", d, ".RData", sep="")
  } else {
    fileName = paste("../../Result/result_", dataName, "_brute_", "m_", m,
                     "_q_", q, "_eig_single_aug_", d, ".RData", sep="")
  }
  
  save(errorABar, errorABarASE, errorPHat, errorPHatASE,
       errorABar0, errorABarASE0, errorPHat0, errorPHatASE0,
       errorABarZG, errorABarUSVT, errorPHatZG, errorPHatUSVT, 
       errorABarZG0, errorABarUSVT0, errorPHatZG0, errorPHatUSVT0, 
       dimZGABar, dimUSVTABar, dimZGPHat, dimUSVTPHat,
       n, M, m, dVec, nIter, file=fileName)
}