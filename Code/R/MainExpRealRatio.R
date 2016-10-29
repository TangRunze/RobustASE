rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")

# mVec <- c(1, 2, 5)
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
    fileName = paste("../../Result/result_", dataName, "_brute_", "m_", m, "_q_", q, "_svd.RData", sep="")
  } else {
    fileName = paste("../../Result/result_", dataName, "_brute_", "m_", m, "_q_", q, "_eig.RData", sep="")
  }
  
  save(errorABar, errorABarASE, errorPHat, errorPHatASE,
       errorABarZG, errorABarUSVT, errorPHatZG, errorPHatUSVT, 
       dimZGABar, dimUSVTABar, dimZGPHat, dimUSVTPHat,
       n, M, m, dVec, nIter, file=fileName)
}