rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")

# mVec <- c(1, 2, 5)
mVec <- 5
q <- 0.9
nIter <- 50
nCores <- 2
# dataName <- "desikan"
dataName <- "CPAC200"

isSVD <- 0

source("function_collection.R")
require(parallel)

tmpList <- ReadDataWeighted(dataName, DA = F)
AListTmp <- tmpList[[1]]
AList <- lapply(1:length(AListTmp), function(ind) {log(AListTmp[[ind]] + 2)})
n <- tmpList[[2]]
M <- tmpList[[3]]
rm(tmpList)
rm(AListTmp)


dVec <- 1:n
nD <- length(dVec)
ASum <- add(AList)

for (m in mVec) {
  print(c(m, isSVD))
  
  errorPHat <- matrix(0, nD, nIter)
  errorABar <- matrix(0, nD, nIter)
  errorPhatASE <- matrix(0, nD, nIter)
  errorABarASE <- matrix(0, nD, nIter)
  
  out <- mclapply(1:nIter, function(x) LognormalAllDim(M, m, dVec, AList, ASum, q, isSVD), 
                  mc.cores=nCores)
  out = array(unlist(out), dim = c(2*nD+6, nIter))
  
  errorABar = out[1,]
  errorABarASE = out[1+(1:nD),]
  errorPHat = out[nD+2,]
  errorPhatASE = out[nD+2+(1:nD),]
  dim_ZG_A_bar = out[2*nD+3,]
  dim_USVT_A_bar = out[2*nD+4,]
  dim_ZG_P_hat = out[2*nD+5,]
  dim_USVT_P_hat = out[2*nD+6,]
  
  errorABar_ZG = rep(0, length(dim_ZG_A_bar))
  errorABar_USVT = rep(0, length(dim_USVT_A_bar))
  errorPHat_ZG = rep(0, length(dim_ZG_P_hat))
  errorPHat_USVT = rep(0, length(dim_USVT_P_hat))
  for (i in 1:length(dim_ZG_A_bar)) {
    errorABar_ZG[i] = errorABarASE[dim_ZG_A_bar[i], i]
    errorABar_USVT[i] = errorABarASE[dim_USVT_A_bar[i], i]
    errorPHat_ZG[i] = errorPhatASE[dim_ZG_P_hat[i], i]
    errorPHat_USVT[i] = errorPhatASE[dim_USVT_P_hat[i], i]
  }
  
  if (isSVD) {
    fileName = paste("../../Result/result_", dataName, "_brute_", "m_", m, "_svd.RData", sep="")
  } else {
    fileName = paste("../../Result/result_", dataName, "_brute_", "m_", m, "_eig.RData", sep="")
  }
  
  save(errorABar, errorABarASE, errorPHat, errorPhatASE,
       errorABar_ZG, errorABar_USVT, errorPHat_ZG, errorPHat_USVT, 
       dim_ZG_A_bar, dim_USVT_A_bar, dim_ZG_P_hat, dim_USVT_P_hat,
       n, M, m, dVec, nIter, file=fileName)
}