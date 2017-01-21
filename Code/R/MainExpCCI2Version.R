rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
# setwd("/cis/home/rtang/RobustASE/Code/R")
source("mylibrary.R")
require(parallel)


###### Parameter Setting ######
# mVec <- c(2, 5, 10)
mVec <- c(2, 5, 10)
q <- 0.9
nIter <- 2
nCores <- 2

isSVD <- 0
isWeighted <- 1

# # dataName2 <- "migrain"
# dataName2 <- "ndmg"
# dataName1 <- "m2g"

dataNameVec <- c("migrain", "ndmg", "m2g")
for (dataName1 in dataNameVec) {
  for (dataName2 in dataNameVec) {
    
    
    ###### Read Data ######
    inputList <- ReadData(dataName1, weighted = isWeighted)
    AList1 <- inputList[[1]]
    n <- inputList[[2]]
    M <- inputList[[3]]
    
    inputList <- ReadData(dataName2, weighted = isWeighted)
    AList2 <- inputList[[1]]
    rm(inputList)
    
    dVec <- 1:n
    nD <- length(dVec)
    
    ###### Main Calculation ######
    for (m in mVec) {
      print(c(m, isSVD))
      
      errorMLE <- matrix(0, nD, nIter)
      errorMLqE <- matrix(0, nD, nIter)
      errorMLEASE <- matrix(0, nD, nIter)
      errorMLqEASE <- matrix(0, nD, nIter)
      out <- mclapply(1:nIter, function(x) ExpRealAllDim(AList1, m, q, isSVD,
                                                         weighted = isWeighted,
                                                         P = add(AList2)/M,
                                                         dVec = dVec),
                      mc.cores=nCores)
      out <- array(unlist(out), dim = c(2*nD+10, nIter))
      
      errorMLE <- out[1, ]
      errorMLEASE <- out[1+(1:nD), ]
      errorMLqE <- out[nD+2, ]
      errorMLqEASE <- out[nD+2+(1:nD), ]
      dimZGMLE <- out[2*nD+3, ]
      dimUSVTMLE <- out[2*nD+4, ]
      dimZGMLqE <- out[2*nD+5, ]
      dimUSVTMLqE <- out[2*nD+6, ]
      errorMLEASEZG <- out[2*nD+7, ]
      errorMLEASEUSVT <- out[2*nD+8, ]
      errorMLqEASEZG <- out[2*nD+9, ]
      errorMLqEASEUSVT <- out[2*nD+10, ]
      
      if (isSVD) {
        strSVD <- "SVD"
      } else {
        strSVD <- "EIG"
      }
      if (isWeighted) {
        strWeighted <- "Weighted"
      } else {
        strWeighted <- "Unweighted"
      }
      fileName = paste0("../../Result/result_brute_", dataName1, "_",
                        dataName2, "_", strWeighted,
                        "_", "m_", m, "_q_", q, "_", strSVD, ".RData")
      
      save(errorMLE, errorMLEASE, errorMLqE, errorMLqEASE,
           errorMLEASEZG, errorMLEASEUSVT, errorMLqEASEZG, errorMLqEASEUSVT, 
           dimZGMLE, dimUSVTMLE, dimZGMLqE, dimUSVTMLqE,
           n, M, m, dVec, isWeighted, nIter, file=fileName)
    }
    
  }
}