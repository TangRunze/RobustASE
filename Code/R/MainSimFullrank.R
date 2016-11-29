
rm(list = ls())

setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")

mVec = c(1,5,10)
isSVD = 1

nIter = 100
nCores = 2

dataName = "CPAC200"
# dataName = "desikan"
# dataName = "JHU"
# dataName = "slab907"
# dataName = "slab1068"
# dataName = "Talairach"

source("function_collection.R")
require(parallel)

tmpList <- ReadDataWeighted(dataName, DA = F, newGraph = T)
A_all = tmpList[[1]]
n = tmpList[[2]]
M = tmpList[[3]]
rm(tmpList)

P = add(A_all)/M
diag(P) = 0

dVec = 1:n
nD = length(dVec)

for (m in mVec) {
  print(c(m, isSVD))
  
  error_P_hat = matrix(0, nD, nIter)
  error_A_bar = matrix(0, nD, nIter)
  
  out <- mclapply(1:nIter, function(x) dim_brute_fullrank(m, dVec, P, isSVD), 
                  mc.cores=nCores)
  out = array(unlist(out), dim = c(nD+1, nIter))
  
  error_A_bar = out[1,]
  error_P_hat = out[2:(nD+1),]
  
  if (isSVD) {
    fileName = paste("../../Result/result_", dataName, "_brute_fullrank1_",
                     "m_", m, "_svd.RData", sep="")
  } else {
    fileName = paste("../../Result/result_", dataName, "_brute_fullrank1_",
                     "m_", m, "_eig.RData", sep="")
  }
  
  save(error_A_bar, error_P_hat, n, m, dVec, nIter, file=fileName)
  
}





nElbow = 3
nIter = 100

source("getElbows.R")
source("USVT.R")

dZG = matrix(0, length(mVec), nIter)
dUSVT = matrix(0, length(mVec), nIter)

for (iM in 1:length(mVec)) {
  m = mVec[iM]
  for (iIter in 1:nIter) {
    print(c(m, iIter))
    A_all = list()
    for (i in 1:m) {
      g = sample_sbm(n, P, rep(1,n), directed=F, loops=F)
      A = as_adj(g, type="both", sparse=FALSE)
      A_all[[i]] = A
    }
    A_bar = add(A_all)/m
    evalVec = ase(A_bar, ceiling(n*3/5), isSVD)[[1]]
    dZG[iM, iIter] = getElbows(evalVec, n=nElbow, plot=F)[[nElbow]]
    dUSVT[iM, iIter] = length(usvt(A_bar, 1, m)$d)
  }
}

dZGMean = rep(0,length(mVec))
dZGL = rep(0,length(mVec))
dZGU = rep(0,length(mVec))
dUSVTMean = rep(0,length(mVec))
dUSVTL = rep(0,length(mVec))
dUSVTU = rep(0,length(mVec))

for (iM in 1:length(mVec)) {
  dZGMean[iM] = mean(dZG[iM,])
  dZGL[iM] = mean(dZG[iM,]) - sqrt(var(dZG[iM,])/nIter)*1.96
  dZGU[iM] = mean(dZG[iM,]) + sqrt(var(dZG[iM,])/nIter)*1.96
  
  dUSVTMean[iM] = mean(dUSVT[iM,])
  dUSVTL[iM] = mean(dUSVT[iM,]) - sqrt(var(dUSVT[iM,])/nIter)*1.96
  dUSVTU[iM] = mean(dUSVT[iM,]) + sqrt(var(dUSVT[iM,])/nIter)*1.96
}

if (isSVD) {
  fileName = paste("../../Result/result_", dataName, "_brute_fullrank1_",
                   "_svd_dim.RData", sep="")
} else {
  fileName = paste("../../Result/result_", dataName, "_brute_fullrank1_",
                   "_eig_dim.RData", sep="")
}

save(dZGMean, dZGL, dZGU, dUSVTMean, dUSVTL, dUSVTU, file=fileName)
