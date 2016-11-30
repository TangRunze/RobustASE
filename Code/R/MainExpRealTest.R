rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")

q <- 0.9
nCores <- 2
# dataName <- "desikan"
dataName <- "CPAC200"
# dataName <- "Talairach"

isSVD <- 0

source("function_collection.R")
require(parallel)

tmpList <- ReadDataWeighted(dataName, DA = F, newGraph = T)
AList <- tmpList[[1]]
n <- tmpList[[2]]
M <- tmpList[[3]]
rm(tmpList)

# dVec <- 1:n
# nD <- length(dVec)
# ASum <- add(AList)

require(Matrix)
pairCrossVec <- combn(M, 2)
pairBetweenVec <- Matrix(1:M, nrow=2)

nCross <- dim(pairCrossVec)[2]
nBetween <- dim(pairBetweenVec)[2]

errorABarCross <- rep(0, nCross)
errorPHatCross <- rep(0, nCross)
errorABarASECross <- rep(0, nCross)
errorPHatASECross <- rep(0, nCross)

errorABarBetween <- rep(0, nBetween)
errorPHatBetween <- rep(0, nBetween)
errorABarASEBetween <- rep(0, nBetween)
errorPHatASEBetween <- rep(0, nBetween)

out <- mclapply(1:nCross, function(iIter) ExpTest(M, AList[pairCrossVec[, iIter]],
                                                  q, isSVD), mc.cores=nCores)
out = array(unlist(out), dim = c(4, nCross))
errorABarCross <- out[1,]
errorPHatCross <- out[2,]
errorABarASECross <- out[3,]
errorPHatASECross <- out[4,]

out <- mclapply(1:nBetween, function(iIter) ExpTest(M, AList[pairBetweenVec[, iIter]],
                                                  q, isSVD), mc.cores=nCores)
out = array(unlist(out), dim = c(4, nBetween))
errorABarBetween <- out[1,]
errorPHatBetween <- out[2,]
errorABarASEBetween <- out[3,]
errorPHatASEBetween <- out[4,]

if (isSVD) {
  fileName = paste("../../Result/result_", dataName, "_test_q_", q, "_svd.RData", sep="")
} else {
  fileName = paste("../../Result/result_", dataName, "_test_q_", q, "_eig.RData", sep="")
}

save(errorABarCross, errorPHatCross, errorABarASECross, errorPHatASECross,
     errorABarBetween, errorPHatBetween, errorABarASEBetween, errorPHatASEBetween,
     n, M, file=fileName)