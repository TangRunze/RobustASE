rm(list = ls())
setwd("E:/GitHub/RobustASE/Code/R")

fileName <- paste("../../Data/fish_4.RData")
load(fileName)

source("function_collection.R")
require(parallel)

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
nCores <- 1


labelScreePlot <- c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8",
                    "N1", "N2", "N3", "D1", "D2", "D3")
classScreePlot <- c(rep(1, 8), rep(2, 3), rep(3, 3))
M <- length(labelScreePlot)

nv <- rep(F, length(dd))
for (i in 1:(length(dd))) {
  if (labelVec[i] %in% labelScreePlot) {
    nv[i] <- T
  }
}

labelVec <- labelVec[nv]
dd <- dd[nv]

set.seed(12345)
nv <- sample(1:length(labelVec), 50)

labelVec <- labelVec[nv]
dd <- dd[nv]

K <- 3

tol <- 1e-5
AList <- dd

require(fossil)
tauStar <- sapply(1:length(labelVec), function(iter) {
  if (substr(labelVec[iter], 1, 1) == "G")
    return(1)
  else if (substr(labelVec[iter], 1, 1) == "N")
    return(2)
  else return(3)})


n <- dim(AList[[1]])[1]
source("getElbows.R")
nElbow <- 2
for (i in 1:length(AList)) {
  diag(AList[[i]]) <- 0
}
muList <- AList[sample(1:length(AList), K)]

dist <- matrix(0, K, length(AList))
# Assignment
for (k in 1:K) {
  dist[k, ] <- sapply(1:(length(AList)), function(iter) {norm(AList[[iter]] - muList[[k]], "F")^2})
}
tau0 <- sapply(1:length(AList), function(iter) {which.min(dist[, iter])})
err0 <- sum(apply(dist, 2, min))
maxTol <- tol*err0

maxIter <- 100

###### MLE ######
tau <- tau0
errOld <- err0
errMin <- errOld
tauMLE <- tau
errDiff <- maxTol + 1
iIter <- 0
while ((abs(errDiff) > maxTol) && (iIter < maxIter)) {
  iIter <- iIter + 1
  print(iIter)
  
  # Average
  for (k in 1:K) {
    nv <- (tau == k)
    AListGroup <- AList[nv]
    muList[[k]] <- add(AListGroup)/length(AListGroup)
  }
  rm(AListGroup)
  
  # Assignment
  for (k in 1:K) {
    dist[k, ] <- sapply(1:(length(AList)), function(iter) {norm(AList[[iter]] - muList[[k]], "F")^2})
  }
  tau <- sapply(1:length(AList), function(iter) {which.min(dist[, iter])})
  errNew <- sum(apply(dist, 2, min))
  if (errNew < errMin) {
    print("Update")
    errMin <- errNew
    tauMLE <- tau
  }
  errDiff <- errNew - errOld
  print(c(errNew, errDiff, maxTol, adj.rand.index(tau, tauStar)))
  errOld <- errNew
}





###### MLE_ASE ######
tau <- tau0
errOld <- err0
errMin <- errOld
tauMLEASE <- tau
errDiff <- maxTol + 1
iIter <- 0
while ((abs(errDiff) > maxTol) && (iIter < maxIter)) {
  iIter <- iIter + 1
  print(iIter)
  
  # Average
  for (k in 1:K) {
    nv <- (tau == k)
    AListGroup <- AList[nv]
    ABar <- add(AListGroup)/length(AListGroup)
    ABarDiagAug <- diag_aug(ABar)
    evalVec <- ase(ABarDiagAug, ceiling(n*3/5), isSVD)[[1]]
    dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
    d <- dZG
    AASE <- ase(ABarDiagAug, dZG, isSVD)
    AHat <- AASE[[3]]%*%diag(AASE[[1]])%*%t(AASE[[2]])
    muList[[k]] <- regularize(AHat)
  }
  
  # Assignment
  for (k in 1:K) {
    dist[k, ] <- sapply(1:(length(AList)), function(iter) {norm(AList[[iter]] - muList[[k]], "F")^2})
  }
  tau <- sapply(1:length(AList), function(iter) {which.min(dist[, iter])})
  errNew <- sum(apply(dist, 2, min))
  if (errNew < errMin) {
    print("Update")
    errMin <- errNew
    tauMLEASE <- tau
  }
  errDiff <- errNew - errOld
  print(c(errNew, errDiff, maxTol, adj.rand.index(tau, tauStar)))
  errOld <- errNew
}




