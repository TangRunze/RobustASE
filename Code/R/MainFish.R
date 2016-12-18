rm(list = ls())
# setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
setwd("E:/GitHub/RobustASE/Code/R")
# setwd("/cis/home/rtang/RobustASE/Code/R")

# # load("../../Data/gg-v2-14_02_01-deltat50-th0.Rbin")
# # load("../../Data/cov-v2-14_02_01-deltat50-th0.RBin")
# load("../../Data/cov-v2-14_02_01-deltat4-th0.RBin")
# load("../../Data/tlab.Rbin")
# fileName <- paste("../../Data/fish_4.RData")
# dd <- lapply(1:(length(dd)), function(i) {abs(dd[[i]])})
# save(dd, tlab, file=fileName)

# Consider the odd id of graphs to make them independent
# dd <- dd[seq(1, length(dd), 2)]

# fileName <- paste("../../Data/fish.RData")
fileName <- paste("../../Data/fish_4.RData")
load(fileName)

source("function_collection.R")
require(parallel)
nCores <- 1

###### Histogram ######
# i <- 1
# j <- 4
# tmp <- sapply(1:length(dd), function(ind) {dd[[ind]][i, j]})
# nv = (tmp <= 1)
# tmp = tmp[nv]
# hist(tmp)

# Make sure every time frame of 4 time belong to the same category
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

# Scree-plot
# labelScreePlot <- c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8",
#                     "N1", "N2", "N3", "D1", "D2", "D3")
# classScreePlot <- c(rep(1, 8), rep(2, 3), rep(3, 3))
# M <- length(labelScreePlot)
# source("getElbows.R")
# require(irlba)
# dMax = floor(n/2)
# elbMat = matrix(0, M, 3)
# evalMat = matrix(0, M, dMax)
# eval3Mat = matrix(0, M, 3)
# plot(1:dMax, type="n", ylim=c(0,100), xlab="embedding dimension", ylab="eigen value")
# for (i in 1:M) {
#   print(labelScreePlot[i])
#   nv <- sapply(1:length(labelVec), function(iter) {labelVec[iter] == labelScreePlot[i]})
#   AList <- dd[nv]
#   A <- add(AList)/length(AList)
#   vecs <- irlba(A, dMax, dMax)$d
#   # A = as.matrix(A_all[[i]])
#   # vecs = eigs_sym(A, dMax, which = "LM")$values
#   evalMat[i,] <- vecs
#   elb <- getElbows(vecs, plot=F)
#   elbMat[i,] <- elb
#   eval3Mat[i,] <- vecs[elb]
#   # points(vecs, type="l", col="grey")
#   points(vecs, type="l", col=classScreePlot[i]+4)
#   points(elb, vecs[elb], pch=19, col=2:4, cex=1)
# }

# plot(1:20, type="n", ylim=c(0,10), xlab="embedding dimension", ylab="eigen value")
# for (i in 1:M) {
#   points(evalMat[i,], type="l", col=(classScreePlot[i]+4))
#   points(elbMat[i,], evalMat[i,elb], pch=19, col=2:4, cex=1)
# }

# ceiling(median(elbMat[,3]))



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

dVec <- seq(1, n, 10)
if (dVec[length(dVec)] != n) {
  dVec <- c(dVec, n)
}
nD <- length(dVec)

out <- ExpAllDimClassify(AListTrain1, labelTrain1, AListTrain2, labelTrain2,
                         AListTest, labelTest, dVec, q, nCores, isSVD)

errorMLEVec <- out[[1]]
errorMLqEVecVec <- out[[2]]
errorMLEASEVec <- out[[3]]
errorMLEASE_ZGVec <- out[[4]]
errorMLqEASEVec <- out[[5]]
errorMLqEASE_ZGVec <- out[[6]]

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

