rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
# setwd("E:/GitHub/RobustASE/Code/R")

source("function_collection.R")
require(parallel)

q <- 0.9
dataName <- "fish"
isSVD <- 0

if (isSVD) {
  fileName = paste("../../Result/result_", dataName, "_dissimilarity", "_q_", q, "_svd.RData", sep="")
} else {
  fileName = paste("../../Result/result_", dataName, "_dissimilarity", "_q_", q, "_eig.RData", sep="")
}

load(fileName)

for (i in 1:length(MLEList)) {
  diag(MLEList[[i]]) <- 0
}

disMatrixMLE <- matrix(0, M, M)
disMatrixMLqE <- matrix(0, M, M)
disMatrixMLEASE <- matrix(0, M, M)
disMatrixMLqEASE <- matrix(0, M, M)
for (i in 1:(M - 1)) {
  for (j in (i + 1):M) {
    disMatrixMLE[i, j] <- norm(MLEList[[i]] - MLEList[[j]], "F")
    disMatrixMLE[j, i] <- disMatrixMLE[i, j]
    disMatrixMLqE[i, j] <- norm(MLqEList[[i]] - MLqEList[[j]], "F")
    disMatrixMLqE[j, i] <- disMatrixMLqE[i, j]
    disMatrixMLEASE[i, j] <- norm(MLEASEList[[i]] - MLEASEList[[j]], "F")
    disMatrixMLEASE[j, i] <- disMatrixMLEASE[i, j]
    disMatrixMLqEASE[i, j] <- norm(MLqEASEList[[i]] - MLqEASEList[[j]], "F")
    disMatrixMLqEASE[j, i] <- disMatrixMLqEASE[i, j]
  }
}

require(Matrix)
image(Matrix(disMatrixMLE))
image(Matrix(disMatrixMLqE))
image(Matrix(disMatrixMLEASE))
image(Matrix(disMatrixMLqEASE))


xMLE <- t(sapply(1:length(MLEList), function(i) {as.vector(MLEList[[i]])}))
xMLqE <- t(sapply(1:length(MLqEList), function(i) {as.vector(MLqEList[[i]])}))
xMLEASE <- t(sapply(1:length(MLEASEList), function(i) {as.vector(MLEASEList[[i]])}))
xMLqEASE <- t(sapply(1:length(MLqEASEList), function(i) {as.vector(MLqEASEList[[i]])}))

# strMethod <- "euclidean"
# strMethod <- "maximum"
strMethod <- "manhattan"
hMLE <- hclust(dist(xMLE, method=strMethod))
plot(hMLE, labels = labelScreePlot, main=paste0("MLE, ", strMethod))
hMLqE <- hclust(dist(xMLqE, method=strMethod))
plot(hMLqE, labels = labelScreePlot, main=paste0("MLqE, ", strMethod))
hMLEASE <- hclust(dist(xMLEASE, method=strMethod))
plot(hMLEASE, labels = labelScreePlot, main=paste0("MLEASE, ", strMethod))
hMLqEASE <- hclust(dist(xMLqEASE, method=strMethod))
plot(hMLqEASE, labels = labelScreePlot, main=paste0("MLqEASE, ", strMethod))


hMLE <- hclust(as.dist(disMatrixMLE))
plot(hMLE, labels = labelScreePlot, main="MLE")
hMLqE <- hclust(as.dist(disMatrixMLqE))
plot(hMLqE, labels = labelScreePlot, main="MLqE")
hMLEASE <- hclust(as.dist(disMatrixMLEASE))
plot(hMLEASE, labels = labelScreePlot, main="MLEASE")
hMLqEASE <- hclust(as.dist(disMatrixMLqEASE))
plot(hMLqEASE, labels = labelScreePlot, main="MLqEASE")

tmp <- MLEList



tmp <- matrix(c(1,2,3,4,5,6,7,8), ncol=2)
dist(tmp, method="euclidean")


MLEList[[1]][1:5,1:5]
MLqEList[[1]][1:5,1:5]
MLEASEList[[1]][1:5,1:5]
MLqEASEList[[1]][1:5,1:5]
