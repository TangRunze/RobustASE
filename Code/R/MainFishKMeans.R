rm(list = ls())
# setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
setwd("E:/GitHub/RobustASE/Code/R")
# setwd("/cis/home/rtang/RobustASE/Code/R")

nSample <- 100
iReplicate <- 1
nCores <- 1

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

# require(Matrix)
# image(Matrix(dd[[1]]))

source("function_collection.R")
require(parallel)

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

set.seed(12345 + iReplicate)
nv <- sample(1:length(labelVec), nSample)

labelVec <- labelVec[nv]
dd <- dd[nv]

K <- 3

result <- ExpKMeans(dd, K, q, labelVec, nCores, isSVD)
tauMLE <- result[[1]]
tauMLqE <- result[[2]]
tauMLEASE <- result[[3]]
tauMLqEASE <- result[[4]]
tauStar <- result[[5]]
muListMLE <- result[[6]]
muListMLqE <- result[[7]]
muListMLEASE <- result[[8]]
muListMLqEASE <- result[[9]]


if (isSVD) {
  fileName = paste("../../Result/result_", dataName, "_kmeans", "_q_", q, "_sample_", nSample,
                   "_K_", K, "_replicate_", iReplicate, "_svd.RData", sep="")
} else {
  fileName = paste("../../Result/result_", dataName, "_kmeans", "_q_", q, "_sample_", nSample,
                   "_K_", K, "_replicate_", iReplicate, "_eig.RData", sep="")
}

indSample <- nv
# save(q, K, indSample, nSample, tauMLE, tauMLEASE, tauStar, muListMLE, muListMLEASE, file=fileName)
save(q, K, labelVec, indSample, nSample, tauMLE, tauMLqE, tauMLEASE, tauMLqEASE, tauStar,
     muListMLE, muListMLqE, muListMLEASE, muListMLqEASE, file=fileName)
