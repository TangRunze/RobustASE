rm(list = ls())
# setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
setwd("E:/GitHub/RobustASE/Code/R")

source("function_collection.R")

nSample <- 50
iReplicate <- 5
K <- 3
q <- 0.9
dataName <- "fish"
isSVD <- 0

if (isSVD) {
  fileName = paste("../../Result/result_", dataName, "_kmeans_threshold", "_q_", q, "_sample_", nSample,
                   "_K_", K, "_replicate_", iReplicate, "_svd.RData", sep="")
} else {
  fileName = paste("../../Result/result_", dataName, "_kmeans_threshold", "_q_", q, "_sample_", nSample,
                   "_K_", K, "_replicate_", iReplicate, "_eig.RData", sep="")
}
load(fileName)

require(fossil)
adj.rand.index(tauMLE, tauStar)
adj.rand.index(tauMLqE, tauStar)
adj.rand.index(tauMLEASE, tauStar)
adj.rand.index(tauMLqEASE, tauStar)

mean(muListMLE[[1]])
mean(muListMLE[[2]])
mean(muListMLE[[3]])

tauMLE



# 
# fileName <- paste("../../Data/fish_4.RData")
# load(fileName)
# labelVec <- rep("0", 1, length(dd))
# for (i in 1:(length(dd))) {
#   if (all(tlab[((i - 1)*length(tlab)/length(dd) + 1):(i*length(tlab)/length(dd))]
#           == tlab[(i*length(tlab)/length(dd))])) {
#     labelVec[i] <- tlab[(i*length(tlab)/length(dd))]
#   }
# }
# # Scree-plot
# labelScreePlot <- c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8",
#                     "N1", "N2", "N3", "D1", "D2", "D3")
# classScreePlot <- c(rep(1, 8), rep(2, 3), rep(3, 3))
# M <- length(labelScreePlot)
# nv <- rep(F, length(dd))
# for (i in 1:(length(dd))) {
#   if (labelVec[i] %in% labelScreePlot) {
#     nv[i] <- T
#   }
# }
# labelVec <- labelVec[nv]
# dd <- dd[nv]
# 
# tmp <- sapply(1:(length(dd)), function(iIter) {mean(dd[[iIter]])})
# hist(tmp[tmp<0.05])
# 
# 
# labelVec <- labelVec[indSample]
# dd <- dd[indSample]
# 
# A1 <- add(dd[tauStar==1])/sum(tauStar==1)
# A2 <- add(dd[tauStar==2])/sum(tauStar==2)
# A3 <- add(dd[tauStar==3])/sum(tauStar==3)
# 
# mean(A1)
# mean(A2)
# mean(A3)
# 
# 
# 
# tauMLE
# mean(muListMLE[[1]])
# mean(muListMLE[[2]])
# mean(muListMLE[[3]])
# 
# 
# indSample[15]
# mean(dd[[15]])
# 
# labelVec[(1:100)[sapply(1:100, function(iIter) {mean(dd[[iIter]])}) > 0.1]]
# (sapply(1:100, function(iIter) {mean(dd[[iIter]])}))[(1:100)[sapply(1:100, function(iIter) {mean(dd[[iIter]])}) > 0.1]]





