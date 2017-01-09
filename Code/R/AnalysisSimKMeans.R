rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
# setwd("E:/GitHub/RobustASE/Code/R")

q <- 0.9
iReplicate <- 1
isSVD <- 0

if (isSVD) {
  fileName = paste("../../Result/result_sim_kmeans", "_q_", q,
                   "_replicate_", iReplicate, "_svd.RData", sep="")
} else {
  fileName = paste("../../Result/result_sim_kmeans", "_q_", q,
                   "_replicate_", iReplicate, "_eig.RData", sep="")
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
