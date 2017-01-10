rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
# setwd("E:/GitHub/RobustASE/Code/R")

m <- 100
n <- 200
eps <- 0.2
q <- 0.9
iReplicate <- 1
isSVD <- 0

if (isSVD) {
  fileName = paste("../../Result/result_sim_kmeans", "_q_", q, "_n_", n, "_m_", m,
                   "_eps_", eps, "_replicate_", iReplicate, "_svd.RData", sep="")
} else {
  fileName = paste("../../Result/result_sim_kmeans", "_q_", q, "_n_", n, "_m_", m,
                   "_eps_", eps, "_replicate_", iReplicate, "_eig.RData", sep="")
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
tauStar
