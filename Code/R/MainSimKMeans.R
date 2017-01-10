rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
# setwd("E:/GitHub/RobustASE/Code/R")
# setwd("/cis/home/rtang/RobustASE/Code/R")

nCores <- 2

K <- 3
m <- 100
n <- 200
rho <- c(0.6, 0.4)
eps <- 0.1

q <- 0.9
isSVD <- 0


for (iReplicate in 1:10) {
  print(iReplicate)
  
  if (isSVD) {
    fileName = paste("../../Result/result_sim_kmeans", "_q_", q, "_n_", n, "_m_", m,
                     "_eps_", eps, "_replicate_", iReplicate, "_svd.RData", sep="")
  } else {
    fileName = paste("../../Result/result_sim_kmeans", "_q_", q, "_n_", n, "_m_", m,
                     "_eps_", eps, "_replicate_", iReplicate, "_eig.RData", sep="")
  }
  
  if (!file.exists(fileName)) {
    
    set.seed(12345+iReplicate)
    
    labelVec <- sample(K, m, replace=T)
    tau <- rep(1:2, round(n*rho))
    
    B <- list()
    C <- list()
    # B[[1]] <- matrix(c(4.2, 2, 2, 7), ncol = 2)
    # C[[1]] <- matrix(c(10, 8, 8, 15), ncol = 2)
    # B[[2]] <- matrix(c(8, 3, 3, 4), ncol = 2)
    # C[[2]] <- matrix(c(14, 5, 5, 7), ncol = 2)
    # B[[3]] <- matrix(c(2, 1, 1, 5), ncol = 2)
    # C[[3]] <- matrix(c(4, 5, 5, 9), ncol = 2)
    B[[1]] <- matrix(c(4.2, 2, 2, 7), ncol = 2)
    C[[1]] <- matrix(c(20, 16, 16, 25), ncol = 2)
    B[[2]] <- matrix(c(8, 3, 3, 4), ncol = 2)
    C[[2]] <- matrix(c(14, 10, 10, 14), ncol = 2)
    B[[3]] <- matrix(c(2, 1, 1, 5), ncol = 2)
    C[[3]] <- matrix(c(7, 5, 5, 12), ncol = 2)
    
    P <- list()
    PC <- list()
    P[[1]] <- B[[1]][tau, tau]
    diag(P[[1]]) <- 0
    P[[2]] <- B[[2]][tau, tau]
    diag(P[[2]]) <- 0
    P[[3]] <- B[[3]][tau, tau]
    diag(P[[3]]) <- 0
    PC[[1]] <- C[[1]][tau, tau]
    diag(C[[1]]) <- 0
    PC[[2]] <- C[[2]][tau, tau]
    diag(C[[2]]) <- 0
    PC[[3]] <- C[[3]][tau, tau]
    diag(C[[3]]) <- 0
    
    source("function_collection.R")
    require(parallel)
    require(igraph)
    
    A_all <- list()
    ind <- lower.tri(matrix(rep(0, n^2), ncol=n), 1)
    for (iGraph in 1:m) {
      contamVec <- (runif(n^2) > eps)
      A <- contamVec*matrix(rexp(n^2, 1/P[[labelVec[iGraph]]]), ncol=n) +
        (1-contamVec)*matrix(rexp(n^2, 1/PC[[labelVec[iGraph]]]), ncol=n)
      A[ind] <- 0
      A <- A + t(A)
      A_all[[iGraph]] <- A
    }
    
    result <- ExpSimKMeans(A_all, K, q, nCores, isSVD)
    tauMLE <- result[[1]]
    tauMLqE <- result[[2]]
    tauMLEASE <- result[[3]]
    tauMLqEASE <- result[[4]]
    muListMLE <- result[[5]]
    muListMLqE <- result[[6]]
    muListMLEASE <- result[[7]]
    muListMLqEASE <- result[[8]]
    
    tauStar <- labelVec
    
    save(B, C, labelVec, tau, K, m, n, rho, eps, q,
         tauMLE, tauMLqE, tauMLEASE, tauMLqEASE, tauStar,
         muListMLE, muListMLqE, muListMLEASE, muListMLqEASE, file=fileName)
  }
}