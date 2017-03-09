rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
source("mylibrary.R")
require(parallel)

###### Parameter Setting ######
m <- 10
n <- 100
qVec = (10:100)/100
isSVD <- 0

nIter <- 200
nCores <- 2

iModel <- 2

if (iModel == 1) {
  B <- matrix(c(4.2, 2, 2, 7), ncol = 2)
  CB <- matrix(c(20, 18, 18, 25), ncol = 2)
  rho <- c(0.5, 0.5)
  K <- length(rho)
  eps <- 0.1
  d <- 2
} else if (iModel == 2) {
  B <- matrix(c(4, 2, 2, 7), ncol = 2)
  CB <- matrix(c(12, 9, 9, 15), ncol = 2)
  rho <- c(0.5, 0.5)
  K <- length(rho)
  eps <- 0.1
  d <- 2
}

for (q in qVec) {
  print(eps)
  
  if (isSVD) {
    fileName <- paste("../../Result/result_sim_", iModel, "_d_", d, "_n_", n, "_m_", m,
                     "_eps_", eps, "_q_", q, "_svd.RData", sep="")
  } else {
    fileName <- paste("../../Result/result_sim_", iModel, "_d_", d, "_n_", n, "_m_", m,
                     "_eps_", eps, "_q_", q, "_eig.RData", sep="")
  }
  
  if (file.exists(fileName) == F) {
    tau <- sample(1:K, n, replace = T, prob = rho)
    P <- B[tau, tau]
    C <- CB[tau, tau]
    diag(P) <- 0
    diag(C) <- 0
    
    error_P_hat <- matrix(0, 2, nIter)
    error_A_bar <- matrix(0, 2, nIter)
    
    out <- mclapply(1:nIter, function(x) SimAll(m, n, tau, B, CB, eps, q, d, isSVD), 
                    mc.cores=nCores)
    out <- array(unlist(out), dim = c(4, nIter))
    
    error_MLE <- out[1,]
    error_MLE_ASE<- out[2,]
    error_MLqE <- out[3,]
    error_MLqE_ASE <- out[4,]
    save(error_MLE, error_MLE_ASE, error_MLqE, error_MLqE_ASE, 
         n, m, rho, tau, B, CB, eps, q, d, nIter, file=fileName)
  }
}
