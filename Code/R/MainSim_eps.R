# Simulation for LLG
rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")

source("function_collection.R")

# ###### Fix m ######
m <- 10
n <- 100
epsVec <- (0:20)/100
isSVD <- 0

nIter <- 100
nCores <- 2

iModel <- 3

if (iModel == 1) {
  B = matrix(c(4.2, 2, 2, 7), ncol = 2)
  CB = matrix(c(20, 18, 18, 25), ncol = 2)
  rho = c(0.5, 0.5)
  K = length(rho)
  eps = 0
  q = 0.8
  d = 2
} else if (iModel == 2) {
  B <- matrix(c(4, 2, 1, 2, 7, 3, 1, 3, 5), ncol = 3)
  CB <- matrix(c(10, 8, 2, 8, 15, 7, 2, 7, 12), ncol = 3)
  rho <- c(0.5, 0.3, 0.2)
  K <- length(rho)
  q <- 0.9
  d <- 3
} else if (iModel == 3) {
  B = matrix(c(4.2, 2, 2, 7), ncol = 2)
  CB = matrix(c(20, 18, 18, 25), ncol = 2)/3
  rho = c(0.5, 0.5)
  K = length(rho)
  eps = 0
  q = 0.8
  d = 2
}

require(parallel)

for (eps in epsVec) {
  print(eps)
  
  if (isSVD) {
    fileName <- paste("../../Result/result_sim_", iModel, "_d_", d, "_n_", n, "_m_", m,
                     "_eps_", eps, "_q_", q, "_svd.RData", sep="")
  } else {
    fileName <- paste("../../Result/result_sim_", iModel, "_d_", d, "_n_", n, "_m_", m,
                     "_eps_", eps, "_q_", q, "_eig.RData", sep="")
  }
  
  if (file.exists(fileName) == F) {
    tau <- rep(1:K, n*rho)
    P = B[tau, tau]
    C = CB[tau, tau]
    diag(P) = 0
    diag(C) = 0
    
    error_P_hat = matrix(0, 2, nIter)
    error_A_bar = matrix(0, 2, nIter)
    
    out <- mclapply(1:nIter, function(x) sim_all(m, n, tau, B, CB, eps, q, d, isSVD), 
                    mc.cores=nCores)
    out = array(unlist(out), dim = c(4, nIter))
    
    error_A_bar = out[1,]
    error_A_bar_ase = out[2,]
    error_P_hat = out[3,]
    error_P_hat_ase = out[4,]
    save(error_A_bar, error_A_bar_ase, error_P_hat, error_P_hat_ase, 
         n, m, rho, tau, B, CB, eps, q, d, nIter, file=fileName)
  }
}


mean(error_A_bar)
mean(error_A_bar_ase)
mean(error_P_hat)
mean(error_P_hat_ase)
