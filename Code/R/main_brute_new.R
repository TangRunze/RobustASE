
# args <- commandArgs(trailingOnly = TRUE)

setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
# setwd("/cis/home/rtang/RobustASE/Code/R")

# m = as.numeric(args[1])
# m = 1
# isSVD = 1

mVec = c(1, 2, 5)
q = 0.8
nIter = 50
nCores = 2
dataName = "desikan"
# dataName = "CPAC200"
source("function_collection.R")

for (m in mVec) {
  for (isSVD in 0) {
    print(c(m, isSVD))

    #     tmpList = read_data_weight(dataName, DA=F)
    tmpList = read_data_weight(dataName, DA=F)
    A_all = tmpList[[1]]
    n = tmpList[[2]]
    M = tmpList[[3]]
    rm(tmpList)
    
    #     for (i in 1:length(A_all)) {
    #       A_all[[i]] = 1*(A_all[[i]]>0)
    #     }
    
    dVec = 1:n
    nD = length(dVec)
    
    A_sum = add(A_all)
    
    error_P_hat = matrix(0, nD, nIter)
    error_A_bar = matrix(0, nD, nIter)
    error_P_hat_ase = matrix(0, nD, nIter)
    error_A_bar_ase = matrix(0, nD, nIter)
    
    require(parallel)
    
    # ptm <- proc.time()
    # proc.time() - ptm
    
    # out <- mclapply(1:nIter, function(x) sapply(dVec, function(d) dim_brute(M, m, d, A_all, A_sum)),
    #                 mc.cores=nCores)
    # out = array(unlist(out), dim = c(2, nD, nIter))
    # error_A_bar = out[1,,]
    # error_P_hat = out[2,,]
    
    out <- mclapply(1:nIter, function(x) dim_brute2_all(M, m, dVec, A_all, A_sum, q, isSVD), 
                    mc.cores=nCores)
    out = array(unlist(out), dim = c(2*nD+6, nIter))
    
    error_A_bar = out[1,]
    error_A_bar_ase = out[1+(1:nD),]
    error_P_hat = out[nD+2,]
    error_P_hat_ase = out[nD+2+(1:nD),]
    dim_ZG_A_bar = out[2*nD+3,]
    dim_USVT_A_bar = out[2*nD+4,]
    dim_ZG_P_hat = out[2*nD+5,]
    dim_USVT_P_hat = out[2*nD+6,]
    
    error_A_bar_ZG = rep(0, length(dim_ZG_A_bar))
    error_A_bar_USVT = rep(0, length(dim_USVT_A_bar))
    error_P_hat_ZG = rep(0, length(dim_ZG_P_hat))
    error_P_hat_USVT = rep(0, length(dim_USVT_P_hat))
    for (i in 1:length(dim_ZG_A_bar)) {
      error_A_bar_ZG[i] = error_A_bar_ase[dim_ZG_A_bar[i], i]
      error_A_bar_USVT[i] = error_A_bar_ase[dim_USVT_A_bar[i], i]
      error_P_hat_ZG[i] = error_P_hat_ase[dim_ZG_P_hat[i], i]
      error_P_hat_USVT[i] = error_P_hat_ase[dim_USVT_P_hat[i], i]
    }
    
    if (isSVD) {
      fileName = paste("../../Result/result_", dataName, "_brute_", "m_", m, "_svd.RData", sep="")
    } else {
      fileName = paste("../../Result/result_", dataName, "_brute_", "m_", m, "_eig.RData", sep="")
    }
    
    save(error_A_bar, error_A_bar_ase, error_P_hat, error_P_hat_ase,
         error_A_bar_ZG, error_A_bar_USVT, error_P_hat_ZG, error_P_hat_USVT, 
         dim_ZG_A_bar, dim_USVT_A_bar, dim_ZG_P_hat, dim_USVT_P_hat,
         n, M, m, dVec, nIter, file=fileName)
    
  }
}