rm(list = ls())
# setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
setwd("E:/GitHub/RobustASE/Code/R")

# # load("../../Data/gg-v2-14_02_01-deltat50-th0.Rbin")
# # load("../../Data/cov-v2-14_02_01-deltat50-th0.RBin")
# load("../../Data/cov-v2-14_02_01-deltat4-th0.RBin")
# load("../../Data/tlab.Rbin")
# fileName <- paste("../../Data/fish_4.RData")
# dd <- lapply(1:(length(dd)), function(i) {abs(dd[[i]])})
# save(dd, tlab, file=fileName)

# Consider the odd id of graphs to make them independent
# dd <- dd[seq(1, length(dd), 2)]

ptm <- proc.time()

# fileName <- paste("../../Data/fish.RData")
fileName <- paste("../../Data/fish_4.RData")
load(fileName)

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


MLEList <- vector("list", M)
MLqEList <- vector("list", M)
MLEASEList <- vector("list", M)
MLqEASEList <- vector("list", M)

source("getElbows.R")
require(irlba)
dMax <- floor(n/2)
elbMat <- matrix(0, M, 3)
elbMat_q <- matrix(0, M, 3)
for (i in 1:M) {
  print(labelScreePlot[i])
  print(proc.time() - ptm)
  nv <- sapply(1:length(labelVec), function(iter) {labelVec[iter] == labelScreePlot[i]})
  AList <- dd[nv]
  MLEList[[i]] <- add(AList)/length(AList)
  
  for (j in 1:length(AList)) {
    AList[[j]][upper.tri(AList[[j]], T)] <- 0
  }
  ATensor <- array(unlist(AList), dim = c(n, n, length(AList)))
  MLqEList[[i]] <- apply(ATensor, c(1, 2), mlqe_exp_solver, q)
  MLqEList[[i]] <- MLqEList[[i]] + t(MLqEList[[i]])
  
  vecs <- irlba(MLEList[[i]], dMax, dMax)$d
  elb <- getElbows(vecs, plot=F)
  elbMat[i,] <- elb
  
  vecs <- irlba(MLqEList[[i]], dMax, dMax)$d
  elb <- getElbows(vecs, plot=F)
  elbMat_q[i,] <- elb
}

# plot(1:20, type="n", ylim=c(0,10), xlab="embedding dimension", ylab="eigen value")
# for (i in 1:M) {
#   points(evalMat[i,], type="l", col=(classScreePlot[i]+4))
#   points(elbMat[i,], evalMat[i,elb], pch=19, col=2:4, cex=1)
# }

dZG <- ceiling(median(elbMat[,3]))
dZG_q <- ceiling(median(elbMat_q[,3]))


for (i in 1:M) {
  MLEDiagAug <- diag_aug(MLEList[[i]])
  AASE = ase(MLEDiagAug, dZG, isSVD)
  if (dZG == 1) {
    AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
  } else {
    AHat <- AASE[[3]][, 1:dZG]%*%diag(AASE[[1]][1:dZG])%*%t(AASE[[2]][ ,1:dZG])
  }
  MLEASEList[[i]] <- regularize(AHat)
  
  MLqEDiagAug <- diag_aug(MLqEList[[i]])
  AASE = ase(MLqEDiagAug, dZG_q, isSVD)
  if (dZG_q == 1) {
    AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
  } else {
    AHat <- AASE[[3]][, 1:dZG_q]%*%diag(AASE[[1]][1:dZG_q])%*%t(AASE[[2]][ ,1:dZG_q])
  }
  MLqEASEList[[i]] <- regularize(AHat)
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

if (isSVD) {
  fileName = paste("../../Result/result_", dataName, "_dissimilarity", "_q_", q, "_svd.RData", sep="")
} else {
  fileName = paste("../../Result/result_", dataName, "_dissimilarity", "_q_", q, "_eig.RData", sep="")
}

save(MLEList, MLqEList, MLEASEList, MLqEASEList, dZG, dZG_q, n, M, labelScreePlot,
     disMatrixMLE, disMatrixMLqE, disMatrixMLEASE, disMatrixMLqEASE, file=fileName)

print(proc.time() - ptm)