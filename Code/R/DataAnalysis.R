rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")

m <- 2
q <- 0.9
d <- 0
nIter <- 200
nCores <- 2
dataName <- "CPAC200"
isSVD <- 0

source("function_collection.R")

tmpList <- ReadDataWeighted(dataName, DA = F)
AList <- tmpList[[1]]
n <- tmpList[[2]]
M <- tmpList[[3]]
rm(tmpList)

# distPair <- matrix(rep(0, M*M), ncol=M)
# for (i in 1:(M - 1)) {
#   for (j in (i + 1):M) {
#     distPair[i, j] <- norm(AList[[i]] - AList[[j]], "F")
#     distPair[j, i] <- distPair[i, j]
#   }
# }

# minDist <- Inf
# for (i in 1:(M - 2)) {
#   print(i)
#   for (j in (i + 1):(M - 1)) {
#     for (k in (j + 1):M) {
#       tmpDist <- distPair[i, j] + distPair[j, k] + distPair[k, i]
#       if (minDist > tmpDist) {
#         minDist = tmpDist
#         i0 <- i
#         j0 <- j
#         k0 <- k
#       }
#     }
#   }
# }

load("distpairmatrix.RData")

i0 <- 96
j0 <- 299
k0 <- 300


# orderVec <- order(distPair)
# x <- rep(0, 3)
# y <- x
# for (i in 1:3) {
#   posTmp <- orderVec[M + i*2]
#   x[i] <- floor((posTmp - 1)/M) + 1
#   y[i] <- posTmp - (x[i] - 1)*M
# }

sampleVec <- c(i0, j0)
sampleVecComplement <- (1:M)[is.na(pmatch(1:M, sampleVec))]
ABar <- add(AList[sampleVec])/m

n <- dim(ABar)[[1]]
ATensor <- array(unlist(AList[sampleVec]), dim = c(n, n, m))
AMLqE <- apply(ATensor, c(1, 2), mlqe_exp_solver, q)

source("getElbows.R")

ABarDiagAug <- diag_aug(ABar)
if (d == 0) {
  # ZG
  nElbow <- 3
  evalVec <- ase(ABarDiagAug, ceiling(n*3/5), isSVD)[[1]]
  d <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
}

AASE <- ase(ABarDiagAug, d, isSVD)
if (d == 1) {
  AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
} else {
  AHat <- AASE[[3]][, 1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][, 1:d])
}
ABarASE <- regularize(AHat)

AMLqEDiagAug <- diag_aug(AMLqE)
if (d == 0) {
  # ZG
  nElbow <- 3
  evalVec <- ase(AMLqEDiagAug, ceiling(n*3/5), isSVD)[[1]]
  d <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
}

AASE <- ase(AMLqEDiagAug, d, isSVD)
if (d == 1) {
  AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
} else {
  AHat <- AASE[[3]][, 1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][, 1:d])
}
PHatASE <- regularize(AHat)



ATest <- AList[[k0]]
norm(ABar - ATest, "F")
norm(ABarASE - ATest, "F")
norm(AMLqE - ATest, "F")
norm(PHatASE - ATest, "F")
