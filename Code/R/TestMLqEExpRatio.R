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

source("getElbows.R")

sampleVec <- sample.int(M, m)
sampleVecComplement <- (1:M)[is.na(pmatch(1:M, sampleVec))]
ABar <- add(AList[sampleVec])/m

n <- dim(ABar)[[1]]
ATensor <- array(unlist(AList[sampleVec]), dim = c(n, n, m))
AMLqE <- apply(ATensor, c(1, 2), mlqe_exp_solver, q)

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


###### AMLqE vs ABar ######
# AMLqE < ABar
ratioList <- lapply(sampleVecComplement, function(i) {
  return(abs(AList[[i]] - AMLqE) < abs(AList[[i]] - ABar))})
ratioAMLqEABar <- add(ratioList)/(M - m)
image(Matrix(ratioAMLqEABar))
sum(ratioAMLqEABar)/n/(n-1)

# AMLqE = ABar
ratioList <- lapply(sampleVecComplement, function(i) {
  return(abs(AList[[i]] - AMLqE) == abs(AList[[i]] - ABar))})
ratioAMLqEABar <- add(ratioList)/(M - m)
image(Matrix(ratioAMLqEABar))
sum(ratioAMLqEABar)/n/(n-1)

# AMLqE <= ABar
ratioList <- lapply(sampleVecComplement, function(i) {
  return(abs(AList[[i]] - AMLqE) <= abs(AList[[i]] - ABar))})
ratioAMLqEABar <- add(ratioList)/(M - m)
image(Matrix(ratioAMLqEABar))
sum(ratioAMLqEABar)/n/(n-1)



###### ABarASE vs ABar ######
# ABarASE < ABar
ratioList <- lapply(sampleVecComplement, function(i) {
  return(abs(AList[[i]] - ABarASE) < abs(AList[[i]] - ABar))})
ratioABarASEABar <- add(ratioList)/(M - m)
image(Matrix(ratioABarASEABar))
sum(ratioABarASEABar)/n/(n-1)

# ABarASE = ABar
ratioList <- lapply(sampleVecComplement, function(i) {
  return(abs(AList[[i]] - ABarASE) == abs(AList[[i]] - ABar))})
ratioABarASEABar <- add(ratioList)/(M - m)
image(Matrix(ratioABarASEABar))
sum(ratioABarASEABar)/n/(n-1)

# ABarASE <= ABar
ratioList <- lapply(sampleVecComplement, function(i) {
  return(abs(AList[[i]] - ABarASE) <= abs(AList[[i]] - ABar))})
ratioABarASEABar <- add(ratioList)/(M - m)
image(Matrix(ratioABarASEABar))
sum(ratioABarASEABar)/n/(n-1)


###### PHatASE vs PHat ######
# PHatASE < PHat
ratioList <- lapply(sampleVecComplement, function(i) {
  return(abs(AList[[i]] - PHatASE) < abs(AList[[i]] - AMLqE))})
ratioPHatASEPHat <- add(ratioList)/(M - m)
image(Matrix(ratioPHatASEPHat))
sum(ratioPHatASEPHat)/n/(n-1)

# PHatASE = PHat
ratioList <- lapply(sampleVecComplement, function(i) {
  return(abs(AList[[i]] - PHatASE) == abs(AList[[i]] - AMLqE))})
ratioPHatASEPHat <- add(ratioList)/(M - m)
image(Matrix(ratioPHatASEPHat))
sum(ratioPHatASEPHat)/n/(n-1)

# PHatASE <= PHat
ratioList <- lapply(sampleVecComplement, function(i) {
  return(abs(AList[[i]] - PHatASE) <= abs(AList[[i]] - AMLqE))})
ratioPHatASEPHat <- add(ratioList)/(M - m)
image(Matrix(ratioPHatASEPHat))
sum(ratioPHatASEPHat)/n/(n-1)



###### PHatASE vs ABarASE ######
# PHatASE < ABarASE
ratioList <- lapply(sampleVecComplement, function(i) {
  return(abs(AList[[i]] - PHatASE) < abs(AList[[i]] - ABarASE))})
ratioPHatASEABarASE <- add(ratioList)/(M - m)
image(Matrix(ratioPHatASEABarASE))
sum(ratioPHatASEABarASE)/n/(n-1)

# PHatASE = ABarASE
ratioList <- lapply(sampleVecComplement, function(i) {
  return(abs(AList[[i]] - PHatASE) == abs(AList[[i]] - ABarASE))})
ratioPHatASEABarASE <- add(ratioList)/(M - m)
image(Matrix(ratioPHatASEABarASE))
sum(ratioPHatASEABarASE)/n/(n-1)

# PHatASE <= ABarASE
ratioList <- lapply(sampleVecComplement, function(i) {
  return(abs(AList[[i]] - PHatASE) <= abs(AList[[i]] - ABarASE))})
ratioPHatASEABarASE <- add(ratioList)/(M - m)
image(Matrix(ratioPHatASEABarASE))
sum(ratioPHatASEABarASE)/n/(n-1)




###### PHatASE vs ABar ######
# PHatASE < ABar
ratioList <- lapply(sampleVecComplement, function(i) {
  return(abs(AList[[i]] - PHatASE) < abs(AList[[i]] - ABar))})
ratioPHatASEABar <- add(ratioList)/(M - m)
image(Matrix(ratioPHatASEABar))
sum(ratioPHatASEABar)/n/(n-1)

# PHatASE = ABar
ratioList <- lapply(sampleVecComplement, function(i) {
  return(abs(AList[[i]] - PHatASE) == abs(AList[[i]] - ABar))})
ratioPHatASEABar <- add(ratioList)/(M - m)
image(Matrix(ratioPHatASEABar))
sum(ratioPHatASEABar)/n/(n-1)

# PHatASE <= ABar
ratioList <- lapply(sampleVecComplement, function(i) {
  return(abs(AList[[i]] - PHatASE) <= abs(AList[[i]] - ABar))})
ratioPHatASEABar <- add(ratioList)/(M - m)
image(Matrix(ratioPHatASEABar))
sum(ratioPHatASEABar)/n/(n-1)