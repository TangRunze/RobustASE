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

ind <- 310
tmp1 <- eigen(diag_aug(AList[[ind]]))
ATmp <- (AList[[ind]] > 0)*1
tmp2 <- eigen(diag_aug(ATmp))
# plot(tmp1$values)
# points(tmp2$values, col="blue")
plot(tmp1$values/sum(abs(tmp1$values)))
points(tmp2$values/sum(abs(tmp2$values)), col="blue")

