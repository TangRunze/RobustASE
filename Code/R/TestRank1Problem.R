rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
# setwd("/cis/home/rtang/RobustASE/Code/R")
source("mylibrary.R")
require(parallel)
library(lattice)

###### Parameter Tunning ######
m <- 2
d <- 1
nMin <- 1
nMax <- 70

###### Parameter Setting ######
set.seed(12345)

q <- 0.9
isSVD <- 0
isWeighted <- 1

# dataName1 <- "m2g"
# dataName2 <- "migrain"

dataName1 <- "migrain"
dataName2 <- "ndmg2"

###### Read Data ######
inputList <- ReadData(dataName1, weighted = isWeighted)
AList1 <- inputList[[1]]
n <- inputList[[2]]
M <- inputList[[3]]

inputList <- ReadData(dataName2, weighted = isWeighted)
AList2 <- inputList[[1]]
rm(inputList)

###### Main Calculation ######

if (isSVD) {
  strSVD <- "SVD"
} else {
  strSVD <- "EIG"
}
if (isWeighted) {
  strWeighted <- "Weighted"
} else {
  strWeighted <- "Unweighted"
}

AList <- AList1
weighted <- isWeighted
P <- add(AList2)/M
valMax <- max(P)

sampleVec <- sample.int(M, m)
A_MLE <- add(AList[sampleVec])/m
valMax <- max(valMax, A_MLE)


myAt <- seq(0, valMax, length.out=20)
myCkey <- list(at = myAt)
new.palette <- colorRampPalette(c("white","black"),space="rgb")
g_P <- levelplot(as.matrix(P[nMin:nMax,nMax:nMin]), col.regions=new.palette(20),
                xlab=list(cex=0), ylab=list(cex=0),
                scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
                main=list(label=paste0("P = Population mean of ", dataName2)),
                at=myAt, colorkey=myCkey, lwd=0)
print(g_P)

g_MLE <- levelplot(as.matrix(A_MLE[nMin:nMax,nMax:nMin]), col.regions=new.palette(20),
                xlab=list(cex=0), ylab=list(cex=0),
                scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
                main=list(label=paste0("MLE of a sample (m = ", m, ") from ", dataName1)),
                at=myAt, colorkey=myCkey, lwd=0)
print(g_MLE)


# MLE_ASE
require(Matrix)
D0 <- Diagonal(n, x = rowSums(A_MLE)/(n-1))
P0 <- LR(A_MLE + D0, d, isSVD)
D1 <- Diagonal(n, x = diag(P0))
P1 <- LR(A_MLE + D1, d)
A_MLE_ASE <- Regularize(P1, weighted)
g_MLE_ASE <- levelplot(as.matrix(A_MLE_ASE[nMin:nMax,nMax:nMin]), col.regions=new.palette(20),
                         xlab=list(cex=0), ylab=list(cex=0),
                         scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
                         main=list(label=paste0("ASE of MLE of a sample (m = ", m, ") from ", dataName1)),
                         at=myAt, colorkey=myCkey, lwd=0)
print(g_MLE_ASE)


(norm(P - A_MLE, "F"))/sqrt(n*(n-1))
(norm(P - A_MLE_ASE, "F"))/sqrt(n*(n-1))
(norm(P, "F"))/sqrt(n*(n-1))
tmp <- Matrix(3000, nrow=n, ncol=n)
diag(tmp) <- 0
(norm(P - tmp, "F"))/sqrt(n*(n-1))


mean(P)
