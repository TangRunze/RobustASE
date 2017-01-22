rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
source("mylibrary.R")
library(lattice)


###### Parameter Setting ######
dataNameVec <- c("migrain", "ndmg", "m2g")
isWeighted <- 1

if (isWeighted) {
  strWeighted <- "Weighted"
} else {
  strWeighted <- "Unweighted"
}

A_MLE <- list()
valMax <- 0
for (iData in 1:length(dataNameVec)) {
  dataName <- dataNameVec[iData]
  inputList <- ReadData(dataName, weighted = isWeighted)
  AList <- inputList[[1]]
  n <- inputList[[2]]
  M <- inputList[[3]]
  A_MLE[[iData]] <- add(AList)/M
  valMax <- max(valMax, max(A_MLE[[iData]]))
}

myAt <- seq(0, valMax, length.out=20)
myCkey <- list(at = myAt)
new.palette <- colorRampPalette(c("white","black"),space="rgb")

for (iData in 1:length(dataNameVec)) {
  dataName <- dataNameVec[iData]
  pdf(paste0("../../Result/Mean_Graph_", dataName, ".pdf"),
      family="Times", width=4, height=4.4)
  g <- levelplot(as.matrix(A_MLE[[iData]][1:n,n:1]), col.regions=new.palette(20),
                 xlab=list(cex=0), ylab=list(cex=0),
                 scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
                 main=list(label=paste0(dataName)),
                 at=myAt, colorkey=myCkey, lwd=0)
  print(g)
  dev.off()
}


xVec <- sapply(1:M, function(iIter) {AList[[iIter]]})
xVec <- c(xVec)
# hist(xVec)
hist(xVec[xVec < 10])




