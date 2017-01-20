rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
source("mylibrary.R")
# library(lattice)


###### Parameter Setting ######
# dataName <- "JHU"
dataName <- "desikan"
# dataName <- "CPAC200"
# dataName <- "Talairach"

isSVD <- 0
isWeighted <- 1

###### Read Data ######
# New Dataset
inputList <- ReadData(dataName, weighted = isWeighted, newGraph = T, DA = F)
AListNew <- inputList[[1]]
n <- inputList[[2]]
M <- inputList[[3]]
# Old Dataset
inputList <- ReadData(dataName, weighted = isWeighted, newGraph = F, DA = F)
AListOld <- inputList[[1]]
rm(inputList)

A_MLE_New <- add(AListNew)/M
# image(Matrix(A_MLE_New), lwd=0, colorkey=TRUE)
A_MLE_Old <- add(AListOld)/M
# image(Matrix(A_MLE_Old), lwd=0, colorkey=TRUE)

# n <- 50
# A_MLE_New <- A_MLE_New[1:n, 1:n]
# A_MLE_Old <- A_MLE_Old[1:n, 1:n]

# n <- 69
# A_MLE_New <- A_MLE_New[2:70, 2:70]
# A_MLE_Old <- A_MLE_Old[1:n, 1:n]



myAt <- seq(0, max(max(A_MLE_New), max(A_MLE_Old)), length.out=20)
myCkey <- list(at=myAt)
new.palette <- colorRampPalette(c("white","black"),space="rgb")

levelplot(as.matrix(A_MLE_Old[1:n,n:1]),col.regions=new.palette(20),xlab=list(cex=0),
          ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
          main=list(label=paste0(dataName, " Old")),
          at=myAt, colorkey=myCkey, lwd=0)

levelplot(as.matrix(A_MLE_New[1:n,n:1]),col.regions=new.palette(20),xlab=list(cex=0),
          ylab=list(cex=0),scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
          main=list(label=paste0(dataName, " New")),
          at=myAt, colorkey=myCkey, lwd=0)


i <- 1
j <- 2
xVec <- sapply(1:M, function(iIter) {AListOld[[iIter]][i, j]})
xVec <- log(xVec + 1)
hist(xVec)



xVec <- sapply(1:M, function(iIter) {AListOld[[iIter]]})
xVec <- c(xVec)
xVec <- log(xVec + 1)
hist(xVec[xVec < 5])




