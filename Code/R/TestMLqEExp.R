rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")

m <- 5
q <- 0.9
# dataName <- "desikan"
dataName <- "CPAC200"

isSVD <- 0

source("function_collection.R")
source("getElbows.R")
source("USVT.R")

###### Read Data ######
tmpList <- ReadDataWeighted(dataName, DA = F)
AList <- tmpList[[1]]
n <- tmpList[[2]]
M <- tmpList[[3]]
rm(tmpList)


###### Visualize Weight Distribution ######
# i <- 3
# j <- 4
# dataTmp <- sapply(1:length(AList), function(x) {AList[[x]][i,j]})
# hist(dataTmp)
# hist(log(dataTmp))
# 
# par(mfrow=c(3,2))
# iVec <- c(1, 2, 3)
# jVec <- c(2, 3, 4)
# for (ind in 1:length(iVec)) {
#   weight <- sapply(1:length(AList), function(x) {AList[[x]][iVec[ind],jVec[ind]]})
#   hist(weight, main=paste0("Histogram of weight for (", iVec[ind], ", ", jVec[ind], ")"))
#   hist(log(weight), main=paste0("Histogram of log of weight for (",
#                                 iVec[ind], ", ", jVec[ind], ")"))
# }



###### Write to Files ######
# pairVec = c(3, 4)
# outputFileName = paste0("../../Data/", dataName, "_pair_", pairVec[1], "_", pairVec[2], ".csv")
# 
# # write.csv(n, file=outputFileName)
# # write.csv(M, file=outputFileName, append=TRUE)
# dataTmp <- sapply(1:length(AList), function(x) {AList[[x]][pairVec[1],pairVec[2]]})
# # write.csv(dataTmp, file=outputFileName, row.names=FALSE, col.names=FALSE)
# write.table(dataTmp, file=outputFileName, sep=",", row.names=FALSE, col.names=FALSE)
# 
# table(dataTmp)



###### Select Element ######
pairVec = c(2, 3)
dataTmp <- sapply(1:length(AList), function(x) {AList[[x]][pairVec[1],pairVec[2]]})
hist(dataTmp)



###### Solve MLqE ######
thetaHat <- mean(dataTmp)
thetaInit <- thetaHat
if (var(dataTmp) > 0) {
  resultOpti <- optim(par = thetaInit, fn = MLqEObjExp, data = dataTmp, q = q)
#   resultOpti <- optim(par = thetaInit, fn = MLqEObjExp, data = dataTmp, q = q,
#                       method = "BFGS")
#   resultOpti <- optim(par = thetaInit, fn = MLqEObjExp, data = dataTmp, q = q,
#                       gr = MLqEGradExp, method = "BFGS")
  thetaHat <- resultOpti$par
  
  
  MLqEGradExp(theta = resultOpti$par, data = dataTmp, q = q)
}


###### Plot Result ######
xVec <- (1:(2*max(dataTmp)*100))/100
fHat <- 1/thetaHat*exp(-xVec/thetaHat)
d <- density(dataTmp)
plot(d)
lines(xVec, fHat)


