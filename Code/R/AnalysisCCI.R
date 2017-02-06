rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
source("mylibrary.R")
library(lattice)


###### Parameter Setting ######
dataNameVec <- c("migrain", "ndmg", "m2g", "ndmg2")
isWeighted <- 1

if (isWeighted) {
  strWeighted <- "Weighted"
} else {
  strWeighted <- "Unweighted"
}

nMin <- 1
# nMax <- 10
nMax <- 70


A_MLE <- list()
A_sd <- list()
valMax <- 0
valMaxSd <- 0
for (iData in 1:length(dataNameVec)) {
  dataName <- dataNameVec[iData]
  inputList <- ReadData(dataName, weighted = isWeighted)
  AList <- inputList[[1]]
  n <- inputList[[2]]
  M <- inputList[[3]]
  A_MLE[[iData]] <- add(AList)/M
  A_MLE[[iData]] <- A_MLE[[iData]][nMin:nMax, nMin:nMax]
  valMax <- max(valMax, max(A_MLE[[iData]]))
  
  A_sd[[iData]] <- matrix(0, ncol=n, nrow=n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      A_sd[[iData]][i, j] <- sd(sapply(1:M, function(iIter) {AList[[iIter]][i, j]}))
    }
  }
  A_sd[[iData]] <- A_sd[[iData]] + t(A_sd[[iData]])
  valMaxSd <- max(valMaxSd, max(A_sd[[iData]]))
}

# m2g < ndmg2 < ndmg < migrain


myAt <- seq(0, valMax, length.out=20)
myCkey <- list(at = myAt)
new.palette <- colorRampPalette(c("white","black"),space="rgb")

for (iData in 1:length(dataNameVec)) {
  dataName <- dataNameVec[iData]
  pdf(paste0("../../Result/Mean_Graph_", dataName, ".pdf"),
      family="Times", width=4, height=4.4)
  g <- levelplot(as.matrix(A_MLE[[iData]][nMin:nMax,nMax:nMin]), col.regions=new.palette(20),
                 xlab=list(cex=0), ylab=list(cex=0),
                 scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
                 main=list(label=paste0(dataName)),
                 at=myAt, colorkey=myCkey, lwd=0)
  print(g)
  dev.off()
}


xVec <- sapply(1:M, function(iIter) {AList[[iIter]]})
xVec <- c(xVec)
hist(xVec)
# hist(xVec[xVec < 10])





myAt <- seq(0, valMaxSd, length.out=20)
myCkey <- list(at = myAt)
new.palette <- colorRampPalette(c("white","black"),space="rgb")

for (iData in 1:length(dataNameVec)) {
  dataName <- dataNameVec[iData]
  pdf(paste0("../../Result/sd_Graph_", dataName, ".pdf"),
      family="Times", width=4, height=4.4)
  g <- levelplot(as.matrix(A_sd[[iData]][nMin:nMax,nMax:nMin]), col.regions=new.palette(20),
                 xlab=list(cex=0), ylab=list(cex=0),
                 scales=list(x=list(draw=FALSE),y=list(draw=FALSE)),
                 main=list(label=paste0(dataName)),
                 at=myAt, colorkey=myCkey, lwd=0)
  print(g)
  dev.off()
}

pp_scree <- list()
eigenResult <- list()
require(Matrix)
# isSVD <- 0
label_y_ub <- 250000
label_y_lb <- -50000
for (iData in 1:length(dataNameVec)) {
  dataName <- dataNameVec[[iData]]
  Abar <- A_MLE[[iData]]
  D0 <- Diagonal(n, x = rowSums(Abar)/(n-1))
  eigenResult[[iData]] <- eigen(Abar + D0)$values
  df <- data.frame(eval=eigenResult[[iData]], k=1:n)
  
  pp_scree[[iData]] <- ggplot(df,aes(x=k,y=eval))+
    geom_line()+
    # scale_linetype_manual(name="",values=c("longdash","dotted","dotdash"))+
    scale_y_continuous(limits = c(label_y_lb, label_y_ub)) + 
    xlab("order in algebraic") + ylab("eigenvalue")+
    theme(panel.grid.major = element_line(colour="grey95"),
          panel.grid.minor = element_blank())+
    theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
    theme(legend.position="none")+
    ggtitle(dataName)
  
  ggsave(paste0("../../Result/screeplot_", dataName, ".pdf"),
         pp_scree[[iData]]+theme(text=element_text(size=10,family="Times")),
         # pp_scree[[iData]]+theme(text=element_text(size=10,family="CM Roman")),
         width=2, height=2)
}




for (iData in 1:length(dataNameVec)) {
  dataName <- dataNameVec[iData]
  print(paste0(dataName, " ", mean(A_MLE[[iData]])))
}


for (iData in 1:length(dataNameVec)) {
  dataName <- dataNameVec[iData]
  print(paste0(dataName, " ", mean(A_sd[[iData]])))
}

