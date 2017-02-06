rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
source("mylibrary.R")
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

###### Parameter Setting ######
iModel <- 7
m <- 10
n <- 70

dMin <- 0
dMax <- n

q <- 0.9
isSVD <- 0
nD <- n

if (isSVD) {
  strSVD <- "SVD"
} else {
  strSVD <- "EIG"
}
strFlipVec <- c("Normal", "Flip")

for (isFlip in 0:1) {
  
  fileName1 <- paste0("../../Result/result_sim_twodata_1_model_", iModel,
                      "_n_", n, "_m_", m, "_q_", q, "_", strSVD, ".RData")
  fileName2 <- paste0("../../Result/result_sim_twodata_2_model_", iModel,
                      "_n_", n, "_m_", m, "_q_", q, "_", strSVD, ".RData")
  
  if (isFlip) {
    fileName <- fileName1
  } else {
    fileName <- fileName2
  }
  
  load(fileName)
  
  if (isFlip == 0) {
    errorMLEMean <- array(0, dim = c(2, 1))
    errorMLELB <- array(0, dim = c(2, 1))
    errorMLEUB <- array(0, dim = c(2, 1))
    errorMLqEMean <- array(0, dim = c(2, 1))
    errorMLqELB <- array(0, dim = c(2, 1))
    errorMLqEUB <- array(0, dim = c(2, 1))
    errorMLEASEMean <- array(0, dim = c(2, n))
    errorMLEASELB <- array(0, dim = c(2, n))
    errorMLEASEUB <- array(0, dim = c(2, n))
    errorMLqEASEMean <- array(0, dim = c(2, n))
    errorMLqEASELB <- array(0, dim = c(2, n))
    errorMLqEASEUB <- array(0, dim = c(2, n))
  }
  
  errorMLEMean[isFlip + 1] <- rep(mean(errorMLE))
  errorMLELB[isFlip + 1] <- errorMLEMean[isFlip + 1] -
    sqrt(var(errorMLE))/sqrt(length(errorMLE))*1.96
  errorMLEUB[isFlip + 1] <- errorMLEMean[isFlip + 1] +
    sqrt(var(errorMLE))/sqrt(length(errorMLE))*1.96
  
  errorMLEASEMean[isFlip + 1, ] <- rowMeans(errorMLEASE)
  errorMLEASELB[isFlip + 1, ] <- errorMLEASEMean[isFlip + 1, ] - 
    sqrt(apply(errorMLEASE, 1, var))/sqrt(dim(errorMLEASE)[2])*1.96
  errorMLEASEUB[isFlip + 1, ] <- errorMLEASEMean[isFlip + 1, ] + 
    sqrt(apply(errorMLEASE, 1, var))/sqrt(dim(errorMLEASE)[2])*1.96
  
  errorMLqEMean[isFlip + 1] <- rep(mean(errorMLqE))
  errorMLqELB[isFlip + 1] <- errorMLqEMean[isFlip + 1] -
    sqrt(var(errorMLqE))/sqrt(length(errorMLqE))*1.96
  errorMLqEUB[isFlip + 1] <- errorMLqEMean[isFlip + 1] +
    sqrt(var(errorMLqE))/sqrt(length(errorMLqE))*1.96
  
  errorMLqEASEMean[isFlip + 1, ] <- rowMeans(errorMLqEASE)
  errorMLqEASELB[isFlip + 1, ] <- errorMLqEASEMean[isFlip + 1, ] - 
    sqrt(apply(errorMLqEASE, 1, var))/sqrt(dim(errorMLqEASE)[2])*1.96
  errorMLqEASEUB[isFlip + 1, ] <- errorMLqEASEMean[isFlip + 1, ] + 
    sqrt(apply(errorMLqEASE, 1, var))/sqrt(dim(errorMLqEASE)[2])*1.96
}

dZGMLEMean <- rep(0, 2)
dZGMLELB <- rep(0, 2)
dZGMLEUB <- rep(0, 2)
dUSVTMLEMean <- rep(0, 2)
dUSVTMLELB <- rep(0, 2)
dUSVTMLEUB <- rep(0, 2)
dZGMLqEMean <- rep(0, 2)
dZGMLqELB <- rep(0, 2)
dZGMLqEUB <- rep(0, 2)
dUSVTMLqEMean <- rep(0, 2)
dUSVTMLqELB <- rep(0, 2)
dUSVTMLqEUB <- rep(0, 2)
errorMLEUSVT <- rep(0, 2)
errorMLEZG <- rep(0, 2)
errorMLqEUSVT <- rep(0, 2)
errorMLqEZG <- rep(0, 2)
for (isFlip in 0:1) {
  dZGMLEMean[isFlip + 1] <- mean(dimZGMLE)
  dZGMLELB[isFlip + 1] <- mean(dimZGMLE) -
    sqrt(var(dimZGMLE))/sqrt(length(dimZGMLE))*1.96
  dZGMLEUB[isFlip + 1] <- mean(dimZGMLE) +
    sqrt(var(dimZGMLE))/sqrt(length(dimZGMLE))*1.96
  
  dUSVTMLEMean[isFlip + 1] <- mean(dimUSVTMLE)
  dUSVTMLELB[isFlip + 1] <- mean(dimUSVTMLE) -
    sqrt(var(dimUSVTMLE))/sqrt(length(dimUSVTMLE))*1.96
  dUSVTMLEUB[isFlip + 1] <- mean(dimUSVTMLE) +
    sqrt(var(dimUSVTMLE))/sqrt(length(dimUSVTMLE))*1.96
  
  dZGMLqEMean[isFlip + 1] <- mean(dimZGMLqE)
  dZGMLqELB[isFlip + 1] <- mean(dimZGMLqE) -
    sqrt(var(dimZGMLqE))/sqrt(length(dimZGMLqE))*1.96
  dZGMLqEUB[isFlip + 1] <- mean(dimZGMLqE) +
    sqrt(var(dimZGMLqE))/sqrt(length(dimZGMLqE))*1.96
  
  dUSVTMLqEMean[isFlip + 1] <- mean(dimUSVTMLqE)
  dUSVTMLqELB[isFlip + 1] <- mean(dimUSVTMLqE) -
    sqrt(var(dimUSVTMLqE))/sqrt(length(dimUSVTMLqE))*1.96
  dUSVTMLqEUB[isFlip + 1] <- mean(dimUSVTMLqE) +
    sqrt(var(dimUSVTMLqE))/sqrt(length(dimUSVTMLqE))*1.96
  
  x <- dUSVTMLEMean[isFlip + 1]
  x1 <- floor(x)
  y1 <- errorMLEASEMean[isFlip + 1, x1]
  x2 <- ceiling(x)
  y2 <- errorMLEASEMean[isFlip + 1, x2]
  errorMLEUSVT[isFlip + 1] <- (y2-y1)/(x2-x1)*(x-x1)+y1
  
  x <- dZGMLEMean[isFlip + 1]
  x1 <- floor(x)
  y1 <- errorMLEASEMean[isFlip + 1, x1]
  x2 <- ceiling(x)
  y2 <- errorMLEASEMean[isFlip + 1, x2]
  errorMLEZG[isFlip + 1] <- (y2-y1)/(x2-x1)*(x-x1)+y1
  
  x <- dUSVTMLqEMean[isFlip + 1]
  x1 <- floor(x)
  y1 <- errorMLqEASEMean[isFlip + 1, x1]
  x2 <- ceiling(x)
  y2 <- errorMLqEASEMean[isFlip + 1, x2]
  errorMLqEUSVT[isFlip + 1] <- (y2-y1)/(x2-x1)*(x-x1)+y1
  
  x <- dZGMLqEMean[isFlip + 1]
  x1 <- floor(x)
  y1 <- errorMLqEASEMean[isFlip + 1, x1]
  x2 <- ceiling(x)
  y2 <- errorMLqEASEMean[isFlip + 1, x2]
  errorMLqEZG[isFlip + 1] <- (y2-y1)/(x2-x1)*(x-x1)+y1
}

# errorByDimDf <- rbind(
#   data.frame(mse = errorMLEMean, lci = errorMLELB, uci = errorMLEUB,
#              which = "ABar", m = mVec, d = 1),
#   data.frame(mse = errorMLEMean, lci = errorMLELB, uci = errorMLEUB,
#              which = "ABar", m = mVec, d = n),
#   data.frame(mse = c(errorMLEASEMean), lci = c(errorMLEASELB), uci = c(errorMLEASEUB),
#              which = "ABarASE", m = rep(mVec,n), d = rep(1:n, each = length(mVec))),
#   data.frame(mse = errorMLqEMean, lci = errorMLqELB, uci = errorMLqEUB,
#              which = "PHat", m = mVec, d = 1),
#   data.frame(mse = errorMLqEMean, lci = errorMLqELB, uci = errorMLqEUB,
#              which = "PHat", m = mVec, d = n),
#   data.frame(mse = c(errorMLqEASEMean), lci = c(errorMLqEASELB), uci = c(errorMLqEASEUB),
#              which = "PHatASE", m = rep(mVec,n), d = rep(1:n, each = length(mVec)))) %>%
#   mutate(m=factor(paste0("m=",m), sapply(mVec, function(m) {paste0("m=", m)})))

errorByDimDf <- rbind(
  data.frame(mse = errorMLEMean, lci = errorMLEMean, uci = errorMLEMean,
             which = "ABar", flip = 0:1, d = max(dMin, 1)),
  data.frame(mse = errorMLEMean, lci = errorMLEMean, uci = errorMLEMean,
             which = "ABar", flip = 0:1, d = min(dMax, n)),
  data.frame(mse = c(errorMLEASEMean), lci = c(errorMLEASEMean), uci = c(errorMLEASEMean),
             which = "ABarASE", flip = rep(0:1,n), d = rep(1:n, each = 2)),
  data.frame(mse = errorMLqEMean, lci = errorMLqEMean, uci = errorMLqEMean,
             which = "PHat", flip = 0:1, d = max(dMin, 1)),
  data.frame(mse = errorMLqEMean, lci = errorMLqEMean, uci = errorMLqEMean,
             which = "PHat", flip = 0:1, d = min(dMax, n)),
  data.frame(mse = c(errorMLqEASEMean), lci = c(errorMLqEASEMean), uci = c(errorMLqEASEMean),
             which = "PHatASE", flip = rep(0:1,n), d = rep(1:n, each = 2))) %>%
  mutate(flip=factor(strFlipVec[flip + 1], strFlipVec))

nv <- ((errorByDimDf$d >= dMin) & (errorByDimDf$d <= dMax))
errorByDimDf <- errorByDimDf[nv, ]

# dimSelectionDf <- rbind(
#   data.frame(mse = errorMLEZG, lci = errorMLEZG, uci = errorMLEZG,
#              which = "ABar ZG 3rd", m = mVec, d = dZGMLEMean),
#   data.frame(mse = errorMLEUSVT, lci = errorMLEUSVT, uci = errorMLEUSVT,
#              which = "ABar USVT c=0.7", m = mVec, d = dUSVTMLEMean),
#   data.frame(mse = errorMLqEZG, lci = errorMLqEZG, uci = errorMLqEZG,
#              which = "PHat ZG 3rd", m = mVec, d = dZGMLqEMean),
#   data.frame(mse = errorMLqEUSVT, lci = errorMLqEUSVT, uci = errorMLqEUSVT,
#              which = "PHat USVT c=0.7", m = mVec, d = dUSVTMLqEMean)) %>%
#   mutate(m=factor(paste0("m=",m), sapply(mVec, function(m) {paste0("m=", m)})))




# nv <- ((dimSelectionDf$d >= dMin) & (dimSelectionDf$d <= dMax))
# dimSelectionDf <- dimSelectionDf[nv, ]

label_y <- with(errorByDimDf, .75*max(mse)+.25*min(mse))

lSize = .8
legendSize = 1.5

gg <- ggplot(errorByDimDf, aes(x = d, y = mse, linetype = factor(which), shape = factor(which))) +
  facet_wrap(~flip) +
  # geom_point(data = dimSelectionDf, size = 3) +
  # scale_linetype_manual(name = "", values = c(1, 0, 0, 2, 3, 0, 0, 4)) +
  # scale_shape_manual(name = "", values = c(-1, 0, 0, -1, -1, 0, 0, -1)) +
  scale_linetype_manual(name = "", values = c(1, 2, 3, 4),
                        labels = c("MLE", "MLE_ASE", "MLqE", "MLqE_ASE")) +
  scale_shape_manual(name = "", values = c(-1, -1, -1, -1),
                     labels = c("MLE", "MLE_ASE", "MLqE", "MLqE_ASE")) +
  geom_line(alpha = 1, size = lSize) +
  geom_linerange(aes(ymin = lci, ymax = uci), alpha = .5, size = 1) +
  xlab("dimension")+ylab("MSE")+
  theme(strip.text.x = element_text(size=20,face="bold"))+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"))+
  theme(panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
  theme(legend.text=element_text(size=20,face="bold"))+
  theme(legend.position="bottom")+
  ggtitle(paste0("Model", iModel, " m=", m, " n=", n, " q=", q, " ", strSVD))+
  theme(legend.key.size=unit(legendSize,"line"))+
  theme(plot.title=element_text(lineheight=.8,size=20,face="bold"))

ggsave(paste0("../../Result/sim_two_data_model_", iModel,
              "_n_", n, "_m_", m, "_q_", q, "_", strSVD, ".pdf"),
       plot=gg+theme(text=element_text(size=10,family="Times")),
       width=8,height=6)

print(gg)


# require(Matrix)
# print(rankMatrix(P)[1])
# print(rankMatrix(C1)[1])
# print(rankMatrix(C2)[1])
# 
# label_y_ub <- max(c(GetEvals(P), GetEvals(C1), GetEvals(C2)))
# label_y_lb <- min(c(GetEvals(P), GetEvals(C1), GetEvals(C2)))
# 
# print(PlotEvals(GetEvals(P), label_y_lb, label_y_ub))
# print(PlotEvals(GetEvals(C1), label_y_lb, label_y_ub))
# print(PlotEvals(GetEvals(C2), label_y_lb, label_y_ub))
# 
# print(PlotEvals(GetEvals((1-eps1)*P+eps1*C1), label_y_lb, label_y_ub))
# print(PlotEvals(GetEvals((1-eps2)*P+eps1*C2), label_y_lb, label_y_ub))

