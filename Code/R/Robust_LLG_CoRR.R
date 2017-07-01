rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/Dissertation/Code/R")
source("mylibrary.R")
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

###### Parameter Setting ######
mVec <- c(2, 5, 10)
# mVec <- c(2, 5)
q <- 0.9

isSVD <- 0
isWeighted <- 1

dMin <- 6
dMax <- Inf

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

dataName1 <- "m2g"
dataName2 <- "ndmg2"
for (iM in 1:length(mVec)) {
  m <- mVec[iM]
  fileName = paste0("../../Result/result_brute_", dataName1, "_",
                    dataName2, "_", strWeighted,
                    "_", "m_", m, "_q_", q, "_", strSVD, ".RData")
  load(fileName)
  
  if (iM == 1) {
    errorMLEMean <- array(0, dim = c(length(mVec), 1))
    errorMLELB <- array(0, dim = c(length(mVec), 1))
    errorMLEUB <- array(0, dim = c(length(mVec), 1))
    errorMLqEMean <- array(0, dim = c(length(mVec), 1))
    errorMLqELB <- array(0, dim = c(length(mVec), 1))
    errorMLqEUB <- array(0, dim = c(length(mVec), 1))
    errorMLEASEMean <- array(0, dim = c(length(mVec), n))
    errorMLEASELB <- array(0, dim = c(length(mVec), n))
    errorMLEASEUB <- array(0, dim = c(length(mVec), n))
    errorMLqEASEMean <- array(0, dim = c(length(mVec), n))
    errorMLqEASELB <- array(0, dim = c(length(mVec), n))
    errorMLqEASEUB <- array(0, dim = c(length(mVec), n))
  }
  
  errorMLEMean[iM] <- rep(mean(errorMLE))
  errorMLELB[iM] <- errorMLEMean[iM] -
    sqrt(var(errorMLE))/sqrt(length(errorMLE))*1.96
  errorMLEUB[iM] <- errorMLEMean[iM] +
    sqrt(var(errorMLE))/sqrt(length(errorMLE))*1.96
  
  errorMLEASEMean[iM, ] <- rowMeans(errorMLEASE)
  errorMLEASELB[iM, ] <- errorMLEASEMean[iM, ] - 
    sqrt(apply(errorMLEASE, 1, var))/sqrt(dim(errorMLEASE)[2])*1.96
  errorMLEASEUB[iM, ] <- errorMLEASEMean[iM, ] + 
    sqrt(apply(errorMLEASE, 1, var))/sqrt(dim(errorMLEASE)[2])*1.96
  
  errorMLqEMean[iM] <- rep(mean(errorMLqE))
  errorMLqELB[iM] <- errorMLqEMean[iM] -
    sqrt(var(errorMLqE))/sqrt(length(errorMLqE))*1.96
  errorMLqEUB[iM] <- errorMLqEMean[iM] +
    sqrt(var(errorMLqE))/sqrt(length(errorMLqE))*1.96
  
  errorMLqEASEMean[iM, ] <- rowMeans(errorMLqEASE)
  errorMLqEASELB[iM, ] <- errorMLqEASEMean[iM, ] - 
    sqrt(apply(errorMLqEASE, 1, var))/sqrt(dim(errorMLqEASE)[2])*1.96
  errorMLqEASEUB[iM, ] <- errorMLqEASEMean[iM, ] + 
    sqrt(apply(errorMLqEASE, 1, var))/sqrt(dim(errorMLqEASE)[2])*1.96
}

dZGMLEMean <- rep(0, length(mVec))
dZGMLELB <- rep(0, length(mVec))
dZGMLEUB <- rep(0, length(mVec))
dUSVTMLEMean <- rep(0, length(mVec))
dUSVTMLELB <- rep(0, length(mVec))
dUSVTMLEUB <- rep(0, length(mVec))
dZGMLqEMean <- rep(0, length(mVec))
dZGMLqELB <- rep(0, length(mVec))
dZGMLqEUB <- rep(0, length(mVec))
dUSVTMLqEMean <- rep(0, length(mVec))
dUSVTMLqELB <- rep(0, length(mVec))
dUSVTMLqEUB <- rep(0, length(mVec))
errorMLEUSVT <- rep(0, length(mVec))
errorMLEZG <- rep(0, length(mVec))
errorMLqEUSVT <- rep(0, length(mVec))
errorMLqEZG <- rep(0, length(mVec))
for (iM in 1:length(mVec)) {
  dZGMLEMean[iM] <- mean(dimZGMLE)
  dZGMLELB[iM] <- mean(dimZGMLE) -
    sqrt(var(dimZGMLE))/sqrt(length(dimZGMLE))*1.96
  dZGMLEUB[iM] <- mean(dimZGMLE) +
    sqrt(var(dimZGMLE))/sqrt(length(dimZGMLE))*1.96
  
  dUSVTMLEMean[iM] <- mean(dimUSVTMLE)
  dUSVTMLELB[iM] <- mean(dimUSVTMLE) -
    sqrt(var(dimUSVTMLE))/sqrt(length(dimUSVTMLE))*1.96
  dUSVTMLEUB[iM] <- mean(dimUSVTMLE) +
    sqrt(var(dimUSVTMLE))/sqrt(length(dimUSVTMLE))*1.96
  
  dZGMLqEMean[iM] <- mean(dimZGMLqE)
  dZGMLqELB[iM] <- mean(dimZGMLqE) -
    sqrt(var(dimZGMLqE))/sqrt(length(dimZGMLqE))*1.96
  dZGMLqEUB[iM] <- mean(dimZGMLqE) +
    sqrt(var(dimZGMLqE))/sqrt(length(dimZGMLqE))*1.96
  
  dUSVTMLqEMean[iM] <- mean(dimUSVTMLqE)
  dUSVTMLqELB[iM] <- mean(dimUSVTMLqE) -
    sqrt(var(dimUSVTMLqE))/sqrt(length(dimUSVTMLqE))*1.96
  dUSVTMLqEUB[iM] <- mean(dimUSVTMLqE) +
    sqrt(var(dimUSVTMLqE))/sqrt(length(dimUSVTMLqE))*1.96
  
  x <- dUSVTMLEMean[iM]
  x1 <- floor(x)
  y1 <- errorMLEASEMean[iM, x1]
  x2 <- ceiling(x)
  y2 <- errorMLEASEMean[iM, x2]
  if (x1 == x2) {
    errorMLEUSVT[iM] <- y1
  } else {
    errorMLEUSVT[iM] <- (y2-y1)/(x2-x1)*(x-x1)+y1
  }
  
  x <- dZGMLEMean[iM]
  x1 <- floor(x)
  y1 <- errorMLEASEMean[iM, x1]
  x2 <- ceiling(x)
  y2 <- errorMLEASEMean[iM, x2]
  if (x1 == x2) {
    errorMLEZG[iM] <- y1
  } else {
    errorMLEZG[iM] <- (y2-y1)/(x2-x1)*(x-x1)+y1
  }
  
  x <- dUSVTMLqEMean[iM]
  x1 <- floor(x)
  y1 <- errorMLqEASEMean[iM, x1]
  x2 <- ceiling(x)
  y2 <- errorMLqEASEMean[iM, x2]
  if (x1 == x2) {
    errorMLqEUSVT[iM] <- y1  
  } else {
    errorMLqEUSVT[iM] <- (y2-y1)/(x2-x1)*(x-x1)+y1
  }
  
  x <- dZGMLqEMean[iM]
  x1 <- floor(x)
  y1 <- errorMLqEASEMean[iM, x1]
  x2 <- ceiling(x)
  y2 <- errorMLqEASEMean[iM, x2]
  if (x1 == x2) {
    errorMLqEZG[iM] <- y1  
  } else {
    errorMLqEZG[iM] <- (y2-y1)/(x2-x1)*(x-x1)+y1
  }
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
             which = "ABar", m = mVec, d = max(dMin, 1)),
  data.frame(mse = errorMLEMean, lci = errorMLEMean, uci = errorMLEMean,
             which = "ABar", m = mVec, d = min(dMax, n)),
  data.frame(mse = c(errorMLEASEMean), lci = c(errorMLEASEMean), uci = c(errorMLEASEMean),
             which = "ABarASE", m = rep(mVec,n), d = rep(1:n, each = length(mVec))),
  data.frame(mse = errorMLqEMean, lci = errorMLqEMean, uci = errorMLqEMean,
             which = "PHat", m = mVec, d = max(dMin, 1)),
  data.frame(mse = errorMLqEMean, lci = errorMLqEMean, uci = errorMLqEMean,
             which = "PHat", m = mVec, d = min(dMax, n)),
  data.frame(mse = c(errorMLqEASEMean), lci = c(errorMLqEASEMean), uci = c(errorMLqEASEMean),
             which = "PHatASE", m = rep(mVec,n), d = rep(1:n, each = length(mVec)))) %>%
  mutate(m=factor(paste0("m=",m), sapply(mVec, function(m) {paste0("m=", m)})))

nv <- ((errorByDimDf$d >= dMin) & (errorByDimDf$d <= dMax))
errorByDimDf <- errorByDimDf[nv, ]

dimSelectionDf <- rbind(
  data.frame(mse = errorMLEZG, lci = errorMLEZG, uci = errorMLEZG,
             which = "ABar ZG 3rd", m = mVec, d = dZGMLEMean),
  # data.frame(mse = errorMLEUSVT, lci = errorMLEUSVT, uci = errorMLEUSVT,
  #            which = "ABar USVT c=0.7", m = mVec, d = dUSVTMLEMean),
  data.frame(mse = errorMLqEZG, lci = errorMLqEZG, uci = errorMLqEZG,
             which = "PHat ZG 3rd", m = mVec, d = dZGMLqEMean)) %>%
  # data.frame(mse = errorMLqEUSVT, lci = errorMLqEUSVT, uci = errorMLqEUSVT,
  #            which = "PHat USVT c=0.7", m = mVec, d = dUSVTMLqEMean)) %>%
  mutate(m=factor(paste0("m=",m), sapply(mVec, function(m) {paste0("m=", m)})))




# nv <- ((dimSelectionDf$d >= dMin) & (dimSelectionDf$d <= dMax))
# dimSelectionDf <- dimSelectionDf[nv, ]

# label_y <- with(errorByDimDf, .75*max(mse)+.25*min(mse))
label_y <- 7e+6

lSize = .8
legendSize = 1.5

gg <- ggplot(errorByDimDf, aes(x = d, y = mse, linetype = factor(which), shape = factor(which))) +
  facet_wrap(~m) +
  geom_point(data = dimSelectionDf, size = 3) +
  # scale_linetype_manual(name = "", values = c(1, 0, 0, 2, 3, 0, 0, 4)) +
  # scale_shape_manual(name = "", values = c(-1, 0, 1, -1, -1, 0, 2, -1)) +
  scale_linetype_manual(name = "", values = c(1, 0, 2, 3, 0, 4),
                        labels = c("MLE", "MLE ZG", "MLE_ASE", "MLqE", "MLqE ZG", "MLqE_ASE")) +
  scale_shape_manual(name = "", values = c(-1, 0, -1, -1, 1, -1),
                     labels = c("MLE", "MLE ZG", "MLE_ASE", "MLqE", "MLqE ZG", "MLqE_ASE")) +
  # scale_linetype_manual(name = "", values = c(1, 2, 3, 4),
  #                       labels = c("MLE", "MLE_ASE", "MLqE", "MLqE_ASE")) +
  # scale_shape_manual(name = "", values = c(-1, -1, -1, -1),
  #                    labels = c("MLE", "MLE_ASE", "MLqE", "MLqE_ASE")) +
  geom_line(alpha = 1, size = lSize) +
  # scale_y_continuous(limits = c(0, label_y)) + 
  # geom_linerange(aes(ymin = lci, ymax = uci), alpha = .5, size = 1) +
  xlab("dimension")+ylab("MSE")+
  theme(strip.text.x = element_text(size=20,face="bold"))+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"))+
  theme(panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
  theme(legend.text=element_text(size=20,face="bold"))+
  theme(legend.position="bottom")+
  ggtitle(paste0(dataName1, " using ", dataName2, " as baseline, n=", n, ", ", M, " graphs"))+
  theme(legend.key.size=unit(legendSize,"line"))+
  theme(plot.title=element_text(lineheight=.8,size=20,face="bold"))

ggsave(paste0("../../Figure/CCI.pdf"),
       plot=gg+theme(text=element_text(size=10,family="Times")),
       # plot=gg+theme(text=element_text(size=10,family="CM Roman")),
       width=8,height=6)




