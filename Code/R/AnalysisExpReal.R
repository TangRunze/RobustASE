rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")

mVec <- c(2, 5, 10)
# mVec <- c(1, 5, 10)
# mVec <- 5
q <- 0.9
nIter <- 100
nCores <- 2
# dataName <- "desikan"
dataName <- "CPAC200"
# dataName <- "JHU"

dMin <- 50
dMax <- Inf

isSVD <- 1

source("function_collection.R")

library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

for (iM in 1:length(mVec)) {
  m <- mVec[iM]
  
  # Eigen-decomposition
  if (isSVD) {
    fileName <- paste("../../Result/result_", dataName, "_brute_", "m_", m, "_q_", q, "_svd.RData", sep="")
  } else {
    fileName <- paste("../../Result/result_", dataName, "_brute_", "m_", m, "_q_", q, "_eig.RData", sep="")
  }
  load(fileName)
  
  if (iM == 1) {
    errorABarMean <- array(0, dim = c(length(mVec), 1))
    errorABarLower <- array(0, dim = c(length(mVec), 1))
    errorABarUpper <- array(0, dim = c(length(mVec), 1))
    errorPHatMean <- array(0, dim = c(length(mVec), 1))
    errorPHatLower <- array(0, dim = c(length(mVec), 1))
    errorPHatUpper <- array(0, dim = c(length(mVec), 1))
    errorABarASEMean <- array(0, dim = c(length(mVec), n))
    errorABarASELower <- array(0, dim = c(length(mVec), n))
    errorABarASEUpper <- array(0, dim = c(length(mVec), n))
    errorPHatASEMean <- array(0, dim = c(length(mVec), n))
    errorPHatASELower <- array(0, dim = c(length(mVec), n))
    errorPHatASEUpper <- array(0, dim = c(length(mVec), n))
  }
  
  errorABarMean[iM] <- rep(mean(errorABar))
  errorABarLower[iM] <- errorABarMean[iM] -
    sqrt(var(errorABar))/sqrt(length(errorABar))*1.96
  errorABarUpper[iM] <- errorABarMean[iM] +
    sqrt(var(errorABar))/sqrt(length(errorABar))*1.96
  
  errorABarASEMean[iM, ] <- rowMeans(errorABarASE)
  errorABarASELower[iM, ] <- errorABarASEMean[iM, ] - 
    sqrt(apply(errorABarASE, 1, var))/sqrt(dim(errorABarASE)[2])*1.96
  errorABarASEUpper[iM, ] <- errorABarASEMean[iM, ] + 
    sqrt(apply(errorABarASE, 1, var))/sqrt(dim(errorABarASE)[2])*1.96
  
  errorPHatMean[iM] <- rep(mean(errorPHat))
  errorPHatLower[iM] <- errorPHatMean[iM] -
    sqrt(var(errorPHat))/sqrt(length(errorPHat))*1.96
  errorPHatUpper[iM] <- errorPHatMean[iM] +
    sqrt(var(errorPHat))/sqrt(length(errorPHat))*1.96
  
  errorPHatASEMean[iM, ] <- rowMeans(errorPHatASE)
  errorPHatASELower[iM, ] <- errorPHatASEMean[iM, ] - 
    sqrt(apply(errorPHatASE, 1, var))/sqrt(dim(errorPHatASE)[2])*1.96
  errorPHatASEUpper[iM, ] <- errorPHatASEMean[iM, ] + 
    sqrt(apply(errorPHatASE, 1, var))/sqrt(dim(errorPHatASE)[2])*1.96
}

dZGABarMean <- rep(0, length(mVec))
dZGABarL <- rep(0, length(mVec))
dZGABarU <- rep(0, length(mVec))
dUSVTABarMean <- rep(0, length(mVec))
dUSVTABarL <- rep(0, length(mVec))
dUSVTABarU <- rep(0, length(mVec))
dZGPHatMean <- rep(0, length(mVec))
dZGPHatL <- rep(0, length(mVec))
dZGPHatU <- rep(0, length(mVec))
dUSVTPHatMean <- rep(0, length(mVec))
dUSVTPHatL <- rep(0, length(mVec))
dUSVTPHatU <- rep(0, length(mVec))
errorABarUSVT <- rep(0, length(mVec))
errorABarZG <- rep(0, length(mVec))
errorPHatUSVT <- rep(0, length(mVec))
errorPHatZG <- rep(0, length(mVec))
for (iM in 1:length(mVec)) {
  
  dZGABarMean[iM] <- mean(dimZGABar)
  dZGABarL[iM] <- mean(dimZGABar) -
    sqrt(var(dimZGABar))/sqrt(length(dimZGABar))*1.96
  dZGABarU[iM] <- mean(dimZGABar) +
    sqrt(var(dimZGABar))/sqrt(length(dimZGABar))*1.96
  
  dUSVTABarMean[iM] <- mean(dimUSVTABar)
  dUSVTABarL[iM] <- mean(dimUSVTABar) -
    sqrt(var(dimUSVTABar))/sqrt(length(dimUSVTABar))*1.96
  dUSVTABarU[iM] <- mean(dimUSVTABar) +
    sqrt(var(dimUSVTABar))/sqrt(length(dimUSVTABar))*1.96
  
  dZGPHatMean[iM] <- mean(dimZGPHat)
  dZGPHatL[iM] <- mean(dimZGPHat) -
    sqrt(var(dimZGPHat))/sqrt(length(dimZGPHat))*1.96
  dZGPHatU[iM] <- mean(dimZGPHat) +
    sqrt(var(dimZGPHat))/sqrt(length(dimZGPHat))*1.96
  
  dUSVTPHatMean[iM] <- mean(dimUSVTPHat)
  dUSVTPHatL[iM] <- mean(dimUSVTPHat) -
    sqrt(var(dimUSVTPHat))/sqrt(length(dimUSVTPHat))*1.96
  dUSVTPHatU[iM] <- mean(dimUSVTPHat) +
    sqrt(var(dimUSVTPHat))/sqrt(length(dimUSVTPHat))*1.96
  
  x <- dUSVTABarMean[iM]
  x1 <- floor(x)
  y1 <- errorABarASEMean[iM, x1]
  x2 <- ceiling(x)
  y2 <- errorABarASEMean[iM, x2]
  errorABarUSVT[iM] <- (y2-y1)/(x2-x1)*(x-x1)+y1
  
  x <- dZGABarMean[iM]
  x1 <- floor(x)
  y1 <- errorABarASEMean[iM, x1]
  x2 <- ceiling(x)
  y2 <- errorABarASEMean[iM, x2]
  errorABarZG[iM] <- (y2-y1)/(x2-x1)*(x-x1)+y1
  
  x <- dUSVTPHatMean[iM]
  x1 <- floor(x)
  y1 <- errorPHatASEMean[iM, x1]
  x2 <- ceiling(x)
  y2 <- errorPHatASEMean[iM, x2]
  errorPHatUSVT[iM] <- (y2-y1)/(x2-x1)*(x-x1)+y1
  
  x <- dZGPHatMean[iM]
  x1 <- floor(x)
  y1 <- errorPHatASEMean[iM, x1]
  x2 <- ceiling(x)
  y2 <- errorPHatASEMean[iM, x2]
  errorPHatZG[iM] <- (y2-y1)/(x2-x1)*(x-x1)+y1
}

# errorByDimDf <- rbind(
#   data.frame(mse = errorABarMean, lci = errorABarLower, uci = errorABarUpper,
#              which = "ABar", m = mVec, d = 1),
#   data.frame(mse = errorABarMean, lci = errorABarLower, uci = errorABarUpper,
#              which = "ABar", m = mVec, d = n),
#   data.frame(mse = c(errorABarASEMean), lci = c(errorABarASELower), uci = c(errorABarASEUpper),
#              which = "ABarASE", m = rep(mVec,n), d = rep(1:n, each = length(mVec))),
#   data.frame(mse = errorPHatMean, lci = errorPHatLower, uci = errorPHatUpper,
#              which = "PHat", m = mVec, d = 1),
#   data.frame(mse = errorPHatMean, lci = errorPHatLower, uci = errorPHatUpper,
#              which = "PHat", m = mVec, d = n),
#   data.frame(mse = c(errorPHatASEMean), lci = c(errorPHatASELower), uci = c(errorPHatASEUpper),
#              which = "PHatASE", m = rep(mVec,n), d = rep(1:n, each = length(mVec)))) %>%
#   mutate(m = factor(paste0("m=", m), c("m=2", "m=5")))

errorByDimDf <- rbind(
  data.frame(mse = errorABarMean, lci = errorABarMean, uci = errorABarMean,
             which = "ABar", m = mVec, d = max(dMin, 1)),
  data.frame(mse = errorABarMean, lci = errorABarMean, uci = errorABarMean,
             which = "ABar", m = mVec, d = min(dMax, n)),
  data.frame(mse = c(errorABarASEMean), lci = c(errorABarASEMean), uci = c(errorABarASEMean),
             which = "ABarASE", m = rep(mVec,n), d = rep(1:n, each = length(mVec))),
  data.frame(mse = errorPHatMean, lci = errorPHatMean, uci = errorPHatMean,
             which = "PHat", m = mVec, d = max(dMin, 1)),
  data.frame(mse = errorPHatMean, lci = errorPHatMean, uci = errorPHatMean,
             which = "PHat", m = mVec, d = min(dMax, n)),
  data.frame(mse = c(errorPHatASEMean), lci = c(errorPHatASEMean), uci = c(errorPHatASEMean),
             which = "PHatASE", m = rep(mVec,n), d = rep(1:n, each = length(mVec)))) %>%
  mutate(m=factor(paste0("m=",m), sapply(mVec, function(m) {paste0("m=", m)})))

nv <- ((errorByDimDf$d >= dMin) & (errorByDimDf$d <= dMax))
errorByDimDf <- errorByDimDf[nv, ]

dimSelectionDf <- rbind(
  data.frame(mse = errorABarZG, lci = errorABarZG, uci = errorABarZG,
             which = "ABar ZG 3rd", m = mVec, d = dZGABarMean),
  data.frame(mse = errorABarUSVT, lci = errorABarUSVT, uci = errorABarUSVT,
             which = "ABar USVT c=0.7", m = mVec, d = dUSVTABarMean),
  data.frame(mse = errorPHatZG, lci = errorPHatZG, uci = errorPHatZG,
             which = "PHat ZG 3rd", m = mVec, d = dZGPHatMean),
  data.frame(mse = errorPHatUSVT, lci = errorPHatUSVT, uci = errorPHatUSVT,
             which = "PHat USVT c=0.7", m = mVec, d = dUSVTPHatMean)) %>%
  mutate(m=factor(paste0("m=",m), sapply(mVec, function(m) {paste0("m=", m)})))

# nv <- ((dimSelectionDf$d >= dMin) & (dimSelectionDf$d <= dMax))
# dimSelectionDf <- dimSelectionDf[nv, ]

label_y <- with(errorByDimDf, .75*max(mse)+.25*min(mse))

lSize = .8
legendSize = 1.5

gg <- ggplot(errorByDimDf, aes(x = d, y = mse, linetype = factor(which), shape = factor(which))) +
  facet_wrap(~m) +
  # geom_point(data = dimSelectionDf, size = 3) +
  # scale_linetype_manual(name = "", values = c(1, 0, 0, 2, 3, 0, 0, 4)) +
  # scale_shape_manual(name = "", values = c(-1, 0, 0, -1, -1, 0, 0, -1)) +
  scale_linetype_manual(name = "", values = c(1, 2, 3, 4)) +
  scale_shape_manual(name = "", values = c(-1, -1, -1, -1)) +
  geom_line(alpha = 1, size = lSize) +
  geom_linerange(aes(ymin = lci, ymax = uci), alpha = .5, size = 1) +
  xlab("Dimension")+ylab("MSE")+
  theme(strip.text.x = element_text(size=20,face="bold"))+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"))+
  theme(panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
  theme(legend.text=element_text(size=20,face="bold"))+
  theme(legend.position="bottom")+
  ggtitle(paste0(dataName, ", N=", n, ", ", M, " graphs"))+
  theme(legend.key.size=unit(legendSize,"line"))+
  theme(plot.title=element_text(lineheight=.8,size=20,face="bold"))

# pp[[3]]=gg

# source("function_collection.R")

# library(gridExtra)
# library(grid)

# grid_arrange_shared_legend2(list(pp[[1]], pp[[2]], pp[[3]]), 3, 1)





gg






