rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")

q <- 0.9
dataName <- "CPAC200"
isSVD <- 0

if (isSVD) {
  fileName = paste("../../Result/result_", dataName, "_test_0_q_", q, "_svd.RData", sep="")
} else {
  fileName = paste("../../Result/result_", dataName, "_test_0_q_", q, "_eig.RData", sep="")
}
load(fileName)
errorABarBetweenAll <- errorABarBetween
errorABarASEBetweenAll <- errorABarASEBetween
errorPHatBetweenAll <- errorPHatBetween
errorPHatASEBetweenAll <- errorPHatASEBetween

nBetween <- length(errorABarBetween)
nGroup <- nBetween*2 - 1
errorABarCrossAll <- rep(0, nGroup*nBetween)
errorABarASECrossAll <- rep(0, nGroup*nBetween)
errorPHatCrossAll <- rep(0, nGroup*nBetween)
errorPHatASECrossAll <- rep(0, nGroup*nBetween)
for (iGroup in 1:nGroup) {
  if (isSVD) {
    fileName = paste("../../Result/result_", dataName, "_test_", iGroup, "_q_", q, "_svd.RData", sep="")
  } else {
    fileName = paste("../../Result/result_", dataName, "_test_", iGroup, "_q_", q, "_eig.RData", sep="")
  }
  load(fileName)
  errorABarCrossAll[((iGroup - 1)*nBetween + 1):(iGroup*nBetween)] <- errorABarCross
  errorABarASECrossAll[((iGroup - 1)*nBetween + 1):(iGroup*nBetween)] <- errorABarASECross
  errorPHatCrossAll[((iGroup - 1)*nBetween + 1):(iGroup*nBetween)] <- errorPHatCross
  errorPHatASECrossAll[((iGroup - 1)*nBetween + 1):(iGroup*nBetween)] <- errorPHatASECross
}

nv <- rep(TRUE, nGroup*nBetween)
ind <- 1
for (i in 453:2) {
  ind <- ind + i
  nv[ind] <- FALSE
}
errorABarCrossAll <- errorABarCrossAll[nv]
errorABarASECrossAll <- errorABarASECrossAll[nv]
errorPHatCrossAll <- errorPHatCrossAll[nv]
errorPHatASECrossAll <- errorPHatASECrossAll[nv]

# library(ggplot2)
# df1 <- data.frame(mse = c(errorABarBetweenAll, errorABarCrossAll),
#                  lines = rep(c("Between", "Cross"),
#                              c(length(errorABarBetweenAll), length(errorABarCrossAll))))
# p1 <- ggplot(df1, aes(x = mse, fill = lines)) +
#   geom_density(alpha = 0.5) + 
#   xlab("MSE") +
#   ggtitle("ABar")
# p1
# 
# df2 <- data.frame(mse = c(errorABarASEBetweenAll, errorABarASECrossAll),
#                  lines = rep(c("Between", "Cross"),
#                              c(length(errorABarASEBetweenAll), length(errorABarASECrossAll))))
# p2 <- ggplot(df2, aes(x = mse, fill = lines)) +
#   geom_density(alpha = 0.5) + 
#   xlab("MSE") +
#   ggtitle("ABarASE")
# p2
# 
# df3 <- data.frame(mse = c(errorPHatBetweenAll, errorPHatCrossAll),
#                   lines = rep(c("Between", "Cross"),
#                               c(length(errorPHatBetweenAll), length(errorPHatCrossAll))))
# p3 <- ggplot(df3, aes(x = mse, fill = lines)) +
#   geom_density(alpha = 0.5) + 
#   xlab("MSE") +
#   ggtitle("PHat")
# p3
# 
# df4 <- data.frame(mse = c(errorPHatASEBetweenAll, errorPHatASECrossAll),
#                   lines = rep(c("Between", "Cross"),
#                               c(length(errorPHatASEBetweenAll), length(errorPHatASECrossAll))))
# p4 <- ggplot(df4, aes(x = mse, fill = lines)) +
#   geom_density(alpha = 0.5) + 
#   xlab("MSE") +
#   ggtitle("PHatASE")
# p4
# 
# library(gridExtra)
# library(grid)
# source("function_collection.R")
# grid_arrange_shared_legend2(list(p1, p2, p3, p4), 2, 2)


plot (density(errorABarBetweenAll))
lines (density(errorABarCrossAll))

lines (density(errorABarASEBetweenAll), col="red")
lines (density(errorABarASECrossAll), col="red")

lines (density(errorPHatBetweenAll), col="blue")
lines (density(errorPHatCrossAll), col="blue")

lines (density(errorPHatASEBetweenAll), col="green")
lines (density(errorPHatASECrossAll), col="green")

legend("topright", col=c("black", "red", "blue", "green"), lty=rep(1, 4),
       legend=c("MLE", "MLE_ASE", "MLqE", "MLqE_ASE"))





