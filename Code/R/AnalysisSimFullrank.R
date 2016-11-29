rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")

lSize = .8
legendSize = 1.5

isSVD <- 1
dataName = "CPAC200"
# dataName = "desikan"
# dataName = "JHU"

mVec = c(1, 5, 10)

library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

for (iM in 1:length(mVec)) {
  m = mVec[iM]
  
  if (isSVD) {
    fileName = paste("../../Result/result_", dataName, "_brute_fullrank1_",
                     "m_", m, "_svd.RData", sep="")
  } else {
    fileName = paste("../../Result/result_", dataName, "_brute_fullrank1_",
                     "m_", m, "_eig.RData", sep="")
  }
  
  load(fileName)
  
  if (iM == 1) {
    errorAbarMean = array(0, dim=c(length(mVec), 1))
    errorAbarLower = array(0, dim=c(length(mVec), 1))
    errorAbarUpper = array(0, dim=c(length(mVec), 1))
    errorPhatEIGMean = array(0, dim=c(length(mVec), n))
    errorPhatEIGLower = array(0, dim=c(length(mVec), n))
    errorPhatEIGUpper = array(0, dim=c(length(mVec), n))
  }
  
  error_A_bar = error_A_bar^2*n*(n-1)
  error_P_hat = error_P_hat^2*n*(n-1)
  
  errorAbarMean[iM] = mean(error_A_bar)
  errorAbarLower[iM] = errorAbarMean[iM] - sqrt(var(error_A_bar))/sqrt(length(error_A_bar))*1.96
  errorAbarUpper[iM] = errorAbarMean[iM] + sqrt(var(error_A_bar))/sqrt(length(error_A_bar))*1.96
  errorPhatEIGMean[iM,] = rowMeans(error_P_hat)
  errorPhatEIGLower[iM,] = errorPhatEIGMean[iM,] - 
    sqrt(apply(error_P_hat, 1, var))/sqrt(dim(error_P_hat)[2])*1.96
  errorPhatEIGUpper[iM,] = errorPhatEIGMean[iM,] + 
    sqrt(apply(error_P_hat, 1, var))/sqrt(dim(error_P_hat)[2])*1.96
}











error_by_dim_df <- rbind(
  data.frame(mse=errorAbarMean,lci=errorAbarLower,uci=errorAbarUpper,
             which="Abar",m=mVec,d=n),
  data.frame(mse=errorAbarMean,lci=errorAbarLower,uci=errorAbarUpper,
             which="Abar",m=mVec,d=1),
  data.frame(mse=c(errorPhatEIGMean),lci=c(errorPhatEIGLower),uci=c(errorPhatEIGUpper),
             which="Phat",m=rep(mVec,n),d=rep(1:n,each=3))) %>%
  mutate(m=factor(paste0("M=",m),c("M=1","M=5","M=10")))

label_y <- with(error_by_dim_df, .75*max(mse)+.25*min(mse))

gg <- ggplot(error_by_dim_df,aes(x=d,y=mse,linetype=factor(which),shape=factor(which)))+
  facet_wrap(~m)+
  scale_linetype_manual(name="",values=c(1,2,0,0))+
  scale_shape_manual(name="",values=c(-1,-1,15,17))+
  geom_line(alpha=1,size=lSize)+
  geom_linerange(aes(ymin=lci,ymax=uci),alpha=.5,size=1)+
  xlab("Dimension")+ylab("MSE")+
  theme(strip.text.x = element_text(size=20,face="bold"))+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"))+
  theme(panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
  theme(legend.text=element_text(size=20,face="bold"))+
  theme(legend.position="bottom")+
  ggtitle(paste0("Simulation based on ", dataName, ", N=", n, ", M=", m))+
  theme(legend.key.size=unit(legendSize,"line"))+
  theme(plot.title=element_text(lineheight=.8,size=20,face="bold"))
print(gg)

