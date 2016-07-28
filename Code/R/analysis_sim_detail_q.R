# Simulation for LLG

rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
source("function_collection.R")


# ###### Fix m ######
m = 10
n = 100
isSVD = 0
iModel = 1

eps = 0.2
qVec = (80:100)/100

B = matrix(c(4.2, 2, 2, 7), ncol = 2)
CB = matrix(c(20, 18, 18, 25), ncol = 2)
rho = c(0.5, 0.5)
K = length(rho)
d = 2

errorAbarVec = rep(0, length(qVec))
errorAbarLBVec = rep(0, length(qVec))
errorAbarUBVec = rep(0, length(qVec))

errorAbarASEVec = rep(0, length(qVec))
errorAbarASELBVec = rep(0, length(qVec))
errorAbarASEUBVec = rep(0, length(qVec))

errorPhatVec = rep(0, length(qVec))
errorPhatLBVec = rep(0, length(qVec))
errorPhatUBVec = rep(0, length(qVec))

errorPhatASEVec = rep(0, length(qVec))
errorPhatASELBVec = rep(0, length(qVec))
errorPhatASEUBVec = rep(0, length(qVec))

for (iQ in 1:length(qVec)) {
  q = qVec[iQ]
  if (isSVD) {
    fileName = paste("../../Result/result_sim_", iModel, "_d_", d, "_n_", n, "_m_", m,
                     "_eps_", eps, "_q_", q, "_svd.RData", sep="")
  } else {
    fileName = paste("../../Result/result_sim_", iModel, "_d_", d, "_n_", n, "_m_", m,
                     "_eps_", eps, "_q_", q, "_eig.RData", sep="")
  }
  
  if (file.exists(fileName) == T) {
    load(fileName)
    errorAbarVec[iQ] = mean(error_A_bar)
    errorAbarLBVec[iQ] = errorAbarVec[iQ] -
      sqrt(var(error_A_bar))/sqrt(length(error_A_bar))*1.96
    errorAbarUBVec[iQ] = errorAbarVec[iQ] +
      sqrt(var(error_A_bar))/sqrt(length(error_A_bar))*1.96
    
    errorAbarASEVec[iQ] = mean(error_A_bar_ase)
    errorAbarASELBVec[iQ] = errorAbarASEVec[iQ] -
      sqrt(var(error_A_bar_ase))/sqrt(length(error_A_bar_ase))*1.96
    errorAbarASEUBVec[iQ] = errorAbarASEVec[iQ] +
      sqrt(var(error_A_bar_ase))/sqrt(length(error_A_bar_ase))*1.96
    
    errorPhatVec[iQ] = mean(error_P_hat)
    errorPhatLBVec[iQ] = errorPhatVec[iQ] -
      sqrt(var(error_P_hat))/sqrt(length(error_P_hat))*1.96
    errorPhatUBVec[iQ] = errorPhatVec[iQ] +
      sqrt(var(error_P_hat))/sqrt(length(error_P_hat))*1.96
    
    errorPhatASEVec[iQ] = mean(error_P_hat_ase)
    errorPhatASELBVec[iQ] = errorPhatASEVec[iQ] -
      sqrt(var(error_P_hat_ase))/sqrt(length(error_P_hat_ase))*1.96
    errorPhatASEUBVec[iQ] = errorPhatASEVec[iQ] +
      sqrt(var(error_P_hat_ase))/sqrt(length(error_P_hat_ase))*1.96
  }
}

dfError <- rbind(
  data.frame(mse=c(errorAbarVec),lci=c(errorAbarLBVec),uci=c(errorAbarUBVec),
             which="MLE",eps=qVec),
  data.frame(mse=c(errorAbarASEVec),lci=c(errorAbarASELBVec),uci=c(errorAbarASEUBVec),
             which="MLE_ASE",eps=qVec),
  data.frame(mse=c(errorPhatVec),lci=c(errorPhatLBVec),uci=c(errorPhatUBVec),
             which="MLqE",eps=qVec),
  data.frame(mse=c(errorPhatASEVec),lci=c(errorPhatASELBVec),uci=c(errorPhatASEUBVec),
             which="MLqE_ASE",eps=qVec))


require(ggplot2)

lSize = .8
legendSize = 1.5

gg <- ggplot(dfError, aes(x=eps, y=mse, group=which, colour=which)) + 
  geom_line(size=2) +
  xlab("q")+ylab("MSE")+
  theme(strip.text.x = element_text(size=20,face="bold"))+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.text=element_text(size=20,face="bold"))+
  theme(legend.position="bottom")+
  ggtitle(paste0("n=", n, ", m=", m, ", epsilon=", eps))+
  theme(legend.key.size=unit(legendSize,"line"))+
  theme(plot.title=element_text(lineheight=.8,size=20,face="bold")) +
  theme(legend.title=element_blank())
gg

