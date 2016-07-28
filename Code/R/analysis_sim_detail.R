# Simulation for LLG

rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
source("function_collection.R")


# ###### Fix m ######
m = 10
# nVec = c(30, 50, 100)
n = 100
isSVD = 0
iModel = 1

epsVec = (0:40)/100

B = matrix(c(4.2, 2, 2, 7), ncol = 2)
CB = matrix(c(20, 18, 18, 25), ncol = 2)
rho = c(0.5, 0.5)
K = length(rho)
q = 0.8
d = 2

errorAbarVec = rep(0, length(epsVec))
errorAbarLBVec = rep(0, length(epsVec))
errorAbarUBVec = rep(0, length(epsVec))

errorAbarASEVec = rep(0, length(epsVec))
errorAbarASELBVec = rep(0, length(epsVec))
errorAbarASEUBVec = rep(0, length(epsVec))

errorPhatVec = rep(0, length(epsVec))
errorPhatLBVec = rep(0, length(epsVec))
errorPhatUBVec = rep(0, length(epsVec))

errorPhatASEVec = rep(0, length(epsVec))
errorPhatASELBVec = rep(0, length(epsVec))
errorPhatASEUBVec = rep(0, length(epsVec))

for (iEps in 1:length(epsVec)) {
  eps = epsVec[iEps]
  if (isSVD) {
    fileName = paste("../../Result/result_sim_", iModel, "_d_", d, "_n_", n, "_m_", m,
                     "_eps_", eps, "_q_", q, "_svd.RData", sep="")
  } else {
    fileName = paste("../../Result/result_sim_", iModel, "_d_", d, "_n_", n, "_m_", m,
                     "_eps_", eps, "_q_", q, "_eig.RData", sep="")
  }
  
  if (file.exists(fileName) == T) {
    load(fileName)
    errorAbarVec[iEps] = mean(error_A_bar)
    errorAbarLBVec[iEps] = errorAbarVec[iEps] -
      sqrt(var(error_A_bar))/sqrt(length(error_A_bar))*1.96
    errorAbarUBVec[iEps] = errorAbarVec[iEps] +
      sqrt(var(error_A_bar))/sqrt(length(error_A_bar))*1.96
    
    errorAbarASEVec[iEps] = mean(error_A_bar_ase)
    errorAbarASELBVec[iEps] = errorAbarASEVec[iEps] -
      sqrt(var(error_A_bar_ase))/sqrt(length(error_A_bar_ase))*1.96
    errorAbarASEUBVec[iEps] = errorAbarASEVec[iEps] +
      sqrt(var(error_A_bar_ase))/sqrt(length(error_A_bar_ase))*1.96
    
    errorPhatVec[iEps] = mean(error_P_hat)
    errorPhatLBVec[iEps] = errorPhatVec[iEps] -
      sqrt(var(error_P_hat))/sqrt(length(error_P_hat))*1.96
    errorPhatUBVec[iEps] = errorPhatVec[iEps] +
      sqrt(var(error_P_hat))/sqrt(length(error_P_hat))*1.96
    
    errorPhatASEVec[iEps] = mean(error_P_hat_ase)
    errorPhatASELBVec[iEps] = errorPhatASEVec[iEps] -
      sqrt(var(error_P_hat_ase))/sqrt(length(error_P_hat_ase))*1.96
    errorPhatASEUBVec[iEps] = errorPhatASEVec[iEps] +
      sqrt(var(error_P_hat_ase))/sqrt(length(error_P_hat_ase))*1.96
  }
}

plot.new()
lines(errorAbarVec, col="black")
lines(errorAbarASEVec, col="red")
lines(errorPhatVec, col="green")
lines(errorPhatASEVec, col="blue")

dfError <- rbind(
  data.frame(mse=c(errorAbarVec),lci=c(errorAbarLBVec),uci=c(errorAbarUBVec),
             which="MLE",eps=epsVec),
  data.frame(mse=c(errorAbarASEVec),lci=c(errorAbarASELBVec),uci=c(errorAbarASEUBVec),
             which="MLE_ASE",eps=epsVec),
  data.frame(mse=c(errorPhatVec),lci=c(errorPhatLBVec),uci=c(errorPhatUBVec),
             which="MLqE",eps=epsVec),
  data.frame(mse=c(errorPhatASEVec),lci=c(errorPhatASELBVec),uci=c(errorPhatASEUBVec),
             which="MLqE_ASE",eps=epsVec))


require(ggplot2)

lSize = .8
legendSize = 1.5

gg <- ggplot(dfError, aes(x=eps, y=mse, group=which, colour=which)) + 
  geom_line(size=2) +
  xlab("epsilon")+ylab("MSE")+
  theme(strip.text.x = element_text(size=20,face="bold"))+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.text=element_text(size=20,face="bold"))+
  theme(legend.position="bottom")+
  ggtitle(paste0("n=", n, ", m=", m, ", q=", q))+
  theme(legend.key.size=unit(legendSize,"line"))+
  theme(plot.title=element_text(lineheight=.8,size=20,face="bold")) +
  theme(legend.title=element_blank())
gg

# gg <- ggplot(dfError,aes(x=eps,y=mse,linetype=factor(which),shape=factor(which)))+
#   #   facet_wrap(~m)+
#   #   geom_point(data=subset(dim_selection_df,which=="ZG 3rd"),size=2,colour="red")+
#   #   geom_point(data=subset(dim_selection_df,which=="USVT c=0.7"),size=2,colour="blue")+
#   #   geom_point(data=dim_selection_df,size=3)+
#   scale_linetype_manual(name="",values=c(1,1,1,1))+
#   scale_shape_manual(name="",values=c(-1,-1,15,17))+
#   #   geom_point(dim_selection_df,aes(shape=which))+
#   geom_line(alpha=1,size=lSize)+
#   geom_linerange(aes(ymin=lci,ymax=uci),alpha=.5,size=1)+
#   #   geom_vline(data=dim_selection_df,
#   #              aes(xintercept=value,color=which,linetype=variable))+
#   #   scale_linetype_manual(name="",values=c(1,2,3,4))+
#   #   geom_text(data=dim_selection_df %>% filter(variable=="mean"),
#   #             aes(x=value+n/30,y=label_y,linetype=variable,label=which,color=which),angle=90)+
#   #   scale_color_discrete(guide=FALSE)+
#   xlab("")+ylab("MSE")+
#   theme(strip.text.x = element_text(size=20,face="bold"))+
#   theme(axis.text=element_text(size=15),
#         axis.title=element_text(size=20,face="bold"))+
#   theme(panel.grid.major = element_line(colour="grey95"),
#         panel.grid.minor = element_blank())+
#   theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
#   theme(legend.text=element_text(size=20,face="bold"))+
#   theme(legend.position="bottom")+
#   ggtitle(paste0(", N=", n, ", ", m, " graphs"))+
#   theme(legend.key.size=unit(legendSize,"line"))+
#   theme(plot.title=element_text(lineheight=.8,size=20,face="bold"))
