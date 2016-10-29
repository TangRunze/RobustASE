# Simulation for LLG

rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
source("function_collection.R")


# ###### Fix m ######
m = 5
isSVD = 0
dataName = "desikan"

if (isSVD) {
  fileName = paste("../../Result/result_", dataName, "_brute_", "m_", m, "_svd.RData", sep="")
} else {
  fileName = paste("../../Result/result_", dataName, "_brute_", "m_", m, "_eig.RData", sep="")
}

if (file.exists(fileName) == T) {
  load(fileName)
  
  errorAbarASEVec = rep(0, length(dVec))
  errorAbarASELBVec = rep(0, length(dVec))
  errorAbarASEUBVec = rep(0, length(dVec))
  
  errorPhatASEVec = rep(0, length(dVec))
  errorPhatASELBVec = rep(0, length(dVec))
  errorPhatASEUBVec = rep(0, length(dVec))
  
  
  errorAbarVec = mean(error_A_bar)
  errorAbarLBVec = errorAbarVec -
    sqrt(var(error_A_bar))/sqrt(length(error_A_bar))*1.96
  errorAbarUBVec = errorAbarVec +
    sqrt(var(error_A_bar))/sqrt(length(error_A_bar))*1.96
  
  for (iD in 1:length(dVec)) {
    errorAbarASEVec[iD] = mean(error_A_bar_ase[iD,])
    errorAbarASELBVec[iD] = errorAbarASEVec[iD] -
      sqrt(var(error_A_bar_ase[iD,]))/sqrt(length(error_A_bar_ase[iD,]))*1.96
    errorAbarASEUBVec[iD] = errorAbarASEVec[iD] +
      sqrt(var(error_A_bar_ase[iD,]))/sqrt(length(error_A_bar_ase[iD,]))*1.96
  }
  
  errorPhatVec = mean(error_P_hat)
  errorPhatLBVec = errorPhatVec -
    sqrt(var(error_P_hat))/sqrt(length(error_P_hat))*1.96
  errorPhatUBVec = errorPhatVec +
    sqrt(var(error_P_hat))/sqrt(length(error_P_hat))*1.96
  
  for (iD in 1:length(dVec)) {
    errorPhatASEVec[iD] = mean(error_P_hat_ase[iD,])
    errorPhatASELBVec[iD] = errorPhatASEVec[iD] -
      sqrt(var(error_P_hat_ase[iD,]))/sqrt(length(error_P_hat_ase[iD,]))*1.96
    errorPhatASEUBVec[iD] = errorPhatASEVec[iD] +
      sqrt(var(error_P_hat_ase[iD,]))/sqrt(length(error_P_hat_ase[iD,]))*1.96
  }
}

dfError <- rbind(
  data.frame(mse=rep(errorAbarVec,2),lci=rep(errorAbarLBVec,2),uci=rep(errorAbarUBVec,2),
             which="MLE",d=dVec),
  data.frame(mse=c(errorAbarASEVec),lci=c(errorAbarASELBVec),uci=c(errorAbarASEUBVec),
             which="MLE_ASE",d=dVec),
  data.frame(mse=rep(errorPhatVec,2),lci=rep(errorPhatLBVec,2),uci=rep(errorPhatUBVec,2),
             which="MLqE",d=dVec),
  data.frame(mse=c(errorPhatASEVec),lci=c(errorPhatASELBVec),uci=c(errorPhatASEUBVec),
             which="MLqE_ASE",d=dVec))


require(ggplot2)

lSize = .8
legendSize = 1.5

gg <- ggplot(dfError, aes(x=d, y=mse, group=which, colour=which)) + 
  geom_line(size=lSize) +
  xlab("d")+ylab("MSE")+
  theme(strip.text.x = element_text(size=20,face="bold"))+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold"))+
  theme(legend.text=element_text(size=20,face="bold"))+
  theme(legend.position="bottom")+
  ggtitle(paste0("n=", n, ", m=", m))+
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
