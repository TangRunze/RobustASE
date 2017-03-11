rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
library(latex2exp)

m <- 20
n <- 100
isSVD<- 0
iModel <- 2

qVec <- (25:50)/50
# qVec <- (5:10)/10

eps <- 0.1
d <- 2

errorMLEVec <- rep(0, length(qVec))
errorMLELBVec <- rep(0, length(qVec))
errorMLEUBVec <- rep(0, length(qVec))

errorMLEASEVec <- rep(0, length(qVec))
errorMLEASELBVec <- rep(0, length(qVec))
errorMLEASEUBVec <- rep(0, length(qVec))

errorMLqEVec <- rep(0, length(qVec))
errorMLqELBVec <- rep(0, length(qVec))
errorMLqEUBVec <- rep(0, length(qVec))

errorMLqEASEVec <- rep(0, length(qVec))
errorMLqEASELBVec <- rep(0, length(qVec))
errorMLqEASEUBVec <- rep(0, length(qVec))

for (iQ in 1:length(qVec)) {
  q <- qVec[iQ]
  if (isSVD) {
    fileName <- paste("../../Result/result_sim_", iModel, "_d_", d, "_n_", n, "_m_", m,
                     "_eps_", eps, "_q_", q, "_svd.RData", sep="")
  } else {
    fileName <- paste("../../Result/result_sim_", iModel, "_d_", d, "_n_", n, "_m_", m,
                     "_eps_", eps, "_q_", q, "_eig.RData", sep="")
  }
  
  if (file.exists(fileName) == T) {
    load(fileName)
    
    errorMLEVec[iQ] <- mean(error_MLE)
    errorMLELBVec[iQ] <- errorMLEVec[iQ] -
      sqrt(var(error_MLE))/sqrt(length(error_MLE))*1.96
    errorMLEUBVec[iQ] <- errorMLEVec[iQ] +
      sqrt(var(error_MLE))/sqrt(length(error_MLE))*1.96
    
    errorMLEASEVec[iQ] <- mean(error_MLE_ASE)
    errorMLEASELBVec[iQ] <- errorMLEASEVec[iQ] -
      sqrt(var(error_MLE_ASE))/sqrt(length(error_MLE_ASE))*1.96
    errorMLEASEUBVec[iQ] <- errorMLEASEVec[iQ] +
      sqrt(var(error_MLE_ASE))/sqrt(length(error_MLE_ASE))*1.96
    
    errorMLqEVec[iQ] <- mean(error_MLqE)
    errorMLqELBVec[iQ] <- errorMLqEVec[iQ] -
      sqrt(var(error_MLqE))/sqrt(length(error_MLqE))*1.96
    errorMLqEUBVec[iQ] <- errorMLqEVec[iQ] +
      sqrt(var(error_MLqE))/sqrt(length(error_MLqE))*1.96
    
    errorMLqEASEVec[iQ] <- mean(error_MLqE_ASE)
    errorMLqEASELBVec[iQ] <- errorMLqEASEVec[iQ] -
      sqrt(var(error_MLqE_ASE))/sqrt(length(error_MLqE_ASE))*1.96
    errorMLqEASEUBVec[iQ] <- errorMLqEASEVec[iQ] +
      sqrt(var(error_MLqE_ASE))/sqrt(length(error_MLqE_ASE))*1.96
  }
}

errorMLEVec <- rep(mean(errorMLEVec), 1, length(errorMLEVec))
errorMLELBVec <- rep(min(errorMLELBVec), 1, length(errorMLEVec))
errorMLEUBVec <- rep(max(errorMLEUBVec), 1, length(errorMLEVec))

errorMLEASEVec <- rep(mean(errorMLEASEVec), 1, length(errorMLEASEVec))
errorMLEASELBVec <- rep(min(errorMLEASELBVec), 1, length(errorMLEASEVec))
errorMLEASEUBVec <- rep(max(errorMLEASEUBVec), 1, length(errorMLEASEVec))



dfError <- rbind(
  data.frame(mse=c(errorMLEVec),lci=c(errorMLELBVec),uci=c(errorMLEUBVec),
             which="MLE",q=qVec),
  data.frame(mse=c(errorMLEASEVec),lci=c(errorMLEASELBVec),uci=c(errorMLEASEUBVec),
             which="MLE_ASE",q=qVec),
  data.frame(mse=c(errorMLqEVec),lci=c(errorMLqELBVec),uci=c(errorMLqEUBVec),
             which="MLqE",q=qVec),
  data.frame(mse=c(errorMLqEASEVec),lci=c(errorMLqEASELBVec),uci=c(errorMLqEASEUBVec),
             which="MLqE_ASE",q=qVec))


require(ggplot2)

lSize = .8
legendSize = 1.5

# gg <- ggplot(dfError, aes(x=q, y=mse, group=which, colour=which)) + 
#   geom_line(size=2) +
#   xlab("q")+ylab("MSE")+
#   theme(strip.text.x = element_text(size=20,face="bold"))+
#   theme(axis.text=element_text(size=15),
#         axis.title=element_text(size=20,face="bold"))+
#   theme(legend.text=element_text(size=20,face="bold"))+
#   theme(legend.position="bottom")+
#   ggtitle(paste0("n=", n, ", m=", m, ", epsilon=", eps))+
#   theme(legend.key.size=unit(legendSize,"line"))+
#   theme(plot.title=element_text(lineheight=.8,size=20,face="bold")) +
#   theme(legend.title=element_blank())


gg <- ggplot(dfError, aes(x=q, y=mse, linetype=which, shape=which)) + 
  scale_linetype_manual(name="", values=c(1, 3, 2, 6),
                        labels=c(expression(widehat(P)^{(1)}),
                                 expression(widetilde(P)^{(1)}),
                                 expression(widehat(P)^{(q)}),
                                 expression(widetilde(P)^{(q)}))) +
  scale_shape_manual(name="", values=c(-1, -1, -1, -1)) +
  # geom_line(size=2) +
  geom_line() +
  geom_linerange(aes(ymin=lci,ymax=uci),linetype=1,alpha=.5,size=.5) +
  xlab(expression(q))+ylab("MSE") +
  theme(panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank()) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey70')) +
  # theme(strip.text.x = element_text(size=20,face="bold")) +
  # theme(axis.text=element_text(size=15),
  #       axis.title=element_text(size=20,face="bold"))+
  # theme(legend.text=element_text(size=15,face="bold"))+
  theme(legend.position="bottom")+
  ggtitle(TeX(sprintf('n = %d, m = %d, $\\epsilon = %.1f$', n, m, eps)))
  theme(legend.key.size=unit(legendSize,"line"))+
  # theme(plot.title=element_text(lineheight=.8,size=20,face="bold")) +
  theme(legend.title=element_blank())
gg
ggsave("../../Result/sim_q.pdf",
       plot=gg+theme(text=element_text(size=15,family="Times")),
       width=5.5,height=4)
