rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")
# source("function_collection.R")

m <- 10
n <- 100
isSVD<- 0
iModel <- 2

epsVec <- (0:10)/50

q <- 0.9
d <- 2

errorMLEVec <- rep(0, length(epsVec))
errorMLELBVec <- rep(0, length(epsVec))
errorMLEUBVec <- rep(0, length(epsVec))

errorMLEASEVec <- rep(0, length(epsVec))
errorMLEASELBVec <- rep(0, length(epsVec))
errorMLEASEUBVec <- rep(0, length(epsVec))

errorMLqEVec <- rep(0, length(epsVec))
errorMLqELBVec <- rep(0, length(epsVec))
errorMLqEUBVec <- rep(0, length(epsVec))

errorMLqEASEVec <- rep(0, length(epsVec))
errorMLqEASELBVec <- rep(0, length(epsVec))
errorMLqEASEUBVec <- rep(0, length(epsVec))

for (iEps in 1:length(epsVec)) {
  eps <- epsVec[iEps]
  if (isSVD) {
    fileName <- paste0("../../Result/result_sim_", iModel, "_d_", d, "_n_", n, "_m_", m,
                       "_eps_", eps, "_q_", q, "_svd.RData")
  } else {
    fileName <- paste0("../../Result/result_sim_", iModel, "_d_", d, "_n_", n, "_m_", m,
                       "_eps_", eps, "_q_", q, "_eig.RData")
  }
  
  if (file.exists(fileName) == T) {
    load(fileName)
    
    errorMLEVec[iEps] <- mean(error_MLE)
    errorMLELBVec[iEps] <- errorMLEVec[iEps] -
      sqrt(var(error_MLE))/sqrt(length(error_MLE))*1.96
    errorMLEUBVec[iEps] <- errorMLEVec[iEps] +
      sqrt(var(error_MLE))/sqrt(length(error_MLE))*1.96
    
    errorMLEASEVec[iEps] <- mean(error_MLE_ASE)
    errorMLEASELBVec[iEps] <- errorMLEASEVec[iEps] -
      sqrt(var(error_MLE_ASE))/sqrt(length(error_MLE_ASE))*1.96
    errorMLEASEUBVec[iEps] <- errorMLEASEVec[iEps] +
      sqrt(var(error_MLE_ASE))/sqrt(length(error_MLE_ASE))*1.96
    
    errorMLqEVec[iEps] <- mean(error_MLqE)
    errorMLqELBVec[iEps] <- errorMLqEVec[iEps] -
      sqrt(var(error_MLqE))/sqrt(length(error_MLqE))*1.96
    errorMLqEUBVec[iEps] <- errorMLqEVec[iEps] +
      sqrt(var(error_MLqE))/sqrt(length(error_MLqE))*1.96
    
    errorMLqEASEVec[iEps] <- mean(error_MLqE_ASE)
    errorMLqEASELBVec[iEps] <- errorMLqEASEVec[iEps] -
      sqrt(var(error_MLqE_ASE))/sqrt(length(error_MLqE_ASE))*1.96
    errorMLqEASEUBVec[iEps] <- errorMLqEASEVec[iEps] +
      sqrt(var(error_MLqE_ASE))/sqrt(length(error_MLqE_ASE))*1.96
  }
}

# plot.new()
# lines(errorMLEVec, col="black")
# lines(errorMLEASEVec, col="red")
# lines(errorMLqEVec, col="green")
# lines(errorMLqEASEVec, col="blue")

dfError <- rbind(
  data.frame(mse=c(errorMLEVec),lci=c(errorMLELBVec),uci=c(errorMLEUBVec),
             which="MLE",eps=epsVec),
  data.frame(mse=c(errorMLEASEVec),lci=c(errorMLEASELBVec),uci=c(errorMLEASEUBVec),
             which="MLE_ASE",eps=epsVec),
  data.frame(mse=c(errorMLqEVec),lci=c(errorMLqELBVec),uci=c(errorMLqEUBVec),
             which="MLqE",eps=epsVec),
  data.frame(mse=c(errorMLqEASEVec),lci=c(errorMLqEASELBVec),uci=c(errorMLqEASEUBVec),
             which="MLqE_ASE",eps=epsVec))

# dfError = dfError[dfError$eps<=0.2,]

require(ggplot2)

lSize = .8
legendSize = 1.5

# gg <- ggplot(dfError, aes(x=eps, y=mse, group=which, colour=which)) + 
#   geom_line(size=2) +
#   xlab("epsilon")+ylab("MSE")+
#   theme(strip.text.x = element_text(size=20,face="bold"))+
#   theme(axis.text=element_text(size=15),
#         axis.title=element_text(size=20,face="bold"))+
#   theme(legend.text=element_text(size=20,face="bold"))+
#   theme(legend.position="bottom")+
#   ggtitle(paste0("n=", n, ", m=", m, ", q=", q))+
#   theme(legend.key.size=unit(legendSize,"line"))+
#   theme(plot.title=element_text(lineheight=.8,size=20,face="bold")) +
#   theme(legend.title=element_blank())

gg <- ggplot(dfError, aes(x=eps, y=mse, linetype=which, shape=which)) + 
  scale_linetype_manual(name="", values=c(1, 3, 2, 6),
                        labels=c(expression(widehat(P)^{(1)}),
                                 expression(widetilde(P)^{(1)}),
                                 expression(widehat(P)^{(q)}),
                                 expression(widetilde(P)^{(q)}))) +
  scale_shape_manual(name="", values=c(-1, -1, -1, -1)) +
  # geom_line(size=2) +
  geom_line() +
  geom_linerange(aes(ymin=lci,ymax=uci),linetype=1,alpha=.5,size=.5) +
  xlab(expression(epsilon))+ylab("MSE") +
  theme(panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank()) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey70')) +
  theme(strip.text.x = element_text(size=20,face="bold")) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20,face="bold")) +
  theme(legend.text=element_text(size=15,face="bold")) +
  theme(legend.position="bottom") +
  ggtitle(TeX(sprintf('n = %d, m = %d, q = %.1f', n, m, q))) +
  theme(legend.key.size=unit(legendSize,"line"))+
  theme(plot.title=element_text(lineheight=.8,size=20,face="bold")) +
  theme(legend.title=element_blank())
gg
ggsave("../../Result/sim_eps.pdf",
       plot=gg+theme(text=element_text(size=10,family="Times")),
       width=5.5,height=4)
