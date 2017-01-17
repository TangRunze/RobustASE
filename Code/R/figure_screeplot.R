rm(list = ls())

###### Entire Population ######

setwd("/Users/Runze/Documents/GitHub/RobustASE/Code/R")

# dataName = "CPAC200"
# dataName = "desikan"
# dataName = "JHU"

dataNameVec = c("JHU", "desikan", "CPAC200")
dataNameDisplayVec = c("JHU", "Desikan", "CPAC200")

# dataNameVec = c("desikan", "CPAC200")
# dataNameDisplayVec = c("Desikan", "CPAC200")

source("function_collection.R")
require(ggplot2)
eigenResult <- list()
pp_scree <- list()
for (iData in 1:length(dataNameVec)) {
  dataName = dataNameVec[iData]
  tmpList = ReadDataWeighted(dataName, DA=F, newGraph=F)
  A_all = tmpList[[1]]
  n = tmpList[[2]]
  M = tmpList[[3]]
  rm(tmpList)
  Abar = add(A_all)/M
  AbarDiagAug = diag_aug(Abar)
  eigenResult[[iData]] = eigen(AbarDiagAug)$values
  yMax = max(eigenResult[[iData]])
  yMin = min(eigenResult[[iData]])
  df <- data.frame(eval=eigenResult[[iData]], k=1:n)
  label_y <- with(df, .75*yMax+.25*yMin)
  
  pp_scree[[iData]] <- ggplot(df,aes(x=k,y=eval))+
    geom_line()+
    scale_linetype_manual(name="",values=c("longdash","dotted","dotdash"))+
    xlab("order in algebraic") + ylab("eigenvalue")+
    theme(panel.grid.major = element_line(colour="grey95"),
          panel.grid.minor = element_blank())+
    theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
    theme(legend.position="none")+
    ggtitle(dataNameDisplayVec[[iData]])
  
  ggsave(paste0("../../Result/screeplot_", dataName, ".pdf"),
         pp_scree[[iData]]+theme(text=element_text(size=10,family="Times")),
         # pp_scree[[iData]]+theme(text=element_text(size=10,family="CM Roman")),
         width=2, height=2)
}





###### Different Sample Size ######

source("function_collection.R")
source("getElbows.R")
source("USVT.R")
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
set.seed(12345)

dataNameVec = c("JHU", "desikan", "CPAC200")
isSVD = 0
pp <- list()

for (iData in 1:length(dataNameVec)) {
  dataName <- dataNameVec[iData]
  
  tmpList = ReadDataWeighted(dataName, DA=F, newGraph=F)
  A_all = tmpList[[1]]
  n = tmpList[[2]]
  M = tmpList[[3]]
  rm(tmpList)
  
  mVec <- c(1, 5, 10, M)
  
  eigenResult <- array(rep(0, n*length(mVec)), dim=c(n, length(mVec)))
  
  for (iM in 1:length(mVec)) {
    m <- mVec[iM]
    if (m == M) {
      sampleVec <- 1:M
    } else {
      sampleVec = sample.int(M, m)
    }
    Abar = add(A_all[sampleVec])/m
    AbarDiagAug = diag_aug(Abar)
    eigenResult[, iM] = eigen(AbarDiagAug)$values
    
  }
  
  error_by_dim_df <- data.frame(eval=c(eigenResult), m=rep(mVec, each=n),
                                d=rep(1:n, times=length(mVec))) %>%
    mutate(m=factor(paste0("M=",m), sapply(mVec, function(m) {paste0("M=", m)})))
  
  label_y <- with(error_by_dim_df, .75*max(eval)+.25*min(eval))
  
  gg <- ggplot(error_by_dim_df, aes(x=d, y=eval))+
    facet_wrap(~m, nrow=1)+
    # scale_linetype_manual(name="",values=c(1,2,0,0))+
    # scale_shape_manual(name="",values=c(-1,-1,15,17))+
    geom_line()+
    xlab("dimension")+ylab("eigenvalue")+
    theme(panel.grid.major = element_line(colour="grey95"),
          panel.grid.minor = element_blank())+
    theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
    theme(legend.position="none")+
    ggtitle(paste0(dataName, ", N=", n))
  
  pp[[iData]]=gg
  
  ggsave(paste0("../../Result/corr_data_eval_", dataName, ".pdf"),
         plot=gg+theme(text=element_text(size=10,family="Times")),
         # plot=gg+theme(text=element_text(size=10,family="CM Roman")),
         width=5.5,height=2.5)
  
  ggsave(paste0("../../Result/corr_data_eval_", dataName, ".png"),
         plot=gg+theme(text=element_text(size=10,family="Times")),
         # plot=gg+theme(text=element_text(size=10,family="CM Roman")),
         width=5.5,height=2.5)
}






###### Different Sample Size Unweighted ######

source("function_collection.R")
source("getElbows.R")
source("USVT.R")
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
set.seed(12345)

dataNameVec = c("JHU", "desikan", "CPAC200")
isSVD = 0
pp <- list()

for (iData in 1:length(dataNameVec)) {
  dataName <- dataNameVec[iData]
  
  tmpList = read_data(dataName, DA=F, newGraph=F)
  A_all = tmpList[[1]]
  n = tmpList[[2]]
  M = tmpList[[3]]
  rm(tmpList)
  
  mVec <- c(1, 5, 10, M)
  
  eigenResult <- array(rep(0, n*length(mVec)), dim=c(n, length(mVec)))
  
  for (iM in 1:length(mVec)) {
    m <- mVec[iM]
    if (m == M) {
      sampleVec <- 1:M
    } else {
      sampleVec = sample.int(M, m)
    }
    Abar = add(A_all[sampleVec])/m
    AbarDiagAug = diag_aug(Abar)
    eigenResult[, iM] = eigen(AbarDiagAug)$values
    
  }
  
  error_by_dim_df <- data.frame(eval=c(eigenResult), m=rep(mVec, each=n),
                                d=rep(1:n, times=length(mVec))) %>%
    mutate(m=factor(paste0("M=",m), sapply(mVec, function(m) {paste0("M=", m)})))
  
  label_y <- with(error_by_dim_df, .75*max(eval)+.25*min(eval))
  
  gg <- ggplot(error_by_dim_df, aes(x=d, y=eval))+
    facet_wrap(~m, nrow=1)+
    # scale_linetype_manual(name="",values=c(1,2,0,0))+
    # scale_shape_manual(name="",values=c(-1,-1,15,17))+
    geom_line()+
    xlab("dimension")+ylab("eigenvalue")+
    theme(panel.grid.major = element_line(colour="grey95"),
          panel.grid.minor = element_blank())+
    theme(panel.background = element_rect(fill = 'white', colour = 'grey70'))+
    theme(legend.position="none")+
    ggtitle(paste0(dataName, ", unweighted, N=", n))
  
  pp[[iData]]=gg
  
  ggsave(paste0("../../Result/corr_data_eval_", dataName, "_unweighted.pdf"),
         plot=gg+theme(text=element_text(size=10,family="Times")),
         # plot=gg+theme(text=element_text(size=10,family="CM Roman")),
         width=5.5,height=2.5)
  
  ggsave(paste0("../../Result/corr_data_eval_", dataName, "_unweighted.png"),
         plot=gg+theme(text=element_text(size=10,family="Times")),
         # plot=gg+theme(text=element_text(size=10,family="CM Roman")),
         width=5.5,height=2.5)
}
