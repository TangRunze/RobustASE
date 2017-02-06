add <- function(x) {
  Reduce("+", x)
}


ReadData <- function(dataName, weighted=T, newGraph=T, DA=F) {
  if (weighted) {
    strWeighted <- "Weighted"
  } else {
    strWeighted <- "Unweighted"
  }
  if (newGraph) {
    strNewGraph <- "New"
  } else {
    strNewGraph <- "Old"
  }
  if (DA) {
    strDA <- "_DA"
  } else {
    strDA <- ""
  }
  if (dataName %in% c("JHU", "desikan", "CPAC200")) {
    fileName <- paste0("../../Data/Data_", dataName, "_", strWeighted, "_", strNewGraph,
                       strDA, ".RData")
  } else {
    fileName <- paste0("../../Data/Data_", dataName, "_", strWeighted, strDA, ".RData")
  }
  if (file.exists(fileName)) {
    load(fileName)
    return(list(AList, n, M))
  } else {
    if (dataName %in% c("JHU", "desikan", "CPAC200")) {
      require(igraph)
      subjectsID <- readLines("../../Data/subnames.txt")
      if (newGraph == F) {
        g <- read_graph(paste0("../../Data/", dataName, "/SWU4_", subjectsID[1], 
                               "_1_", dataName, "_sg.graphml"), format="graphml")
      } else {
        g <- read_graph(paste0("../../Data/", dataName, "_new/SWU4_", subjectsID[1], 
                               "_1_DTI_", dataName, ".graphml"), format="graphml")      
      }
      n <- vcount(g)
      M <- length(subjectsID)*2;
      AList <- list()
      for (iSub in 1:length(subjectsID)) {
        for (iSession in 1:2) {
          if (newGraph == F) {
            g <- read_graph(paste0("../../Data/", dataName, "/SWU4_", subjectsID[iSub], 
                                   "_", iSession, "_", dataName, "_sg.graphml"),
                            format = "graphml")
          } else {
            g <- read_graph(paste0("../../Data/", dataName, "_new/SWU4_", subjectsID[iSub], 
                                   "_", iSession, "_DTI_", dataName, ".graphml"),
                            format = "graphml")
          }
          if (weighted) {
            A <- as_adj(g, attr="weight", type="both", sparse=FALSE)
          } else {
            A <- as_adj(g, type="both", sparse=FALSE)
          }
          #CHANGE HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          # if (DA) {
          #   A <- diag_aug(A)
          # }
          AList[[(iSub-1)*2 + iSession]] <- A;
        }
      }
    } else {
      load(paste0("../../Data/", dataName, ".Rbin"))
      n <- nrow(fibergraph.list[[1]])
      M <- length(fibergraph.list)
      AList <- list()
      for (i in 1:length(fibergraph.list)) {
        if (weighted) {
          AList[[i]] <- as.matrix(fibergraph.list[[i]])
        } else {
          AList[[i]] <- as.matrix((fibergraph.list[[i]] > 0)*1)
        }
        if (dataName == "migrain") {
          AList[[i]] <- AList[[i]] + t(AList[[i]])
        }
      }
    }
    save(AList, n, M, file=fileName)
    return(list(AList, n, M))
  }
}



DiagAug <- function(A, isSVD = 0, method = 1, param = 3) {
  require(Matrix)
  n <- dim(A)[1]
  D0 <- Diagonal(n, x = rowSums(A)/(n-1))
  d <- DimSelect(A + D0, isSVD, method, param)
  P0 <- LR(A + D0, d, isSVD)
  D1 <- Diagonal(n, x = diag(P0))
  return(A + D1)
}


LR_Estimate <- function(A, weighted = 1, isSVD = 0, method = 1, param = 3) {
  require(Matrix)
  n <- dim(A)[1]
  D0 <- Diagonal(n, x = rowSums(A)/(n-1))
  d <- DimSelect(A + D0, isSVD, method, param)
  P0 <- LR(A + D0, d, isSVD)
  D1 <- Diagonal(n, x = diag(P0))
  P1 <- LR(A + D1, d)
  return(Regularize(P1, weighted))
}




GetElbows <- function(dat, n = 3, threshold = FALSE, plot = TRUE) {
  if (is.matrix(dat))
    d <- sort(apply(dat,2,sd), decreasing=TRUE)
  else
    d <- sort(dat,decreasing=TRUE)
  
  if (!is.logical(threshold))
    d <- d[d > threshold]
  
  p <- length(d)
  if (p == 0)
    stop(paste("d must have elements that are larger than the threshold ",
               threshold), "!", sep="")
  
  lq <- rep(0.0, p)                     # log likelihood, function of q
  for (q in 1:p) {
    mu1 <- mean(d[1:q])
    mu2 <- mean(d[-(1:q)])              # = NaN when q = p
    sigma2 <- (sum((d[1:q] - mu1)^2) + sum((d[-(1:q)] - mu2)^2)) /
      (p - 1 - (q < p))
    lq[q] <- sum( dnorm(  d[1:q ], mu1, sqrt(sigma2), log=TRUE) ) +
      sum( dnorm(d[-(1:q)], mu2, sqrt(sigma2), log=TRUE) )
  }
  
  q <- which.max(lq)
  if (n > 1 && q < p) {
    q <- c(q, q + GetElbows(d[(q+1):p], n-1, plot=FALSE))
  }
  
  if (plot==TRUE) {
    if (is.matrix(dat)) {
      sdv <- d # apply(dat,2,sd)
      plot(sdv,type="b",xlab="dim",ylab="stdev")
      points(q,sdv[q],col=2,pch=19)
    } else {
      plot(dat, type="b")
      points(q,dat[q],col=2,pch=19)
    }
  }
  return(q)
}


USVT <- function(A, m, minS = 1, c = 0.7){
  if(class(A)=="igraph"){
    A <- get.adjacency(A)
  }
  n <- nrow(A)
  # the threshold
  tau <- c*sqrt(n/m)
  usv <- svd(A)
  s <- sum(usv$d >= tau)
  usv$d <- usv$d[1:s]
  return(usv)
}


DimSelect <- function(A, isSVD = 0, method = 1, param = 3) {
  if (method == 1) {
    # Method 1: Zhu & Ghodsi
    nElbow <- param
    n <- dim(A)[1]
    evalVec <- ASE(A, ceiling(n*3/5), isSVD)[[1]]
    d <- GetElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
  } else {
    # Method 2: USVT
    m <- param
    d <- length(USVT(A, m)$d)
  }
}



ASE <- function(A, dim = NA, isSVD = 0) {
  if (is.na(dim)) {
    dim <- dim(A)[1]
  }
  if (isSVD) {
    if (nrow(A) >= 400) {
      require(irlba)
      A.svd <- irlba(A, nu = dim, nv = dim)
      A.values <- A.svd$d
      A.lvectors <- A.svd$u
      A.rvectors <- A.svd$v
    } else {
      A.svd <- svd(A)
      A.values <- A.svd$d[1:dim]
      A.lvectors <- A.svd$u[, 1:dim]
      A.rvectors <- A.svd$v[, 1:dim]
    }
  } else {
    if (nrow(A) >= 400) {
      require(rARPACK)
      A.eig <- eigs_sym(matrix(A, ncol=dim(A)[1]), dim, which = "LA")
      A.values <- A.eig$values
      A.lvectors <- A.eig$vectors
      A.rvectors <- A.lvectors
    } else {
      A.eig <- eigen(A, symmetric = T)
      A.values <- A.eig$values[1:dim]
      A.lvectors <- A.eig$vectors[, 1:dim]
      A.rvectors <- A.lvectors
    }
  }
  return(list(A.values, A.rvectors, A.lvectors))
}



LR <- function(A, d, isSVD = 0) {
  resultASE <- ASE(A, d, isSVD)
  if (d == 1) {
    LR <- resultASE[[1]]*resultASE[[3]] %*% t(resultASE[[2]])
  } else {
    LR <- resultASE[[3]] %*% diag(resultASE[[1]]) %*% t(resultASE[[2]])
  }
}



Regularize <- function(A, weighted = 1) {
  diag(A) <- 0
  A[A < 0] <- 0
  if (!weighted) {
    A[A > 1] <- 1
  }
  return(A)
}




MLqE_Exp_Fun <- function(x, theta, q) {
  # return ((x - theta)/theta^2*(exp(-x/theta)/theta)^(1-q))
  return ((x - theta)*(exp(-(1-q)*x/theta)))
}

MLqE_Exp_Solver <- function(xVec, q, tol=1e-6) {
  if (q == 1) {
    return(mean(xVec))
  }
  thetaMin <- min(xVec)
  thetaMax <- mean(xVec)
  if (thetaMin == thetaMax) {
    return(thetaMin)
  }
  sumBase <- sum(sapply(xVec, MLqE_Exp_Fun, mean(xVec), q))
  sumTmpMin <- abs(sumBase)
  thetaFinal <- mean(xVec)
  theta <- (thetaMin + thetaMax)/2
  sumTmp <- sum(sapply(xVec, MLqE_Exp_Fun, theta, q))
  if (abs(sumTmp) < sumTmpMin) {
    sumTmpMin <- abs(sumTmp)
    thetaFinal <- theta
  }
  maxIter <- 100
  iter <- 0
  while ((abs(sumTmp) > tol*sumBase) && (iter <= maxIter)) {
    iter <- iter + 1
    if (sumTmp > 0) {
      thetaMin <- theta
    } else {
      thetaMax <- theta
    }
    theta <- (thetaMin + thetaMax)/2
    sumTmp <- sum(sapply(xVec, MLqE_Exp_Fun, theta, q))
    if (is.nan(sumTmp)) {
      return(thetaFinal)
    }
    if (abs(sumTmp) < sumTmpMin) {
      sumTmpMin <- abs(sumTmp)
      thetaFinal <- theta
    }
  }
  return(thetaFinal)
}





ExpRealAllDim <- function(AList, m, q, isSVD = 0, weighted = 1, P = NA, dVec = NA) {
  
  n <- dim(AList[[1]])[1]
  M <- length(AList)
  if (any(is.na(dVec))) {
    dVec <- 1:n
  }
  nD <- length(dVec)
  dMax <- max(dVec)
  
  result <- rep(NaN, 2*nD+10)
  
  if (m < M) {
    sampleVec <- sample.int(M, m)
  } else {
    sampleVec <- 1:M
  }
  A_MLE <- add(AList[sampleVec])
  if (any(is.na(P))) {
    P <- (add(AList) - A_MLE)/(M - m)
    # P <- ASum/M
  }
  
  # MLE
  A_MLE <- A_MLE/m
  result[1] <- (norm(P - A_MLE, "F"))^2/(n*(n-1))
  
  
  
  # MLqE  
  AListTmp <- AList[sampleVec]
  nv <- lower.tri(AListTmp[[1]], T)
  for (i in 1:length(AListTmp)) {
    AListTmp[[i]][nv] <- 0
  }
  ATensor <- array(unlist(AListTmp), dim = c(n, n, m))
  A_MLqE <- apply(ATensor, c(1, 2), MLqE_Exp_Solver, q)
  A_MLqE <- A_MLqE + t(A_MLqE)
  # ATensor <- array(unlist(AList[sampleVec]), dim = c(n, n, m))
  # A_MLqE <- apply(ATensor, c(1, 2), MLqE_Exp_Solver, q)
  result[nD + 2] <- (norm(P - A_MLqE, "F"))^2/(n*(n-1))
  
  
  # MLE_ASE
  require(Matrix)
  D0 <- Diagonal(n, x = rowSums(A_MLE)/(n-1))
  dZG <- DimSelect(A_MLE + D0, isSVD, method = 1)
  dUSVT <- DimSelect(A_MLE + D0, isSVD = 1, method = 2, param = m)
  result[nD*2 + 3] <- dZG
  result[nD*2 + 4] <- dUSVT
  for (iD in 1:nD) {
    d <- dVec[iD]
    P0 <- LR(A_MLE + D0, d, isSVD)
    D1 <- Diagonal(n, x = diag(P0))
    P1 <- LR(A_MLE + D1, d)
    A_MLE_ASE <- Regularize(P1, weighted)
    result[1 + iD] <- (norm(P - A_MLE_ASE, "F"))^2/(n*(n-1))
  }
  
  A_MLE_ASE_ZG <- LR_Estimate(A_MLE, weighted, isSVD, method = 1)
  result[2*nD+7] <- (norm(P - A_MLE_ASE_ZG, "F"))^2/(n*(n-1))
  A_MLE_ASE_USVT <- LR_Estimate(A_MLE, weighted, isSVD = 1, method = 2, param = m)
  result[2*nD+8] <- (norm(P - A_MLE_ASE_USVT, "F"))^2/(n*(n-1))
  
  # MLqE_ASE
  D0 <- Diagonal(n, x = rowSums(A_MLqE)/(n-1))
  dZG <- DimSelect(A_MLqE + D0, isSVD, method = 1)
  dUSVT <- DimSelect(A_MLqE + D0, isSVD = 1, method = 2, param = m)
  result[nD*2 + 5] <- dZG
  result[nD*2 + 6] <- dUSVT
  for (iD in 1:nD) {
    d <- dVec[iD]
    P0 <- LR(A_MLqE + D0, d, isSVD)
    D1 <- Diagonal(n, x = diag(P0))
    P1 <- LR(A_MLqE + D1, d)
    A_MLqE_ASE <- Regularize(P1, weighted)
    result[nD + 2 + iD] <- (norm(P - A_MLqE_ASE, "F"))^2/(n*(n-1))
  }
  
  A_MLqE_ASE_ZG <- LR_Estimate(A_MLqE, weighted, isSVD, method = 1)
  result[2*nD+9] <- (norm(P - A_MLqE_ASE_ZG, "F"))^2/(n*(n-1))
  A_MLqE_ASE_USVT <- LR_Estimate(A_MLqE, weighted, isSVD = 1, method = 2, param = m)
  result[2*nD+10] <- (norm(P - A_MLqE_ASE_USVT, "F"))^2/(n*(n-1))
  
  return(result)
}

















SimTwoData <- function(m, n, tau, P, C1, C2, eps1, eps2, q, d, isSVD = 0) {
  # result = rep(NaN, 4)
  
  # P <- B[tau, tau] + rnorm(n^2, mean = 0, sd = 0)*(B[tau, tau])/5
  
  # ind <- lower.tri(P, 1)
  # diag(P) <- 0
  
  # P1 <- (1-eps1)*P + eps1*C1
  # D0 <- Diagonal(n, x = rowSums(P1)/(n-1))
  # eigenResult <- eigen(P1 + D0)$values
  # plot(eigenResult)
  # 
  # P2 <- (1-eps2)*P + eps2*C2
  # D0 <- Diagonal(n, x = rowSums(P2)/(n-1))
  # eigenResult <- eigen(P2 + D0)$values
  # plot(eigenResult)
  
  
  AList1 <- list()
  AList2 <- list()
  
  dVec <- 1:n
  nD <- length(dVec)
  
  ind <- lower.tri(P, 1)
  for (i in 1:m) {
    contamVec1 <- (runif(n^2) > eps1)
    A <- contamVec1*matrix(rexp(n^2, 1/P), ncol=n) +
      (1 - contamVec1)*matrix(rexp(n^2, 1/C1), ncol=n)
    A[ind] <- 0
    A <- A + t(A)
    AList1[[i]] <- A
    
    contamVec2 <- (runif(n^2) > eps2)
    # contamVec2 <- contamVec1
    A <- contamVec2*matrix(rexp(n^2, 1/P), ncol=n) +
      (1 - contamVec2)*matrix(rexp(n^2, 1/C2), ncol=n)
    A[ind] <- 0
    A <- A + t(A)
    AList2[[i]] <- A
  }
  
  out1 <- ExpRealAllDim(AList1, m, q, isSVD, P = add(AList2)/m)
  out2 <- ExpRealAllDim(AList2, m, q, isSVD, P = add(AList1)/m)
  
  return(c(out1, out2))
}




ExpAllDimCompareLR <- function(AList, m, q, isSVD = 0, dVec = NA) {
  
  weighted <- 1
  
  n <- dim(AList[[1]])[1]
  M <- length(AList)
  if (any(is.na(dVec))) {
    dVec <- 1:n
  }
  nD <- length(dVec)
  dMax <- max(dVec)
  
  result <- rep(NaN, 2*nD+11)
  
  if (m < M) {
    sampleVec <- sample.int(M, m)
  } else {
    sampleVec <- 1:M
  }
  A_MLE <- add(AList[sampleVec])
  P <- (add(AList) - A_MLE)/(M - m)
  
  require(Matrix)
  D0 <- Diagonal(n, x = rowSums(P)/(n-1))
  dZG <- DimSelect(P + D0, isSVD, method = 1)
  d <- dZG
  P0 <- LR(P + D0, d, isSVD)
  D1 <- Diagonal(n, x = diag(P0))
  P1 <- LR(P + D1, d)
  P <- Regularize(P1, weighted)
  result[2*nD+11] <- d
  
  # MLE
  A_MLE <- A_MLE/m
  result[1] <- (norm(P - A_MLE, "F"))^2/(n*(n-1))
  
  
  
  # MLqE  
  AListTmp <- AList[sampleVec]
  nv <- lower.tri(AListTmp[[1]], T)
  for (i in 1:length(AListTmp)) {
    AListTmp[[i]][nv] <- 0
  }
  ATensor <- array(unlist(AListTmp), dim = c(n, n, m))
  A_MLqE <- apply(ATensor, c(1, 2), MLqE_Exp_Solver, q)
  A_MLqE <- A_MLqE + t(A_MLqE)
  # ATensor <- array(unlist(AList[sampleVec]), dim = c(n, n, m))
  # A_MLqE <- apply(ATensor, c(1, 2), MLqE_Exp_Solver, q)
  result[nD + 2] <- (norm(P - A_MLqE, "F"))^2/(n*(n-1))
  
  
  # MLE_ASE
  require(Matrix)
  D0 <- Diagonal(n, x = rowSums(A_MLE)/(n-1))
  dZG <- DimSelect(A_MLE + D0, isSVD, method = 1)
  dUSVT <- DimSelect(A_MLE + D0, isSVD = 1, method = 2, param = m)
  result[nD*2 + 3] <- dZG
  result[nD*2 + 4] <- dUSVT
  for (iD in 1:nD) {
    d <- dVec[iD]
    P0 <- LR(A_MLE + D0, d, isSVD)
    D1 <- Diagonal(n, x = diag(P0))
    P1 <- LR(A_MLE + D1, d)
    A_MLE_ASE <- Regularize(P1, weighted)
    result[1 + iD] <- (norm(P - A_MLE_ASE, "F"))^2/(n*(n-1))
  }
  
  A_MLE_ASE_ZG <- LR_Estimate(A_MLE, weighted, isSVD, method = 1)
  result[2*nD+7] <- (norm(P - A_MLE_ASE_ZG, "F"))^2/(n*(n-1))
  A_MLE_ASE_USVT <- LR_Estimate(A_MLE, weighted, isSVD = 1, method = 2, param = m)
  result[2*nD+8] <- (norm(P - A_MLE_ASE_USVT, "F"))^2/(n*(n-1))
  
  # MLqE_ASE
  D0 <- Diagonal(n, x = rowSums(A_MLqE)/(n-1))
  dZG <- DimSelect(A_MLqE + D0, isSVD, method = 1)
  dUSVT <- DimSelect(A_MLqE + D0, isSVD = 1, method = 2, param = m)
  result[nD*2 + 5] <- dZG
  result[nD*2 + 6] <- dUSVT
  for (iD in 1:nD) {
    d <- dVec[iD]
    P0 <- LR(A_MLqE + D0, d, isSVD)
    D1 <- Diagonal(n, x = diag(P0))
    P1 <- LR(A_MLqE + D1, d)
    A_MLqE_ASE <- Regularize(P1, weighted)
    result[nD + 2 + iD] <- (norm(P - A_MLqE_ASE, "F"))^2/(n*(n-1))
  }
  
  A_MLqE_ASE_ZG <- LR_Estimate(A_MLqE, weighted, isSVD, method = 1)
  result[2*nD+9] <- (norm(P - A_MLqE_ASE_ZG, "F"))^2/(n*(n-1))
  A_MLqE_ASE_USVT <- LR_Estimate(A_MLqE, weighted, isSVD = 1, method = 2, param = m)
  result[2*nD+10] <- (norm(P - A_MLqE_ASE_USVT, "F"))^2/(n*(n-1))
  
  return(result) 
}




GetEvals <- function(A) {
  require(Matrix)
  D0 <- Diagonal(n, x = rowSums(A)/(n-1))
  return(eigen(A + D0)$values)
}

PlotEvals <- function(xVec, label_y_lb=NA, label_y_ub=NA) {
  df <- data.frame(eval=xVec, k=1:n)
  if ((is.na(label_y_lb)) || (is.na(label_y_ub))) {
    gg <- ggplot(df,aes(x=k,y=eval))+
      geom_line()+
      xlab("order in algebraic") + ylab("eigenvalue")
  } else {
    gg <- ggplot(df,aes(x=k,y=eval))+
      geom_line()+
      scale_y_continuous(limits = c(label_y_lb, label_y_ub)) +
      xlab("order in algebraic") + ylab("eigenvalue")
  }
  return(gg)
}
















# grid_arrange_shared_legend2 <- function(plots, nrows = 1, ncols = 2) {
#   library(gridExtra)
#   g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
#   legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
#   lheight <- sum(legend$height)
#   pl  <- lapply(plots, function(x) x + theme(legend.position="none"))
#   tmp <- do.call(arrangeGrob, c(pl, list(ncol=ncols, nrow=nrows)))
#   grid.arrange(tmp, legend, ncol=1, heights = unit.c(unit(1, "npc") - lheight, lheight))
# }
# 
# 
# # Diagonal Augmentation
# diag_aug <- function(A, d=0) {
#   if (d == 0) {
#     require(Matrix)
#     n = dim(A)[1]
#     return(A + Diagonal(n, x=rowSums(A))/(n-1))
#   } else {
#     for (iIter in 1:1) {
#       tmp = ase(A, d, 0)
#       if (d == 1)
#         diag(A) = diag(tmp[[1]] * tmp[[3]] %*% t(tmp[[2]]))
#       else
#         diag(A) = diag(tmp[[3]][,1:d] %*% diag(tmp[[1]][1:d]) %*% t(tmp[[2]][,1:d]))
#     }
#     return(A)
#   }
# }
# 
# 
# 
# 
# # Regularize probability matrix
# regularize <- function(A) {
#   diag(A) = 0
#   #   A[A > 1] = 1
#   A[A < 0] = 0
#   return(A)
# }
# 
# # Regularize probability matrix with truncation
# regularizetruncate <- function(A, A0) {
#   A[A < 0] = 0
#   diag(A) = 0
#   A <- pmin(A, A0)
#   return(A)
# }
# 
# 
# add <- function(x) Reduce("+", x)
# 
# dim_brute <- function(m, n, rho, tau, B, dVec, isSVD=1) {
#   result = rep(NaN, nD+1)
#   
#   require(igraph)
#   A_all = list()
#   for (i in 1:m) {
#     g = sample_sbm(n, B, n*rho, directed=F, loops=F)
#     A = as_adj(g, type="both", sparse=FALSE)
#     A_all[[i]] = A
#   }
#   
#   tau = rep(1:K,n*rho)
#   P = B[tau,tau]
#   diag(P) = 0
#   
#   A_bar = add(A_all)/m
#   result[1] = norm(P - A_bar, "F")/n/(n-1)
#   
#   dMax = max(dVec)
#   nD = length(dVec)
#   
#   A.ase = ase(diag_aug(A_bar), dMax, isSVD)
#   for (iD in 1:nD) {
#     d = dVec[iD]
#     if (d == 1)
#       Ahat = A.ase[[1]][1] * A.ase[[3]][,1:d] %*% t(A.ase[[2]][,1:d])
#     else
#       Ahat <- A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
#     P_hat = regularize(Ahat)
#     result[iD+1] = norm(P - P_hat, "F")/n/(n-1)
#   }
#   
#   return(result)
# }
# 
# 
# dim_brute1 <- function(M, m, dVec, A_all, A_sum, isSVD=1) {
#   result = rep(NaN, nD+1)
#   
#   sampleVec = sample.int(M, m)
#   A_bar = add(A_all[sampleVec])
#   P_bar = (A_sum - A_bar)/(M - m)
#   #   P_bar = A_sum/M
#   A_bar = A_bar/m
#   result[1] = norm(P_bar - A_bar, "F")/n/(n-1)
#   
#   dMax = max(dVec)
#   nD = length(dVec)
#   
#   A.ase = ase(diag_aug(A_bar), dMax, isSVD)
#   for (iD in 1:nD) {
#     d = dVec[iD]
#     if (d == 1)
#       Ahat = A.ase[[1]] * A.ase[[3]] %*% t(A.ase[[2]])
#     else
#       Ahat <- A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
#     P_hat = regularize(Ahat)
#     result[iD+1] = norm(P_bar - P_hat, "F")/n/(n-1)
#   }
#   
#   return(result)
# }
# 
# 
# dim_brute2 <- function(M, m, dVec, A_all, A_sum, isSVD=1) {
#   source("getElbows.R")
#   source("USVT.R")
#   
#   result = rep(NaN, nD+1)
#   
#   sampleVec = sample.int(M, m)
#   A_bar = add(A_all[sampleVec])
#   P_bar = (A_sum - A_bar)/(M - m)
#   #   P_bar = A_sum/M
#   A_bar = A_bar/m
#   result[1] = norm(P_bar - A_bar, "F")/n/(n-1)
#   
#   dMax = max(dVec)
#   nD = length(dVec)
#   
#   A_bar_diag_aug = diag_aug(A_bar)
#   
#   # ZG
#   nElbow = 3
#   evalVec = ase(A_bar_diag_aug, ceiling(n*3/5), isSVD)[[1]]
#   dZG = getElbows(evalVec, n=nElbow, plot=F)[[nElbow]]
#   
#   # USVT
#   dUSVT = length(usvt(A_bar_diag_aug, 1, m)$d)
#   
#   A.ase = ase(diag_aug(A_bar), dMax, isSVD)
#   for (iD in 1:nD) {
#     d = dVec[iD]
#     if (d == 1)
#       Ahat = A.ase[[1]] * A.ase[[3]] %*% t(A.ase[[2]])
#     else
#       Ahat <- A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
#     P_hat = regularize(Ahat)
#     result[iD+1] = norm(P_bar - P_hat, "F")/n/(n-1)
#   }
#   
#   result[nD+2] = dZG
#   result[nD+3] = dUSVT
#   
#   return(result)
# }
# 
# 
# dim_brute2_all <- function(M, m, dVec, A_all, A_sum, q, isSVD=1) {
#   source("getElbows.R")
#   source("USVT.R")
#   
#   nD = length(dVec)
#   dMax = max(dVec)
#   result = rep(NaN, 2*nD+6)
#   
#   sampleVec = sample.int(M, m)
#   A_bar = add(A_all[sampleVec])
#   P_bar = (A_sum - A_bar)/(M - m)
#   #   P_bar = A_sum/M
#   A_bar = A_bar/m
#   result[1] = norm(P_bar - A_bar, "F")/n/(n-1)
#   
#   n = dim(A_bar)[[1]]
#   A_all_unlist = array(unlist(A_all[sampleVec]), dim=c(n, n, m))
#   time0 = proc.time()
#   A_mlqe = apply(A_all_unlist, c(1, 2), mlqe_exp_solver, q)
#   time1 = proc.time() - time0
#   result[nD+2] = norm(P_bar - A_mlqe, "F")/n/(n-1)
#   
#   A_bar_diag_aug = diag_aug(A_bar)
#   # ZG
#   nElbow = 3
#   evalVec = ase(A_bar_diag_aug, ceiling(n*3/5), isSVD)[[1]]
#   dZG = getElbows(evalVec, n=nElbow, plot=F)[[nElbow]]
#   # USVT
#   dUSVT = length(usvt(A_bar_diag_aug, 1, m)$d)
#   result[nD*2+3] = dZG
#   result[nD*2+4] = dUSVT
#   
#   A.ase = ase(A_bar_diag_aug, dMax, isSVD)
#   for (iD in 1:nD) {
#     d = dVec[iD]
#     if (d == 1)
#       Ahat = A.ase[[1]] * A.ase[[3]] %*% t(A.ase[[2]])
#     else
#       Ahat <- A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
#     P_hat = regularize(Ahat)
#     result[1+iD] = norm(P_bar - P_hat, "F")/n/(n-1)
#   }
#   #   plot(result[[2]])
#   
#   A_mlqe_diag_aug = diag_aug(A_mlqe)
#   # ZG
#   nElbow = 3
#   evalVec = ase(A_mlqe_diag_aug, ceiling(n*3/5), isSVD)[[1]]
#   dZG = getElbows(evalVec, n=nElbow, plot=F)[[nElbow]]
#   # USVT
#   dUSVT = length(usvt(A_mlqe_diag_aug, 1, m)$d)
#   result[nD*2+5] = dZG
#   result[nD*2+6] = dUSVT
#   
#   A.ase = ase(A_mlqe_diag_aug, dMax, isSVD)
#   for (iD in 1:nD) {
#     d = dVec[iD]
#     if (d == 1)
#       Ahat = A.ase[[1]] * A.ase[[3]] %*% t(A.ase[[2]])
#     else
#       Ahat <- A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
#     P_hat_ase = regularize(Ahat)
#     result[nD+2+iD] = norm(P_bar - P_hat_ase, "F")/n/(n-1)
#   }
#   #   plot(result[[4]])
#   return(result)
# }
# 
# 
# dim_brute_robust <- function(M, m, dVec, A_all, A_all_unlist, isSVD=1) {
#   result = rep(NaN, nD+1)
#   
#   sampleVec = sample.int(M, m)
#   A_bar_robust = apply(A_all_unlist[,,sampleVec], 1:2, median)
#   A_bar_robust = medianlist(A_all[sampleVec])
#   P_bar = (A_sum - A_bar)/(M - m)
#   #   P_bar = A_sum/M
#   A_bar = A_bar/m
#   result[1] = norm(P_bar - A_bar, "F")/n/(n-1)
#   
#   dMax = max(dVec)
#   nD = length(dVec)
#   
#   A.ase = ase(diag_aug(A_bar), dMax, isSVD)
#   for (iD in 1:nD) {
#     d = dVec[iD]
#     if (d == 1)
#       Ahat = A.ase[[1]] * A.ase[[3]] %*% t(A.ase[[2]])
#     else
#       Ahat <- A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
#     P_hat = regularize(Ahat)
#     result[iD+1] = norm(P_bar - P_hat, "F")/n/(n-1)
#   }
#   
#   return(result)
# }
# 
# 
# 
# # Using ZG to choose dimension
# llg_ZG <- function(M, m, A_all, A_sum, isSVD=1) {
#   
#   sampleVec = sample.int(M, m)
#   A_bar = add(A_all[sampleVec])
#   P_bar = (A_sum - A_bar)/(M - m);
#   A_bar = A_bar/m;
#   
#   P_hat = regularize(ase.Ahat(diag_aug(A_bar), d, isSVD))
#   
#   return(c(norm(P_bar - A_bar, "F")/n/(n-1), norm(P_bar - P_hat, "F")/n/(n-1)), d)
# }
# 
# 
# 
# 
# llg_d <- function(M, m, A_all, A_sum, d, isSVD=1) {
#   
#   sampleVec = sample.int(M, m)
#   A_bar = add(A_all[sampleVec])
#   P_bar = (A_sum - A_bar)/(M - m);
#   A_bar = A_bar/m;
#   
#   P_hat = regularize(ase.Ahat(diag_aug(A_bar), d, isSVD))
#   
#   return(c(norm(P_bar - A_bar, "F")/n/(n-1), norm(P_bar - P_hat, "F")/n/(n-1), d))
# }
# 
# 
# 
# # ASE using SVD or eigen-decomposition.
# ase <- function(A, dim, isSVD=1){
#   if (isSVD) {
#     if(nrow(A) >= 400){
#       require(irlba)
#       A.svd = irlba(A, nu = dim, nv = dim)
#       A.values = A.svd$d
#       A.lvectors = A.svd$u
#       A.rvectors = A.svd$v
#     } else{
#       A.svd = svd(A)
#       A.values = A.svd$d[1:dim]
#       A.lvectors = A.svd$u[,1:dim]
#       A.rvectors = A.svd$v[,1:dim]
#     }
#   } else {
#     if(nrow(A) >= 400){
#       require(rARPACK)
#       A.eig = eigs_sym(matrix(A, ncol=dim(A)[1]), dim, which = "LA")
#       A.values = A.eig$values
#       A.lvectors = A.eig$vectors
#       A.rvectors = A.lvectors
#     } else{
#       A.eig = eigen(A, symmetric = T)
#       A.values = A.eig$values[1:dim]
#       A.lvectors = A.eig$vectors[,1:dim]
#       A.rvectors = A.lvectors
#     }
#   }
#   return(list(A.values, A.rvectors, A.lvectors))
# }
# 
# 
# # ASE return xhat
# ase.x <- function(A, dim, isSVD=1){
#   A.ase = ase(A, dim, isSVD)
#   if(dim == 1)
#     A.x = sqrt(A.ase[[1]]) * A.ase[[2]]
#   else
#     A.x <- A.ase[[2]] %*% diag(sqrt(A.ase[[1]]))
#   return(A.x)
# }
# 
# 
# # ASE return Ahat
# ase.Ahat <- function(A, dim, isSVD=1){
#   A.ase = ase(A, dim, isSVD)
#   if(dim == 1)
#     Ahat = A.ase[[1]] * A.ase[[3]] %*% t(A.ase[[2]])
#   else
#     Ahat <- A.ase[[3]] %*% diag(A.ase[[1]]) %*% t(A.ase[[2]])
#   return(Ahat)
# }
# 
# 
# 
# multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
#   library(grid)
#   
#   # Make a list from the ... arguments and plotlist
#   plots <- c(list(...), plotlist)
#   
#   numPlots = length(plots)
#   
#   # If layout is NULL, then use 'cols' to determine layout
#   if (is.null(layout)) {
#     # Make the panel
#     # ncol: Number of columns of plots
#     # nrow: Number of rows needed, calculated from # of cols
#     layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
#                      ncol = cols, nrow = ceiling(numPlots/cols))
#   }
#   
#   if (numPlots==1) {
#     print(plots[[1]])
#     
#   } else {
#     # Set up the page
#     grid.newpage()
#     pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
#     
#     # Make each plot, in the correct location
#     for (i in 1:numPlots) {
#       # Get the i,j matrix positions of the regions that contain this subplot
#       matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
#       
#       print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
#                                       layout.pos.col = matchidx$col))
#     }
#   }
# }
# 
# 
# mlqe_exp_fun <- function(x, theta, q) {
#   return ((x - theta)/theta^2*(exp(-x/theta)/theta)^(1-q))
# }
# 
# mlqe_exp_solver <- function(xVec, q, tol=1e-6) {
#   if (q == 1) {
#     return(mean(xVec))
#   }
#   thetaMin = min(xVec)
#   #   thetaMax = max(xVec)
#   thetaMax = mean(xVec)
#   if (thetaMin == thetaMax) {
#     return(thetaMin)
#   }
#   theta = (thetaMin + thetaMax)/2
#   sumTmp = sum(sapply(xVec, mlqe_exp_fun, theta, q))
#   sumTmp0 = sumTmp
#   sumBase = abs(sumTmp)
#   maxIter = 100
#   iter = 0
#   while (abs(sumTmp) > tol*sumBase) {
#     iter = iter + 1
#     if (sumTmp > 0) {
#       thetaMin = theta
#     } else {
#       thetaMax = theta
#     }
#     theta = (thetaMin + thetaMax)/2
#     sumTmp = sum(sapply(xVec, mlqe_exp_fun, theta, q))
#     if (is.nan(sumTmp)) {
#       return(mean(xVec))
#     }
#     if ((iter >= maxIter) && (abs(sumTmp) > sumTmp0)) {
#       return(mean(xVec))
#     }
#     if (iter >= maxIter) {
#       return(theta)
#     }
#     #     sumTmp1 = sum(sapply(xVec, mlqe_exp_fun, (thetaMin + theta)/2, q))
#     #     sumTmp2 = sum(sapply(xVec, mlqe_exp_fun, (thetaMax + theta)/2, q))
#     #     if ((abs(sumTmp1) < abs(sumTmp)) && (abs(sumTmp2) < abs(sumTmp))) {
#     #       if (sumTmp > 0) {
#     #         thetaMin = theta
#     #       } else {
#     #         thetaMax = theta
#     #       }
#     #     } else if (abs(sumTmp1) < abs(sumTmp)) {
#     #       thetaMax = theta
#     #     } else {
#     #       thetaMin = theta
#     #     }
#     #     theta = (thetaMin + thetaMax)/2
#     #     sumTmp = sum(sapply(xVec, mlqe_exp_fun, theta, q))
#   }
#   return(theta)
# }
# 
# sim_all <- function(m, n, tau, B, CB, eps, q, d, isSVD=1) {
#   result = rep(NaN, 4)
#   
#   P = B[tau, tau]
#   C = CB[tau, tau]
#   diag(P) = 0
#   diag(C) = 0
#   
#   A_all = array(0, dim=c(m, n, n))
#   ind = lower.tri(A_all[1,,], 1)
#   #   A_all1 = list()
#   for (i in 1:m) {
#     contamVec = (runif(n^2) > eps)
#     A = contamVec*matrix(rexp(n^2, 1/P), ncol=n) +
#       (1-contamVec)*matrix(rexp(n^2, 1/C), ncol=n)
#     A[ind] = 0
#     A = A + t(A)
#     A_all[i,,] = A
#     #     A_all1[[i]] = A
#   }
#   
#   A_bar = apply(A_all, MARGIN=c(2, 3), sum)/m
#   #   A_bar = add(A_all1)/m
#   A.ase = ase(diag_aug(A_bar, d), d, isSVD)
#   if (d == 1) {
#     A_bar_ase = A.ase[[1]][1] * A.ase[[3]][,1:d] %*% t(A.ase[[2]][,1:d])
#   } else {
#     A_bar_ase = A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
#   }
#   A_bar_ase = regularize(A_bar_ase)
#   
#   #   tmp = apply(A_all, 1, mlqe_exp_solver, q)
#   #   A_mlqe = array(0, dim=c(n, n))
#   #   for (i in 1:(n-1)) {
#   #     print(i)
#   #     for (j in (i+1):n) {
#   #       A_mlqe[i,j] = mlqe_exp_solver(A_all[,i,j], q)
#   #       A_mlqe[j,i] = A_mlqe[i,j]
#   #     }
#   #   }
#   A_mlqe = apply(A_all, c(2, 3), mlqe_exp_solver, q)
#   
#   A.ase = ase(diag_aug(A_mlqe, d), d, isSVD)
#   if (d == 1) {
#     A_mlqe_ase = A.ase[[1]][1] * A.ase[[3]][,1:d] %*% t(A.ase[[2]][,1:d])
#   } else {
#     A_mlqe_ase = A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
#   }
#   A_mlqe_ase = regularize(A_mlqe_ase)
#   
#   result[1] = norm(A_bar - P, "F")/n/(n-1)
#   result[2] = norm(A_bar_ase - P, "F")/n/(n-1)
#   result[3] = norm(A_mlqe - P, "F")/n/(n-1)
#   result[4] = norm(A_mlqe_ase - P, "F")/n/(n-1)
#   
#   return(result)
# }
# 
# 
# 
# dim_brute_fullrank <- function(m, dVec, P, isSVD=1, contamVec=1) {
#   
#   dMax = max(dVec)
#   nD = length(dVec)
#   result = rep(NaN, nD+1)
#   
#   require(igraph)
#   A_all = list()
#   for (i in 1:m) {
#     A = contamVec*matrix(rexp(n^2, 1/P), ncol=n) +
#       (1-contamVec)*matrix(rexp(n^2, 1/P/10), ncol=n)
#     ind = lower.tri(A, 1)
#     A[ind] = 0
#     A = A + t(A)
#     A_all[[i]] = A
#   }
#   
#   A_bar = add(A_all)/m
#   result[1] = norm(P - A_bar, "F")/n/(n-1)
#   
#   A.ase = ase(diag_aug(A_bar), dMax, isSVD)
#   
#   for (iD in 1:nD) {
#     d = dVec[iD]
#     #     A.ase = ase(diag_aug(A_bar, d), d, isSVD)
#     if (d == 1)
#       Ahat = A.ase[[1]][1] * A.ase[[3]] %*% t(A.ase[[2]])
#     else
#       Ahat <- A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
#     P_hat = regularize(Ahat)
#     result[iD+1] = norm(P - P_hat, "F")/n/(n-1)
#   }
#   
#   return(result)
# }
# 
# 
# 
# test_rank1 <- function(m, dVec, P, isSVD=1) {
#   
#   dMax = max(dVec)
#   nD = length(dVec)
#   result = rep(NaN, nD)
#   
#   require(igraph)
#   A_all = list()
#   for (i in 1:m) {
#     g = sample_sbm(n, P, rep(1,n), directed=F, loops=F)
#     A = as_adj(g, type="both", sparse=FALSE)
#     A_all[[i]] = A
#   }
#   
#   A_bar = add(A_all)/m
#   
#   A.ase = eigen(diag_aug(A_bar), dMax)
#   P.ase = eigen(diag_aug(P), dMax)
#   
#   for (iD in 1:nD) {
#     d = dVec[iD]
#     A_hat = A.ase$values[d] * A.ase$vectors[,d] %*% t(A.ase$vectors[,d])
#     P_hat = P.ase$values[d] * P.ase$vectors[,d] %*% t(P.ase$vectors[,d])
#     result[iD] = norm(A_hat - P_hat, "F")/n/(n-1)
#   }
#   
#   return(result)
# }
# 
# 
# 
# 
# 
# MLqEObjLognormal <- function(theta, data, q) {
#   mu <- theta[1]
#   sigma <- theta[2]^2
#   fVal <- 1/sigma/sqrt(2*pi)/data*exp(-(log(data) - mu)^2/2/(sigma^2))
#   fqVal <- fVal^(1 - q)
#   val <- sum(fqVal*(log(data) - mu))^2 +
#     sum(fqVal*(sigma - (log(data) - mu)^2/(sigma)))^2
#   return(val)
# }
# 
# MLqEObjLognormal1 <- function(theta, data, q) {
#   mu <- theta[1]
#   sigma <- theta[2]^2
#   fVal <- 1/sigma/sqrt(2*pi)/data*exp(-(log(data) - mu)^2/2/(sigma^2))
#   if (q == 1) {
#     val <- sum(log(fVal))
#   } else {
#     val <- sum((fVal^(1 - q) - 1)/(1 - q))
#   }
#   return(-val)
# }
# 
# 
# MLqESolverLognormal <- function(data, q) {
#   muHat <- mean(log(data))
#   sigmaHat <- sqrt(mean((log(data) - muHat)^2))
#   # To avoid the constraint that sigma > 0, we optimize over sqrt(sigma)
#   thetaInit <- c(muHat, sqrt(sigmaHat))
#   if (sigmaHat > 0) {
#     resultOpti <- optim(par = thetaInit, fn = MLqEObjLognormal1, data = data, q = q,
#                         method = "BFGS")
#     muHat <- resultOpti$par[1]
#     sigmaHat <- resultOpti$par[2]^2
#   }
#   val <- exp(muHat + sigmaHat^2/2)
#   return(val)
# }
# 
# 
# 
# LognormalAllDim <- function(M, m, dVec, AList, ASum, q, isSVD=1) {
#   source("getElbows.R")
#   source("USVT.R")
#   
#   nD <- length(dVec)
#   dMax <- max(dVec)
#   result <- rep(NaN, 2*nD+6)
#   
#   sampleVec <- sample.int(M, m)
#   ABar <- add(AList[sampleVec])
#   PBar <- (ASum - ABar)/(M - m)
#   #   PBar <- ASum/M
#   ABar <- ABar/m
#   result[1] <- (norm(PBar - ABar, "F"))^2/n/(n-1)
#   
#   n <- dim(ABar)[[1]]
#   ATensor <- array(unlist(AList[sampleVec]), dim = c(n, n, m))
#   AMLqE <- apply(ATensor, c(1, 2), MLqESolverLognormal, q)
#   result[nD + 2] <- (norm(PBar - AMLqE, "F"))^2/n/(n-1)
#   
#   A_bar_diag_aug = diag_aug(ABar)
#   # ZG
#   nElbow = 3
#   evalVec = ase(A_bar_diag_aug, ceiling(n*3/5), isSVD)[[1]]
#   dZG = getElbows(evalVec, n=nElbow, plot=F)[[nElbow]]
#   # USVT
#   dUSVT = length(usvt(A_bar_diag_aug, 1, m)$d)
#   result[nD*2+3] = dZG
#   result[nD*2+4] = dUSVT
#   
#   A.ase = ase(A_bar_diag_aug, dMax, isSVD)
#   for (iD in 1:nD) {
#     d = dVec[iD]
#     if (d == 1)
#       Ahat = A.ase[[1]] * A.ase[[3]] %*% t(A.ase[[2]])
#     else
#       Ahat <- A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
#     P_hat = regularize(Ahat)
#     result[1+iD] = norm(PBar - P_hat, "F")/n/(n-1)
#   }
#   #   plot(result[[2]])
#   
#   A_mlqe_diag_aug = diag_aug(AMLqE)
#   # ZG
#   nElbow = 3
#   evalVec = ase(A_mlqe_diag_aug, ceiling(n*3/5), isSVD)[[1]]
#   dZG = getElbows(evalVec, n=nElbow, plot=F)[[nElbow]]
#   # USVT
#   dUSVT = length(usvt(A_mlqe_diag_aug, 1, m)$d)
#   result[nD*2+5] = dZG
#   result[nD*2+6] = dUSVT
#   
#   A.ase = ase(A_mlqe_diag_aug, dMax, isSVD)
#   for (iD in 1:nD) {
#     d = dVec[iD]
#     if (d == 1)
#       Ahat = A.ase[[1]] * A.ase[[3]] %*% t(A.ase[[2]])
#     else
#       Ahat <- A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
#     P_hat_ase = regularize(Ahat)
#     result[nD+2+iD] = norm(PBar - P_hat_ase, "F")/n/(n-1)
#   }
#   #   plot(result[[4]])
#   return(result)
# }
# 
# 
# 
# MLqEObjExp <- function(theta, data, q) {
#   fVal <- 1/theta*exp(-data/theta)
#   if (q == 1) {
#     val <- sum(log(fVal))
#   } else {
#     val <- sum((fVal^(1 - q) - 1)/(1 - q))
#   }
#   return(-val)
# }
# 
# MLqEGradExp <- function(theta, data, q) {
#   if (q == 1) {
#     val <- sum(-1/theta + data/(theta^2))
#   } else {
#     fVal <- 1/theta*exp(-data/theta)
#     val <- sum((fVal^(1 - q)*(-1/theta + data/(theta^2))))
#   }
#   return(val)
# }
# 
# 
# MLqESolverExp <- function(data, q) {
#   thetaHat <- mean(data)
#   thetaInit <- thetaHat
#   if (var(data) > 0) {
#     #     resultOpti <- optim(par = thetaInit, fn = MLqEObjExp, data = data, q = q,
#     #                         gr = MLqEGradExp, method = "BFGS")
#     #     resultOpti <- optim(par = thetaInit, fn = MLqEObjExp, data = data, q = q,
#     #                         method = "BFGS")
#     resultOpti <- optim(par = thetaInit, fn = MLqEObjExp, data = data, q = q)
#     thetaHat <- resultOpti$par
#   }
#   return(thetaHat)
# }
# 
# 
# ExpAllDim <- function(M, m, dVec, AList, ASum, q, isSVD=1, PBar=NA) {
#   source("getElbows.R")
#   source("USVT.R")
#   
#   nD <- length(dVec)
#   dMax <- max(dVec)
#   result <- rep(NaN, 2*nD+6)
#   
#   sampleVec <- sample.int(M, m)
#   ABar <- add(AList[sampleVec])
#   if (any(is.na(PBar))) {
#     PBar <- (ASum - ABar)/(M - m)
#   }
#   #   PBar <- ASum/M
#   ABar <- ABar/m
#   result[1] <- (norm(PBar - ABar, "F"))^2/n/(n-1)
#   
#   n <- dim(ABar)[[1]]
#   ATensor <- array(unlist(AList[sampleVec]), dim = c(n, n, m))
#   #   AMLqE <- apply(ATensor, c(1, 2), MLqESolverExp, q)  
#   AMLqE <- apply(ATensor, c(1, 2), mlqe_exp_solver, q)
#   result[nD + 2] <- (norm(PBar - AMLqE, "F"))^2/n/(n-1)
#   
#   ABarDiagAug <- diag_aug(ABar)
#   # ABarDiagAug <- ABar
#   # ZG
#   nElbow <- 3
#   evalVec <- ase(ABarDiagAug, ceiling(n*3/5), isSVD)[[1]]
#   dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#   # USVT
#   dUSVT <- length(usvt(ABarDiagAug, 1, m)$d)
#   result[nD*2 + 3] <- dZG
#   result[nD*2 + 4] <- dUSVT
#   
#   AASE = ase(ABarDiagAug, dMax, isSVD)
#   for (iD in 1:nD) {
#     d <- dVec[iD]
#     if (d == 1) {
#       AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
#     } else {
#       AHat <- AASE[[3]][, 1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][ ,1:d])
#     }
#     PHat <- regularize(AHat)
#     result[1 + iD] <- (norm(PBar - PHat, "F"))^2/n/(n-1)
#   }
#   
#   # AASE1 = ase(ABarDiagAug, dMax, 1)
#   # AASE0 = ase(ABarDiagAug, dMax, 0)
#   # for (iD in 1:nD) {
#   #   d <- dVec[iD]
#   #   if (d == 1) {
#   #     AHat1 <- AASE1[[1]]*AASE1[[3]]%*%t(AASE1[[2]])
#   #     AHat0 <- AASE0[[1]]*AASE0[[3]]%*%t(AASE0[[2]])
#   #   } else {
#   #     AHat1 <- AASE1[[3]][, 1:d]%*%diag(AASE1[[1]][1:d])%*%t(AASE1[[2]][ ,1:d])
#   #     AHat0 <- AASE0[[3]][, 1:d]%*%diag(AASE0[[1]][1:d])%*%t(AASE0[[2]][ ,1:d])
#   #   }
#   #   
#   #   (norm(ABarDiagAug - AHat1, "F"))^2/n/(n-1)
#   #   (norm(ABarDiagAug - AHat0, "F"))^2/n/(n-1)
#   #   
#   #   PHat1 <- regularize(AHat1)
#   #   PHat0 <- regularize(AHat0)
#   #   
#   #   (norm(ABarDiagAug - PHat1, "F"))^2/n/(n-1)
#   #   (norm(ABarDiagAug - PHat0, "F"))^2/n/(n-1)
#   #   
#   #   (norm(PBar - PHat1, "F"))^2/n/(n-1)
#   #   (norm(PBar - PHat0, "F"))^2/n/(n-1)
#   #   (norm(PBar, "F"))^2/n/(n-1)
#   #   (norm(ABarDiagAug, "F"))^2/n/(n-1)
#   # }
#   
#   
#   # (norm(PBar, "F"))^2/n/(n-1)
#   # (norm(ABar, "F"))^2/n/(n-1)
#   # (norm(AHat, "F"))^2/n/(n-1)
#   # AHat1 <- AHat
#   # diag(AHat1) <- 0
#   # (norm(AHat, "F"))^2/n/(n-1)
#   # (norm(PHat, "F"))^2/n/(n-1)
#   
#   
#   
#   AMLqEDiagAug <- diag_aug(AMLqE)
#   # AMLqEDiagAug <- AMLqE
#   # ZG
#   nElbow <- 3
#   evalVec <- ase(AMLqEDiagAug, ceiling(n*3/5), isSVD)[[1]]
#   dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#   # USVT
#   dUSVT <- length(usvt(AMLqEDiagAug, 1, m)$d)
#   result[nD*2 + 5] <- dZG
#   result[nD*2 + 6] <- dUSVT
#   
#   AASE <- ase(AMLqEDiagAug, dMax, isSVD)
#   for (iD in 1:nD) {
#     d <- dVec[iD]
#     if (d == 1) {
#       AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
#     } else {
#       AHat <- AASE[[3]][ ,1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][ ,1:d])
#     }
#     PHatASE <- regularize(AHat)
#     result[nD + 2 + iD] <- (norm(PBar - PHatASE, "F"))^2/n/(n-1)
#   }
#   return(result)
# }
# 
# 
# 
# 
# 
# ExpAllDimAug <- function(M, m, dVec, AList, ASum, q, isSVD=1) {
#   source("getElbows.R")
#   source("USVT.R")
#   
#   nD <- length(dVec)
#   dMax <- max(dVec)
#   result <- rep(NaN, 2*nD+6)
#   
#   sampleVec <- sample.int(M, m)
#   ABar <- add(AList[sampleVec])
#   PBar <- (ASum - ABar)/(M - m)
#   #   PBar <- ASum/M
#   ABar <- ABar/m
#   tmpDiff <- PBar - ABar
#   diag(tmpDiff) <- 0
#   result[1] <- (norm(tmpDiff, "F"))^2/n/(n-1)
#   
#   n <- dim(ABar)[[1]]
#   ATensor <- array(unlist(AList[sampleVec]), dim = c(n, n, m))
#   #   AMLqE <- apply(ATensor, c(1, 2), MLqESolverExp, q)  
#   AMLqE <- apply(ATensor, c(1, 2), mlqe_exp_solver, q)
#   tmpDiff <- PBar - AMLqE
#   diag(tmpDiff) <- 0
#   result[nD + 2] <- (norm(tmpDiff, "F"))^2/n/(n-1)
#   
#   #   ABarDiagAug <- diag_aug(ABar)
#   ABarDiagAug <- ABar
#   # ZG
#   nElbow <- 3
#   evalVec <- ase(ABarDiagAug, ceiling(n*3/5), isSVD)[[1]]
#   dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#   # USVT
#   dUSVT <- length(usvt(ABarDiagAug, 1, m)$d)
#   result[nD*2 + 3] <- dZG
#   result[nD*2 + 4] <- dUSVT
#   
#   AASE = ase(ABarDiagAug, dMax, isSVD)
#   for (iD in 1:nD) {
#     d <- dVec[iD]
#     if (d == 1) {
#       AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
#     } else {
#       AHat <- AASE[[3]][, 1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][ ,1:d])
#     }
#     PHat <- regularize(AHat)
#     tmpDiff <- PBar - PHat
#     diag(tmpDiff) <- 0
#     result[1 + iD] <- (norm(tmpDiff, "F"))^2/n/(n-1)
#   }
#   
#   #   AMLqEDiagAug <- diag_aug(AMLqE)
#   AMLqEDiagAug <- AMLqE
#   # ZG
#   nElbow <- 3
#   evalVec <- ase(AMLqEDiagAug, ceiling(n*3/5), isSVD)[[1]]
#   dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#   # USVT
#   dUSVT <- length(usvt(AMLqEDiagAug, 1, m)$d)
#   result[nD*2 + 5] <- dZG
#   result[nD*2 + 6] <- dUSVT
#   
#   AASE <- ase(AMLqEDiagAug, dMax, isSVD)
#   for (iD in 1:nD) {
#     d <- dVec[iD]
#     if (d == 1) {
#       AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
#     } else {
#       AHat <- AASE[[3]][ ,1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][ ,1:d])
#     }
#     PHatASE <- regularize(AHat)
#     tmpDiff <- PBar - PHatASE
#     diag(tmpDiff) <- 0
#     result[nD + 2 + iD] <- (norm(tmpDiff, "F"))^2/n/(n-1)
#   }
#   return(result)
# }
# 
# 
# 
# 
# 
# 
# 
# 
# ExpAllDimRatio <- function(M, m, d = 0, AList, ASum, q, isSVD = 1) {
#   source("getElbows.R")
#   
#   dMax <- max(dVec)
#   result <- rep(NaN, 2*nD+6)
#   
#   sampleVec <- sample.int(M, m)
#   sampleVecComplement <- (1:M)[is.na(pmatch(1:M, sampleVec))]
#   ABar <- add(AList[sampleVec])/m
#   
#   n <- dim(ABar)[[1]]
#   ATensor <- array(unlist(AList[sampleVec]), dim = c(n, n, m))
#   #   AMLqE <- apply(ATensor, c(1, 2), MLqESolverExp, q)  
#   AMLqE <- apply(ATensor, c(1, 2), mlqe_exp_solver, q)
#   
#   ABarDiagAug <- diag_aug(ABar)
#   if (d == 0) {
#     # ZG
#     nElbow <- 3
#     evalVec <- ase(ABarDiagAug, ceiling(n*3/5), isSVD)[[1]]
#     d <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#   }
#   
#   AASE <- ase(ABarDiagAug, d, isSVD)
#   if (d == 1) {
#     AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
#   } else {
#     AHat <- AASE[[3]][, 1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][, 1:d])
#   }
#   ABarASE <- regularize(AHat)
#   
#   AMLqEDiagAug <- diag_aug(AMLqE)
#   if (d == 0) {
#     # ZG
#     nElbow <- 3
#     evalVec <- ase(AMLqEDiagAug, ceiling(n*3/5), isSVD)[[1]]
#     d <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#   }
#   
#   AASE <- ase(AMLqEDiagAug, d, isSVD)
#   if (d == 1) {
#     AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
#   } else {
#     AHat <- AASE[[3]][, 1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][, 1:d])
#   }
#   PHatASE <- regularize(AHat)
#   
#   ratioList <- lapply(sampleVecComplement, function(i) {
#     return(abs(AList[[i]] - ABar) == abs(AList[[i]] - AMLqE))})
#   ratioAMLqEABar <- add(ratioList)/(M - m)
#   image(Matrix(ratioAMLqEABar))
#   sum(ratioAMLqEABar)/n/(n-1)
#   
#   return(result)
# }
# 
# 
# 
# 
# 
# 
# ExpAllDimSingleAug <- function(M, m, dVec, AList, ASum, AList0, ASum0, q, isSVD=1, PBar=NA) {
#   source("getElbows.R")
#   source("USVT.R")
#   
#   nD <- length(dVec)
#   dMax <- max(dVec)
#   result <- rep(NaN, 4*nD+8)
#   
#   sampleVec <- sample.int(M, m)
#   ABar <- add(AList[sampleVec])
#   ABar0 <- add(AList0[sampleVec])
#   if (any(is.na(PBar))) {
#     PBar <- (ASum - ABar)/(M - m)
#     PBar0 <- (ASum0 - ABar0)/(M - m)
#   }
#   #   PBar <- ASum/M
#   ABar <- ABar/m
#   result[1] <- (norm(PBar - ABar, "F"))^2/n/(n-1)
#   result[2*nD+7] <- (norm(PBar0 - ABar, "F"))^2/n/(n-1)
#   
#   n <- dim(ABar)[[1]]
#   
#   # ATensor <- array(unlist(AList[sampleVec]), dim = c(n, n, m))
#   # AMLqE1 <- apply(ATensor, c(1, 2), mlqe_exp_solver, q)
#   
#   AList1 <- AList[sampleVec]
#   for (i in 1:m) {
#     AList1[[i]][upper.tri(AList1[[i]])] <- 0
#   }
#   ATensor <- array(unlist(AList1), dim = c(n, n, m))
#   AMLqE <- apply(ATensor, c(1, 2), mlqe_exp_solver, q)
#   AMLqE <- AMLqE + t(AMLqE)
#   
#   result[nD + 2] <- (norm(PBar - AMLqE, "F"))^2/n/(n-1)
#   result[3*nD+8] <- (norm(PBar0 - AMLqE, "F"))^2/n/(n-1)
#   
#   
#   ### MLqE ###
#   # ATensor <- array(unlist(AList), dim = c(n, n, m))
#   # AMLqE <- apply(ATensor, c(1, 2), mlqe_exp_solver, q)
#   AList1 <- AList
#   AList1[[1]][upper.tri(AList1[[1]])] <- 0
#   AList1[[2]][upper.tri(AList1[[2]])] <- 0
#   ATensor <- array(unlist(AList1), dim = c(n, n, m))
#   AMLqE <- apply(ATensor, c(1, 2), mlqe_exp_solver, q)
#   AMLqE <- AMLqE + t(AMLqE)
#   
#   
#   
#   
#   
#   ABarDiagAug <- diag_aug(ABar)
#   # ABarDiagAug <- ABar
#   # ZG
#   nElbow <- 3
#   evalVec <- ase(ABarDiagAug, ceiling(n*3/5), isSVD)[[1]]
#   dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#   # USVT
#   dUSVT <- length(usvt(ABarDiagAug, 1, m)$d)
#   result[nD*2 + 3] <- dZG
#   result[nD*2 + 4] <- dUSVT
#   
#   AASE = ase(ABarDiagAug, dMax, isSVD)
#   for (iD in 1:nD) {
#     d <- dVec[iD]
#     if (d == 1) {
#       AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
#     } else {
#       AHat <- AASE[[3]][, 1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][ ,1:d])
#     }
#     PHat <- regularize(AHat)
#     result[1 + iD] <- (norm(PBar - PHat, "F"))^2/n/(n-1)
#     result[2*nD + 7 + iD] <- (norm(PBar0 - PHat, "F"))^2/n/(n-1)
#   }
#   
#   AMLqEDiagAug <- diag_aug(AMLqE)
#   # AMLqEDiagAug <- AMLqE
#   # ZG
#   nElbow <- 3
#   evalVec <- ase(AMLqEDiagAug, ceiling(n*3/5), isSVD)[[1]]
#   dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#   # USVT
#   dUSVT <- length(usvt(AMLqEDiagAug, 1, m)$d)
#   result[nD*2 + 5] <- dZG
#   result[nD*2 + 6] <- dUSVT
#   
#   AASE <- ase(AMLqEDiagAug, dMax, isSVD)
#   for (iD in 1:nD) {
#     d <- dVec[iD]
#     if (d == 1) {
#       AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
#     } else {
#       AHat <- AASE[[3]][ ,1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][ ,1:d])
#     }
#     PHatASE <- regularize(AHat)
#     result[nD + 2 + iD] <- (norm(PBar - PHatASE, "F"))^2/n/(n-1)
#     result[3*nD + 8 + iD] <- (norm(PBar0 - PHatASE, "F"))^2/n/(n-1)
#   }
#   return(result)
# }
# 
# 
# 
# ExpTest <- function(M, AList, q, isSVD=1) {
#   source("getElbows.R")
#   source("USVT.R")
#   m <- 2
#   
#   result <- rep(NaN, 4)
#   
#   ### MLE ###
#   ABar <- add(AList)
#   ABar <- ABar/m
#   # result[1] <- (sum((ABar - AList[[1]])^2) + sum((ABar - AList[[2]])^2))/n/(n-1)
#   result[1] <- (2*sum((ABar - AList[[1]])^2))/n/(n-1)
#   
#   ### MLqE ###
#   n <- dim(ABar)[[1]]
#   # ATensor <- array(unlist(AList), dim = c(n, n, m))
#   # AMLqE <- apply(ATensor, c(1, 2), mlqe_exp_solver, q)
#   AList1 <- AList
#   AList1[[1]][upper.tri(AList1[[1]])] <- 0
#   AList1[[2]][upper.tri(AList1[[2]])] <- 0
#   ATensor <- array(unlist(AList1), dim = c(n, n, m))
#   AMLqE <- apply(ATensor, c(1, 2), mlqe_exp_solver, q)
#   AMLqE <- AMLqE + t(AMLqE)
#   result[2] <- (sum((AMLqE - AList[[1]])^2) + sum((AMLqE - AList[[2]])^2))/n/(n-1)
#   
#   ### ASE o MLE ###
#   ABarDiagAug <- diag_aug(ABar)
#   # ZG
#   nElbow <- 3
#   evalVec <- ase(ABarDiagAug, ceiling(n*3/5), isSVD)[[1]]
#   dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#   
#   AASE = ase(ABarDiagAug, dZG, isSVD)
#   if (dZG == 1) {
#     AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
#   } else {
#     AHat <- AASE[[3]]%*%diag(AASE[[1]])%*%t(AASE[[2]])
#   }
#   PHat <- regularize(AHat)
#   result[3] <- (sum((PHat - AList[[1]])^2) + sum((PHat - AList[[2]])^2))/n/(n-1)
#   
#   ### ASE o MLqE ###
#   AMLqEDiagAug <- diag_aug(AMLqE)
#   # ZG
#   nElbow <- 3
#   evalVec <- ase(AMLqEDiagAug, ceiling(n*3/5), isSVD)[[1]]
#   dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#   
#   AASE <- ase(AMLqEDiagAug, dZG, isSVD)
#   if (dZG == 1) {
#     AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
#   } else {
#     AHat <- AASE[[3]]%*%diag(AASE[[1]])%*%t(AASE[[2]])
#   }
#   PHatASE <- regularize(AHat)
#   result[4] <- (sum((PHatASE - AList[[1]])^2) + sum((PHatASE - AList[[2]])^2))/n/(n-1)
#   
#   return(result)
# }
# 
# 
# 
# ExpAllDimTruncate <- function(M, m, dVec, AList, ASum, q, isSVD=1, PBar=NA) {
#   source("getElbows.R")
#   source("USVT.R")
#   
#   nD <- length(dVec)
#   dMax <- max(dVec)
#   result <- rep(NaN, 2*nD+6)
#   
#   sampleVec <- sample.int(M, m)
#   ABar <- add(AList[sampleVec])
#   if (any(is.na(PBar))) {
#     PBar <- (ASum - ABar)/(M - m)
#   }
#   #   PBar <- ASum/M
#   ABar <- ABar/m
#   result[1] <- (norm(PBar - ABar, "F"))^2/n/(n-1)
#   
#   n <- dim(ABar)[[1]]
#   ATensor <- array(unlist(AList[sampleVec]), dim = c(n, n, m))
#   #   AMLqE <- apply(ATensor, c(1, 2), MLqESolverExp, q)  
#   AMLqE <- apply(ATensor, c(1, 2), mlqe_exp_solver, q)
#   result[nD + 2] <- (norm(PBar - AMLqE, "F"))^2/n/(n-1)
#   
#   ABarDiagAug <- diag_aug(ABar)
#   # ABarDiagAug <- ABar
#   # ZG
#   nElbow <- 3
#   evalVec <- ase(ABarDiagAug, ceiling(n*3/5), isSVD)[[1]]
#   dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#   # USVT
#   dUSVT <- length(usvt(ABarDiagAug, 1, m)$d)
#   result[nD*2 + 3] <- dZG
#   result[nD*2 + 4] <- dUSVT
#   
#   AASE = ase(ABarDiagAug, dMax, isSVD)
#   ABarDiagAugMatrix <- as.matrix(ABarDiagAug)
#   for (iD in 1:nD) {
#     d <- dVec[iD]
#     if (d == 1) {
#       AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
#     } else {
#       AHat <- AASE[[3]][, 1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][ ,1:d])
#     }
#     PHat <- regularizetruncate(AHat, ABarDiagAugMatrix)
#     result[1 + iD] <- (norm(PBar - PHat, "F"))^2/n/(n-1)
#   }
#   
#   AMLqEDiagAug <- diag_aug(AMLqE)
#   # AMLqEDiagAug <- AMLqE
#   # ZG
#   nElbow <- 3
#   evalVec <- ase(AMLqEDiagAug, ceiling(n*3/5), isSVD)[[1]]
#   dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#   # USVT
#   dUSVT <- length(usvt(AMLqEDiagAug, 1, m)$d)
#   result[nD*2 + 5] <- dZG
#   result[nD*2 + 6] <- dUSVT
#   
#   AASE <- ase(AMLqEDiagAug, dMax, isSVD)
#   AMLqEDiagAugMatrix <- as.matrix(AMLqEDiagAug)
#   for (iD in 1:nD) {
#     d <- dVec[iD]
#     if (d == 1) {
#       AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
#     } else {
#       AHat <- AASE[[3]][ ,1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][ ,1:d])
#     }
#     PHatASE <- regularizetruncate(AHat, AMLqEDiagAugMatrix)
#     result[nD + 2 + iD] <- (norm(PBar - PHatASE, "F"))^2/n/(n-1)
#   }
#   return(result)
# }
# 
# 
# ExpAllDimClassify <- function(AListTrain1, labelTrain1, AListTrain2, labelTrain2,
#                               AListTest, labelTest, dVec, q, ncores=1, isSVD=1) {
#   source("getElbows.R")
#   source("USVT.R")
#   
#   nD <- length(dVec)
#   dMax <- max(dVec)
#   
#   # MLE
#   m1 <- length(AListTrain1)
#   m2 <- length(AListTrain2)
#   ABar1 <- add(AListTrain1)/m1
#   ABar2 <- add(AListTrain2)/m2
#   
#   errorMLEVec <- rep(0, length(AListTest))
#   for (iTest in 1:length(AListTest)) {
#     if (substr(labelTest[iTest], 1, 1) == substr(labelTrain1, 1, 1)) {
#       errorMLEVec[iTest] <- (norm(AListTest[[iTest]] - ABar1, "F") > norm(AListTest[[iTest]] - ABar2, "F"))
#     } else if (substr(labelTest[iTest], 1, 1) == substr(labelTrain2, 1, 1)) {
#       errorMLEVec[iTest] <- (norm(AListTest[[iTest]] - ABar1, "F") < norm(AListTest[[iTest]] - ABar2, "F"))
#     } else {
#       errorMLEVec[iTest] <- errorMLEVec[iTest] + "error"
#     }
#   }
#   
#   # MLqE
#   n <- dim(ABar1)[[1]]
#   AListTmp <- AListTrain1
#   for (i in 1:length(AListTmp)) {
#     AListTmp[[i]][lower.tri(AListTmp[[i]], T)] <- 0
#   }
#   ATensor <- array(unlist(AListTmp), dim = c(n, n, length(AListTmp)))
#   # AMLqE1 <- apply(ATensor, c(1, 2), mlqe_exp_solver, q)
#   
#   out <- mclapply(1:(n*(n-1)/2), function(k) {
#     i <- n - 2 - floor(sqrt(-8*(k - 1) + 4*n*(n-1)-7)/2.0 - 0.5);
#     j <- k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
#     i <- i + 1;
#     mlqe_exp_solver(ATensor[i, j, ], q)}, mc.cores = nCores)
#   
#   AMLqE1 <- matrix(0, n, n)
#   AMLqE1[lower.tri(AMLqE1, diag=FALSE)] <- unlist(out)
#   AMLqE1 <- AMLqE1 + t(AMLqE1)
#   
#   AListTmp <- AListTrain2
#   for (i in 1:length(AListTmp)) {
#     AListTmp[[i]][lower.tri(AListTmp[[i]], T)] <- 0
#   }
#   ATensor <- array(unlist(AListTmp), dim = c(n, n, length(AListTmp)))
#   # AMLqE2 <- apply(ATensor, c(1, 2), mlqe_exp_solver, q)
#   
#   out <- mclapply(1:(n*(n-1)/2), function(k) {
#     i <- n - 2 - floor(sqrt(-8*(k - 1) + 4*n*(n-1)-7)/2.0 - 0.5);
#     j <- k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
#     i <- i + 1;
#     mlqe_exp_solver(ATensor[i, j, ], q)}, mc.cores = nCores)
#   
#   AMLqE2 <- matrix(0, n, n)
#   AMLqE2[lower.tri(AMLqE2, diag=FALSE)] <- unlist(out)
#   AMLqE2 <- AMLqE2 + t(AMLqE2)
#   
#   errorMLqEVec <- rep(0, length(AListTest))
#   for (iTest in 1:length(AListTest)) {
#     if (substr(labelTest[iTest], 1, 1) == substr(labelTrain1, 1, 1)) {
#       errorMLqEVec[iTest] <- (norm(AListTest[[iTest]] - AMLqE1, "F") > norm(AListTest[[iTest]] - AMLqE2, "F"))
#     } else if (substr(labelTest[iTest], 1, 1) == substr(labelTrain2, 1, 1)) {
#       errorMLqEVec[iTest] <- (norm(AListTest[[iTest]] - AMLqE1, "F") < norm(AListTest[[iTest]] - AMLqE2, "F"))
#     } else {
#       errorMLqEVec[iTest] <- errorMLqEVec[iTest] + "error"
#     }
#   }
#   
#   # MLE_ASE
#   errorMLEASEVec <- matrix(0, nD, length(AListTest))
#   
#   ABarDiagAug1 <- diag_aug(ABar1)
#   ABarDiagAug2 <- diag_aug(ABar2)
#   
#   AASE1 = ase(ABarDiagAug1, dMax, isSVD)
#   AASE2 = ase(ABarDiagAug2, dMax, isSVD)
#   
#   # AListTestZG <- AListTest
#   nElbow <- 3
#   AASEList <- AListTest
#   for (iTest in 1:length(AListTest)) {
#     AASEList[[iTest]] = ase(diag_aug(AListTest[[iTest]]), dMax, isSVD)
#     # evalVec <- ase(diag_aug(AListTest[[iTest]]), ceiling(n*3/5), isSVD)[[1]]
#     # dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#     # AHat <- AASEList[[iTest]][[3]][, 1:dZG]%*%diag(AASEList[[iTest]][[1]][1:dZG])%*%t(AASEList[[iTest]][[2]][ ,1:dZG])
#     # AListTestZG[[iTest]] <- regularize(AHat)
#   }
#   
#   for (iD in 1:nD) {
#     d <- dVec[iD]
#     if (d == 1) {
#       AHat1 <- AASE1[[1]]*AASE1[[3]]%*%t(AASE1[[2]])
#       AHat2 <- AASE2[[1]]*AASE2[[3]]%*%t(AASE2[[2]])
#     } else {
#       AHat1 <- AASE1[[3]][, 1:d]%*%diag(AASE1[[1]][1:d])%*%t(AASE1[[2]][ ,1:d])
#       AHat2 <- AASE2[[3]][, 1:d]%*%diag(AASE2[[1]][1:d])%*%t(AASE2[[2]][ ,1:d])
#     }
#     PHat1 <- regularize(AHat1)
#     PHat2 <- regularize(AHat2)
#     for (iTest in 1:length(AListTest)) {
#       if (d == 1) {
#         AHat <- AASEList[[iTest]][[1]]*AASEList[[iTest]][[3]]%*%t(AASEList[[iTest]][[2]])
#       } else {
#         AHat <- AASEList[[iTest]][[3]][, 1:d]%*%diag(AASEList[[iTest]][[1]][1:d])%*%t(AASEList[[iTest]][[2]][ ,1:d])
#       }
#       PHat <- regularize(AHat)
#       if (substr(labelTest[iTest], 1, 1) == substr(labelTrain1, 1, 1)) {
#         errorMLEASEVec[iD, iTest] <- (norm(PHat - PHat1, "F") > norm(PHat - PHat2, "F"))
#       } else if (substr(labelTest[iTest], 1, 1) == substr(labelTrain2, 1, 1)) {
#         errorMLEASEVec[iD, iTest] <- (norm(PHat - PHat1, "F") < norm(PHat - PHat2, "F"))
#       } else {
#         errorMLEASEVec[iD, iTest] <- errorMLEASEVec[iD, iTest] + "error"
#       }
#     }
#   }
#   
#   # ZG
#   nElbow <- 3
#   evalVec <- ase(ABarDiagAug1, ceiling(n*3/5), isSVD)[[1]]
#   dZG1 <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#   evalVec <- ase(ABarDiagAug2, ceiling(n*3/5), isSVD)[[1]]
#   dZG2 <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#   d <- max(dZG1, dZG2)
#   AHat1 <- AASE1[[3]][, 1:d]%*%diag(AASE1[[1]][1:d])%*%t(AASE1[[2]][ ,1:d])
#   AHat2 <- AASE2[[3]][, 1:d]%*%diag(AASE2[[1]][1:d])%*%t(AASE2[[2]][ ,1:d])
#   PHat1 <- regularize(AHat1)
#   PHat2 <- regularize(AHat2)
#   errorMLEASE_ZGVec <- rep(0, length(AListTest))
#   for (iTest in 1:length(AListTest)) {
#     if (d == 1) {
#       AHat <- AASEList[[iTest]][[1]]*AASEList[[iTest]][[3]]%*%t(AASEList[[iTest]][[2]])
#     } else {
#       AHat <- AASEList[[iTest]][[3]][, 1:d]%*%diag(AASEList[[iTest]][[1]][1:d])%*%t(AASEList[[iTest]][[2]][ ,1:d])
#     }
#     PHat <- regularize(AHat)
#     if (substr(labelTest[iTest], 1, 1) == substr(labelTrain1, 1, 1)) {
#       errorMLEASE_ZGVec[iTest] <- (norm(PHat - PHat1, "F") > norm(PHat - PHat2, "F"))
#     } else if (substr(labelTest[iTest], 1, 1) == substr(labelTrain2, 1, 1)) {
#       errorMLEASE_ZGVec[iTest] <- (norm(PHat - PHat1, "F") < norm(PHat - PHat2, "F"))
#     } else {
#       errorMLEASE_ZGVec[iTest] <- errorMLEASE_ZGVec[iTest] + "error"
#     }
#   }
#   
#   # MLqE_ASE
#   errorMLqEASEVec <- matrix(0, nD, length(AListTest))
#   
#   AMLqEDiagAug1 <- diag_aug(AMLqE1)
#   AMLqEDiagAug2 <- diag_aug(AMLqE2)
#   
#   AASE1 <- ase(AMLqEDiagAug1, dMax, isSVD)
#   AASE2 <- ase(AMLqEDiagAug2, dMax, isSVD)
#   
#   for (iD in 1:nD) {
#     d <- dVec[iD]
#     if (d == 1) {
#       AHat1 <- AASE1[[1]]*AASE1[[3]]%*%t(AASE1[[2]])
#       AHat2 <- AASE2[[1]]*AASE2[[3]]%*%t(AASE2[[2]])
#     } else {
#       AHat1 <- AASE1[[3]][ ,1:d]%*%diag(AASE1[[1]][1:d])%*%t(AASE1[[2]][ ,1:d])
#       AHat2 <- AASE2[[3]][ ,1:d]%*%diag(AASE2[[1]][1:d])%*%t(AASE2[[2]][ ,1:d])
#     }
#     PHatASE1 <- regularize(AHat1)
#     PHatASE2 <- regularize(AHat2)
#     for (iTest in 1:length(AListTest)) {
#       if (d == 1) {
#         AHat <- AASEList[[iTest]][[1]]*AASEList[[iTest]][[3]]%*%t(AASEList[[iTest]][[2]])
#       } else {
#         AHat <- AASEList[[iTest]][[3]][, 1:d]%*%diag(AASEList[[iTest]][[1]][1:d])%*%t(AASEList[[iTest]][[2]][ ,1:d])
#       }
#       PHat <- regularize(AHat)
#       if (substr(labelTest[iTest], 1, 1) == substr(labelTrain1, 1, 1)) {
#         errorMLqEASEVec[iD, iTest] <- (norm(PHat - PHatASE1, "F") > norm(PHat - PHatASE2, "F"))
#       } else if (substr(labelTest[iTest], 1, 1) == substr(labelTrain2, 1, 1)) {
#         errorMLqEASEVec[iD, iTest] <- (norm(PHat - PHatASE1, "F") < norm(PHat - PHatASE2, "F"))
#       } else {
#         errorMLqEASEVec[iD, iTest] <- errorMLqEASEVec[iD, iTest] + "error"
#       }
#     }
#   }
#   
#   # ZG
#   nElbow <- 3
#   evalVec <- ase(AMLqEDiagAug1, ceiling(n*3/5), isSVD)[[1]]
#   dZG1 <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#   evalVec <- ase(AMLqEDiagAug2, ceiling(n*3/5), isSVD)[[1]]
#   dZG2 <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#   d <- max(dZG1, dZG2)
#   AHat1 <- AASE1[[3]][, 1:dZG1]%*%diag(AASE1[[1]][1:dZG1])%*%t(AASE1[[2]][ ,1:dZG1])
#   AHat2 <- AASE2[[3]][, 1:dZG2]%*%diag(AASE2[[1]][1:dZG2])%*%t(AASE2[[2]][ ,1:dZG2])
#   PHatASE1 <- regularize(AHat1)
#   PHatASE2 <- regularize(AHat2)
#   errorMLqEASE_ZGVec <- rep(0, length(AListTest))
#   for (iTest in 1:length(AListTest)) {
#     if (d == 1) {
#       AHat <- AASEList[[iTest]][[1]]*AASEList[[iTest]][[3]]%*%t(AASEList[[iTest]][[2]])
#     } else {
#       AHat <- AASEList[[iTest]][[3]][, 1:d]%*%diag(AASEList[[iTest]][[1]][1:d])%*%t(AASEList[[iTest]][[2]][ ,1:d])
#     }
#     PHat <- regularize(AHat)
#     if (substr(labelTest[iTest], 1, 1) == substr(labelTrain1, 1, 1)) {
#       errorMLqEASE_ZGVec[iTest] <- (norm(PHat - PHatASE1, "F") > norm(PHat - PHatASE2, "F"))
#     } else if (substr(labelTest[iTest], 1, 1) == substr(labelTrain2, 1, 1)) {
#       errorMLqEASE_ZGVec[iTest] <- (norm(PHat - PHatASE1, "F") < norm(PHat - PHatASE2, "F"))
#     } else {
#       errorMLqEASE_ZGVec[iTest] <- errorMLqEASE_ZGVec[iTest] + "error"
#     }
#   }
#   
#   result <- list(errorMLEVec, errorMLqEVec, errorMLEASEVec, errorMLEASE_ZGVec,
#                  errorMLqEASEVec, errorMLqEASE_ZGVec)
# }
# 
# 
# 
# 
# BuildMcNemarMatrix <- function(xVec, yVec) {
#   testMatrix <- matrix(0, 2, 2)
#   testMatrix[1, 1] <- sum(!(xVec | yVec))
#   testMatrix[1, 2] <- sum(!(xVec | !yVec))
#   testMatrix[2, 1] <- sum(!(!xVec | yVec))
#   testMatrix[2, 2] <- sum(!(!xVec | !yVec))
#   return(testMatrix)
# }
# 
# ExpKMeans <- function(AList, K, q, labelVec, nCores=1, isSVD=1) {
#   # require(fossil)
#   tauStar <- sapply(1:length(labelVec), function(iter) {
#     if (substr(labelVec[iter], 1, 1) == "G")
#       return(1)
#     else if (substr(labelVec[iter], 1, 1) == "N")
#       return(2)
#     else return(3)})
#   
#   
#   n <- dim(AList[[1]])[1]
#   source("getElbows.R")
#   nElbow <- 2
#   for (i in 1:length(AList)) {
#     diag(AList[[i]]) <- 0
#   }
#   muList <- AList[sample(1:length(AList), K)]
#   
#   dist <- matrix(0, K, length(AList))
#   # Assignment
#   for (k in 1:K) {
#     dist[k, ] <- sapply(1:(length(AList)), function(iter) {norm(AList[[iter]] - muList[[k]], "F")^2})
#   }
#   tau0 <- sapply(1:length(AList), function(iter) {which.min(dist[, iter])})
#   err0 <- sum(apply(dist, 2, min))
#   maxTol <- (1e-6)*err0
#   
#   maxIter <- 100
#   
#   ###### MLE ######
#   tau <- tau0
#   errOld <- err0
#   errMin <- errOld
#   tauMLE <- tau
#   errDiff <- maxTol + 1
#   iIter <- 0
#   while ((abs(errDiff) > maxTol) && (iIter < maxIter)) {
#     iIter <- iIter + 1
#     print(iIter)
#     
#     # Average
#     for (k in 1:K) {
#       nv <- (tau == k)
#       AListGroup <- AList[nv]
#       muList[[k]] <- add(AListGroup)/length(AListGroup)
#     }
#     rm(AListGroup)
#     
#     # Assignment
#     for (k in 1:K) {
#       dist[k, ] <- sapply(1:(length(AList)), function(iter) {norm(AList[[iter]] - muList[[k]], "F")^2})
#     }
#     tau <- sapply(1:length(AList), function(iter) {which.min(dist[, iter])})
#     errNew <- sum(apply(dist, 2, min))
#     if (errNew < errMin) {
#       print("Update")
#       errMin <- errNew
#       tauMLE <- tau
#       muListMLE <- muList
#     }
#     errDiff <- errNew - errOld
#     # print(c(errNew, errDiff, maxTol, adj.rand.index(tau, tauStar)))
#     print(c(errNew, errDiff, maxTol))
#     errOld <- errNew
#   }
#   
#   
#   
#   
#   ###### MLqE ######
#   tau <- tau0
#   errOld <- err0
#   errMin <- errOld
#   tauMLqE <- tau
#   errDiff <- maxTol + 1
#   iIter <- 0
#   while ((abs(errDiff) > maxTol) && (iIter < maxIter)) {
#     iIter <- iIter + 1
#     print(iIter)
#     
#     # Average
#     for (k in 1:K) {
#       nv <- (tau == k)
#       AListTmp <- AList[nv]
#       for (i in 1:length(AListTmp)) {
#         AListTmp[[i]][lower.tri(AListTmp[[i]], T)] <- 0
#       }
#       ATensor <- array(unlist(AListTmp), dim = c(n, n, length(AListTmp)))
#       
#       out <- mclapply(1:(n*(n-1)/2), function(k) {
#         i <- n - 2 - floor(sqrt(-8*(k - 1) + 4*n*(n-1)-7)/2.0 - 0.5);
#         j <- k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
#         i <- i + 1;
#         mlqe_exp_solver(ATensor[i, j, ], q, tol = 1e-3)}, mc.cores = nCores)
#       
#       AMLqE <- matrix(0, n, n)
#       AMLqE[lower.tri(AMLqE, diag=FALSE)] <- unlist(out)
#       AMLqE <- AMLqE + t(AMLqE)
#       muList[[k]] <- AMLqE
#     }
#     
#     
#     # Assignment
#     for (k in 1:K) {
#       dist[k, ] <- sapply(1:(length(AList)), function(iter) {norm(AList[[iter]] - muList[[k]], "F")^2})
#     }
#     tau <- sapply(1:length(AList), function(iter) {which.min(dist[, iter])})
#     errNew <- sum(apply(dist, 2, min))
#     if (errNew < errMin) {
#       print("Update")
#       errMin <- errNew
#       tauMLqE <- tau
#       muListMLqE <- muList
#     }
#     errDiff <- errNew - errOld
#     # print(c(errNew, errDiff, maxTol, adj.rand.index(tau, tauStar)))
#     print(c(errNew, errDiff, maxTol))
#     errOld <- errNew
#   }
#   rm(AListTmp)
#   
#   ###### MLE_ASE ######
#   tau <- tau0
#   errOld <- err0
#   errMin <- errOld
#   tauMLEASE <- tau
#   errDiff <- maxTol + 1
#   iIter <- 0
#   while ((abs(errDiff) > maxTol) && (iIter < maxIter)) {
#     iIter <- iIter + 1
#     print(iIter)
#     
#     # Average
#     for (k in 1:K) {
#       nv <- (tau == k)
#       AListGroup <- AList[nv]
#       ABar <- add(AListGroup)/length(AListGroup)
#       ABarDiagAug <- diag_aug(ABar)
#       evalVec <- ase(ABarDiagAug, ceiling(n*3/5), isSVD)[[1]]
#       dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#       d <- dZG
#       AASE <- ase(ABarDiagAug, dZG, isSVD)
#       AHat <- AASE[[3]]%*%diag(AASE[[1]])%*%t(AASE[[2]])
#       muList[[k]] <- regularize(AHat)
#     }
#     
#     # Assignment
#     for (k in 1:K) {
#       dist[k, ] <- sapply(1:(length(AList)), function(iter) {norm(AList[[iter]] - muList[[k]], "F")^2})
#     }
#     tau <- sapply(1:length(AList), function(iter) {which.min(dist[, iter])})
#     errNew <- sum(apply(dist, 2, min))
#     if (errNew < errMin) {
#       print("Update")
#       errMin <- errNew
#       tauMLEASE <- tau
#       muListMLEASE <- muList
#     }
#     errDiff <- errNew - errOld
#     # print(c(errNew, errDiff, maxTol, adj.rand.index(tau, tauStar)))
#     print(c(errNew, errDiff, maxTol))
#     errOld <- errNew
#   }
#   
#   
#   
#   
#   ###### MLqE_ASE ######
#   tau <- tau0
#   errOld <- err0
#   errMin <- errOld
#   tauMLqEASE <- tau
#   errDiff <- maxTol + 1
#   iIter <- 0
#   while ((abs(errDiff) > maxTol) && (iIter < maxIter)) {
#     iIter <- iIter + 1
#     print(iIter)
#     
#     # Average
#     for (k in 1:K) {
#       nv <- (tau == k)
#       AListTmp <- AList[nv]
#       for (i in 1:length(AListTmp)) {
#         AListTmp[[i]][lower.tri(AListTmp[[i]], T)] <- 0
#       }
#       ATensor <- array(unlist(AListTmp), dim = c(n, n, length(AListTmp)))
#       
#       out <- mclapply(1:(n*(n-1)/2), function(k) {
#         i <- n - 2 - floor(sqrt(-8*(k - 1) + 4*n*(n-1)-7)/2.0 - 0.5);
#         j <- k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
#         i <- i + 1;
#         mlqe_exp_solver(ATensor[i, j, ], q, tol = 1e-3)}, mc.cores = nCores)
#       
#       AMLqE <- matrix(0, n, n)
#       AMLqE[lower.tri(AMLqE, diag=FALSE)] <- unlist(out)
#       AMLqE <- AMLqE + t(AMLqE)
#       
#       AMLqEDiagAug <- diag_aug(AMLqE)
#       evalVec <- ase(AMLqEDiagAug, ceiling(n*3/5), isSVD)[[1]]
#       dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#       d <- dZG
#       AASE <- ase(AMLqEDiagAug, dZG, isSVD)
#       AHat <- AASE[[3]]%*%diag(AASE[[1]])%*%t(AASE[[2]])
#       muList[[k]] <- regularize(AHat)
#     }
#     
#     # Assignment
#     for (k in 1:K) {
#       dist[k, ] <- sapply(1:(length(AList)), function(iter) {norm(AList[[iter]] - muList[[k]], "F")^2})
#     }
#     tau <- sapply(1:length(AList), function(iter) {which.min(dist[, iter])})
#     errNew <- sum(apply(dist, 2, min))
#     if (errNew < errMin) {
#       print("Update")
#       errMin <- errNew
#       tauMLqEASE <- tau
#       muListMLqEASE <- muList
#     }
#     errDiff <- errNew - errOld
#     # print(c(errNew, errDiff, maxTol, adj.rand.index(tau, tauStar)))
#     print(c(errNew, errDiff, maxTol))
#     errOld <- errNew
#   }
#   
#   
#   result <- list(tauMLE, tauMLqE, tauMLEASE, tauMLqEASE, tauStar, muListMLE, muListMLqE, muListMLEASE, muListMLqEASE)
#   # result <- list(tauMLE, tauMLEASE, tauStar, muListMLE, muListMLEASE)
# }
# 
# 
# ExpSimKMeans <- function(AList, K, q, nCores=1, isSVD=1) {
#   n <- dim(AList[[1]])[1]
#   source("getElbows.R")
#   nElbow <- 2
#   for (i in 1:length(AList)) {
#     diag(AList[[i]]) <- 0
#   }
#   muList <- AList[sample(1:length(AList), K)]
#   
#   dist <- matrix(0, K, length(AList))
#   # Assignment
#   for (k in 1:K) {
#     dist[k, ] <- sapply(1:(length(AList)), function(iter) {norm(AList[[iter]] - muList[[k]], "F")^2})
#   }
#   tau0 <- sapply(1:length(AList), function(iter) {which.min(dist[, iter])})
#   err0 <- sum(apply(dist, 2, min))
#   maxTol <- (1e-6)*err0
#   
#   maxIter <- 100
#   
#   ###### MLE ######
#   tau <- tau0
#   errOld <- err0
#   errMin <- errOld
#   tauMLE <- tau
#   errDiff <- maxTol + 1
#   iIter <- 0
#   while ((abs(errDiff) > maxTol) && (iIter < maxIter)) {
#     iIter <- iIter + 1
#     print(iIter)
#     
#     # Average
#     for (k in 1:K) {
#       nv <- (tau == k)
#       AListGroup <- AList[nv]
#       muList[[k]] <- add(AListGroup)/length(AListGroup)
#     }
#     rm(AListGroup)
#     
#     # Assignment
#     for (k in 1:K) {
#       dist[k, ] <- sapply(1:(length(AList)), function(iter) {norm(AList[[iter]] - muList[[k]], "F")^2})
#     }
#     tau <- sapply(1:length(AList), function(iter) {which.min(dist[, iter])})
#     errNew <- sum(apply(dist, 2, min))
#     if (errNew < errMin) {
#       print("Update")
#       errMin <- errNew
#       tauMLE <- tau
#       muListMLE <- muList
#     }
#     errDiff <- errNew - errOld
#     # print(c(errNew, errDiff, maxTol, adj.rand.index(tau, tauStar)))
#     print(c(errNew, errDiff, maxTol))
#     errOld <- errNew
#   }
#   
#   
#   
#   
#   ###### MLqE ######
#   tau <- tau0
#   errOld <- err0
#   errMin <- errOld
#   tauMLqE <- tau
#   errDiff <- maxTol + 1
#   iIter <- 0
#   while ((abs(errDiff) > maxTol) && (iIter < maxIter)) {
#     iIter <- iIter + 1
#     print(iIter)
#     
#     # Average
#     for (k in 1:K) {
#       nv <- (tau == k)
#       AListTmp <- AList[nv]
#       for (i in 1:length(AListTmp)) {
#         AListTmp[[i]][lower.tri(AListTmp[[i]], T)] <- 0
#       }
#       ATensor <- array(unlist(AListTmp), dim = c(n, n, length(AListTmp)))
#       
#       out <- mclapply(1:(n*(n-1)/2), function(k) {
#         i <- n - 2 - floor(sqrt(-8*(k - 1) + 4*n*(n-1)-7)/2.0 - 0.5);
#         j <- k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
#         i <- i + 1;
#         mlqe_exp_solver(ATensor[i, j, ], q, tol = 1e-3)}, mc.cores = nCores)
#       
#       AMLqE <- matrix(0, n, n)
#       AMLqE[lower.tri(AMLqE, diag=FALSE)] <- unlist(out)
#       AMLqE <- AMLqE + t(AMLqE)
#       muList[[k]] <- AMLqE
#     }
#     
#     
#     # Assignment
#     for (k in 1:K) {
#       dist[k, ] <- sapply(1:(length(AList)), function(iter) {norm(AList[[iter]] - muList[[k]], "F")^2})
#     }
#     tau <- sapply(1:length(AList), function(iter) {which.min(dist[, iter])})
#     errNew <- sum(apply(dist, 2, min))
#     if (errNew < errMin) {
#       print("Update")
#       errMin <- errNew
#       tauMLqE <- tau
#       muListMLqE <- muList
#     }
#     errDiff <- errNew - errOld
#     # print(c(errNew, errDiff, maxTol, adj.rand.index(tau, tauStar)))
#     print(c(errNew, errDiff, maxTol))
#     errOld <- errNew
#   }
#   rm(AListTmp)
#   
#   ###### MLE_ASE ######
#   AList0 <- AList
#   for (i in 1:length(AList)) {
#     # d <- 55
#     d <- 10
#     AASE <- ase(AList[[i]], d, isSVD)
#     if (d == 1) {
#       AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
#     } else {
#       AHat <- AASE[[3]][ ,1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][ ,1:d])
#     }
#     AList[[i]] <- regularize(AHat)
#   }
#   
#   
#   
#   tau <- tau0
#   errOld <- err0
#   errMin <- errOld
#   tauMLEASE <- tau
#   errDiff <- maxTol + 1
#   iIter <- 0
#   while ((abs(errDiff) > maxTol) && (iIter < maxIter)) {
#     iIter <- iIter + 1
#     print(iIter)
#     
#     # Average
#     for (k in 1:K) {
#       nv <- (tau == k)
#       AListGroup <- AList[nv]
#       ABar <- add(AListGroup)/length(AListGroup)
#       ABarDiagAug <- diag_aug(ABar)
#       evalVec <- ase(ABarDiagAug, ceiling(n*3/5), isSVD)[[1]]
#       dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#       d <- dZG
#       AASE <- ase(ABarDiagAug, dZG, isSVD)
#       AHat <- AASE[[3]]%*%diag(AASE[[1]])%*%t(AASE[[2]])
#       muList[[k]] <- regularize(AHat)
#     }
#     
#     # Assignment
#     for (k in 1:K) {
#       dist[k, ] <- sapply(1:(length(AList)), function(iter) {norm(AList[[iter]] - muList[[k]], "F")^2})
#     }
#     tau <- sapply(1:length(AList), function(iter) {which.min(dist[, iter])})
#     errNew <- sum(apply(dist, 2, min))
#     if (errNew < errMin) {
#       print("Update")
#       errMin <- errNew
#       tauMLEASE <- tau
#       muListMLEASE <- muList
#     }
#     errDiff <- errNew - errOld
#     # print(c(errNew, errDiff, maxTol, adj.rand.index(tau, tauStar)))
#     print(c(errNew, errDiff, maxTol))
#     errOld <- errNew
#   }
#   
#   
#   ###### MLqE_ASE ######
#   tau <- tau0
#   errOld <- err0
#   errMin <- errOld
#   tauMLqEASE <- tau
#   errDiff <- maxTol + 1
#   iIter <- 0
#   while ((abs(errDiff) > maxTol) && (iIter < maxIter)) {
#     iIter <- iIter + 1
#     print(iIter)
#     
#     # Average
#     for (k in 1:K) {
#       nv <- (tau == k)
#       AListTmp <- AList[nv]
#       for (i in 1:length(AListTmp)) {
#         AListTmp[[i]][lower.tri(AListTmp[[i]], T)] <- 0
#       }
#       ATensor <- array(unlist(AListTmp), dim = c(n, n, length(AListTmp)))
#       
#       out <- mclapply(1:(n*(n-1)/2), function(k) {
#         i <- n - 2 - floor(sqrt(-8*(k - 1) + 4*n*(n-1)-7)/2.0 - 0.5);
#         j <- k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
#         i <- i + 1;
#         mlqe_exp_solver(ATensor[i, j, ], q, tol = 1e-3)}, mc.cores = nCores)
#       
#       AMLqE <- matrix(0, n, n)
#       AMLqE[lower.tri(AMLqE, diag=FALSE)] <- unlist(out)
#       AMLqE <- AMLqE + t(AMLqE)
#       
#       AMLqEDiagAug <- diag_aug(AMLqE)
#       evalVec <- ase(AMLqEDiagAug, ceiling(n*3/5), isSVD)[[1]]
#       dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
#       d <- dZG
#       AASE <- ase(AMLqEDiagAug, dZG, isSVD)
#       AHat <- AASE[[3]]%*%diag(AASE[[1]])%*%t(AASE[[2]])
#       muList[[k]] <- regularize(AHat)
#     }
#     
#     # Assignment
#     for (k in 1:K) {
#       dist[k, ] <- sapply(1:(length(AList)), function(iter) {norm(AList[[iter]] - muList[[k]], "F")^2})
#     }
#     tau <- sapply(1:length(AList), function(iter) {which.min(dist[, iter])})
#     errNew <- sum(apply(dist, 2, min))
#     if (errNew < errMin) {
#       print("Update")
#       errMin <- errNew
#       tauMLqEASE <- tau
#       muListMLqEASE <- muList
#     }
#     errDiff <- errNew - errOld
#     # print(c(errNew, errDiff, maxTol, adj.rand.index(tau, tauStar)))
#     print(c(errNew, errDiff, maxTol))
#     errOld <- errNew
#   }
#   
#   
#   result <- list(tauMLE, tauMLqE, tauMLEASE, tauMLqEASE, muListMLE, muListMLqE, muListMLEASE, muListMLqEASE)
#   # result <- list(tauMLE, tauMLEASE, tauStar, muListMLE, muListMLEASE)
# }