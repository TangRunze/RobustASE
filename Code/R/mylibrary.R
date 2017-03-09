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






SimAll <- function(m, n, tau, B, CB, eps, q, d, isSVD=1) {
  result = rep(NaN, 4)
  
  P <- B[tau, tau]
  C <- CB[tau, tau]
  diag(P) <- 0
  diag(C) <- 0
  
  AList <- list()
  array(0, dim = c(m, n, n))
  ind <- lower.tri(array(0, dim=c(n, n)), 1)
  for (i in 1:m) {
    contamVec <- (runif(n^2) > eps)
    A <- contamVec*matrix(rexp(n^2, 1/P), ncol = n) +
      (1 - contamVec)*matrix(rexp(n^2, 1/C), ncol = n)
    A[ind] <- 0
    A <- A + t(A)
    AList[[i]] <- A
  }
  
  # MLE
  A_MLE <- add(AList)/m
  result[1] <- (norm(P - A_MLE, "F"))^2/(n*(n-1))
  
  # MLE_ASE
  require(Matrix)
  D0 <- Diagonal(n, x = rowSums(A_MLE)/(n-1))
  # dZG <- DimSelect(A_MLE + D0, isSVD, method = 1)
  P0 <- LR(A_MLE + D0, d, isSVD)
  D1 <- Diagonal(n, x = diag(P0))
  P1 <- LR(A_MLE + D1, d)
  A_MLE_ASE <- Regularize(P1, T)
  result[2] <- (norm(P - A_MLE_ASE, "F"))^2/(n*(n-1))
  
  # MLqE  
  AListTmp <- AList
  nv <- lower.tri(AListTmp[[1]], T)
  for (i in 1:length(AListTmp)) {
    AListTmp[[i]][nv] <- 0
  }
  ATensor <- array(unlist(AListTmp), dim = c(n, n, m))
  A_MLqE <- apply(ATensor, c(1, 2), MLqE_Exp_Solver, q)
  A_MLqE <- A_MLqE + t(A_MLqE)
  result[3] = (norm(P - A_MLqE, "F"))^2/(n*(n-1))
  
  # MLqE_ASE
  D0 <- Diagonal(n, x = rowSums(A_MLqE)/(n-1))
  # dZG <- DimSelect(A_MLqE + D0, isSVD, method = 1)
  P0 <- LR(A_MLqE + D0, d, isSVD)
  D1 <- Diagonal(n, x = diag(P0))
  P1 <- LR(A_MLqE + D1, d)
  A_MLqE_ASE <- Regularize(P1, T)
  result[4] <- (norm(P - A_MLqE_ASE, "F"))^2/(n*(n-1))
  
  return(result)
}







