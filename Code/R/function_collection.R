grid_arrange_shared_legend2 <- function(plots, nrows = 1, ncols = 2) {
  library(gridExtra)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  pl  <- lapply(plots, function(x) x + theme(legend.position="none"))
  tmp <- do.call(arrangeGrob, c(pl, list(ncol=ncols, nrow=nrows)))
  grid.arrange(tmp, legend, ncol=1, heights = unit.c(unit(1, "npc") - lheight, lheight))
}

# Read M graphs
read_data <- function(dataName, DA=T) {
  if (DA) {
    fileName = paste("../../Data/data_", dataName, "_DA.RData", sep="")
  } else {
    fileName = paste("../../Data/data_", dataName, ".RData", sep="")
  }
  if (file.exists(fileName)) {
    load(fileName)
    return(list(A_all, n, M))
  } else {
    require(igraph)
    subjectsID = readLines("../../Data/subnames.txt")
    if (dataName == "DS01876") {
      g = read_graph(paste("../../Data/", dataName, "/SWU4_", subjectsID[1], 
                           "_1_DTI_", dataName, ".graphml", sep =""), format="graphml")      
    } else {
      g = read_graph(paste("../../Data/", dataName, "/SWU4_", subjectsID[1], 
                           "_1_", dataName, "_sg.graphml", sep =""), format="graphml")
    }
    n = vcount(g)
    
    M = 227*2;
    A_all = list()
    for (sub in 1:227) {
      for (session in 1:2) {
        if (dataName == "DS01876") {
          g = read_graph(paste("../../Data/", dataName, "/SWU4_", subjectsID[sub], 
                               "_", session, "_DTI_", dataName, ".graphml",sep=""),
                         format = "graphml")          
        } else {
          g = read_graph(paste("../../Data/", dataName, "/SWU4_", subjectsID[sub], 
                               "_", session, "_", dataName, "_sg.graphml",sep=""),
                         format = "graphml")
        }
        A = as_adj(g, type="both", sparse=FALSE)
        if (DA) {
          A = diag_aug(A)
        }
        A_all[[(sub-1)*2 + session]] = A;
      }
    }
    
    save(A_all, n, M, file=fileName)
    return(list(A_all, n, M))
  }  
}


ReadDataWeighted <- function(dataName, DA=T, newGraph=T) {
  if (DA) {
    if (newGraph == F) {
      fileName = paste("../../Data/data_", dataName, "_DA.RData", sep="")
    } else {
      fileName = paste("../../Data/data_", dataName, "_new_DA.RData", sep="")
    }
  } else {
    if (newGraph == F) {
      fileName = paste("../../Data/data_", dataName, ".RData", sep="")
    } else {
      fileName = paste("../../Data/data_", dataName, "_new.RData", sep="")
    }
  }
  if (file.exists(fileName)) {
    load(fileName)
    return(list(A_all, n, M))
  } else {
    require(igraph)
    subjectsID = readLines("../../Data/subnames.txt")
    if (dataName == "DS01876") {
      g = read_graph(paste("../../Data/", dataName, "/SWU4_", subjectsID[1], 
                           "_1_DTI_", dataName, ".graphml", sep =""), format="graphml")
    } else {
      if (newGraph == F) {
        g = read_graph(paste("../../Data/", dataName, "/SWU4_", subjectsID[1], 
                             "_1_", dataName, "_sg.graphml", sep =""), format="graphml")
      } else {
        g = read_graph(paste("../../Data/", dataName, "_new/SWU4_", subjectsID[1], 
                             "_1_DTI_", dataName, ".graphml", sep =""), format="graphml")      
      }
    }
    n = vcount(g)
    
    M = 227*2;
    A_all = list()
    for (sub in 1:227) {
      for (session in 1:2) {
        if (dataName == "DS01876") {
          if (newGraph == F) {
            g = read_graph(paste("../../Data/", dataName, "/SWU4_", subjectsID[sub], 
                                 "_", session, "_DTI_", dataName, ".graphml",sep=""),
                           format = "graphml")
          } else {
            g = read_graph(paste("../../Data/", dataName, "_new/SWU4_", subjectsID[sub], 
                                 "_", session, "_DTI_", dataName, ".graphml",sep=""),
                           format = "graphml")
          }
        } else {
          if (newGraph == F) {
            g = read_graph(paste("../../Data/", dataName, "/SWU4_", subjectsID[sub], 
                                 "_", session, "_", dataName, "_sg.graphml",sep=""),
                           format = "graphml")
          } else {
            g = read_graph(paste("../../Data/", dataName, "_new/SWU4_", subjectsID[sub], 
                                 "_", session, "_DTI_", dataName, ".graphml",sep=""),
                           format = "graphml")
          }
        }
        A = as_adj(g, attr="weight", type="both", sparse=FALSE)
        if (DA) {
          A = diag_aug(A)
        }
        A_all[[(sub-1)*2 + session]] = A;
      }
    }
    
    save(A_all, n, M, file=fileName)
    return(list(A_all, n, M))
  }  
}


# Diagonal Augmentation
diag_aug <- function(A, d=0) {
  if (d == 0) {
    require(Matrix)
    n = dim(A)[1]
    return(A + Diagonal(n, x=rowSums(A))/(n-1))
  } else {
    for (iIter in 1:1) {
      tmp = ase(A, d, 0)
      if (d == 1)
        diag(A) = diag(tmp[[1]] * tmp[[3]] %*% t(tmp[[2]]))
      else
        diag(A) = diag(tmp[[3]][,1:d] %*% diag(tmp[[1]][1:d]) %*% t(tmp[[2]][,1:d]))
    }
    return(A)
  }
}




# Regularize probability matrix
regularize <- function(A) {
  diag(A) = 0
  #   A[A > 1] = 1
  A[A < 0] = 0
  return(A)
}


add <- function(x) Reduce("+", x)

dim_brute <- function(m, n, rho, tau, B, dVec, isSVD=1) {
  result = rep(NaN, nD+1)
  
  require(igraph)
  A_all = list()
  for (i in 1:m) {
    g = sample_sbm(n, B, n*rho, directed=F, loops=F)
    A = as_adj(g, type="both", sparse=FALSE)
    A_all[[i]] = A
  }
  
  tau = rep(1:K,n*rho)
  P = B[tau,tau]
  diag(P) = 0
  
  A_bar = add(A_all)/m
  result[1] = norm(P - A_bar, "F")/n/(n-1)
  
  dMax = max(dVec)
  nD = length(dVec)
  
  A.ase = ase(diag_aug(A_bar), dMax, isSVD)
  for (iD in 1:nD) {
    d = dVec[iD]
    if (d == 1)
      Ahat = A.ase[[1]][1] * A.ase[[3]][,1:d] %*% t(A.ase[[2]][,1:d])
    else
      Ahat <- A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
    P_hat = regularize(Ahat)
    result[iD+1] = norm(P - P_hat, "F")/n/(n-1)
  }
  
  return(result)
}


dim_brute1 <- function(M, m, dVec, A_all, A_sum, isSVD=1) {
  result = rep(NaN, nD+1)
  
  sampleVec = sample.int(M, m)
  A_bar = add(A_all[sampleVec])
  P_bar = (A_sum - A_bar)/(M - m)
  #   P_bar = A_sum/M
  A_bar = A_bar/m
  result[1] = norm(P_bar - A_bar, "F")/n/(n-1)
  
  dMax = max(dVec)
  nD = length(dVec)
  
  A.ase = ase(diag_aug(A_bar), dMax, isSVD)
  for (iD in 1:nD) {
    d = dVec[iD]
    if (d == 1)
      Ahat = A.ase[[1]] * A.ase[[3]] %*% t(A.ase[[2]])
    else
      Ahat <- A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
    P_hat = regularize(Ahat)
    result[iD+1] = norm(P_bar - P_hat, "F")/n/(n-1)
  }
  
  return(result)
}


dim_brute2 <- function(M, m, dVec, A_all, A_sum, isSVD=1) {
  source("getElbows.R")
  source("USVT.R")
  
  result = rep(NaN, nD+1)
  
  sampleVec = sample.int(M, m)
  A_bar = add(A_all[sampleVec])
  P_bar = (A_sum - A_bar)/(M - m)
  #   P_bar = A_sum/M
  A_bar = A_bar/m
  result[1] = norm(P_bar - A_bar, "F")/n/(n-1)
  
  dMax = max(dVec)
  nD = length(dVec)
  
  A_bar_diag_aug = diag_aug(A_bar)
  
  # ZG
  nElbow = 3
  evalVec = ase(A_bar_diag_aug, ceiling(n*3/5), isSVD)[[1]]
  dZG = getElbows(evalVec, n=nElbow, plot=F)[[nElbow]]
  
  # USVT
  dUSVT = length(usvt(A_bar_diag_aug, 1, m)$d)
  
  A.ase = ase(diag_aug(A_bar), dMax, isSVD)
  for (iD in 1:nD) {
    d = dVec[iD]
    if (d == 1)
      Ahat = A.ase[[1]] * A.ase[[3]] %*% t(A.ase[[2]])
    else
      Ahat <- A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
    P_hat = regularize(Ahat)
    result[iD+1] = norm(P_bar - P_hat, "F")/n/(n-1)
  }
  
  result[nD+2] = dZG
  result[nD+3] = dUSVT
  
  return(result)
}


dim_brute2_all <- function(M, m, dVec, A_all, A_sum, q, isSVD=1) {
  source("getElbows.R")
  source("USVT.R")
  
  nD = length(dVec)
  dMax = max(dVec)
  result = rep(NaN, 2*nD+6)
  
  sampleVec = sample.int(M, m)
  A_bar = add(A_all[sampleVec])
  P_bar = (A_sum - A_bar)/(M - m)
  #   P_bar = A_sum/M
  A_bar = A_bar/m
  result[1] = norm(P_bar - A_bar, "F")/n/(n-1)
  
  n = dim(A_bar)[[1]]
  A_all_unlist = array(unlist(A_all[sampleVec]), dim=c(n, n, m))
  time0 = proc.time()
  A_mlqe = apply(A_all_unlist, c(1, 2), mlqe_exp_solver, q)
  time1 = proc.time() - time0
  result[nD+2] = norm(P_bar - A_mlqe, "F")/n/(n-1)
  
  A_bar_diag_aug = diag_aug(A_bar)
  # ZG
  nElbow = 3
  evalVec = ase(A_bar_diag_aug, ceiling(n*3/5), isSVD)[[1]]
  dZG = getElbows(evalVec, n=nElbow, plot=F)[[nElbow]]
  # USVT
  dUSVT = length(usvt(A_bar_diag_aug, 1, m)$d)
  result[nD*2+3] = dZG
  result[nD*2+4] = dUSVT
  
  A.ase = ase(A_bar_diag_aug, dMax, isSVD)
  for (iD in 1:nD) {
    d = dVec[iD]
    if (d == 1)
      Ahat = A.ase[[1]] * A.ase[[3]] %*% t(A.ase[[2]])
    else
      Ahat <- A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
    P_hat = regularize(Ahat)
    result[1+iD] = norm(P_bar - P_hat, "F")/n/(n-1)
  }
  #   plot(result[[2]])
  
  A_mlqe_diag_aug = diag_aug(A_mlqe)
  # ZG
  nElbow = 3
  evalVec = ase(A_mlqe_diag_aug, ceiling(n*3/5), isSVD)[[1]]
  dZG = getElbows(evalVec, n=nElbow, plot=F)[[nElbow]]
  # USVT
  dUSVT = length(usvt(A_mlqe_diag_aug, 1, m)$d)
  result[nD*2+5] = dZG
  result[nD*2+6] = dUSVT
  
  A.ase = ase(A_mlqe_diag_aug, dMax, isSVD)
  for (iD in 1:nD) {
    d = dVec[iD]
    if (d == 1)
      Ahat = A.ase[[1]] * A.ase[[3]] %*% t(A.ase[[2]])
    else
      Ahat <- A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
    P_hat_ase = regularize(Ahat)
    result[nD+2+iD] = norm(P_bar - P_hat_ase, "F")/n/(n-1)
  }
  #   plot(result[[4]])
  return(result)
}


dim_brute_robust <- function(M, m, dVec, A_all, A_all_unlist, isSVD=1) {
  result = rep(NaN, nD+1)
  
  sampleVec = sample.int(M, m)
  A_bar_robust = apply(A_all_unlist[,,sampleVec], 1:2, median)
  A_bar_robust = medianlist(A_all[sampleVec])
  P_bar = (A_sum - A_bar)/(M - m)
  #   P_bar = A_sum/M
  A_bar = A_bar/m
  result[1] = norm(P_bar - A_bar, "F")/n/(n-1)
  
  dMax = max(dVec)
  nD = length(dVec)
  
  A.ase = ase(diag_aug(A_bar), dMax, isSVD)
  for (iD in 1:nD) {
    d = dVec[iD]
    if (d == 1)
      Ahat = A.ase[[1]] * A.ase[[3]] %*% t(A.ase[[2]])
    else
      Ahat <- A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
    P_hat = regularize(Ahat)
    result[iD+1] = norm(P_bar - P_hat, "F")/n/(n-1)
  }
  
  return(result)
}



# Using ZG to choose dimension
llg_ZG <- function(M, m, A_all, A_sum, isSVD=1) {
  
  sampleVec = sample.int(M, m)
  A_bar = add(A_all[sampleVec])
  P_bar = (A_sum - A_bar)/(M - m);
  A_bar = A_bar/m;
  
  P_hat = regularize(ase.Ahat(diag_aug(A_bar), d, isSVD))
  
  return(c(norm(P_bar - A_bar, "F")/n/(n-1), norm(P_bar - P_hat, "F")/n/(n-1)), d)
}




llg_d <- function(M, m, A_all, A_sum, d, isSVD=1) {
  
  sampleVec = sample.int(M, m)
  A_bar = add(A_all[sampleVec])
  P_bar = (A_sum - A_bar)/(M - m);
  A_bar = A_bar/m;
  
  P_hat = regularize(ase.Ahat(diag_aug(A_bar), d, isSVD))
  
  return(c(norm(P_bar - A_bar, "F")/n/(n-1), norm(P_bar - P_hat, "F")/n/(n-1), d))
}



# ASE using SVD or eigen-decomposition.
ase <- function(A, dim, isSVD=1){
  if (isSVD) {
    if(nrow(A) >= 400){
      require(irlba)
      A.svd = irlba(A, nu = dim, nv = dim)
      A.values = A.svd$d
      A.lvectors = A.svd$u
      A.rvectors = A.svd$v
    } else{
      A.svd = svd(A)
      A.values = A.svd$d[1:dim]
      A.lvectors = A.svd$u[,1:dim]
      A.rvectors = A.svd$v[,1:dim]
    }
  } else {
    if(nrow(A) >= 400){
      require(rARPACK)
      A.eig = eigs_sym(matrix(A, ncol=dim(A)[1]), dim, which = "LA")
      A.values = A.eig$values
      A.lvectors = A.eig$vectors
      A.rvectors = A.lvectors
    } else{
      A.eig = eigen(A, symmetric = T)
      A.values = A.eig$values[1:dim]
      A.lvectors = A.eig$vectors[,1:dim]
      A.rvectors = A.lvectors
    }
  }
  return(list(A.values, A.rvectors, A.lvectors))
}


# ASE return xhat
ase.x <- function(A, dim, isSVD=1){
  A.ase = ase(A, dim, isSVD)
  if(dim == 1)
    A.x = sqrt(A.ase[[1]]) * A.ase[[2]]
  else
    A.x <- A.ase[[2]] %*% diag(sqrt(A.ase[[1]]))
  return(A.x)
}


# ASE return Ahat
ase.Ahat <- function(A, dim, isSVD=1){
  A.ase = ase(A, dim, isSVD)
  if(dim == 1)
    Ahat = A.ase[[1]] * A.ase[[3]] %*% t(A.ase[[2]])
  else
    Ahat <- A.ase[[3]] %*% diag(A.ase[[1]]) %*% t(A.ase[[2]])
  return(Ahat)
}



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


mlqe_exp_fun <- function(x, theta, q) {
  return ((x - theta)/theta^2*(exp(-x/theta)/theta)^(1-q))
}

mlqe_exp_solver <- function(xVec, q, tol=1e-6) {
  if (q == 1) {
    return(mean(xVec))
  }
  thetaMin = min(xVec)
#   thetaMax = max(xVec)
  thetaMax = mean(xVec)
  if (thetaMin == thetaMax) {
    return(thetaMin)
  }
  theta = (thetaMin + thetaMax)/2
  sumTmp = sum(sapply(xVec, mlqe_exp_fun, theta, q))
  sumTmp0 = sumTmp
  sumBase = abs(sumTmp)
  maxIter = 100
  iter = 0
  while (abs(sumTmp) > tol*sumBase) {
    iter = iter + 1
    if (sumTmp > 0) {
      thetaMin = theta
    } else {
      thetaMax = theta
    }
    theta = (thetaMin + thetaMax)/2
    sumTmp = sum(sapply(xVec, mlqe_exp_fun, theta, q))
    if (is.nan(sumTmp)) {
      return(mean(xVec))
    }
    if ((iter >= maxIter) && (abs(sumTmp) > sumTmp0)) {
      return(mean(xVec))
    }
    if (iter >= maxIter) {
      return(theta)
    }
    #     sumTmp1 = sum(sapply(xVec, mlqe_exp_fun, (thetaMin + theta)/2, q))
    #     sumTmp2 = sum(sapply(xVec, mlqe_exp_fun, (thetaMax + theta)/2, q))
    #     if ((abs(sumTmp1) < abs(sumTmp)) && (abs(sumTmp2) < abs(sumTmp))) {
    #       if (sumTmp > 0) {
    #         thetaMin = theta
    #       } else {
    #         thetaMax = theta
    #       }
    #     } else if (abs(sumTmp1) < abs(sumTmp)) {
    #       thetaMax = theta
    #     } else {
    #       thetaMin = theta
    #     }
    #     theta = (thetaMin + thetaMax)/2
    #     sumTmp = sum(sapply(xVec, mlqe_exp_fun, theta, q))
  }
  return(theta)
}

sim_all <- function(m, n, tau, B, CB, eps, q, d, isSVD=1) {
  result = rep(NaN, 4)
  
  P = B[tau, tau]
  C = CB[tau, tau]
  diag(P) = 0
  diag(C) = 0
  
  A_all = array(0, dim=c(m, n, n))
  ind = lower.tri(A_all[1,,], 1)
  #   A_all1 = list()
  for (i in 1:m) {
    contamVec = (runif(n^2) > eps)
    A = contamVec*matrix(rexp(n^2, 1/P), ncol=n) +
      (1-contamVec)*matrix(rexp(n^2, 1/C), ncol=n)
    A[ind] = 0
    A = A + t(A)
    A_all[i,,] = A
    #     A_all1[[i]] = A
  }
  
  A_bar = apply(A_all, MARGIN=c(2, 3), sum)/m
  #   A_bar = add(A_all1)/m
  A.ase = ase(diag_aug(A_bar, d), d, isSVD)
  if (d == 1) {
    A_bar_ase = A.ase[[1]][1] * A.ase[[3]][,1:d] %*% t(A.ase[[2]][,1:d])
  } else {
    A_bar_ase = A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
  }
  A_bar_ase = regularize(A_bar_ase)
  
  #   tmp = apply(A_all, 1, mlqe_exp_solver, q)
  #   A_mlqe = array(0, dim=c(n, n))
  #   for (i in 1:(n-1)) {
  #     print(i)
  #     for (j in (i+1):n) {
  #       A_mlqe[i,j] = mlqe_exp_solver(A_all[,i,j], q)
  #       A_mlqe[j,i] = A_mlqe[i,j]
  #     }
  #   }
  A_mlqe = apply(A_all, c(2, 3), mlqe_exp_solver, q)
  
  A.ase = ase(diag_aug(A_mlqe, d), d, isSVD)
  if (d == 1) {
    A_mlqe_ase = A.ase[[1]][1] * A.ase[[3]][,1:d] %*% t(A.ase[[2]][,1:d])
  } else {
    A_mlqe_ase = A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
  }
  A_mlqe_ase = regularize(A_mlqe_ase)
  
  result[1] = norm(A_bar - P, "F")/n/(n-1)
  result[2] = norm(A_bar_ase - P, "F")/n/(n-1)
  result[3] = norm(A_mlqe - P, "F")/n/(n-1)
  result[4] = norm(A_mlqe_ase - P, "F")/n/(n-1)
  
  return(result)
}



dim_brute_fullrank <- function(m, dVec, P, isSVD=1, contamVec=1) {
  
  dMax = max(dVec)
  nD = length(dVec)
  result = rep(NaN, nD+1)
  
  require(igraph)
  A_all = list()
  for (i in 1:m) {
    A = contamVec*matrix(rexp(n^2, 1/P), ncol=n) +
      (1-contamVec)*matrix(rexp(n^2, 1/P/10), ncol=n)
    ind = lower.tri(A, 1)
    A[ind] = 0
    A = A + t(A)
    A_all[[i]] = A
  }
  
  A_bar = add(A_all)/m
  result[1] = norm(P - A_bar, "F")/n/(n-1)
  
  A.ase = ase(diag_aug(A_bar), dMax, isSVD)
  
  for (iD in 1:nD) {
    d = dVec[iD]
    #     A.ase = ase(diag_aug(A_bar, d), d, isSVD)
    if (d == 1)
      Ahat = A.ase[[1]][1] * A.ase[[3]] %*% t(A.ase[[2]])
    else
      Ahat <- A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
    P_hat = regularize(Ahat)
    result[iD+1] = norm(P - P_hat, "F")/n/(n-1)
  }
  
  return(result)
}



test_rank1 <- function(m, dVec, P, isSVD=1) {
  
  dMax = max(dVec)
  nD = length(dVec)
  result = rep(NaN, nD)
  
  require(igraph)
  A_all = list()
  for (i in 1:m) {
    g = sample_sbm(n, P, rep(1,n), directed=F, loops=F)
    A = as_adj(g, type="both", sparse=FALSE)
    A_all[[i]] = A
  }
  
  A_bar = add(A_all)/m
  
  A.ase = eigen(diag_aug(A_bar), dMax)
  P.ase = eigen(diag_aug(P), dMax)
  
  for (iD in 1:nD) {
    d = dVec[iD]
    A_hat = A.ase$values[d] * A.ase$vectors[,d] %*% t(A.ase$vectors[,d])
    P_hat = P.ase$values[d] * P.ase$vectors[,d] %*% t(P.ase$vectors[,d])
    result[iD] = norm(A_hat - P_hat, "F")/n/(n-1)
  }
  
  return(result)
}





MLqEObjLognormal <- function(theta, data, q) {
  mu <- theta[1]
  sigma <- theta[2]^2
  fVal <- 1/sigma/sqrt(2*pi)/data*exp(-(log(data) - mu)^2/2/(sigma^2))
  fqVal <- fVal^(1 - q)
  val <- sum(fqVal*(log(data) - mu))^2 +
    sum(fqVal*(sigma - (log(data) - mu)^2/(sigma)))^2
  return(val)
}

MLqEObjLognormal1 <- function(theta, data, q) {
  mu <- theta[1]
  sigma <- theta[2]^2
  fVal <- 1/sigma/sqrt(2*pi)/data*exp(-(log(data) - mu)^2/2/(sigma^2))
  if (q == 1) {
    val <- sum(log(fVal))
  } else {
    val <- sum((fVal^(1 - q) - 1)/(1 - q))
  }
  return(-val)
}


MLqESolverLognormal <- function(data, q) {
  muHat <- mean(log(data))
  sigmaHat <- sqrt(mean((log(data) - muHat)^2))
  # To avoid the constraint that sigma > 0, we optimize over sqrt(sigma)
  thetaInit <- c(muHat, sqrt(sigmaHat))
  if (sigmaHat > 0) {
    resultOpti <- optim(par = thetaInit, fn = MLqEObjLognormal1, data = data, q = q,
                        method = "BFGS")
    muHat <- resultOpti$par[1]
    sigmaHat <- resultOpti$par[2]^2
  }
  val <- exp(muHat + sigmaHat^2/2)
  return(val)
}



LognormalAllDim <- function(M, m, dVec, AList, ASum, q, isSVD=1) {
  source("getElbows.R")
  source("USVT.R")
  
  nD <- length(dVec)
  dMax <- max(dVec)
  result <- rep(NaN, 2*nD+6)
  
  sampleVec <- sample.int(M, m)
  ABar <- add(AList[sampleVec])
  PBar <- (ASum - ABar)/(M - m)
  #   PBar <- ASum/M
  ABar <- ABar/m
  result[1] <- (norm(PBar - ABar, "F"))^2/n/(n-1)
  
  n <- dim(ABar)[[1]]
  ATensor <- array(unlist(AList[sampleVec]), dim = c(n, n, m))
  AMLqE <- apply(ATensor, c(1, 2), MLqESolverLognormal, q)
  result[nD + 2] <- (norm(PBar - AMLqE, "F"))^2/n/(n-1)
  
  A_bar_diag_aug = diag_aug(ABar)
  # ZG
  nElbow = 3
  evalVec = ase(A_bar_diag_aug, ceiling(n*3/5), isSVD)[[1]]
  dZG = getElbows(evalVec, n=nElbow, plot=F)[[nElbow]]
  # USVT
  dUSVT = length(usvt(A_bar_diag_aug, 1, m)$d)
  result[nD*2+3] = dZG
  result[nD*2+4] = dUSVT
  
  A.ase = ase(A_bar_diag_aug, dMax, isSVD)
  for (iD in 1:nD) {
    d = dVec[iD]
    if (d == 1)
      Ahat = A.ase[[1]] * A.ase[[3]] %*% t(A.ase[[2]])
    else
      Ahat <- A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
    P_hat = regularize(Ahat)
    result[1+iD] = norm(PBar - P_hat, "F")/n/(n-1)
  }
  #   plot(result[[2]])
  
  A_mlqe_diag_aug = diag_aug(AMLqE)
  # ZG
  nElbow = 3
  evalVec = ase(A_mlqe_diag_aug, ceiling(n*3/5), isSVD)[[1]]
  dZG = getElbows(evalVec, n=nElbow, plot=F)[[nElbow]]
  # USVT
  dUSVT = length(usvt(A_mlqe_diag_aug, 1, m)$d)
  result[nD*2+5] = dZG
  result[nD*2+6] = dUSVT
  
  A.ase = ase(A_mlqe_diag_aug, dMax, isSVD)
  for (iD in 1:nD) {
    d = dVec[iD]
    if (d == 1)
      Ahat = A.ase[[1]] * A.ase[[3]] %*% t(A.ase[[2]])
    else
      Ahat <- A.ase[[3]][,1:d] %*% diag(A.ase[[1]][1:d]) %*% t(A.ase[[2]][,1:d])
    P_hat_ase = regularize(Ahat)
    result[nD+2+iD] = norm(PBar - P_hat_ase, "F")/n/(n-1)
  }
  #   plot(result[[4]])
  return(result)
}



MLqEObjExp <- function(theta, data, q) {
  fVal <- 1/theta*exp(-data/theta)
  if (q == 1) {
    val <- sum(log(fVal))
  } else {
    val <- sum((fVal^(1 - q) - 1)/(1 - q))
  }
  return(-val)
}

MLqEGradExp <- function(theta, data, q) {
  if (q == 1) {
    val <- sum(-1/theta + data/(theta^2))
  } else {
    fVal <- 1/theta*exp(-data/theta)
    val <- sum((fVal^(1 - q)*(-1/theta + data/(theta^2))))
  }
  return(val)
}


MLqESolverExp <- function(data, q) {
  thetaHat <- mean(data)
  thetaInit <- thetaHat
  if (var(data) > 0) {
#     resultOpti <- optim(par = thetaInit, fn = MLqEObjExp, data = data, q = q,
#                         gr = MLqEGradExp, method = "BFGS")
#     resultOpti <- optim(par = thetaInit, fn = MLqEObjExp, data = data, q = q,
#                         method = "BFGS")
    resultOpti <- optim(par = thetaInit, fn = MLqEObjExp, data = data, q = q)
    thetaHat <- resultOpti$par
  }
  return(thetaHat)
}


ExpAllDim <- function(M, m, dVec, AList, ASum, q, isSVD=1, PBar=NA) {
  source("getElbows.R")
  source("USVT.R")
  
  nD <- length(dVec)
  dMax <- max(dVec)
  result <- rep(NaN, 2*nD+6)
  
  sampleVec <- sample.int(M, m)
  ABar <- add(AList[sampleVec])
  if (any(is.na(PBar))) {
    PBar <- (ASum - ABar)/(M - m)
  }
  #   PBar <- ASum/M
  ABar <- ABar/m
  result[1] <- (norm(PBar - ABar, "F"))^2/n/(n-1)
  
  n <- dim(ABar)[[1]]
  ATensor <- array(unlist(AList[sampleVec]), dim = c(n, n, m))
#   AMLqE <- apply(ATensor, c(1, 2), MLqESolverExp, q)  
  AMLqE <- apply(ATensor, c(1, 2), mlqe_exp_solver, q)
  result[nD + 2] <- (norm(PBar - AMLqE, "F"))^2/n/(n-1)
  
  ABarDiagAug <- diag_aug(ABar)
  # ABarDiagAug <- ABar
  # ZG
  nElbow <- 3
  evalVec <- ase(ABarDiagAug, ceiling(n*3/5), isSVD)[[1]]
  dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
  # USVT
  dUSVT <- length(usvt(ABarDiagAug, 1, m)$d)
  result[nD*2 + 3] <- dZG
  result[nD*2 + 4] <- dUSVT
  
  AASE = ase(ABarDiagAug, dMax, isSVD)
  for (iD in 1:nD) {
    d <- dVec[iD]
    if (d == 1) {
      AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
    } else {
      AHat <- AASE[[3]][, 1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][ ,1:d])
    }
    PHat <- regularize(AHat)
    result[1 + iD] <- (norm(PBar - PHat, "F"))^2/n/(n-1)
  }

  # AASE1 = ase(ABarDiagAug, dMax, 1)
  # AASE0 = ase(ABarDiagAug, dMax, 0)
  # for (iD in 1:nD) {
  #   d <- dVec[iD]
  #   if (d == 1) {
  #     AHat1 <- AASE1[[1]]*AASE1[[3]]%*%t(AASE1[[2]])
  #     AHat0 <- AASE0[[1]]*AASE0[[3]]%*%t(AASE0[[2]])
  #   } else {
  #     AHat1 <- AASE1[[3]][, 1:d]%*%diag(AASE1[[1]][1:d])%*%t(AASE1[[2]][ ,1:d])
  #     AHat0 <- AASE0[[3]][, 1:d]%*%diag(AASE0[[1]][1:d])%*%t(AASE0[[2]][ ,1:d])
  #   }
  #   
  #   (norm(ABarDiagAug - AHat1, "F"))^2/n/(n-1)
  #   (norm(ABarDiagAug - AHat0, "F"))^2/n/(n-1)
  #   
  #   PHat1 <- regularize(AHat1)
  #   PHat0 <- regularize(AHat0)
  #   
  #   (norm(ABarDiagAug - PHat1, "F"))^2/n/(n-1)
  #   (norm(ABarDiagAug - PHat0, "F"))^2/n/(n-1)
  #   
  #   (norm(PBar - PHat1, "F"))^2/n/(n-1)
  #   (norm(PBar - PHat0, "F"))^2/n/(n-1)
  #   (norm(PBar, "F"))^2/n/(n-1)
  #   (norm(ABarDiagAug, "F"))^2/n/(n-1)
  # }
  
    
  # (norm(PBar, "F"))^2/n/(n-1)
  # (norm(ABar, "F"))^2/n/(n-1)
  # (norm(AHat, "F"))^2/n/(n-1)
  # AHat1 <- AHat
  # diag(AHat1) <- 0
  # (norm(AHat, "F"))^2/n/(n-1)
  # (norm(PHat, "F"))^2/n/(n-1)
  
  
  
  AMLqEDiagAug <- diag_aug(AMLqE)
  # AMLqEDiagAug <- AMLqE
  # ZG
  nElbow <- 3
  evalVec <- ase(AMLqEDiagAug, ceiling(n*3/5), isSVD)[[1]]
  dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
  # USVT
  dUSVT <- length(usvt(AMLqEDiagAug, 1, m)$d)
  result[nD*2 + 5] <- dZG
  result[nD*2 + 6] <- dUSVT
  
  AASE <- ase(AMLqEDiagAug, dMax, isSVD)
  for (iD in 1:nD) {
    d <- dVec[iD]
    if (d == 1) {
      AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
    } else {
      AHat <- AASE[[3]][ ,1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][ ,1:d])
    }
    PHatASE <- regularize(AHat)
    result[nD + 2 + iD] <- (norm(PBar - PHatASE, "F"))^2/n/(n-1)
  }
  return(result)
}





ExpAllDimAug <- function(M, m, dVec, AList, ASum, q, isSVD=1) {
  source("getElbows.R")
  source("USVT.R")
  
  nD <- length(dVec)
  dMax <- max(dVec)
  result <- rep(NaN, 2*nD+6)
  
  sampleVec <- sample.int(M, m)
  ABar <- add(AList[sampleVec])
  PBar <- (ASum - ABar)/(M - m)
  #   PBar <- ASum/M
  ABar <- ABar/m
  tmpDiff <- PBar - ABar
  diag(tmpDiff) <- 0
  result[1] <- (norm(tmpDiff, "F"))^2/n/(n-1)
  
  n <- dim(ABar)[[1]]
  ATensor <- array(unlist(AList[sampleVec]), dim = c(n, n, m))
  #   AMLqE <- apply(ATensor, c(1, 2), MLqESolverExp, q)  
  AMLqE <- apply(ATensor, c(1, 2), mlqe_exp_solver, q)
  tmpDiff <- PBar - AMLqE
  diag(tmpDiff) <- 0
  result[nD + 2] <- (norm(tmpDiff, "F"))^2/n/(n-1)
  
  #   ABarDiagAug <- diag_aug(ABar)
  ABarDiagAug <- ABar
  # ZG
  nElbow <- 3
  evalVec <- ase(ABarDiagAug, ceiling(n*3/5), isSVD)[[1]]
  dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
  # USVT
  dUSVT <- length(usvt(ABarDiagAug, 1, m)$d)
  result[nD*2 + 3] <- dZG
  result[nD*2 + 4] <- dUSVT
  
  AASE = ase(ABarDiagAug, dMax, isSVD)
  for (iD in 1:nD) {
    d <- dVec[iD]
    if (d == 1) {
      AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
    } else {
      AHat <- AASE[[3]][, 1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][ ,1:d])
    }
    PHat <- regularize(AHat)
    tmpDiff <- PBar - PHat
    diag(tmpDiff) <- 0
    result[1 + iD] <- (norm(tmpDiff, "F"))^2/n/(n-1)
  }
  
  #   AMLqEDiagAug <- diag_aug(AMLqE)
  AMLqEDiagAug <- AMLqE
  # ZG
  nElbow <- 3
  evalVec <- ase(AMLqEDiagAug, ceiling(n*3/5), isSVD)[[1]]
  dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
  # USVT
  dUSVT <- length(usvt(AMLqEDiagAug, 1, m)$d)
  result[nD*2 + 5] <- dZG
  result[nD*2 + 6] <- dUSVT
  
  AASE <- ase(AMLqEDiagAug, dMax, isSVD)
  for (iD in 1:nD) {
    d <- dVec[iD]
    if (d == 1) {
      AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
    } else {
      AHat <- AASE[[3]][ ,1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][ ,1:d])
    }
    PHatASE <- regularize(AHat)
    tmpDiff <- PBar - PHatASE
    diag(tmpDiff) <- 0
    result[nD + 2 + iD] <- (norm(tmpDiff, "F"))^2/n/(n-1)
  }
  return(result)
}








ExpAllDimRatio <- function(M, m, d = 0, AList, ASum, q, isSVD = 1) {
  source("getElbows.R")
  
  dMax <- max(dVec)
  result <- rep(NaN, 2*nD+6)
  
  sampleVec <- sample.int(M, m)
  sampleVecComplement <- (1:M)[is.na(pmatch(1:M, sampleVec))]
  ABar <- add(AList[sampleVec])/m
  
  n <- dim(ABar)[[1]]
  ATensor <- array(unlist(AList[sampleVec]), dim = c(n, n, m))
  #   AMLqE <- apply(ATensor, c(1, 2), MLqESolverExp, q)  
  AMLqE <- apply(ATensor, c(1, 2), mlqe_exp_solver, q)
  
  ABarDiagAug <- diag_aug(ABar)
  if (d == 0) {
    # ZG
    nElbow <- 3
    evalVec <- ase(ABarDiagAug, ceiling(n*3/5), isSVD)[[1]]
    d <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
  }
  
  AASE <- ase(ABarDiagAug, d, isSVD)
  if (d == 1) {
    AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
  } else {
    AHat <- AASE[[3]][, 1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][, 1:d])
  }
  ABarASE <- regularize(AHat)
  
  AMLqEDiagAug <- diag_aug(AMLqE)
  if (d == 0) {
    # ZG
    nElbow <- 3
    evalVec <- ase(AMLqEDiagAug, ceiling(n*3/5), isSVD)[[1]]
    d <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
  }
  
  AASE <- ase(AMLqEDiagAug, d, isSVD)
  if (d == 1) {
    AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
  } else {
    AHat <- AASE[[3]][, 1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][, 1:d])
  }
  PHatASE <- regularize(AHat)
  
  ratioList <- lapply(sampleVecComplement, function(i) {
    return(abs(AList[[i]] - ABar) == abs(AList[[i]] - AMLqE))})
  ratioAMLqEABar <- add(ratioList)/(M - m)
  image(Matrix(ratioAMLqEABar))
  sum(ratioAMLqEABar)/n/(n-1)
  
  return(result)
}






ExpAllDimSingleAug <- function(M, m, dVec, AList, ASum, AList0, ASum0, q, isSVD=1, PBar=NA) {
  source("getElbows.R")
  source("USVT.R")
  
  nD <- length(dVec)
  dMax <- max(dVec)
  result <- rep(NaN, 4*nD+8)
  
  sampleVec <- sample.int(M, m)
  ABar <- add(AList[sampleVec])
  ABar0 <- add(AList0[sampleVec])
  if (any(is.na(PBar))) {
    PBar <- (ASum - ABar)/(M - m)
    PBar0 <- (ASum0 - ABar0)/(M - m)
  }
  #   PBar <- ASum/M
  ABar <- ABar/m
  result[1] <- (norm(PBar - ABar, "F"))^2/n/(n-1)
  result[2*nD+7] <- (norm(PBar0 - ABar, "F"))^2/n/(n-1)
  
  n <- dim(ABar)[[1]]
  ATensor <- array(unlist(AList[sampleVec]), dim = c(n, n, m))
  #   AMLqE <- apply(ATensor, c(1, 2), MLqESolverExp, q)  
  AMLqE <- apply(ATensor, c(1, 2), mlqe_exp_solver, q)
  result[nD + 2] <- (norm(PBar - AMLqE, "F"))^2/n/(n-1)
  result[3*nD+8] <- (norm(PBar0 - AMLqE, "F"))^2/n/(n-1)
  
  ABarDiagAug <- diag_aug(ABar)
  # ABarDiagAug <- ABar
  # ZG
  nElbow <- 3
  evalVec <- ase(ABarDiagAug, ceiling(n*3/5), isSVD)[[1]]
  dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
  # USVT
  dUSVT <- length(usvt(ABarDiagAug, 1, m)$d)
  result[nD*2 + 3] <- dZG
  result[nD*2 + 4] <- dUSVT
  
  AASE = ase(ABarDiagAug, dMax, isSVD)
  for (iD in 1:nD) {
    d <- dVec[iD]
    if (d == 1) {
      AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
    } else {
      AHat <- AASE[[3]][, 1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][ ,1:d])
    }
    PHat <- regularize(AHat)
    result[1 + iD] <- (norm(PBar - PHat, "F"))^2/n/(n-1)
    result[2*nD + 7 + iD] <- (norm(PBar0 - PHat, "F"))^2/n/(n-1)
  }
  
  AMLqEDiagAug <- diag_aug(AMLqE)
  # AMLqEDiagAug <- AMLqE
  # ZG
  nElbow <- 3
  evalVec <- ase(AMLqEDiagAug, ceiling(n*3/5), isSVD)[[1]]
  dZG <- getElbows(evalVec, n = nElbow, plot = F)[[nElbow]]
  # USVT
  dUSVT <- length(usvt(AMLqEDiagAug, 1, m)$d)
  result[nD*2 + 5] <- dZG
  result[nD*2 + 6] <- dUSVT
  
  AASE <- ase(AMLqEDiagAug, dMax, isSVD)
  for (iD in 1:nD) {
    d <- dVec[iD]
    if (d == 1) {
      AHat <- AASE[[1]]*AASE[[3]]%*%t(AASE[[2]])
    } else {
      AHat <- AASE[[3]][ ,1:d]%*%diag(AASE[[1]][1:d])%*%t(AASE[[2]][ ,1:d])
    }
    PHatASE <- regularize(AHat)
    result[nD + 2 + iD] <- (norm(PBar - PHatASE, "F"))^2/n/(n-1)
    result[3*nD + 8 + iD] <- (norm(PBar0 - PHatASE, "F"))^2/n/(n-1)
  }
  return(result)
}
