#     dataName = "CPAC200"
dataName = "desikan"
#     dataName = "JHU"

source("function_collection.R")
DA = F
fileName = paste("../../Data/data_", dataName, ".RData", sep="")
# if (file.exists(fileName)) {
#   load(fileName)
#   return(list(A_all, n, M))
# } else {
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
      A = as_adj(g, attr="weight", type="both", sparse=FALSE)
      if (DA) {
        A = diag_aug(A)
      }
      A_all[[(sub-1)*2 + session]] = A;
    }
  }
  
  #   save(A_all, n, M, file=fileName)
  #   return(list(A_all, n, M))
# }  

i = 1
j = 2
data = c()
for (iter in 1:length(A_all)) {
  data[iter] = A_all[[iter]][i,j]
}

hist(data)

mlqe_poisson_function = function(x, theta, q) {
  tmp = 1
  for (i in 1:x) {
    tmp = tmp*theta/i
    if (i <= theta) {
      tmp = tmp*exp(-1)
    }
  }
  tmp = tmp*exp(-(theta - floor(theta)))
  tmp = tmp^(1-q)*(x - theta)
  return(tmp)
}

# mlqe_gaussian_function = function(x, theta, q) {
#   tmp = -1/sqrt(2*pi)*exp(-(x-theta)^2/2/theta)*
#     (theta^(-3/2)/2 - theta^(-5/2)/2*(x^2-theta^2))
#   return(tmp)
# }

mlqe_gaussian_function = function(x, theta, q) {
  tmp = -1/2/theta + (x - theta)/theta + (x-theta)^2/2/theta^2
  tmp = tmp*(dnorm(x, theta, sqrt(theta)))^(1-q)
  return(tmp)
}

mlqe_poisson_solver = function(xVec, q) {
  thetaMin = min(xVec)
  thetaMax = max(xVec)
  theta = mean(xVec)
  
  fTmp = sum(sapply(xVec, mlqe_gaussian_function, theta, q))
  while (abs(fTmp) > tol) {
    if (fTmp < 0) {
      
    } else {
      
    }
    fTmp = sum(sapply(xVec, mlqe_gaussian_function, theta, q))
  }
}

tmp = c()
for (theta in 1:10000) {
  print(theta)
  tmp[theta] = sum(sapply(xVec, mlqe_gaussian_function, theta, q))
}