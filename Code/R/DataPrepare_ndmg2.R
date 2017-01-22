rm(list = ls())
setwd("/Users/Runze/Documents/GitHub/RobustASE/Data/desikan_ndmg2/")
require(igraph)

dataName <- "ndmg2"
fileName <- paste0("../", dataName, ".Rbin")

files <- list.files(pattern = "\\.graphml$")
fibergraph.list <- list()
for (iFile in 1:length(files)) {
  fibergraph.list[[iFile]] <- as_adj(read_graph(files[iFile], format="graphml"),
                                     type = "both", sparse = 1)
}

save(fibergraph.list, file=fileName)