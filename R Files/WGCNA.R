#### SECTION 1 ####

#https://kateto.net/network-visualization

rm(list=ls()) 
setwd("M:/w2k/BigDataBio/Brassica Report")  
load("brassicaData.rda")
load("candidateAnalysis.rda")
library(igraph)

#### SECTION 2 ####

## 2.1 Constructing corrrelation matrix ##
covariance <- cor(t(rpkm[candidates,]), method="pearson")
correlation <- cov2cor(covariance)

## 2.2 Forming adjacency matrix ##
threshold <- 0.5 
adjacency <- abs(covariance)
adjacency[adjacency<threshold] <- 0

## 2.3 Constructing co-expression network ##
network <- graph_from_adjacency_matrix(adjacency,
                                       diag = F, "undirected", weighted=T)

## 2.4 Finding optimal clusters within candidate genes ##
clust <- cluster_optimal(network)

#### SECTION 3 ####

## 3.1 Setting node and edge sizes ##
V(network)$size <- 10+(1.2*degree(network))
E(network)$width <- 10*(E(network)$weight*E(network)$weight)

## 3.2 Partitioning genes into clusters ##
clusters <- clust$membership
list_clust <- list()
for(i in 1:max(clusters)){
  list_clust[i] <- list((which(clusters==i)))
}

## 3.3 Node, edge and cluster colors ##
V(network)$color <- rep("#42f72a", vcount(network))
V(network)[c(suppressor1, suppressor2)]$color <- "#ff737e"
  col.clusters <- c("#c5e5e7", "#e8d3d4", "#ead8f2", "#ecd89a")

## 3.4 Plotting network graph ##
plot(network,
     edge.color = "gray40",
     vertex.label.family = "Helvetica",
     vertex.label.cex = 0.7,
     vertex.label = c(1:21),
     mark.groups = list_clust[1:4],
     mark.col = col.clusters,
     mark.border = NA)
legend(x= -1, y=-1.1, c("Candidates", "Candidate 'suppressors'"), pch=21,
       col="#777777", pt.bg=c("#42f72a", "#ff737e"), pt.cex=2, cex=.8, bty="n", ncol=1)

#### SECTION 4 ####

save(adjacency, network, list_clust, clust, threshold, file="WGCNA.rda")
load("WGCNA.rda")
