library(igraph)
library(corpustools)


x=  cor(as.matrix(t(gex.M)))
x[1:10, 1:10]

vec= sort(x["Angptl4",],decreasing = T)
vec[1:150]
hist(vec, breaks = 100)
saveRDS(x, "output/gene.corr.matrix.rds")
x= readRDS("output/gene.corr.matrix.rds")

d= graph_from_adjacency_matrix(x, mode = "undirected", weighted= T)
rm(x)

saveRDS(d, "output/gene.corr.net.rds")

hpo_net= backbone_filter(d, alpha =0.05 )
