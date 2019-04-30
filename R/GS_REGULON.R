## Load package
library(igraph)


regulon <- read.table(file="data/network_tf_gene.txt", sep="\t", header=FALSE)

#View(regulon)

edge_list <- regulon[,1:2]
edge_list <- unique(edge_list)

genes <- c(as.character(regulon[,1]), as.character(regulon[,2]))
genes <- unique(genes)

tfs <- as.character(regulon[,1])
tfs <- unique(tfs)
length(tfs)

first_genes <- genes[1:length(tfs)]

length(intersect(tfs, first_genes))

# search_string <- "-"
# coord <- grep(search_string,tfs)
# tfs[coord]
# 
# search_string2 <- "H-NS"
# coord2 <- grep(search_string2,tfs)
# tfs[coord2]
# 
# coord <- setdiff(coord,coord2)
# 
# coord <- grep(search_string,edge_list[,1])
# coord2 <- grep(search_string2,edge_list[,1])
# coord <- setdiff(coord,coord2)
# edge_list[coord,]
# 
# length(coord)
# length(coord2)
# length(edge_list[,1])


#######FIND INTERSECTION IN SETS OF GENES IN EXPRESSION DATA AND TRANSCRIPTIONAL NETWORK
genes_data <- rownames(count_data)
genes_trn <- genes

length(genes_data)
length(genes_trn)

genes_intersect <- intersect(genes_data, genes_trn)
tfs_intersect <- intersect(genes_data, tfs)


length(genes_intersect)
length(tfs_intersect)






net <- graph_from_data_frame(d=edge_list, directed=T)

plot(net, layout=layout_with_kk, edge.arrow.size=.3, vertex.size = 2, vertex.label=NA)

adjMat <- get.adjacency(net)

genenames <- colnames(adjMat)
length(genenames)
sum(adjMat)

genenumbers <- 0:(length(genenames)-1)

rownames(adjMat) <- genenumbers
colnames(adjMat) <- genenumbers

net2 <-graph_from_adjacency_matrix(adjMat)
head(net2)

numbernet <- get.edgelist(net2)
class(numbernet)

numbernet <- as.data.frame(numbernet)

write.csv(numbernet, file = "EcoliGS.csv")

