# rm(list=ls())
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("cluster", version = "3.8")
library(edgeR)
library(DESeq2)
library(stringr)
library(ggrepel)
 
 library(org.EcK12.eg.db)
 library(DOSE)
 library(ggraph)
 library(clusterProfiler)
 # library(MetaboSignal)

library(Rtsne) # Load package

#GINCKO DATA
count_data <- read.csv(file="data/count_data_gincko.csv", header=T,  sep=",")

sum(rowSums(is.na(count_data)))

data2 <- count_data[,-1]
rownames(data2) <- count_data[,1]

count_data <- data2
rm(data2)

dim(count_data)
tail(count_data)

#Filtering genes that have too many zeros
checkzeros <- count_data == 0
nullcounts <- rowSums(checkzeros)
allcounts <- rowSums(count_data)

plot(nullcounts, allcounts, log = "xy")

count_data <- count_data[nullcounts < 35 & allcounts > 8800,]
dim(count_data)

hist(nullcounts, labels=F, col="gray", right=F, breaks=seq(min(nullcounts),max(nullcounts),l=max(nullcounts)-min(nullcounts)))


RM <- colSums(count_data)
hist(RM, breaks=seq(min(RM),max(RM),l=50), freq=F)
RM <- colSums(cpm_data)
hist(RM, breaks=seq(min(RM),max(RM),l=50), freq=F)

meta_data <- read.csv(file="data/meta_data_gincko.csv", header=T,  sep=",")
meta_data <- meta_data[,2:8]
meta_data_original <- meta_data
#View(meta_data_original)

meta_data[,2] <- str_sub(meta_data[,2], 8, 8)
meta_data[,3] <- str_sub(meta_data[,3], 2, str_length(meta_data[,3])-8)
meta_data[,5] <- str_sub(meta_data[,5], 1, str_length(meta_data[,5])-3)
meta_data[,6] <- str_sub(meta_data[,6], 1, str_length(meta_data[,6])-3)
meta_data[,7] <- str_sub(meta_data[,7], 2, str_length(meta_data[,7])-5)

newIDs <- paste("S", meta_data[,2], sep = "")
newIDs <- paste(newIDs, meta_data[,3], sep = "_T")
newIDs <- paste(newIDs, meta_data[,7], sep = "_t")
newIDs <- paste(newIDs, meta_data[,5], sep = "_I")
newIDs <- paste(newIDs, meta_data[,6], sep = "_A")
newIDs <- paste(newIDs, meta_data[,4], sep = "_r")

rownames(meta_data) <- newIDs
colnames(count_data) <- newIDs

#View(meta_data)

dim(count_data)
#####SELECT ACTIVATED time = 5 and  time = 18 separately

search_string <- "t5"
straincoord <- grep(search_string, colnames(count_data))
count_data_t5 <- count_data[,straincoord]
dim(count_data_t5)

search_string <- "t18"
straincoord <- grep(search_string, colnames(count_data))
count_data_t18 <- count_data[,straincoord]
dim(count_data_t18)


cpm_data_t5 <- cpm(count_data_t5, log=TRUE)
cpm_data_t18 <- cpm(count_data_t18, log=TRUE)

dim(count_data_t5)+dim(count_data_t18)
dim(cpm_data_t5)+dim(cpm_data_t18)




# pca <- prcomp(t(count_data_t18), scale. = TRUE)
# plot(pca$x[, 1], pca$x[, 2], pch = 20, xlab = "PC1", ylab = "PC2", main="PCA")
# 

# ### PLots t-SNE
# library(tsne) # Load package
# 
# tsne_model_1 <- tsne(count_data_t5, initial_config = NULL, k=2, initial_dims=30, perplexity=30, max_iter = 200, min_cost=0, epoch_callback=NULL,whiten=TRUE, epoch=100 )
# plot(tsne_model_1, asp=1) # Plot the result
# 
# ### PLots t-SNE###################################
# set.seed(42) # Sets seed for reproducibility
# tsne_out <- Rtsne(as.matrix(t(count_data_t5)), dims = 2, perplexity = 1) # Run TSNE
# plot(tsne_out$Y, asp=1) # Plot the result
# text(tsne_out$Y, row.names(t(count_data_t5)), cex=0.4, pos=4, col="blue")
# 
# tsne <- tsne_out$Y
# rownames(tsne) <- rownames(t(count_data_t5))
# 
# ggplot(as.data.frame(tsne), aes(tsne[,1], tsne[,2], label = rownames(tsne))) +
#   geom_text_repel() +
#   geom_point(color = 'red') +
#   theme_classic(base_size = 18) +
#   ggtitle("t-sne embedding of E. coli gene expression")
# 
# 
# ############################################################
# 
