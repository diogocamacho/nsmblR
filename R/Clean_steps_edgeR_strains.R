rm(list=ls())

# load required libraries
library(edgeR)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.EcK12.eg.db)
library(DOSE)
library(pathview)
library(purrr)
library(ggraph)
library(clusterProfiler)
library(annotables)
library(MetaboSignal)


####SELECTING DATA TO CONTRAST S2 VS S1

count_data_strain <- cbind(count_data_strain1, count_data_strain2)
dim(count_data_strain)

groups <- colnames(count_data_strain)
group <- str_sub(groups, 1, str_length(groups)-6)
group

y <- DGEList(counts = count_data_strain, group = factor(group))
y <- calcNormFactors(y)
y <- estimateDisp(y)

y$samples

#######################STRAIN EFFECTS#############################
S2vsS1_T30_t18 <- exactTest(y, pair=c(1,5))
S2vsS1_T30_t5 <- exactTest(y, pair=c(2,6))
S2vsS1_T37_t18  <- exactTest(y, pair=c(3,7))
S2vsS1_T37_t5 <- exactTest(y, pair=c(4,8))

deS2vsS1_T30_t18 <- decideTestsDGE(S2vsS1_T30_t18, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS2vsS1_T30_t18)
deS2vsS1_T30_t5 <- decideTestsDGE(S2vsS1_T30_t5, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS2vsS1_T30_t5)
deS2vsS1_T37_t18 <- decideTestsDGE(S2vsS1_T37_t18, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS2vsS1_T37_t18)
deS2vsS1_T37_t5 <- decideTestsDGE(S2vsS1_T37_t5, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS2vsS1_T37_t5)
##################################################################

#####4-way Venn
ytags1 <- rownames(deS2vsS1_T30_t18)[as.logical(deS2vsS1_T30_t18)] 
ytags2 <- rownames(deS2vsS1_T30_t5)[as.logical(deS2vsS1_T30_t5)]
ytags3 <- rownames(deS2vsS1_T37_t18)[as.logical(deS2vsS1_T37_t18)]
ytags4 <- rownames(deS2vsS1_T37_t5)[as.logical(deS2vsS1_T37_t5)]

# Create Venn-Diagram
A <- ytags1
B <- ytags2
C <- ytags3
D <- ytags4

AB <- list(t(A), t(B), t(C), t(D))

library(VennDiagram)
pdf("venn_diagram.pdf")
venn.plot <- venn.diagram(AB, NULL, fill=c("red", "green", "blue", "yellow"), alpha=c(0.5,0.5, 0.5, 0.5), cex = 2, cat.fontface=4, category.names=c("T30_t18", "T30_t5", "T37_t18", "T37_t5"))
grid.draw(venn.plot)
dev.off()


###########################FUNCTIONAL ENRICHMENTS#############################

####BACKGROUND TO RE-USE, RUN ONLY ONCE
background <- mapIds(org.EcK12.eg.db, rownames(count_data_strain), 'ENTREZID', 'SYMBOL')
background <- background[!is.na(background)]
background <- matrix(background)
dim(background)
background_kegg <- MS_convertGene(background, "eco", output = "vector", orthology = F)
background_kegg <- matrix(background_kegg)
dim(background_kegg)
background_kegg <- str_sub(background_kegg, 5, 9)

########################
intersect_list <- intersect(ytags1, ytags2)
intersect_list <- intersect(intersect_list, ytags3)
intersect_list <- intersect(intersect_list, ytags4)
dim(as.matrix(intersect_list))

union_list <- union(ytags1, ytags2)
union_list <- union(union_list, ytags3)
union_list <- union(union_list, ytags4)
union_list <- unique(union_list)
dim(as.matrix(union_list))

####GENE LIST
symbols <- as.matrix(intersect_list)

querylist_entrez <- mapIds(org.EcK12.eg.db, symbols, 'ENTREZID', 'SYMBOL')
querylist_entrez <- querylist_entrez[!is.na(querylist_entrez)]
querylist_entrez <- unique(querylist_entrez)
background <- unique(background)
querylist_entrez <- matrix(querylist_entrez)
dim(querylist_entrez)
write.csv(querylist_entrez, "results/Ecoli.csv")
querylist_kegg <- MS_convertGene(querylist_entrez, "eco", output = "vector", orthology = F)
querylist_kegg <- str_sub(querylist_kegg, 5, 9)
querylist_kegg <- matrix(querylist_kegg)
dim(querylist_kegg)
####################

#Run KEGG enrichment analysis
KEGGEResults <- enrichKEGG(querylist_kegg, universe = background_kegg, organism="eco", pvalueCutoff=0.05)
## Output results
write.csv(KEGGEResults, "results/NEWKEGGResults.csv", quote=F)


## Run GO enrichment analysis
egoResults <- enrichGO(gene = querylist_entrez, keyType = "ENTREZID", universe = background, OrgDb = org.EcK12.eg.db, ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
## Output results from GO analysis to a table
write.csv(egoResults, "results/NEWGOresults.csv")

##########################################################################
####SELECTING DATA TO CONTRAST S3 VS S1

count_data_strain <- cbind(count_data_strain1, count_data_strain3)
dim(count_data_strain)

groups <- colnames(count_data_strain)
group <- str_sub(groups, 1, str_length(groups)-6)
group

y <- DGEList(counts = count_data_strain, group = factor(group))
y <- calcNormFactors(y)
y <- estimateDisp(y)

y$samples

#######################STRAIN EFFECTS#############################
S3vsS1_T30_t18 <- exactTest(y, pair=c(1,5))
S3vsS1_T30_t5 <- exactTest(y, pair=c(2,6))
S3vsS1_T37_t18  <- exactTest(y, pair=c(3,7))
S3vsS1_T37_t5 <- exactTest(y, pair=c(4,8))

deS3vsS1_T30_t18 <- decideTestsDGE(S3vsS1_T30_t18, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS3vsS1_T30_t18)
deS3vsS1_T30_t5 <- decideTestsDGE(S3vsS1_T30_t5, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS3vsS1_T30_t5)
deS3vsS1_T37_t18 <- decideTestsDGE(S3vsS1_T37_t18, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS3vsS1_T37_t18)
deS3vsS1_T37_t5 <- decideTestsDGE(S3vsS1_T37_t5, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS3vsS1_T37_t5)
##################################################################

#####4-way Venn
ytags1 <- rownames(deS3vsS1_T30_t18)[as.logical(deS3vsS1_T30_t18)] 
ytags2 <- rownames(deS3vsS1_T30_t5)[as.logical(deS3vsS1_T30_t5)]
ytags3 <- rownames(deS3vsS1_T37_t18)[as.logical(deS3vsS1_T37_t18)]
ytags4 <- rownames(deS3vsS1_T37_t5)[as.logical(deS3vsS1_T37_t5)]

# Create Venn-Diagram
A <- ytags1
B <- ytags2
C <- ytags3
D <- ytags4

AB <- list(t(A), t(B), t(C), t(D))

library(VennDiagram)
pdf("venn_diagram.pdf")
venn.plot <- venn.diagram(AB, NULL, fill=c("red", "green", "blue", "yellow"), alpha=c(0.5,0.5, 0.5, 0.5), cex = 2, cat.fontface=4, category.names=c("T30_t18", "T30_t5", "T37_t18", "T37_t5"))
grid.draw(venn.plot)
dev.off()


###########################FUNCTIONAL ENRICHMENTS#############################

########################
intersect_list <- intersect(ytags1, ytags2)
intersect_list <- intersect(intersect_list, ytags3)
intersect_list <- intersect(intersect_list, ytags4)
dim(as.matrix(intersect_list))

union_list <- union(ytags1, ytags2)
union_list <- union(union_list, ytags3)
union_list <- union(union_list, ytags4)
union_list <- unique(union_list)
dim(as.matrix(union_list))

####GENE LIST
symbols <- as.matrix(intersect_list)

querylist_entrez <- mapIds(org.EcK12.eg.db, symbols, 'ENTREZID', 'SYMBOL')
querylist_entrez <- querylist_entrez[!is.na(querylist_entrez)]
querylist_entrez <- unique(querylist_entrez)
background <- unique(background)
querylist_entrez <- matrix(querylist_entrez)
dim(querylist_entrez)
write.csv(querylist_entrez, "results/Ecoli.csv")
querylist_kegg <- MS_convertGene(querylist_entrez, "eco", output = "vector", orthology = F)
querylist_kegg <- str_sub(querylist_kegg, 5, 9)
querylist_kegg <- matrix(querylist_kegg)
dim(querylist_kegg)
####################

#Run KEGG enrichment analysis
KEGGEResults <- enrichKEGG(querylist_kegg, universe = background_kegg, organism="eco", pvalueCutoff=0.05)
## Output results
write.csv(KEGGEResults, "results/NEWKEGGResults.csv", quote=F)


## Run GO enrichment analysis
egoResults <- enrichGO(gene = querylist_entrez, keyType = "ENTREZID", universe = background, OrgDb = org.EcK12.eg.db, ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
## Output results from GO analysis to a table
write.csv(egoResults, "results/NEWGOresults.csv")



##########################################################################
####SELECTING DATA TO CONTRAST S4 VS S1

count_data_strain <- cbind(count_data_strain1, count_data_strain4)
dim(count_data_strain)

groups <- colnames(count_data_strain)
group <- str_sub(groups, 1, str_length(groups)-6)
group

y <- DGEList(counts = count_data_strain, group = factor(group))
y <- calcNormFactors(y)
y <- estimateDisp(y)

y$samples

#######################STRAIN EFFECTS#############################
S4vsS1_T30_t18 <- exactTest(y, pair=c(1,5))
S4vsS1_T30_t5 <- exactTest(y, pair=c(2,6))
S4vsS1_T37_t18  <- exactTest(y, pair=c(3,7))
S4vsS1_T37_t5 <- exactTest(y, pair=c(4,8))

deS4vsS1_T30_t18 <- decideTestsDGE(S4vsS1_T30_t18, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS4vsS1_T30_t18)
deS4vsS1_T30_t5 <- decideTestsDGE(S4vsS1_T30_t5, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS4vsS1_T30_t5)
deS4vsS1_T37_t18 <- decideTestsDGE(S4vsS1_T37_t18, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS4vsS1_T37_t18)
deS4vsS1_T37_t5 <- decideTestsDGE(S4vsS1_T37_t5, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS4vsS1_T37_t5)
##################################################################

#####4-way Venn
ytags1 <- rownames(deS4vsS1_T30_t18)[as.logical(deS4vsS1_T30_t18)] 
ytags2 <- rownames(deS4vsS1_T30_t5)[as.logical(deS4vsS1_T30_t5)]
ytags3 <- rownames(deS4vsS1_T37_t18)[as.logical(deS4vsS1_T37_t18)]
ytags4 <- rownames(deS4vsS1_T37_t5)[as.logical(deS4vsS1_T37_t5)]

# Create Venn-Diagram
A <- ytags1
B <- ytags2
C <- ytags3
D <- ytags4

AB <- list(t(A), t(B), t(C), t(D))

library(VennDiagram)
pdf("venn_diagram.pdf")
venn.plot <- venn.diagram(AB, NULL, fill=c("red", "green", "blue", "yellow"), alpha=c(0.5,0.5, 0.5, 0.5), cex = 2, cat.fontface=4, category.names=c("T30_t18", "T30_t5", "T37_t18", "T37_t5"))
grid.draw(venn.plot)
dev.off()


###########################FUNCTIONAL ENRICHMENTS#############################

########################
intersect_list <- intersect(ytags1, ytags2)
intersect_list <- intersect(intersect_list, ytags3)
intersect_list <- intersect(intersect_list, ytags4)
dim(as.matrix(intersect_list))

union_list <- union(ytags1, ytags2)
union_list <- union(union_list, ytags3)
union_list <- union(union_list, ytags4)
union_list <- unique(union_list)
dim(as.matrix(union_list))

####GENE LIST
symbols <- as.matrix(intersect_list)

querylist_entrez <- mapIds(org.EcK12.eg.db, symbols, 'ENTREZID', 'SYMBOL')
querylist_entrez <- querylist_entrez[!is.na(querylist_entrez)]
querylist_entrez <- unique(querylist_entrez)
background <- unique(background)
querylist_entrez <- matrix(querylist_entrez)
dim(querylist_entrez)
write.csv(querylist_entrez, "results/Ecoli.csv")
querylist_kegg <- MS_convertGene(querylist_entrez, "eco", output = "vector", orthology = F)
querylist_kegg <- str_sub(querylist_kegg, 5, 9)
querylist_kegg <- matrix(querylist_kegg)
dim(querylist_kegg)
####################

#Run KEGG enrichment analysis
KEGGEResults <- enrichKEGG(querylist_kegg, universe = background_kegg, organism="eco", pvalueCutoff=0.05)
## Output results
write.csv(KEGGEResults, "results/NEWKEGGResults.csv", quote=F)


## Run GO enrichment analysis
egoResults <- enrichGO(gene = querylist_entrez, keyType = "ENTREZID", universe = background, OrgDb = org.EcK12.eg.db, ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
## Output results from GO analysis to a table
write.csv(egoResults, "results/NEWGOresults.csv")


##########################################################################
####SELECTING DATA TO CONTRAST S5 VS S1

count_data_strain <- cbind(count_data_strain1, count_data_strain5)
dim(count_data_strain)

groups <- colnames(count_data_strain)
group <- str_sub(groups, 1, str_length(groups)-6)
group

y <- DGEList(counts = count_data_strain, group = factor(group))
y <- calcNormFactors(y)
y <- estimateDisp(y)

y$samples

#######################STRAIN EFFECTS#############################
S5vsS1_T30_t18 <- exactTest(y, pair=c(1,5))
S5vsS1_T30_t5 <- exactTest(y, pair=c(2,6))
S5vsS1_T37_t18  <- exactTest(y, pair=c(3,7))
S5vsS1_T37_t5 <- exactTest(y, pair=c(4,8))

deS5vsS1_T30_t18 <- decideTestsDGE(S5vsS1_T30_t18, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS5vsS1_T30_t18)
deS5vsS1_T30_t5 <- decideTestsDGE(S5vsS1_T30_t5, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS5vsS1_T30_t5)
deS5vsS1_T37_t18 <- decideTestsDGE(S5vsS1_T37_t18, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS5vsS1_T37_t18)
deS5vsS1_T37_t5 <- decideTestsDGE(S5vsS1_T37_t5, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS5vsS1_T37_t5)
##################################################################

#####4-way Venn
ytags1 <- rownames(deS5vsS1_T30_t18)[as.logical(deS5vsS1_T30_t18)] 
ytags2 <- rownames(deS5vsS1_T30_t5)[as.logical(deS5vsS1_T30_t5)]
ytags3 <- rownames(deS5vsS1_T37_t18)[as.logical(deS5vsS1_T37_t18)]
ytags4 <- rownames(deS5vsS1_T37_t5)[as.logical(deS5vsS1_T37_t5)]

# Create Venn-Diagram
A <- ytags1
B <- ytags2
C <- ytags3
D <- ytags4

AB <- list(t(A), t(B), t(C), t(D))

library(VennDiagram)
pdf("venn_diagram.pdf")
venn.plot <- venn.diagram(AB, NULL, fill=c("red", "green", "blue", "yellow"), alpha=c(0.5,0.5, 0.5, 0.5), cex = 2, cat.fontface=4, category.names=c("T30_t18", "T30_t5", "T37_t18", "T37_t5"))
grid.draw(venn.plot)
dev.off()


###########################FUNCTIONAL ENRICHMENTS#############################

########################
intersect_list <- intersect(ytags1, ytags2)
intersect_list <- intersect(intersect_list, ytags3)
intersect_list <- intersect(intersect_list, ytags4)
dim(as.matrix(intersect_list))

union_list <- union(ytags1, ytags2)
union_list <- union(union_list, ytags3)
union_list <- union(union_list, ytags4)
union_list <- unique(union_list)
dim(as.matrix(union_list))

####GENE LIST
symbols <- as.matrix(intersect_list)

querylist_entrez <- mapIds(org.EcK12.eg.db, symbols, 'ENTREZID', 'SYMBOL')
querylist_entrez <- querylist_entrez[!is.na(querylist_entrez)]
querylist_entrez <- unique(querylist_entrez)
background <- unique(background)
querylist_entrez <- matrix(querylist_entrez)
dim(querylist_entrez)
write.csv(querylist_entrez, "results/Ecoli.csv")
querylist_kegg <- MS_convertGene(querylist_entrez, "eco", output = "vector", orthology = F)
querylist_kegg <- str_sub(querylist_kegg, 5, 9)
querylist_kegg <- matrix(querylist_kegg)
dim(querylist_kegg)
####################

#Run KEGG enrichment analysis
KEGGEResults <- enrichKEGG(querylist_kegg, universe = background_kegg, organism="eco", pvalueCutoff=0.05)
## Output results
write.csv(KEGGEResults, "results/NEWKEGGResults.csv", quote=F)


## Run GO enrichment analysis
egoResults <- enrichGO(gene = querylist_entrez, keyType = "ENTREZID", universe = background, OrgDb = org.EcK12.eg.db, ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
## Output results from GO analysis to a table
write.csv(egoResults, "results/NEWGOresults.csv")


##########################################################################
####SELECTING DATA TO CONTRAST S6 VS S1

count_data_strain <- cbind(count_data_strain1, count_data_strain6)
dim(count_data_strain)

groups <- colnames(count_data_strain)
group <- str_sub(groups, 1, str_length(groups)-6)
group

y <- DGEList(counts = count_data_strain, group = factor(group))
y <- calcNormFactors(y)
y <- estimateDisp(y)

y$samples

#######################STRAIN EFFECTS#############################
S6vsS1_T30_t18 <- exactTest(y, pair=c(1,5))
S6vsS1_T30_t5 <- exactTest(y, pair=c(2,6))
S6vsS1_T37_t18  <- exactTest(y, pair=c(3,7))
S6vsS1_T37_t5 <- exactTest(y, pair=c(4,8))

deS6vsS1_T30_t18 <- decideTestsDGE(S6vsS1_T30_t18, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS6vsS1_T30_t18)
deS6vsS1_T30_t5 <- decideTestsDGE(S6vsS1_T30_t5, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS6vsS1_T30_t5)
deS6vsS1_T37_t18 <- decideTestsDGE(S6vsS1_T37_t18, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS6vsS1_T37_t18)
deS6vsS1_T37_t5 <- decideTestsDGE(S6vsS1_T37_t5, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS6vsS1_T37_t5)
##################################################################

#####4-way Venn
ytags1 <- rownames(deS6vsS1_T30_t18)[as.logical(deS6vsS1_T30_t18)] 
ytags2 <- rownames(deS6vsS1_T30_t5)[as.logical(deS6vsS1_T30_t5)]
ytags3 <- rownames(deS6vsS1_T37_t18)[as.logical(deS6vsS1_T37_t18)]
ytags4 <- rownames(deS6vsS1_T37_t5)[as.logical(deS6vsS1_T37_t5)]

# Create Venn-Diagram
A <- ytags1
B <- ytags2
C <- ytags3
D <- ytags4

AB <- list(t(A), t(B), t(C), t(D))

library(VennDiagram)
pdf("venn_diagram.pdf")
venn.plot <- venn.diagram(AB, NULL, fill=c("red", "green", "blue", "yellow"), alpha=c(0.5,0.5, 0.5, 0.5), cex = 2, cat.fontface=4, category.names=c("T30_t18", "T30_t5", "T37_t18", "T37_t5"))
grid.draw(venn.plot)
dev.off()


###########################FUNCTIONAL ENRICHMENTS#############################

########################
intersect_list <- intersect(ytags1, ytags2)
intersect_list <- intersect(intersect_list, ytags3)
intersect_list <- intersect(intersect_list, ytags4)
dim(as.matrix(intersect_list))
####GENE LIST
symbols <- as.matrix(intersect_list)

querylist_entrez <- mapIds(org.EcK12.eg.db, symbols, 'ENTREZID', 'SYMBOL')
querylist_entrez <- querylist_entrez[!is.na(querylist_entrez)]
querylist_entrez <- unique(querylist_entrez)
background <- unique(background)
querylist_entrez <- matrix(querylist_entrez)
dim(querylist_entrez)
write.csv(querylist_entrez, "results/Ecoli.csv")
querylist_kegg <- MS_convertGene(querylist_entrez, "eco", output = "vector", orthology = F)
querylist_kegg <- str_sub(querylist_kegg, 5, 9)
querylist_kegg <- matrix(querylist_kegg)
dim(querylist_kegg)
####################

#Run KEGG enrichment analysis
KEGGEResults <- enrichKEGG(querylist_kegg, universe = background_kegg, organism="eco", pvalueCutoff=0.05)
## Output results
write.csv(KEGGEResults, "results/NEWKEGGResults.csv", quote=F)


## Run GO enrichment analysis
egoResults <- enrichGO(gene = querylist_entrez, keyType = "ENTREZID", universe = background, OrgDb = org.EcK12.eg.db, ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
## Output results from GO analysis to a table
write.csv(egoResults, "results/NEWGOresults.csv")





#############################ALL CONTRASTS#################

degS2vsS1_T30_t18 <- rownames(deS2vsS1_T30_t18)[as.logical(deS2vsS1_T30_t18)] 
degS2vsS1_T30_t5 <- rownames(deS2vsS1_T30_t5)[as.logical(deS2vsS1_T30_t5)]
degS2vsS1_T37_t18 <- rownames(deS2vsS1_T37_t18)[as.logical(deS2vsS1_T37_t18)]
degS2vsS1_T37_t5 <- rownames(deS2vsS1_T37_t5)[as.logical(deS2vsS1_T37_t5)]

degS3vsS1_T30_t18 <- rownames(deS3vsS1_T30_t18)[as.logical(deS3vsS1_T30_t18)] 
degS3vsS1_T30_t5 <- rownames(deS3vsS1_T30_t5)[as.logical(deS3vsS1_T30_t5)]
degS3vsS1_T37_t18 <- rownames(deS3vsS1_T37_t18)[as.logical(deS3vsS1_T37_t18)]
degS3vsS1_T37_t5 <- rownames(deS3vsS1_T37_t5)[as.logical(deS3vsS1_T37_t5)]

degS4vsS1_T30_t18 <- rownames(deS4vsS1_T30_t18)[as.logical(deS4vsS1_T30_t18)] 
degS4vsS1_T30_t5 <- rownames(deS4vsS1_T30_t5)[as.logical(deS4vsS1_T30_t5)]
degS4vsS1_T37_t18 <- rownames(deS4vsS1_T37_t18)[as.logical(deS4vsS1_T37_t18)]
degS4vsS1_T37_t5 <- rownames(deS4vsS1_T37_t5)[as.logical(deS4vsS1_T37_t5)]

degS5vsS1_T30_t18 <- rownames(deS5vsS1_T30_t18)[as.logical(deS5vsS1_T30_t18)] 
degS5vsS1_T30_t5 <- rownames(deS5vsS1_T30_t5)[as.logical(deS5vsS1_T30_t5)]
degS5vsS1_T37_t18 <- rownames(deS5vsS1_T37_t18)[as.logical(deS5vsS1_T37_t18)]
degS5vsS1_T37_t5 <- rownames(deS5vsS1_T37_t5)[as.logical(deS5vsS1_T37_t5)]

degS6vsS1_T30_t18 <- rownames(deS6vsS1_T30_t18)[as.logical(deS6vsS1_T30_t18)] 
degS6vsS1_T30_t5 <- rownames(deS6vsS1_T30_t5)[as.logical(deS6vsS1_T30_t5)]
degS6vsS1_T37_t18 <- rownames(deS6vsS1_T37_t18)[as.logical(deS6vsS1_T37_t18)]
degS6vsS1_T37_t5 <- rownames(deS6vsS1_T37_t5)[as.logical(deS6vsS1_T37_t5)]


##################################################T30 t18#####################
#####5-way Venn
ytags1 <- degS2vsS1_T30_t18 
ytags2 <- degS3vsS1_T30_t18
ytags3 <- degS4vsS1_T30_t18
ytags4 <- degS5vsS1_T30_t18
ytags5 <- degS6vsS1_T30_t18

# Create Venn-Diagram
A <- ytags1
B <- ytags2
C <- ytags3
D <- ytags4
E <- ytags5

AB <- list(t(A), t(B), t(C), t(D), t(E))

library(VennDiagram)
pdf("venn_diagram.pdf")
venn.plot <- venn.diagram(AB, NULL, fill=c("red", "green", "blue", "yellow", "purple"), alpha=c(0.5,0.5, 0.5, 0.5, 0.5), cex = 2, cat.fontface=4, category.names=c("S2_T30_t18", "S3_T30_t18", "S4_T30_t18", "S5_T30_t18", "S6_T30_t18"))
grid.draw(venn.plot)
dev.off()


########################
intersect_list <- intersect(ytags1, ytags2)
intersect_list <- intersect(intersect_list, ytags3)
intersect_list <- intersect(intersect_list, ytags4)
intersect_list <- intersect(intersect_list, ytags5)
dim(as.matrix(intersect_list))
intersection_T30_t18 <- intersect_list
####GENE LIST
symbols <- as.matrix(intersect_list)

querylist_entrez <- mapIds(org.EcK12.eg.db, symbols, 'ENTREZID', 'SYMBOL')
querylist_entrez <- querylist_entrez[!is.na(querylist_entrez)]
querylist_entrez <- unique(querylist_entrez)
background <- unique(background)
querylist_entrez <- matrix(querylist_entrez)
dim(querylist_entrez)
write.csv(querylist_entrez, "results/Ecoli.csv")
querylist_kegg <- MS_convertGene(querylist_entrez, "eco", output = "vector", orthology = F)
querylist_kegg <- str_sub(querylist_kegg, 5, 9)
querylist_kegg <- matrix(querylist_kegg)
dim(querylist_kegg)
####################

#Run KEGG enrichment analysis
KEGGEResults <- enrichKEGG(querylist_kegg, universe = background_kegg, organism="eco", pvalueCutoff=0.05)
## Output results
write.csv(KEGGEResults, "results/NEWKEGGResults.csv", quote=F)


## Run GO enrichment analysis
egoResults <- enrichGO(gene = querylist_entrez, keyType = "ENTREZID", universe = background, OrgDb = org.EcK12.eg.db, ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
## Output results from GO analysis to a table
write.csv(egoResults, "results/NEWGOresults.csv")

#####################################T30 t5###########################################
#####5-way Venn
ytags1 <- degS2vsS1_T30_t5
ytags2 <- degS3vsS1_T30_t5
ytags3 <- degS4vsS1_T30_t5
ytags4 <- degS5vsS1_T30_t5
ytags5 <- degS6vsS1_T30_t5

# Create Venn-Diagram
A <- ytags1
B <- ytags2
C <- ytags3
D <- ytags4
E <- ytags5

AB <- list(t(A), t(B), t(C), t(D), t(E))

library(VennDiagram)
pdf("venn_diagram.pdf")
venn.plot <- venn.diagram(AB, NULL, fill=c("red", "green", "blue", "yellow", "purple"), alpha=c(0.5,0.5, 0.5, 0.5, 0.5), cex = 2, cat.fontface=4, category.names=c("S2_T30_t5", "S3_T30_t5", "S4_T30_t5", "S5_T30_t5", "S6_T30_t5"))
grid.draw(venn.plot)
dev.off()


########################
intersect_list <- intersect(ytags1, ytags2)
intersect_list <- intersect(intersect_list, ytags3)
intersect_list <- intersect(intersect_list, ytags4)
intersect_list <- intersect(intersect_list, ytags5)
dim(as.matrix(intersect_list))
intersection_T30_t5 <- intersect_list
####GENE LIST
symbols <- as.matrix(intersect_list)

querylist_entrez <- mapIds(org.EcK12.eg.db, symbols, 'ENTREZID', 'SYMBOL')
querylist_entrez <- querylist_entrez[!is.na(querylist_entrez)]
querylist_entrez <- unique(querylist_entrez)
background <- unique(background)
querylist_entrez <- matrix(querylist_entrez)
dim(querylist_entrez)
write.csv(querylist_entrez, "results/Ecoli.csv")
querylist_kegg <- MS_convertGene(querylist_entrez, "eco", output = "vector", orthology = F)
querylist_kegg <- str_sub(querylist_kegg, 5, 9)
querylist_kegg <- matrix(querylist_kegg)
dim(querylist_kegg)
####################

#Run KEGG enrichment analysis
KEGGEResults <- enrichKEGG(querylist_kegg, universe = background_kegg, organism="eco", pvalueCutoff=0.05)
## Output results
write.csv(KEGGEResults, "results/NEWKEGGResults.csv", quote=F)


## Run GO enrichment analysis
egoResults <- enrichGO(gene = querylist_entrez, keyType = "ENTREZID", universe = background, OrgDb = org.EcK12.eg.db, ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
## Output results from GO analysis to a table
write.csv(egoResults, "results/NEWGOresults.csv")


#####################################T37 t18###########################################
#####5-way Venn
ytags1 <- degS2vsS1_T37_t18
ytags2 <- degS3vsS1_T37_t18
ytags3 <- degS4vsS1_T37_t18
ytags4 <- degS5vsS1_T37_t18
ytags5 <- degS6vsS1_T37_t18

# Create Venn-Diagram
A <- ytags1
B <- ytags2
C <- ytags3
D <- ytags4
E <- ytags5

AB <- list(t(A), t(B), t(C), t(D), t(E))

library(VennDiagram)
pdf("venn_diagram.pdf")
venn.plot <- venn.diagram(AB, NULL, fill=c("red", "green", "blue", "yellow", "purple"), alpha=c(0.5,0.5, 0.5, 0.5, 0.5), cex = 2, cat.fontface=4, category.names=c("S2_T37_t18", "S3_T37_t18", "S4_T37_t18", "S5_T37_t18", "S6_T37_t18"))
grid.draw(venn.plot)
dev.off()


########################
intersect_list <- intersect(ytags1, ytags2)
intersect_list <- intersect(intersect_list, ytags3)
intersect_list <- intersect(intersect_list, ytags4)
intersect_list <- intersect(intersect_list, ytags5)
dim(as.matrix(intersect_list))
intersection_T37_t18 <- intersect_list
####GENE LIST
symbols <- as.matrix(intersect_list)

querylist_entrez <- mapIds(org.EcK12.eg.db, symbols, 'ENTREZID', 'SYMBOL')
querylist_entrez <- querylist_entrez[!is.na(querylist_entrez)]
querylist_entrez <- unique(querylist_entrez)
background <- unique(background)
querylist_entrez <- matrix(querylist_entrez)
dim(querylist_entrez)
write.csv(querylist_entrez, "results/Ecoli.csv")
querylist_kegg <- MS_convertGene(querylist_entrez, "eco", output = "vector", orthology = F)
querylist_kegg <- str_sub(querylist_kegg, 5, 9)
querylist_kegg <- matrix(querylist_kegg)
dim(querylist_kegg)
####################

#Run KEGG enrichment analysis
KEGGEResults <- enrichKEGG(querylist_kegg, universe = background_kegg, organism="eco", pvalueCutoff=0.05)
## Output results
write.csv(KEGGEResults, "results/NEWKEGGResults.csv", quote=F)


## Run GO enrichment analysis
egoResults <- enrichGO(gene = querylist_entrez, keyType = "ENTREZID", universe = background, OrgDb = org.EcK12.eg.db, ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
## Output results from GO analysis to a table
write.csv(egoResults, "results/NEWGOresults.csv")


#####################################T37 t5###########################################
#####5-way Venn
ytags1 <- degS2vsS1_T37_t5
ytags2 <- degS3vsS1_T37_t5
ytags3 <- degS4vsS1_T37_t5
ytags4 <- degS5vsS1_T37_t5
ytags5 <- degS6vsS1_T37_t5

# Create Venn-Diagram
A <- ytags1
B <- ytags2
C <- ytags3
D <- ytags4
E <- ytags5

ABCDE <- list(t(A), t(B), t(C), t(D), t(E))

library(VennDiagram)
pdf("venn_diagram.pdf")
venn.plot <- venn.diagram(ABCDE, NULL, fill=c("red", "green", "blue", "yellow", "purple"), alpha=c(0.5,0.5, 0.5, 0.5, 0.5), cex = 2, cat.fontface=4, category.names=c("S2_T37_t5", "S3_T37_t5", "S4_T37_t5", "S5_T37_t5", "S6_T37_t5"))
grid.draw(venn.plot)
dev.off()


########################
intersect_list <- intersect(ytags1, ytags2)
intersect_list <- intersect(intersect_list, ytags3)
intersect_list <- intersect(intersect_list, ytags4)
intersect_list <- intersect(intersect_list, ytags5)
dim(as.matrix(intersect_list))
intersection_T37_t5 <- intersect_list
####GENE LIST
symbols <- as.matrix(intersect_list)

querylist_entrez <- mapIds(org.EcK12.eg.db, symbols, 'ENTREZID', 'SYMBOL')
querylist_entrez <- querylist_entrez[!is.na(querylist_entrez)]
querylist_entrez <- unique(querylist_entrez)
background <- unique(background)
querylist_entrez <- matrix(querylist_entrez)
dim(querylist_entrez)
write.csv(querylist_entrez, "results/Ecoli.csv")
querylist_kegg <- MS_convertGene(querylist_entrez, "eco", output = "vector", orthology = F)
querylist_kegg <- str_sub(querylist_kegg, 5, 9)
querylist_kegg <- matrix(querylist_kegg)
dim(querylist_kegg)
####################

#Run KEGG enrichment analysis
KEGGEResults <- enrichKEGG(querylist_kegg, universe = background_kegg, organism="eco", pvalueCutoff=0.05)
## Output results
write.csv(KEGGEResults, "results/NEWKEGGResults.csv", quote=F)


## Run GO enrichment analysis
egoResults <- enrichGO(gene = querylist_entrez, keyType = "ENTREZID", universe = background, OrgDb = org.EcK12.eg.db, ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
## Output results from GO analysis to a table
write.csv(egoResults, "results/NEWGOresults.csv")


#####INTERSECTION OF INTERSECTIONS#####################################################################

#####5-way Venn
ytags1 <- intersection_T30_t18
ytags2 <- intersection_T30_t5
ytags3 <- intersection_T37_t18
ytags4 <- intersection_T37_t5

# Create Venn-Diagram
A <- ytags1
B <- ytags2
C <- ytags3
D <- ytags4

ABCD <- list(t(A), t(B), t(C), t(D))

library(VennDiagram)
pdf("venn_diagram.pdf")
venn.plot <- venn.diagram(ABCD, NULL, fill=c("red", "green", "blue", "yellow"), alpha=c(0.5,0.5, 0.5, 0.5), cex = 2, cat.fontface=4, category.names=c("T37_t5", "T37_t5", "T37_t5", "T37_t5"))
grid.draw(venn.plot)
dev.off()


########################
intersect_list <- intersect(ytags1, ytags2)
intersect_list <- intersect(intersect_list, ytags3)
intersect_list <- intersect(intersect_list, ytags4)
intersect_list <- intersect(intersect_list, ytags5)
dim(as.matrix(intersect_list))
intersection_T37_t5 <- intersect_list
####GENE LIST
symbols <- as.matrix(intersect_list)

querylist_entrez <- mapIds(org.EcK12.eg.db, symbols, 'ENTREZID', 'SYMBOL')
querylist_entrez <- querylist_entrez[!is.na(querylist_entrez)]
querylist_entrez <- unique(querylist_entrez)
background <- unique(background)
querylist_entrez <- matrix(querylist_entrez)
dim(querylist_entrez)
write.csv(querylist_entrez, "results/Ecoli.csv")
querylist_kegg <- MS_convertGene(querylist_entrez, "eco", output = "vector", orthology = F)
querylist_kegg <- str_sub(querylist_kegg, 5, 9)
querylist_kegg <- matrix(querylist_kegg)
dim(querylist_kegg)
####################

#Run KEGG enrichment analysis
KEGGEResults <- enrichKEGG(querylist_kegg, universe = background_kegg, organism="eco", pvalueCutoff=0.05)
## Output results
write.csv(KEGGEResults, "results/NEWKEGGResults.csv", quote=F)


## Run GO enrichment analysis
egoResults <- enrichGO(gene = querylist_entrez, keyType = "ENTREZID", universe = background, OrgDb = org.EcK12.eg.db, ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
## Output results from GO analysis to a table
write.csv(egoResults, "results/NEWGOresults.csv")