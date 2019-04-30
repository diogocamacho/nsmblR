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


####SELECTING DATA TO CONTRAST S1_I38 VS S1_I0

groups <- colnames(count_data_strain1)
group <- str_sub(groups, 1, str_length(groups)-6)
group

y <- DGEList(counts = count_data_strain1, group = factor(group))
y <- calcNormFactors(y)
y <- estimateDisp(y)

y$samples

#######################STRAIN EFFECTS#############################
S1_T30_t18 <- exactTest(y, pair=c(1,2))
S1_T30_t5 <- exactTest(y, pair=c(3,4))
S1_T37_t18  <- exactTest(y, pair=c(5,6))
S1_T37_t5 <- exactTest(y, pair=c(7,8))

deS1_T30_t18 <- decideTestsDGE(S1_T30_t18, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS1_T30_t18)
deS1_T30_t5 <- decideTestsDGE(S1_T30_t5, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS1_T30_t5)
deS1_T37_t18 <- decideTestsDGE(S1_T37_t18, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS1_T37_t18)
deS1_T37_t5 <- decideTestsDGE(S1_T37_t5, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS1_T37_t5)
##################################################################

#####4-way Venn
ytags1 <- rownames(deS1_T30_t18)[as.logical(deS1_T30_t18)] 
ytags2 <- rownames(deS1_T30_t5)[as.logical(deS1_T30_t5)]
ytags3 <- rownames(deS1_T37_t18)[as.logical(deS1_T37_t18)]
ytags4 <- rownames(deS1_T37_t5)[as.logical(deS1_T37_t5)]

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
####SELECTING DATA TO CONTRAST S2_I38 VS S2_I0

groups <- colnames(count_data_strain2)
group <- str_sub(groups, 1, str_length(groups)-6)
group

y <- DGEList(counts = count_data_strain2, group = factor(group))
y <- calcNormFactors(y)
y <- estimateDisp(y)

y$samples

#######################STRAIN EFFECTS#############################
S2_T30_t18 <- exactTest(y, pair=c(1,2))
S2_T30_t5 <- exactTest(y, pair=c(3,4))
S2_T37_t18  <- exactTest(y, pair=c(5,6))
S2_T37_t5 <- exactTest(y, pair=c(7,8))

deS2_T30_t18 <- decideTestsDGE(S2_T30_t18, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS2_T30_t18)
deS2_T30_t5 <- decideTestsDGE(S2_T30_t5, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS2_T30_t5)
deS2_T37_t18 <- decideTestsDGE(S2_T37_t18, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS2_T37_t18)
deS2_T37_t5 <- decideTestsDGE(S2_T37_t5, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS2_T37_t5)
##################################################################

#####4-way Venn
ytags1 <- rownames(deS2_T30_t18)[as.logical(deS2_T30_t18)] 
ytags2 <- rownames(deS2_T30_t5)[as.logical(deS2_T30_t5)]
ytags3 <- rownames(deS2_T37_t18)[as.logical(deS2_T37_t18)]
ytags4 <- rownames(deS2_T37_t5)[as.logical(deS2_T37_t5)]

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
####SELECTING DATA TO CONTRAST S3_I38 VS S3_I0

groups <- colnames(count_data_strain3)
group <- str_sub(groups, 1, str_length(groups)-6)
group

y <- DGEList(counts = count_data_strain3, group = factor(group))
y <- calcNormFactors(y)
y <- estimateDisp(y)

y$samples

#######################STRAIN EFFECTS#############################
S3_T30_t18 <- exactTest(y, pair=c(1,2))
S3_T30_t5 <- exactTest(y, pair=c(3,4))
S3_T37_t18  <- exactTest(y, pair=c(5,6))
S3_T37_t5 <- exactTest(y, pair=c(7,8))

deS3_T30_t18 <- decideTestsDGE(S3_T30_t18, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS3_T30_t18)
deS3_T30_t5 <- decideTestsDGE(S3_T30_t5, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS3_T30_t5)
deS3_T37_t18 <- decideTestsDGE(S3_T37_t18, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS3_T37_t18)
deS3_T37_t5 <- decideTestsDGE(S3_T37_t5, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS3_T37_t5)
##################################################################

#####4-way Venn
ytags1 <- rownames(deS3_T30_t18)[as.logical(deS3_T30_t18)] 
ytags2 <- rownames(deS3_T30_t5)[as.logical(deS3_T30_t5)]
ytags3 <- rownames(deS3_T37_t18)[as.logical(deS3_T37_t18)]
ytags4 <- rownames(deS3_T37_t5)[as.logical(deS3_T37_t5)]

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
####SELECTING DATA TO CONTRAST S4_I38 VS S4_I0

groups <- colnames(count_data_strain4)
group <- str_sub(groups, 1, str_length(groups)-6)
group

y <- DGEList(counts = count_data_strain4, group = factor(group))
y <- calcNormFactors(y)
y <- estimateDisp(y)

y$samples

#######################STRAIN EFFECTS#############################
S4_T30_t18 <- exactTest(y, pair=c(1,2))
S4_T30_t5 <- exactTest(y, pair=c(3,4))
S4_T37_t18  <- exactTest(y, pair=c(5,6))
S4_T37_t5 <- exactTest(y, pair=c(7,8))

deS4_T30_t18 <- decideTestsDGE(S4_T30_t18, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS4_T30_t18)
deS4_T30_t5 <- decideTestsDGE(S4_T30_t5, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS4_T30_t5)
deS4_T37_t18 <- decideTestsDGE(S4_T37_t18, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS4_T37_t18)
deS4_T37_t5 <- decideTestsDGE(S4_T37_t5, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS4_T37_t5)
##################################################################

#####4-way Venn
ytags1 <- rownames(deS4_T30_t18)[as.logical(deS4_T30_t18)] 
ytags2 <- rownames(deS4_T30_t5)[as.logical(deS4_T30_t5)]
ytags3 <- rownames(deS4_T37_t18)[as.logical(deS4_T37_t18)]
ytags4 <- rownames(deS4_T37_t5)[as.logical(deS4_T37_t5)]

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
####SELECTING DATA TO CONTRAST S6_I38 VS S6_I0

groups <- colnames(count_data_strain6)
group <- str_sub(groups, 1, str_length(groups)-6)
group

y <- DGEList(counts = count_data_strain6, group = factor(group))
y <- calcNormFactors(y)
y <- estimateDisp(y)

y$samples

#######################STRAIN EFFECTS#############################
S6_T30_t18 <- exactTest(y, pair=c(1,2))
S6_T30_t5 <- exactTest(y, pair=c(3,4))
S6_T37_t18  <- exactTest(y, pair=c(5,6))
S6_T37_t5 <- exactTest(y, pair=c(7,8))

deS6_T30_t18 <- decideTestsDGE(S6_T30_t18, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS6_T30_t18)
deS6_T30_t5 <- decideTestsDGE(S6_T30_t5, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS6_T30_t5)
deS6_T37_t18 <- decideTestsDGE(S6_T37_t18, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS6_T37_t18)
deS6_T37_t5 <- decideTestsDGE(S6_T37_t5, adjust.method="BH", p.value=0.05, lfc=0.5849625)
summary(deS6_T37_t5)
##################################################################

#####4-way Venn
ytags1 <- rownames(deS6_T30_t18)[as.logical(deS6_T30_t18)] 
ytags2 <- rownames(deS6_T30_t5)[as.logical(deS6_T30_t5)]
ytags3 <- rownames(deS6_T37_t18)[as.logical(deS6_T37_t18)]
ytags4 <- rownames(deS6_T37_t5)[as.logical(deS6_T37_t5)]

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

degS1_T30_t18 <- rownames(deS1_T30_t18)[as.logical(deS1_T30_t18)] 
degS1_T30_t5 <- rownames(deS1_T30_t5)[as.logical(deS1_T30_t5)]
degS1_T37_t18 <- rownames(deS1_T37_t18)[as.logical(deS1_T37_t18)]
degS1_T37_t5 <- rownames(deS1_T37_t5)[as.logical(deS1_T37_t5)]

degS2_T30_t18 <- rownames(deS2_T30_t18)[as.logical(deS2_T30_t18)] 
degS2_T30_t5 <- rownames(deS2_T30_t5)[as.logical(deS2_T30_t5)]
degS2_T37_t18 <- rownames(deS2_T37_t18)[as.logical(deS2_T37_t18)]
degS2_T37_t5 <- rownames(deS2_T37_t5)[as.logical(deS2_T37_t5)]

degS3_T30_t18 <- rownames(deS3_T30_t18)[as.logical(deS3_T30_t18)] 
degS3_T30_t5 <- rownames(deS3_T30_t5)[as.logical(deS3_T30_t5)]
degS3_T37_t18 <- rownames(deS3_T37_t18)[as.logical(deS3_T37_t18)]
degS3_T37_t5 <- rownames(deS3_T37_t5)[as.logical(deS3_T37_t5)]

degS4_T30_t18 <- rownames(deS4_T30_t18)[as.logical(deS4_T30_t18)] 
degS4_T30_t5 <- rownames(deS4_T30_t5)[as.logical(deS4_T30_t5)]
degS4_T37_t18 <- rownames(deS4_T37_t18)[as.logical(deS4_T37_t18)]
degS4_T37_t5 <- rownames(deS4_T37_t5)[as.logical(deS4_T37_t5)]

degS6_T30_t18 <- rownames(deS6_T30_t18)[as.logical(deS6_T30_t18)]
degS6_T30_t5 <- rownames(deS6_T30_t5)[as.logical(deS6_T30_t5)]
degS6_T37_t18 <- rownames(deS6_T37_t18)[as.logical(deS6_T37_t18)]
degS6_T37_t5 <- rownames(deS6_T37_t5)[as.logical(deS6_T37_t5)]


##################################################T30 t18#####################
#####5-way Venn
ytags1 <- degS1_T30_t18 
ytags2 <- degS2_T30_t18
ytags3 <- degS3_T30_t18
ytags4 <- degS4_T30_t18
ytags5 <- degS6_T30_t18

# Create Venn-Diagram
A <- ytags1
B <- ytags2
C <- ytags3
D <- ytags4
E <- ytags5

AB <- list(t(A), t(B), t(C), t(D), t(E))

library(VennDiagram)
pdf("venn_diagram.pdf")
venn.plot <- venn.diagram(AB, NULL, fill=c("red", "green", "blue", "yellow", "purple"), alpha=c(0.5,0.5, 0.5, 0.5, 0.5), cex = 2, cat.fontface=4, category.names=c("S1_T30_t18", "S2_T30_t18", "S3_T30_t18", "S4_T30_t18", "S6_T30_t18"))
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
ytags1 <- degS1_T30_t5
ytags2 <- degS2_T30_t5
ytags3 <- degS3_T30_t5
ytags4 <- degS4_T30_t5
ytags5 <- degS6_T30_t5

# Create Venn-Diagram
A <- ytags1
B <- ytags2
C <- ytags3
D <- ytags4
E <- ytags5

AB <- list(t(A), t(B), t(C), t(D), t(E))

library(VennDiagram)
pdf("venn_diagram.pdf")
venn.plot <- venn.diagram(AB, NULL, fill=c("red", "green", "blue", "yellow", "purple"), alpha=c(0.5,0.5, 0.5, 0.5, 0.5), cex = 2, cat.fontface=4, category.names=c("S1_T30_t5", "S2_T30_t5", "S3_T30_t5", "S4_T30_t5", "S6_T30_t5"))
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
ytags1 <- degS1_T37_t18
ytags2 <- degS2_T37_t18
ytags3 <- degS3_T37_t18
ytags4 <- degS4_T37_t18
ytags5 <- degS6_T37_t18

# Create Venn-Diagram
A <- ytags1
B <- ytags2
C <- ytags3
D <- ytags4
E <- ytags5

AB <- list(t(A), t(B), t(C), t(D), t(E))

library(VennDiagram)
pdf("venn_diagram.pdf")
venn.plot <- venn.diagram(AB, NULL, fill=c("red", "green", "blue", "yellow", "purple"), alpha=c(0.5,0.5, 0.5, 0.5, 0.5), cex = 2, cat.fontface=4, category.names=c("S1_T37_t18", "S2_T37_t18", "S3_T37_t18", "S4_T37_t18", "S6_T37_t18"))
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
ytags1 <- degS1_T37_t5
ytags2 <- degS2_T37_t5
ytags3 <- degS3_T37_t5
ytags4 <- degS4_T37_t5
ytags5 <- degS6_T37_t5

# Create Venn-Diagram
A <- ytags1
B <- ytags2
C <- ytags3
D <- ytags4
E <- ytags5

ABCDE <- list(t(A), t(B), t(C), t(D), t(E))

library(VennDiagram)
pdf("venn_diagram.pdf")
venn.plot <- venn.diagram(ABCDE, NULL, fill=c("red", "green", "blue", "yellow", "purple"), alpha=c(0.5,0.5, 0.5, 0.5, 0.5), cex = 2, cat.fontface=4, category.names=c("S1_T37_t5", "S2_T37_t5", "S3_T37_t5", "S4_T37_t5", "S4_T37_t5"))
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

#####4-way Venn
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
venn.plot <- venn.diagram(ABCD, NULL, fill=c("red", "green", "blue", "yellow"), alpha=c(0.5,0.5, 0.5, 0.5), cex = 2, cat.fontface=4, category.names=c("T30_t18", "T30_t5", "T37_t18", "T37_t5"))
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


####UNIONS
degS1_ALL <- union(degS1_T30_t18, degS1_T30_t5)
degS1_ALL <- union(degS1_ALL, degS1_T37_t18)
degS1_ALL <- union(degS1_ALL, degS1_T37_t5)
dim(as.matrix(degS1_ALL))

degS2_ALL <- union(degS2_T30_t18, degS2_T30_t5)
degS2_ALL <- union(degS2_ALL, degS2_T37_t18)
degS2_ALL <- union(degS2_ALL, degS2_T37_t5)
dim(as.matrix(degS2_ALL))

degS3_ALL <- union(degS3_T30_t18, degS3_T30_t5)
degS3_ALL <- union(degS3_ALL, degS3_T37_t18)
degS3_ALL <- union(degS3_ALL, degS3_T37_t5)
dim(as.matrix(degS3_ALL))

degS4_ALL <- union(degS4_T30_t18, degS4_T30_t5)
degS4_ALL <- union(degS4_ALL, degS4_T37_t18)
degS4_ALL <- union(degS4_ALL, degS4_T37_t5)
dim(as.matrix(degS4_ALL))

degS6_ALL <- union(degS6_T30_t18, degS6_T30_t5)
degS6_ALL <- union(degS6_ALL, degS6_T37_t18)
degS6_ALL <- union(degS6_ALL, degS6_T37_t5)
dim(as.matrix(degS6_ALL))

#####5-way Venn
ytags1 <- degS1_ALL
ytags2 <- degS2_ALL
ytags3 <- degS3_ALL
ytags4 <- degS4_ALL
ytags5 <- degS6_ALL

# Create Venn-Diagram
A <- ytags1
B <- ytags2
C <- ytags3
D <- ytags4
E <- ytags5

ABCDE <- list(t(A), t(B), t(C), t(D), t(E))

library(VennDiagram)
pdf("venn_diagram.pdf")
venn.plot <- venn.diagram(ABCDE, NULL, fill=c("red", "green", "blue", "yellow", "purple"), alpha=c(0.5,0.5, 0.5, 0.5, 0.5), cex = 2, cat.fontface=4, category.names=c("S1", "S2", "S3", "S4", "S6"))
grid.draw(venn.plot)
dev.off()

####UNIONS MINUS S1 CHANGES
degS2_ALL_noS1 <- setdiff(degS2_ALL, degS1_ALL)
degS3_ALL_noS1 <- setdiff(degS3_ALL, degS1_ALL)
degS4_ALL_noS1 <- setdiff(degS4_ALL, degS1_ALL)
degS6_ALL_noS1 <- setdiff(degS6_ALL, degS1_ALL)

#####4-way Venn
ytags1 <- degS2_ALL_noS1
ytags2 <- degS3_ALL_noS1
ytags3 <- degS4_ALL_noS1
ytags4 <- degS6_ALL_noS1

# Create Venn-Diagram
A <- ytags1
B <- ytags2
C <- ytags3
D <- ytags4

ABCD <- list(t(A), t(B), t(C), t(D))

library(VennDiagram)
pdf("venn_diagram.pdf")
venn.plot <- venn.diagram(ABCD, NULL, fill=c("red", "green", "blue", "yellow"), alpha=c(0.5,0.5, 0.5, 0.5), cex = 2, cat.fontface=4, category.names=c("S2", "S3", "S4", "S6"))
grid.draw(venn.plot)
dev.off()

####GENE LIST
symbols <- as.matrix(ytags2)

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


