# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("minet", version = "3.8")
#install.packages("Matrix")

library(minet)
library(Matrix)
# 
# library(help = minet)
# 
# vignette("minet")
# 
# data(syn.net)
# data(syn.data)
# mim <- build.mim(syn.data,estimator="spearman")
# table <- validate(mim,syn.net)
# auc.roc(table)
# show.roc(table, col="red",type="b")
# show.pr(table, col="red",type="b")
# 
# aracne.net <- aracne(mim, eps=0)
# table <- validate(aracne.net,syn.net)
# show.roc(table, col="red",type="b")
# show.pr(table, col="red",type="b")



#####E. coli data analysis###########
#######SPEARMAN######################
####################################
start_time <- Sys.time()
mim1 <- build.mim(t(count_data_t5), estimator = "spearman", disc = "none", nbins = sqrt(NROW(t(count_data_t5))))
mim2 <- build.mim(t(count_data_t18), estimator = "spearman", disc = "none", nbins = sqrt(NROW(t(count_data_t18))))
mim3 <- build.mim(t(cpm_data_t5), estimator = "spearman", disc = "none", nbins = sqrt(NROW(t(cpm_data_t5))))
mim4 <- build.mim(t(cpm_data_t18), estimator = "spearman", disc = "none", nbins = sqrt(NROW(t(cpm_data_t18))))
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
clr.net1 <- clr(mim1)
clr.net2 <- clr(mim2)
clr.net3 <- clr(mim3)
clr.net4 <- clr(mim4)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
aracne.net1 <- aracne(mim1, eps=0)
aracne.net2 <- aracne(mim2, eps=0)
aracne.net3 <- aracne(mim3, eps=0)
aracne.net4 <- aracne(mim4, eps=0)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
mrnet.net1 <- mrnet(mim1)
mrnet.net2 <- mrnet(mim2)
mrnet.net3 <- mrnet(mim3)
mrnet.net4 <- mrnet(mim4)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
mrnetb.net1 <- mrnetb(mim1)
mrnetb.net2 <- mrnetb(mim2)
mrnetb.net3 <- mrnetb(mim3)
mrnetb.net4 <- mrnetb(mim4)
end_time <- Sys.time()
end_time - start_time

######CREATE ENSEMBLE MIM#######################
mask <- mim1[upper.tri(mim1)]
#length(mask)
#nrow(mim1)*(nrow(mim1)-1)/2
ranked_mask <- rank(mask, na.last = TRUE)
mim1[upper.tri(mim1, diag=FALSE)] <- ranked_mask

mask <- mim2[upper.tri(mim2)]
ranked_mask <- rank(mask, na.last = TRUE)
mim2[upper.tri(mim2, diag=FALSE)] <- ranked_mask

mask <- mim3[upper.tri(mim3)]
ranked_mask <- rank(mask, na.last = TRUE)
mim3[upper.tri(mim3, diag=FALSE)] <- ranked_mask

mask <- mim4[upper.tri(mim4)]
ranked_mask <- rank(mask, na.last = TRUE)
mim4[upper.tri(mim4, diag=FALSE)] <- ranked_mask

ensemble_mim <- (mim1 + mim2 + mim3 + mim4)/4
rm(mim1, mim2, mim3, mim4)
##################################################

######CREATE ENSEMBLE CLR#########################
mask <- clr.net1[upper.tri(clr.net1)]
ranked_mask <- rank(mask, na.last = TRUE)
clr.net1[upper.tri(clr.net1, diag=FALSE)] <- ranked_mask

mask <- clr.net2[upper.tri(clr.net2)]
ranked_mask <- rank(mask, na.last = TRUE)
clr.net2[upper.tri(clr.net2, diag=FALSE)] <- ranked_mask

mask <- clr.net3[upper.tri(clr.net3)]
ranked_mask <- rank(mask, na.last = TRUE)
clr.net3[upper.tri(clr.net3, diag=FALSE)] <- ranked_mask

mask <- clr.net4[upper.tri(clr.net4)]
ranked_mask <- rank(mask, na.last = TRUE)
clr.net4[upper.tri(clr.net4, diag=FALSE)] <- ranked_mask

ensemble_clr.net <- (clr.net1 + clr.net2 + clr.net3 + clr.net4)/4
rm(clr.net1, clr.net2, clr.net3, clr.net4)
##################################################

######CREATE ENSEMBLE ARACNE#########################
mask <- aracne.net1[upper.tri(aracne.net1)]
ranked_mask <- rank(mask, na.last = TRUE)
aracne.net1[upper.tri(aracne.net1, diag=FALSE)] <- ranked_mask

mask <- aracne.net2[upper.tri(aracne.net2)]
ranked_mask <- rank(mask, na.last = TRUE)
aracne.net2[upper.tri(aracne.net2, diag=FALSE)] <- ranked_mask

mask <- aracne.net3[upper.tri(aracne.net3)]
ranked_mask <- rank(mask, na.last = TRUE)
aracne.net3[upper.tri(aracne.net3, diag=FALSE)] <- ranked_mask

mask <- aracne.net4[upper.tri(aracne.net4)]
ranked_mask <- rank(mask, na.last = TRUE)
aracne.net4[upper.tri(aracne.net4, diag=FALSE)] <- ranked_mask

ensemble_aracne.net <- (aracne.net1 + aracne.net2 + aracne.net3 + aracne.net4)/4
rm(aracne.net1, aracne.net2, aracne.net3, aracne.net4)
##################################################

######CREATE ENSEMBLE MRNET#########################
mask <- mrnet.net1[upper.tri(mrnet.net1)]
ranked_mask <- rank(mask, na.last = TRUE)
mrnet.net1[upper.tri(mrnet.net1, diag=FALSE)] <- ranked_mask

mask <- mrnet.net2[upper.tri(mrnet.net2)]
ranked_mask <- rank(mask, na.last = TRUE)
mrnet.net2[upper.tri(mrnet.net2, diag=FALSE)] <- ranked_mask

mask <- mrnet.net3[upper.tri(mrnet.net3)]
ranked_mask <- rank(mask, na.last = TRUE)
mrnet.net3[upper.tri(mrnet.net3, diag=FALSE)] <- ranked_mask

mask <- mrnet.net4[upper.tri(mrnet.net4)]
ranked_mask <- rank(mask, na.last = TRUE)
mrnet.net4[upper.tri(mrnet.net4, diag=FALSE)] <- ranked_mask

ensemble_mrnet.net <- (mrnet.net1 + mrnet.net2 + mrnet.net3 + mrnet.net4)/4
rm(mrnet.net1, mrnet.net2, mrnet.net3, mrnet.net4)
##################################################

######CREATE ENSEMBLE MRNETB#########################
mask <- mrnetb.net1[upper.tri(mrnetb.net1)]
ranked_mask <- rank(mask, na.last = TRUE)
mrnetb.net1[upper.tri(mrnetb.net1, diag=FALSE)] <- ranked_mask

mask <- mrnetb.net2[upper.tri(mrnetb.net2)]
ranked_mask <- rank(mask, na.last = TRUE)
mrnetb.net2[upper.tri(mrnetb.net2, diag=FALSE)] <- ranked_mask

mask <- mrnetb.net3[upper.tri(mrnetb.net3)]
ranked_mask <- rank(mask, na.last = TRUE)
mrnetb.net3[upper.tri(mrnetb.net3, diag=FALSE)] <- ranked_mask

mask <- mrnetb.net4[upper.tri(mrnetb.net4)]
ranked_mask <- rank(mask, na.last = TRUE)
mrnetb.net4[upper.tri(mrnetb.net4, diag=FALSE)] <- ranked_mask

ensemble_mrnetb.net <- (mrnetb.net1 + mrnetb.net2 + mrnetb.net3 + mrnetb.net4)/4
rm(mrnetb.net1, mrnetb.net2, mrnetb.net3, mrnetb.net4)
##################################################
rm(mask, ranked_mask)
#################################################
#############FINAL###############################
#################################################
super_ensemble <- (ensemble_mim + ensemble_clr.net + ensemble_aracne.net + ensemble_mrnet.net + ensemble_mrnetb.net)/5
#################################################
#################################################
#################################################
max(ensemble_mim[upper.tri(ensemble_mim)])
min(ensemble_mim[upper.tri(ensemble_mim)])
max(ensemble_mim[lower.tri(ensemble_mim)])
min(ensemble_mim[lower.tri(ensemble_mim)])

max(ensemble_clr.net[upper.tri(ensemble_clr.net)])
min(ensemble_clr.net[upper.tri(ensemble_clr.net)])
max(ensemble_clr.net[lower.tri(ensemble_clr.net)])
min(ensemble_clr.net[lower.tri(ensemble_clr.net)])

max(ensemble_aracne.net[upper.tri(ensemble_aracne.net)])
min(ensemble_aracne.net[upper.tri(ensemble_aracne.net)])
max(ensemble_aracne.net[lower.tri(ensemble_aracne.net)])
min(ensemble_aracne.net[lower.tri(ensemble_aracne.net)])

max(ensemble_mrnet.net[upper.tri(ensemble_mrnet.net)])
min(ensemble_mrnet.net[upper.tri(ensemble_mrnet.net)])
max(ensemble_mrnet.net[lower.tri(ensemble_mrnet.net)])
min(ensemble_mrnet.net[lower.tri(ensemble_mrnet.net)])

max(ensemble_mrnetb.net[upper.tri(ensemble_mrnetb.net)])
min(ensemble_mrnetb.net[upper.tri(ensemble_mrnetb.net)])
max(ensemble_mrnetb.net[lower.tri(ensemble_mrnetb.net)])
min(ensemble_mrnetb.net[lower.tri(ensemble_mrnetb.net)])

max(super_ensemble[upper.tri(super_ensemble)])
min(super_ensemble[upper.tri(super_ensemble)])
max(super_ensemble[lower.tri(super_ensemble)])
min(super_ensemble[lower.tri(super_ensemble)])

cor(ensemble_mim[upper.tri(ensemble_mim)], ensemble_clr.net[upper.tri(ensemble_clr.net)], method="spearman")
cor(ensemble_mim[upper.tri(ensemble_mim)], ensemble_aracne.net[upper.tri(ensemble_aracne.net)], method="spearman")
cor(ensemble_mim[upper.tri(ensemble_mim)], ensemble_mrnet.net[upper.tri(ensemble_mrnet.net)], method="spearman")
cor(ensemble_mim[upper.tri(ensemble_mim)], ensemble_mrnetb.net[upper.tri(ensemble_mrnetb.net)], method="spearman")

cor(ensemble_mim[lower.tri(ensemble_mim)], ensemble_clr.net[lower.tri(ensemble_clr.net)], method="spearman")
cor(ensemble_mim[lower.tri(ensemble_mim)], ensemble_aracne.net[lower.tri(ensemble_aracne.net)], method="spearman")
cor(ensemble_mim[lower.tri(ensemble_mim)], ensemble_mrnet.net[lower.tri(ensemble_mrnet.net)], method="spearman")
cor(ensemble_mim[lower.tri(ensemble_mim)], ensemble_mrnetb.net[lower.tri(ensemble_mrnetb.net)], method="spearman")

cor(ensemble_mim[upper.tri(ensemble_mim)], t(ensemble_mim)[upper.tri(t(ensemble_mim))], method="spearman")
cor(ensemble_clr.net[upper.tri(ensemble_clr.net)], t(ensemble_clr.net)[upper.tri(t(ensemble_clr.net))], method="spearman")
cor(ensemble_aracne.net[upper.tri(ensemble_aracne.net)], t(ensemble_aracne.net)[upper.tri(t(ensemble_aracne.net))], method="spearman")
cor(ensemble_mrnet.net[upper.tri(ensemble_mrnet.net)], t(ensemble_mrnet.net)[upper.tri(t(ensemble_mrnet.net))], method="spearman")
cor(ensemble_mrnetb.net[upper.tri(ensemble_mrnetb.net)], t(ensemble_mrnetb.net)[upper.tri(t(ensemble_mrnetb.net))], method="spearman")


######CREATE EXAMPLE GOLD STANDARD BASED ON MIM ENSEMBLE PREDICTION

edge_num <- 10000
thresh <- nrow(ensemble_mim)*(nrow(ensemble_mim)-1)/2 - edge_num

mask <- ensemble_mim[upper.tri(ensemble_mim)]
ranked_mask <- rank(mask, na.last = TRUE)
ensemble_mim[upper.tri(ensemble_mim, diag=FALSE)] <- ranked_mask

sig <- ensemble_mim > thresh
dim(sig)
GOLD <- sig * 1
sum(GOLD)

###################################################################

table <- validate(ensemble_clr.net, GOLD)
auc.roc(table)
auc.pr(table)
table2 <- table[1:10e5,]
show.pr(table2, col="red",type="b")
#show.roc(table2, col="red",type="b")
rm(table,table2)

table <- validate(ensemble_aracne.net, GOLD)
auc.roc(table)
auc.pr(table)
table2 <- table[1:10e5,]
show.pr(table2, col="red",type="b")
#show.roc(table2, col="red",type="b")

table <- validate(ensemble_mrnet.net, GOLD)
auc.roc(table)
auc.pr(table)
table2 <- table[1:10e5,]
show.pr(table2, col="red",type="b")
#show.roc(table2, col="red",type="b")

table <- validate(ensemble_mrnetb.net, GOLD)
auc.roc(table)
auc.pr(table)
table2 <- table[1:10e5,]
show.pr(table2, col="red",type="b")
#show.roc(table2, col="red",type="b")

k