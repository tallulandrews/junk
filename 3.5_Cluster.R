source("~/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")
source("~/Collaborations/LiverOrganoids/Laura_Pipeline/99_Functions.R")

args <- commandArgs(trailingOnly=TRUE) # rds file for data, prefix for output, number of cores
n.cores = args[3];
outprefix = args[2];
expr_type <- "lognorm"
k_min = 4
if (length(args) > 3) {
	expr_type = args[4];
}

require("scater")
SCE <- readRDS(args[1]);
#SCE <- readRDS("Test_QCed_BatchCor.rds")
#SCE <- SCE[,SCE$Type=="CCA1"]

set.seed(123)

max_ks <- ceiling(dim(SCE)[2]/25)

# SC3 - farm job 4 cores
require("SC3")

SCE_orig <- SCE

# Remove CC genes
CC_genes <- load_CC("cycling")
CC_genes <- c(as.character(CC_genes$Whitfield[,2]), as.character(CC_genes$Tirosh[,1]))
SCE <- SCE[!( fData(SCE)$feature_symbol %in% CC_genes ),]

# SCE Clustering
exprs(SCE) <- get_exprs(SCE, expr_type)
SCE <- sc3_prepare(SCE, ks=2:max_ks)
SCE <- sc3_estimate_k(SCE)
estimated_k <- SCE@sc3$k_estimation # Does this work?
if( estimated_k > max_ks ) {
	max_ks <- estimated_k
}
SCE <- sc3(SCE, ks=2:max_ks, n_cores=as.numeric(n.cores), biology=FALSE)

print("SC3 - finished")
## Redo clustering ##
# use consensus agreement as distance directly for clustering
# Calculate silhouette & stability & internal-vs-external consistency
# Algorithmicly select k
# save best k clusters

# Overall cluster scores at each k

reclustered_optim <- get_optimal_k(SCE)
optimal_k <- reclustered_optim$optim
fine_k <- reclustered_optim$fine

png(paste(outprefix, "k", optimal_k, "SC3_noCC_consensus.png", sep="_"), width=7, height=7, units="in", res=300)
out <- my_clustering(SCE, as.character(optimal_k), suppress.plot=FALSE)
dev.off()

pData(SCE_orig)$clusters_noCC_coarse <- factor(out$Cs); # Save optimal K clusters

png(paste(outprefix, "k", fine_k, "SC3_noCC_consensus.png", sep="_"), width=7, height=7, units="in", res=300)
out <- my_clustering(SCE, as.character(fine_k), suppress.plot=FALSE)
dev.off()

pData(SCE_orig)$clusters_noCC_fine <- factor(out$Cs); # Save optimal K clusters


## Refine Clusters ##
# calculate markers & refine clusters
# save refined clusters

refined_clusters <- refine_clusters(SCE, expr_type)
newCs <- refined_clusters$newCs
markers <- refined_clusters$markers

pData(SCE_orig)$clusters_clean <- factor(newCs);

identical(rownames(markers), rownames(fData(SCE_orig)))
colnames(markers) <- paste("fine_noCC_marker", colnames(markers), sep="_");
fData(SCE_orig) <- cbind(fData(SCE_orig), markers)

saveRDS(SCE_orig, file=paste(outprefix,"SC3_noCC.rds", sep="_"))
