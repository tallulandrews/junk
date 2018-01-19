source("~/Collaborations/LiverOrganoids/Laura_Pipeline/0_ColourScheme.R")

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

my_clustering <- function(SCE, k, suppress.plot=TRUE) {
	require("CellTypeProfiles")
	
	consensus_mat <- SCE@sc3$consensus[[k]]
        consensus_mat <- consensus_mat[[1]]
        distances <- as.dist(1-consensus_mat)
        Cs <- cutree(hclust(distances, method="complete"), k=as.numeric(k))
        Cs2 <- cutree(hclust(distances, method="complete"), h=1-10^-10)
	if (max(Cs2) > max(Cs)) {Cs <-Cs2}

	Cns <- CellTypeProfiles::factor_counts(factor(Cs))

	if (!suppress.plot) {
		palette <- cluster_col(max(Cs))
		require("gplots")
		heatmap.2(consensus_mat, distfun=function(x){as.dist(1-x)}, hclustfun=function(x){hclust(x, method="complete")},
			trace="none", ColSideColors=palette[Cs])
	}

	require("cluster")
	sil <- cluster::silhouette(Cs, distances)
	#m <- mean(sil[,3])
	#s <- sd(sil[,3])
	#p <- p.adjust(pnorm(abs(sil[,3]-m)/s, lower.tail=FALSE), method="fdr")
	silhouettes <- mean(sil[,3])
	neighbours <- table(sil[,1], sil[,2])/Cns

	sets <- split(seq(length(Cs)), factor(Cs))
        score <- sapply(sets, function(a) {mean(consensus_mat[a,a]) - mean(consensus_mat[a,-a])})
        scores <- sum(score*Cns)/ncol(consensus_mat)

	return(list(sil_nn=neighbours, c_scores=score, overallSil=silhouettes, overallScore=scores, Cs=Cs));
}

# Overall cluster scores at each k
silhouettes <- vector()
scores <- vector()

for (i in names(SCE@sc3$consensus)) {
	clust_name <- paste("sc3", i, "clusters", sep="_");
	stuff <- my_clustering(SCE, i)

	pData(SCE)[,clust_name] <- stuff$Cs;

	silhouettes <- c(silhouettes, stuff$overallSil)
	scores <- c(scores, stuff$overallScore)
}

composite_score <- apply(cbind(silhouettes, scores), 1, min);
optimal_k <- which(composite_score == max(composite_score))+1 # Select optimal K
fine_k <- optimal_k
if (optimal_k < k_min) {
	composite_score[optimal_k-1] = 0
	fine_k <- which(composite_score == max(composite_score))+1
}

print("optimal k found")

png(paste(outprefix, "k", optimal_k, "SC3_consensus.png", sep="_"), width=7, height=7, units="in", res=300)
out <- my_clustering(SCE, as.character(optimal_k), suppress.plot=FALSE)
dev.off()

pData(SCE_orig)$clusters_coarse <- factor(out$Cs); # Save optimal K clusters

png(paste(outprefix, "k", fine_k, "SC3_consensus.png", sep="_"), width=7, height=7, units="in", res=300)
out <- my_clustering(SCE, as.character(fine_k), suppress.plot=FALSE)
dev.off()

pData(SCE_orig)$clusters_fine <- factor(out$Cs); # Save optimal K clusters


## Refine Clusters ##
# calculate markers & refine clusters
# save refined clusters

# Calculate Markers
require("CellTypeProfiles")

markers <- complex_markers(get_exprs(SCE, expr_type), pData(SCE_orig)$clusters_fine)
sig <- markers[markers$q.value < 0.05 & markers$q.value > 0,]
good <- sig[sig$AUC > 0.7,]

assigned <- good[,-c(1, ncol(good), ncol(good)-1)]
a_names <- colnames(assigned)
assigned <- t(apply(assigned, 1, function(a){
		a <- as.numeric(a)
		if (mean(a) > 0.5){return (1-a)}
		else {return(a)}
		}))
colnames(assigned) <- a_names

print("markers calculated")

# Unique vs Shared Markers
key_markers <- assigned[rowSums(assigned) == 1,]
shared_markers <- assigned[rowSums(assigned) > 1,]
# How evenly are they spread?
expected <- qbinom(0.05/ncol(assigned), size=nrow(key_markers), prob=1/ncol(key_markers))
# How big are the clusters?
nCs <- factor_counts(pData(SCE_orig)$clusters_fine)
min_C_size <- ncol(SCE)*0.05

keepC <- colnames(assigned)[colSums(key_markers) > expected] # lots of unique markers == real

# Do the Refining
allC <- as.character(levels(pData(SCE_orig)$clusters_fine))
rawCs <- as.character(pData(SCE_orig)$clusters_fine)
newCs <- rawCs
for(i in allC) {
	new <- i
	if (nCs[as.numeric(i)] == 1) {
		new <- "Outliers"
		newCs[rawCs==i] <- new[1];
		next;
	}
	# Small & cohesive = outliers
	if (min_C_size > nCs[as.numeric(i)] & out$c_scores[as.numeric(i)] >= min(out$c_scores[as.numeric(keepC)])) {
		new <- "Outliers"
	} 
	if ( !(i %in% keepC) ) {
		# if more cohesive than a keptC keep it too
		if (out$c_scores[as.numeric(i)] > min(out$c_scores[as.numeric(keepC)])) {
                        new <- i
                } else {
			# Proportion of shared markers that are shared with cluster X
			this_C_markers <- shared_markers[,colnames(shared_markers) == i] == 1
			m_score <- colSums(shared_markers[this_C_markers,])/sum(this_C_markers)
			m_score[colnames(shared_markers) == i] <- 0
			# Proportion of cells where the next closest cluster is cluster X
			sil_score <- out$sil_nn[as.numeric(i),]

			m_closest <- names(m_score)[m_score == max(m_score)]
			sil_closest <- names(sil_score)[sil_score == max(sil_score)]

			if ( length(intersect(m_closest, sil_closest)) == 1) {
				# Aggreement on closest cluster
				if (max(m_score) > 0.5 & max(sil_score) > 0.75) {
					# similar enough to be equivalent
					new <- unique(newCs[rawCs==intersect(m_closest, sil_closest)])
				}
			} else {
				if (min_C_size > nCs[as.numeric(i)]) {
					new <- "Outliers"
				}
			}
		}
	}
	newCs[rawCs==i] <- new[1];
}

print("clusters refined")

pData(SCE_orig)$clusters_clean <- factor(newCs);

identical(rownames(markers), rownames(fData(SCE_orig)))
colnames(markers) <- paste("fine_marker", colnames(markers), sep="_");
fData(SCE_orig) <- cbind(fData(SCE_orig), markers)

saveRDS(SCE_orig, file=paste(outprefix,"SC3.rds", sep="_"))
