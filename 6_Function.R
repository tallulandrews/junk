# GO enrichments of markers
# External markers
# diff CC phase of same functional cell-type?

args <- commandArgs(trailingOnly=TRUE) # SCE RDSs
nSCEs <- length(args)
expr_type <- "lognorm"

SCE_list <- list();
keep_genes <- c()
background_genes <- c()
max_ngenes <- 0;
for (f in args) {
	require("scater")
	obj <- readRDS(f);
	if (nrow(obj) > max_ngenes) {max_ngenes <- nrow(obj)}
#	keep_genes <- c(keep_genes, rownames(obj)[fData(obj)$KW_is.Feature & fData(obj)$pct_dropout < 90]);
	keep_genes <- c(keep_genes, as.character(fData(obj)[ fData(obj)$pct_dropout < 90, "Symbol"]));
	background_genes <- c(background_genes, rownames(obj)[ fData(obj)$pct_dropout <= 99]);
	tmp <- unlist(strsplit(f, "\\."))
	SCE_list[[tmp[1]]] <- obj;
}
keep_genes <- unique(keep_genes);
background_genes <- unique(background_genes);


# Cell-Type Markers
require("CellTypeProfiles")
source("/nfs/users/nfs_t/ta6/NetworkInferencePipeline/Dropouts/My_R_packages/CellTypeProfiles/R/Markers.R")
require("gProfileR")
require("proxy")
n_top_rich = 30;
for (i in 1:nSCEs) {
	obj <- SCE_list[[i]]
	obj_orig <- obj;

	obj <- obj[rownames(obj) %in% background_genes,]
#	markers1 <- complex_markers(get_exprs(obj, expr_type), pData(obj)$clusters_raw, n_max=1)
#	markers1 <- markers1[order(markers1$AUC, descending=TRUE),]
	obj <- obj[, pData(obj)$clusters_clean != "Outliers"]
	markers1 <- complex_markers(get_exprs(obj, expr_type), factor(pData(obj)$clusters_clean), n_max=1)
	markers1 <- markers1[order(markers1$AUC, decreasing=TRUE),]
	markers1 <- markers1[ markers1$q.value < 0.05 & markers1$q.value >= 0 , ]
	m <- markers1[match(rownames(obj), rownames(markers1)),]
	colnames(m) <- paste("Markers_N1", colnames(m), sep="_")
	fData(obj_orig) <- cbind(fData(obj_orig), m)
	SCE_list[[i]] <- obj_orig;
	
	# get gene list

	get_top_markers <- function(group) {
		m <- markers1[markers1$Group == group,]
		m <- m[1:10,]
		return(m)
	}
	get_richments <- function(group) {
		gene_list <- rownames(markers1)[markers1$Group == group]
		   enrichments <- gprofiler(gene_list, organism="hsapiens", 
		   ordered_query=T, significant=T, custom_bg=background_genes, 
		   hier_filtering="moderate", max_set_size=10000,
		   src_filter=c("GO:BP", "KEGG", "REAC", "HPA"), 
		   correction_method="fdr", min_isect_size=3, min_set_size=10)
		enrichments<-enrichments[order(enrichments$p.value),]
		# enrichment filtering
		# remove hpa low
		exclude <- enrichments$domain == "hpa" & grepl("Low", enrichments$term.name)
		exclude <- exclude | enrichments$domain == "hpa" & grepl("Not detected", enrichments$term.name)
		enrichments <- enrichments[!exclude,]
		enrichments <- enrichments[1:min(nrow(enrichments), n_top_rich),]
		enrichments$GroupID = rep(group, times=nrow(enrichments))
		return(enrichments)
## ID similar terms
#	riched_gene_lists<-lapply(enrichments[,14], function(x){strsplit(x, ",")})
#	d <- proxy::dist(riched_gene_lists, function(x, y){length(intersect(unlist(x), unlist(y)))/length(union(unlist(x), unlist(y)))})
#	dendro <- hclust(1-d, method="single")
#	merge <- cutree(dendro, h=0.5)
	}
	top_marker_table <- factor()
	Rich_table <- factor()
	for(g in levels(factor(markers1$Group))) {
		top_marker_table <- rbind(top_marker_table, get_top_markers(g))
		Rich_table <- rbind(Rich_table, get_richments(g))
	}	
	write.table( Rich_table, file=paste(names(SCE_list)[i], "marker_GOrich.txt", sep="_"))
	write.table(top_marker_table, file=paste(names(SCE_list)[i], "TopMarkerTable.txt", sep="_"))
}


## External Markers ##

# Read in Putative Markers
Lit_Markers <- read.delim("~/Collaborations/LiverOrganoids/Literature_Markers.txt", sep="\t", header=TRUE)


CC_Stem <- read.delim("~/Collaborations/LiverOrganoids/Markers_CC_Stemness.txt", sep=" ", header=TRUE)
HCC_Stem <- read.delim("~/Collaborations/LiverOrganoids/Markers_HCC_Stemness.txt", sep=" ", header=FALSE)
CHC_Stem <- read.delim("~/Collaborations/LiverOrganoids/Markers_CHC_Stemness.txt", sep=" ", header=TRUE)


Hepato_MSigdb <- read.delim("~/Collaborations/LiverOrganoids/Markers_HCC_Differentiation.txt")

HPA <- readRDS("/lustre/scratch117/cellgen/team218/TA/OtherDownloadedData/Good_HPA_Hep_Chol_Markers.rds")

Samp <- readRDS("/lustre/scratch117/cellgen/team218/TA/LiverOrganoids/ExternalData/Sampaziotis_MarkersV2.rds")

Hepato_Camp <- read.table("~/Collaborations/LiverOrganoids/FromLaura_Camp_Hepatocyte_Markers.txt")
Stem_Camp <- read.table("~/Collaborations/LiverOrganoids/FromLaura_Camp_Hepatoblast_Markers.txt")
MSC_Camp <- read.table("~/Collaborations/LiverOrganoids/FromLaura_Camp_Mesenchymal_Markers.txt")

CC_1 <- read.table("~/Data/Whitfield_CC.txt")
CC_2 <- read.table("~/Collaborations/LiverOrganoids/New_CC_171117.txt", header=FALSE)
CC_3 <- read.delim("/lustre/scratch117/cellgen/team218/TA/Mirrored_Annotations/GO/hsapiens_80_GO_Annoations_Emsembl.out", sep="\t", header=FALSE)
CC_3 <- CC_3[CC_3[,3] == "cell cycle",1]
source("~/R-Scripts/Ensembl_Stuff.R")
CC_3 <- unique(map_symbol_ensg(as.character(CC_3), is.org="Hsap", is.name="ensg"))

All_CC <- unique(c(as.character(CC_1[,2]), as.character(CC_2[,1]), as.character(CC_3)))

All_Hepato <-c(unique(as.character(Hepato_MSigdb[,1])),
	unique(as.character(Hepato_Camp[,1])),
	unique(as.character(names(Samp$HB))),
	#as.character(HPA[HPA[,8] < HPA[,16],2]),
	unique(as.character(Lit_Markers[Lit_Markers[,2]=="hepatocyte",1]))
		)
Good_Hepato <- All_Hepato[duplicated(All_Hepato)]
All_Chol <- c(unique(as.character(HPA[HPA[,8] > HPA[,16],2])),
	unique(as.character(Lit_Markers[Lit_Markers[,2] %in% c("cholangiocyte", "mature cholangiocyte"),1])),
	unique(as.character(names(Samp$mChol)))
	)
Good_Chol <- All_Chol[duplicated(All_Chol)]
All_Stem <- c(unique(as.character(CC_Stem[CC_Stem[,2]==1,1])),
	unique(as.character(CHC_Stem[CHC_Stem[,2]==1,1])),
	unique(as.character(HCC_Stem[HCC_Stem[,2]==1,1])),
	unique(as.character(Stem_Camp[,1])), 
	unique(names(Samp$Stem)))
Good_Stem <- All_Stem[duplicated(All_Stem)]

All_MSC <- c(as.character(MSC_Camp[,1]))

All_Hepato <- unique(All_Hepato); All_Hepato <- All_Hepato[!is.na(All_Hepato)]
All_Chol <- unique(All_Chol); All_Chol <- All_Chol[!is.na(All_Chol)]
All_Stem <- unique(All_Stem); All_Stem <- All_Stem[!is.na(All_Stem)]
All_MSC <- unique(All_MSC); All_MSC <- All_MSC[!is.na(All_MSC)]

# Remove any gene that appears on conflicting lists
All_Markers <- c(All_Hepato, All_Chol, All_Stem, All_MSC, All_CC)
exclude <- All_Markers[duplicated(All_Markers)]
All_Markers <- unique(All_Markers); All_Markers <- All_Markers[!All_Markers %in% exclude];
All_Hepato <- All_Hepato[!All_Hepato %in% exclude]
All_Chol <- All_Chol[!All_Chol %in% exclude]
All_Stem <- All_Stem[!All_Stem %in% exclude]
All_MSC <- All_MSC[!All_MSC %in% exclude]

for (i in 1:nSCEs) {
	obj <- SCE_list[[i]]
	fData(obj)$MType <- rep("", times=nrow(obj));
	fData(obj)$MType[fData(obj)$Symbol %in% All_Hepato] <- "Hepato"
	fData(obj)$MType[fData(obj)$Symbol %in% All_Chol] <- "Chol"
	fData(obj)$MType[fData(obj)$Symbol %in% All_Stem] <- "Stem"
	fData(obj)$MType[fData(obj)$Symbol %in% All_MSC] <- "MSC"
	fData(obj)$MType[fData(obj)$Symbol %in% All_CC] <- "CC"
	SCE_list[[i]] <- obj;
}

# Test if overrepresented in kept genes:
kept_test <- function(x) {
	q <- sum(x %in% keep_genes)
	m <- length(keep_genes)
	n <- max_ngenes-m
	k <- length(x)
	return(phyper(q, m, n, k, lower.tail=FALSE))
}
if (kept_test(All_Hepato) > 10^-5) {print("Warning Hepato markers not enriched among kept genes")}
if (kept_test(All_Chol) > 10^-5) {print("Warning Chol markers not enriched among kept genes")}
if (kept_test(All_Stem) > 10^-5) {print("Warning Stem markers not enriched among kept genes")}
if (kept_test(All_MSC) > 10^-5) {print("Warning MSC markers not enriched among kept genes")}

### Move the below to a separate script I think ###
# Should I look for GO enrichments for each cell-type? then weed out genes not associated with an enriched function?


# Remove low & non-DE genes from the lists
All_Hepato <- All_Hepato[All_Hepato %in% keep_genes]
All_Chol <- All_Chol[All_Chol %in% keep_genes]
All_Stem <- All_Stem[All_Stem %in% keep_genes]
All_MSC <- All_MSC[All_MSC %in% keep_genes]
All_Markers <- All_Markers[All_Markers %in% keep_genes]
All_Markers <- All_Markers[!All_Markers %in% All_CC]

require("CellTypeProfiles");

	
# gene-gene correlations of markers -> cluster gene-gene network and ID clusters enriched for markers with a particular annotation. 
# Clustering = complete, theshold = ?

mark_cols = c("forestgreen", "navy", "purple", "salmon")

Score <- rep(0, times=length(All_Markers))
names(Score) <- All_Markers
for (i in 1:nSCEs) {
	obj <- SCE_list[[i]]
	obj <- calculateQCMetrics(obj)
	obj <- obj[fData(obj)$Symbol %in% All_Markers & fData(obj)$pct_dropout < 99,]
	expr_mat <- get_exprs(obj, expr_type)
	rownames(expr_mat) <- fData(obj)$Symbol
	c_mat <- cor(t(expr_mat), method="spearman")
#	heatout <- heatmap.3(c_mat, distfun=function(x){as.dist(1-x)}, 
#		hclustfun=function(x){hclust(x, method="ward.D")},
#		ColSideColors=matrix(mark_cols[factor(fData(obj)$MType)], nrow=nrow(c_mat)), ColSideColorsSize=1)
	score <- apply(c_mat, 1, min)
	score2 <- apply(c_mat, 1, function(x){max(x[x < 1])})
	overall <- apply(cbind(rank(-1*score),rank(score2)), 1, min)
	overall <- overall[match(names(Score), names(overall))]
	overall[is.na(overall)] <- 0
	Score <- Score+overall
}

Score <- rev(sort(Score))

all_plot <- function(gene) {
	require("RColorBrewer")
	get_splits <- function(gene) {
		vals = c();
		for (obj in SCE_list) {
			mat <- get_exprs(obj, expr_type)
			vals <- c(vals, mat[fData(obj)$Symbol == gene, ]) 
		}
		splits <- seq(from=min(vals), to=max(vals), length=5)
		return(splits);
	}

	splits <- get_splits(gene);
	a <- ceiling(sqrt(length(SCE_list)))
	par(mfrow=c(a,a))
	par(mar=c(1,4,3,1))
	
	for( i in 1:length(SCE_list) ) {
		obj <- SCE_list[[i]]
		binned <- cut( get_exprs(obj, expr_type)[fData(obj)$Symbol == gene,], breaks=splits, include.lowest=TRUE )
		my_cols <- colorRampPalette(brewer.pal(6, "Blues"))(length(splits)-1)

		my_title <- names(SCE_list)[i]
		my_title <- unlist(strsplit(my_title, "_"))
		my_title <- my_title[1]
		pt_col <- my_cols[binned]
		plot(obj$DM1, obj$DM2, pch=16, main=my_title, xaxt="n", yaxt="n", ylab=gene, bty="n", col=pt_col)
	} 	
}

# Or is this more informative:
i = 1
obj <- SCE_list[[i]]
Hep_score <- colSums(get_exprs(obj, expr_type)[fData(obj)$Symbol %in% Good_Hepato,])
Chol_score <- colSums(get_exprs(obj, expr_type)[fData(obj)$Symbol %in% Good_Chol,])
Stem_score <- colSums(get_exprs(obj, expr_type)[fData(obj)$Symbol %in% Good_Stem,])
i = 2
obj <- SCE_list[[i]]
Hep_score<- c(Hep_score, colSums(get_exprs(obj, expr_type)[fData(obj)$Symbol %in% Good_Hepato,]))
Chol_score<- c(Chol_score, colSums(get_exprs(obj, expr_type)[fData(obj)$Symbol %in% Good_Chol,]))
Stem_score<- c(Stem_score, colSums(get_exprs(obj, expr_type)[fData(obj)$Symbol %in% Good_Stem,]))
i = 3
obj <- SCE_list[[i]]
Hep_score<- c(Hep_score, colSums(get_exprs(obj, expr_type)[fData(obj)$Symbol %in% Good_Hepato,]))
Chol_score<- c(Chol_score, colSums(get_exprs(obj, expr_type)[fData(obj)$Symbol %in% Good_Chol,]))
Stem_score<- c(Stem_score, colSums(get_exprs(obj, expr_type)[fData(obj)$Symbol %in% Good_Stem,]))
donor_id <- rep(names(SCE_list),  times=c(ncol(SCE_list[[1]]), ncol(SCE_list[[2]]), ncol(SCE_list[[3]])))
plot(Stem_score, Hep_score, col=c("forestgreen", "purple", "black")[factor(donor_id)], pch=16)
plot(Chol_score, Hep_score, col=c("forestgreen", "purple", "black")[factor(donor_id)], pch=16)
plot(Chol_score, Stem_score, col=c("forestgreen", "purple", "black")[factor(donor_id)], pch=16)
plot(Hep_score-Chol_score, Stem_score, col=c("forestgreen", "purple", "black")[factor(donor_id)], pch=16, xlab="chol - hepa", ylab="stem-ness")


