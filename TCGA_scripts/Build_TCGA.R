ann = read.delim("transcriptomes_metadata.txt", sep="\t", header=F)
data = read.delim("ExprMat.txt", sep="\t", header=F, stringsAsFactor=FALSE)
data[1,] <- c("000000_Gene", data[1,1:(length(data[1,])-1)])
data <- data[,order(data[1,])]
ann <- ann[order(ann[,1]),]
ann <- ann[ann[,1] %in% data[1,],]
colnames(data) <- data[1,]
rownames(data) <- data[,1]
data<- data[-1,-1]
data <- data.frame(data)
tmp <- apply(data, 2, as.numeric)
rownames(tmp) <- rownames(data)
obj <- list(data=tmp, ann=ann)
saveRDS(obj, "TCGA_Transcriptome.rds")

ann2 = read.delim("Chol_metadata.txt", sep="\t", header=F)
data2 = read.delim("Chol_ExprMat.txt", sep="\t", header=F, stringsAsFactor=FALSE)
data2[1,] <- c("000000_Gene", data2[1,1:(length(data2[1,])-1)])
data2 <- data2[,order(data2[1,])]
ann2 <- ann2[order(ann2[,1]),]
ann2 <- ann2[ann2[,1] %in% data2[1,],]
colnames(data2) <- data2[1,]
rownames(data2) <- data2[,1]
data2 <- data2[-1,-1]
data2 <- data.frame(data2)
tmp2 <- apply(data2, 2, as.numeric)
rownames(tmp2) <- rownames(data2)
obj <- list(data=tmp2, ann=ann2)
saveRDS(obj, "TCGA_Chol_Transcriptome.rds")

comb_tmp <- cbind(tmp, tmp2)
comb_ann <- rbind(ann, ann2)

obj <- list(data=comb_tmp, ann=comb_ann)
saveRDS(obj, "All_TCGA_Transcriptome.rds")
