library(WGCNA)
library(Seurat)
library(tidyverse)
options(stringsAsFactors = FALSE)
#args = commandArgs(T)
enableWGCNAThreads(8)
#exprMat <- args[1]
calculate_fpkm <- function(counts, lengths) {
  total_counts <- colSums(counts)
  fpkm <- t(t(counts) / total_counts) * 1e6  # 转换为每百万映射的reads
  fpkm <- fpkm / (lengths / 1e3)  # 转换为每千碱基
  return(fpkm)
}
calculate_tpm <- function(fpkm) {
  scaling_factor <- colSums(fpkm)
  tpm <- t(t(fpkm) / scaling_factor) * 1e6
  return(tpm)
}

seurat_obj <- readRDS("/data/wangmeiyue/cooperation/WK_wheat_regenerate/03divGene/L7.rds")
avg_exp <- as.data.frame(AggregateExpression(seurat_obj,group.by = "seurat_clusters")$Spatial)
exp <- read.table("/data/wangmeiyue/CS_multi_tissue_and_treat/lncRNA/all_sample_lnc_count.GeneV1.1.20200901.txt",sep = "\t",stringsAsFactors = F,row.names = 1,header = T)
len <- exp[rownames(avg_exp),"Length"]
fpkm <- calculate_fpkm(avg_exp,len)
tpm <- calculate_tpm(na.omit(fpkm))
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)
dataExpr <- tpm
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > 0),]
dataExpr <- as.data.frame(t(dataExprVar))
#save(dataExpr,file = "dataExpr.RData")
gsg = goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5)
power = sft$powerEstimate
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=F,
                       #saveTOMFileBase = paste0("WGCNA_net",".tom"),
                       verbose = 3)
save(net,file = "WGCNA_net.RData")
GM <- data.frame(net$colors)
GM$id <- rownames(GM)
write.table(GM,"GMs.out",row.names=F,col.names=T,quote=F,sep="\t")

adjacency <- adjacency(dataExpr, power = power)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM
save(TOM,file="TOM.RData")

