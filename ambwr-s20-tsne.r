### AMBWR - Automated Metagenomic Binning (With Receipts)

#metadata
metadata <- read.table("metagenomics_metadata_ids.txt",sep="\t",head=T,stringsAsFactors=F)
# rownames(metadata) <- metadata$Sample.name
rownames(metadata) <- metadata$plate_uniq_id


### dimensional reduction
#first with Fast PCA
# covtab <- unique(input.num) #may not want to do it this way, could lose contigs, instead pull !duplicated(input.meta$node)
covtab <- input.num[!duplicated(input.meta$node),]
cov.pca <- FastPCA(covtab,50) #centers, doesn't scale automatically
cov.pca.uniq <- cov.pca$x
rownames(cov.pca.uniq) <- rownames(covtab)
colnames(covtab) <- gsub("-","_",x=gsub(".mapped_sorted.bam","",x=colnames(covtab)))


# bh-tsne 2D
# cov.tsne.list <- Rtsne(cov.pca.uniq, dims=2, max_iter=2000, verbose=TRUE, theta=0.5, pca=FALSE, check_duplicates=F)
# saveRDS(cov.tsne.list,file=file.path(outdir,"tsne-embed-2d.rds"))
cov.tsne.list <- readRDS(file=file.path(outdir,"tsne-embed-2d.rds"))
cov.tsne <- cov.tsne.list$Y
rownames(cov.tsne) <- rownames(covtab)
colnames(cov.tsne) <- c("x","y")

#---


#### extract single-linkage clusters for re-assembly
#for testing only

library(fastcluster) #for fast hclust()
library(dendextend) #for cutree()
library(ggplot2)
library(ape) #for cophenetic()
#library(parallelDist)

tsne.dist <- dist(cov.tsne) # slow, could use parDist
tsne.hclust <- fastcluster::hclust(tsne.dist, method="single")
# it's so fast I'm crying, and works on long vectors

tsne.td <- cophenetic(tsne.hclust) #slow
tsne.tmat <- as.matrix(tsne.td)

#test correlation between raw distance matrix and clustered distances
# cor(tsne.tmat[1:2^30], tsne.dmat[1:2^30]) # = 0.4, logarithmic better fit
# aa<-sample(1:length(tsne.td),1e5)
# plot(tsne.tmat[aa],sqrt(tsne.dmat[aa]),pch='.')

rm(tsne.td)
rm(tsne.dist)

covprop <- prop.table(as.matrix(covtab),mar=2)
covprop[is.nan(covprop)] <- 0

covmeans <- rowMeans(covprop)
covsd <- apply(covprop,1,sd)

# save.image(paste0(prefix,".Rdata"))

source("bin/ambwr-s30-loadmoredata.r",echo=T)