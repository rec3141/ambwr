# AMBWR - Automated Metagenomic Binning (With Receipts)

#install.packages(c("gatepoints","Rtsne","Matrix","rgl","scales"))
#install_github("zarquon42b/Morpho", local=FALSE)
# had to module load some things to get these to compile
# library(umap)
# library(NbClust)
# library(mclust)
# library(cluster)
# library(MCL)
# library("FactoMineR")
# library("factoextra")
# library(klaR)
# library(cba)
# library(parallelDist)


### IDEAS
# paths to better bins
# do reassembly with semiperfect reads or longer kmers with 5% error
# add tnf clustering?
# add scaffolds < 2000 bp as 2nd class, don't use them for first stage single link clustering step

# cleanup step to remove isolated contigs that may act as linkers? machine learning
# remove clusters with 5 or fewer contigs at cutoff=1.0/2
# which(as.numeric(bintable$contig_count) <6) #282 bins
# add to bin_0

#add graph linkages to small tsne plot
#bin accumulation plots -- as you decrease threshold, how many bins do you get at each threshold?
#what is the size distribution of those bins?
#how many connected components per bin? 

#ternery diagram animation of % of marker genes + trnas + viral genes as a function of cutoff; color by redundancy (green, yellow, red); size by ratio of nonhypo
#plot of contig length versus % mapped to genomes
#plot of cumulative length in contigs sorted by contig length
### IDEAS

library(gatepoints)
library(Rtsne)
# library(Morpho)
library(Matrix)
library(data.table)
library(gplots)
# library(rgl)
library(scales)
library(parallel)
library(igraph)
source("bin/FastPCA.R")

options(stringsAsFactors = FALSE)

outdir = "annotation/kaiju/spades_all_norm3_combined20190220_2019-02-20-13-34-25/scaffolds-gt2000.assembly"
prefix = "scaffolds-gt2000-2"
headerfile = "mapped/spades_all_norm3_combined20190220_2019-02-20-13-34-25/combined20190220_1/combined20190220.headers.txt"
header.in = read.table(headerfile,stringsAsFactors=F)
kaijudir = "annotation/kaiju"
bindir = file.path("annotation/bins",prefix)
dir.create(bindir)


# eventually get rid of the bash crap and do it all in R

inputfile = "annotation/kaiju/spades_all_norm3_combined20190220_2019-02-20-13-34-25/scaffolds-gt2000.assembly/combined20190220.node_orf_bin_ann_tax_depths.tsv"
inputDT <- fread(inputfile,sep="\t",head=F,stringsAsFactors=F)

#headers = c("orf","taxonomy","node","metabat","description","contiglen","totalAvgDepth", as.character(header.in[4:length(header.in)]))
headers=c("node","orf","metabat","description","taxonomy","contiglen","totalAvgDepth",as.character(header.in[4:length(header.in)]))
colnames(inputDT) <- headers

rownames(inputDT) <- make.unique(inputDT$node)
input.meta <- as.data.frame(inputDT[,1:7,with=F])
rownames(input.meta) <- rownames(inputDT)

input.num <- as.data.frame(inputDT[,8:ncol(inputDT),with=F])
input.num <- input.num[,grep(".var",colnames(input.num),invert=T)]
rownames(input.num) <- rownames(inputDT)

rm(inputDT)

source("bin/ambwr-s20-tsne.r",echo=T)