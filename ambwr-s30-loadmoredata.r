
###### LOAD MORE INPUT DATA

gfafile = "assembly/spades/spades_all_norm3_combined20190220_2019-02-20-13-34-25/assembly_graph_with_scaffolds.gfa"
scaffile = "assembly/spades/spades_all_norm3_combined20190220_2019-02-20-13-34-25/scaffolds-gt2000.assembly.fasta"

virfinderfile = "annotation/virfinder/virfinder.pred.tsv"
virsorterfile = "annotation/virsorter/VIRSorter_global-phage-signal.csv"
virsortercircfile = "annotation/virsorter/fasta/VIRSorter_circu.list"

# oldmarkerfile = "bacteria_marker_genes.txt"
bactinfofile = "annotation/gtdbtk/gtdbtk.bac120.marker_info.tsv"
archinfofile = "annotation/gtdbtk/gtdbtk.ar122.marker_info.tsv"
eukinfofile = "annotation/busco/genes.txt"
pfammarkfile = "annotation/gtdbtk/marker_genes/scaffolds-gt2000/scaffolds-gt2000_pfam_tophit.tsv"
tigrmarkfile = "annotation/gtdbtk/marker_genes/scaffolds-gt2000/scaffolds-gt2000_tigrfam_tophit.tsv"
eukmarkfile = "annotation/busco/scaffolds-gt2000/scaffolds-gt2000.tsv"

bacrrnafile = "annotation/barrnap/scaffolds-gt2000/bac.rrna.gff"
arcrrnafile = "annotation/barrnap/scaffolds-gt2000/arc.rrna.gff"
eukrrnafile = "annotation/barrnap/scaffolds-gt2000/euk.rrna.gff"
mitorrnafile = "annotation/barrnap/scaffolds-gt2000/mito.rrna.gff"



### FIND MARKER CONTIGS

input.meta$orfcount <- 1
input.meta$contigcount <- rownames(input.meta) %in% input.meta$node + 0
input.meta$contigsize <- input.meta$contigcount * input.meta$contiglen
input.meta$contigcov <- input.meta$contigcount * input.meta$totalAvgDepth

archmarkers <- cbind("Archaea",read.table(archinfofile,stringsAsFactors=F,sep="\t",quote="",head=1))
bactmarkers <- cbind("Bacteria",read.table(bactinfofile,stringsAsFactors=F,sep="\t",quote="",head=1))
colnames(archmarkers) <- c("domain","hmm","name","description","length")
colnames(bactmarkers) <- c("domain","hmm","name","description","length")
markerinfo <- rbind(archmarkers,bactmarkers,stringsAsFactors=F)
markerinfo$domain[duplicated(markerinfo$hmm) | duplicated(markerinfo$hmm,fromLast=T)] <- "Archaea+Bacteria"

pfam.in <- read.table(pfammarkfile,sep="\t",quote="",head=T,row.names=1)
pfam.in[,1] <- gsub(";.*$","",pfam.in[,1]) #just drop extra hits per orf
pfam.marks <- do.call(rbind,strsplit(pfam.in[,1],","))
pfam.marks <- as.data.frame(pfam.marks,stringsAsFactors=F)
rownames(pfam.marks) <- gsub("_\\d+$","",rownames(pfam.in))
colnames(pfam.marks) <- c("hmm","evalue","bitscore")
pfam.marks$evalue <- as.numeric(pfam.marks$evalue)
pfam.marks$bitscore <- as.numeric(pfam.marks$bitscore)

tigr.in <- read.table(tigrmarkfile,sep="\t",quote="",head=T,row.names=1)
tigr.in[,1] <- gsub(";.*$","",tigr.in[,1]) #just drop extra hits per orf
tigr.marks <- do.call(rbind,strsplit(tigr.in[,1],","))
tigr.marks <- as.data.frame(tigr.marks,stringsAsFactors=F)
rownames(tigr.marks) <- gsub("_\\d+$","",rownames(tigr.in))
colnames(tigr.marks) <- c("hmm","evalue","bitscore")
tigr.marks$evalue <- as.numeric(tigr.marks$evalue)
tigr.marks$bitscore <- as.numeric(tigr.marks$bitscore)

no_col <- max(count.fields(eukmarkfile, sep = ""))
euk.in <- read.table(eukmarkfile,header=F,sep="",quote="",comment="#",fill=T, col.names=paste("V",1:no_col,sep=""))
euk.marks <- euk.in[,c(1,3,5,6)]
colnames(euk.marks) <- c("orf","hmm","evalue","bitscore")
euk.marks <- euk.marks[order(euk.marks$evalue),]
euk.marks <- euk.marks[!duplicated(euk.marks$orf),] #remove any duplicates
rownames(euk.marks) <- euk.marks$orf
euk.marks$orf <- NULL

eukmarkers <- read.table(eukinfofile,head=T)
eukmarkers <- cbind("Eukarya",eukmarkers$gene,"","","")
colnames(eukmarkers) <- colnames(markerinfo)
markerinfo <- rbind(markerinfo,eukmarkers)

#two orfs weren't found, lost somewhere? let's get rid of them here to avoid problems downstream
# do it for all of them
euk.marks <- euk.marks[rownames(euk.marks) %in% input.meta$orf,]
tigr.marks <- tigr.marks[rownames(tigr.marks) %in% input.meta$orf,]
pfam.marks <- pfam.marks[rownames(pfam.marks) %in% input.meta$orf,]

input.meta$tigr.mark <- input.meta$orf %in% rownames(tigr.marks) + 0
input.meta$pfam.mark <- input.meta$orf %in% rownames(pfam.marks) + 0
input.meta$euk.mark <- input.meta$orf %in% rownames(euk.marks) + 0 

input.meta$marker <- (input.meta$tigr.mark | input.meta$pfam.mark | input.meta$euk.mark) + 0
input.meta$mark.hmm <- NA
input.meta$mark.hmm[input.meta$tigr.mark > 0] <- tigr.marks[input.meta$orf[input.meta$tigr.mark > 0],"hmm"]
input.meta$mark.hmm[input.meta$pfam.mark > 0] <- pfam.marks[input.meta$orf[input.meta$pfam.mark > 0],"hmm"]
input.meta$mark.hmm[input.meta$euk.mark > 0] <- euk.marks[input.meta$orf[input.meta$euk.mark > 0],"hmm"]

#for each mark.hmm, get matching domain name from markerinfo
input.meta$bact.mark <- (grepl("Bacteria",markerinfo$domain[match(input.meta$mark.hmm, markerinfo$hmm)]))+0
input.meta$arch.mark <- (grepl("Archaea",markerinfo$domain[match(input.meta$mark.hmm, markerinfo$hmm)]))+0
input.meta$euk.mark <- (grepl("Eukarya",markerinfo$domain[match(input.meta$mark.hmm, markerinfo$hmm)]))+0
input.meta$bact.mark[is.na(input.meta$bact.mark)] <- 0
input.meta$arch.mark[is.na(input.meta$arch.mark)] <- 0
input.meta$euk.mark[is.na(input.meta$euk.mark)] <- 0

#some of the markers are never found (wrong evalues cutoffs?)
# bactmarkcount <- sum(grepl("Bacteria",markerinfo$domain))
# archmarkcount <- sum(grepl("Archaea",markerinfo$domain))
# eukmarkcount <- sum(grepl("Eukarya",markerinfo$domain))

bactmarkcount <- length(unique(input.meta$mark.hmm[input.meta$bact.mark > 0]))
archmarkcount <- length(unique(input.meta$mark.hmm[input.meta$arch.mark > 0]))
eukmarkcount <- length(unique(input.meta$mark.hmm[input.meta$euk.mark > 0]))


### FIND VIRUS CONTIGS

#virfinder / doesn't work in screen
# source("https://bioconductor.org/biocLite.R")
# biocLite("qvalue")
# source("parVF.pred")
# vir.pred <- parVF.pred(scaffile)
# write.table(vir.pred,file="virfinder.pred.tsv",sep="\t", quote=F)
# run again on ORFs instead of contigs

vf.pred <- read.table(file=virfinderfile,sep="\t",quote="")
vf.nodes <- vf.pred$name[vf.pred$pvalue < 0.01]
vflike.nodes <- vf.pred$name[vf.pred$pvalue < 0.05 & vf.pred$pvalue >= 0.01]

input.meta$vf.phage1 <- rownames(input.meta) %in% vf.nodes + 0
input.meta$vf.phage2 <- rownames(input.meta) %in% vflike.nodes + 0

#virsorter / very slow (overnight)

vs.pred <- read.csv(virsorterfile,quote="",head=F,stringsAsFactors=F)
vs.head <- read.table(virsorterfile,sep=",",quote="",head=T,comment="",skip=1,nrows=1)
colnames(vs.pred) <- colnames(vs.head)
colnames(vs.pred)[1] <- "vs.id"
vs.cats <- do.call(rbind,strsplit(x=as.character(vs.pred$vs.id[grep("category",vs.pred$vs.id)]),split=" - ",fixed=T))[,2]
vs.num <- grep("category",vs.pred$vs.id)
vs.pred$Category <- paste(c("",rep.int(vs.cats, c(vs.num[-1],nrow(vs.pred)) - vs.num)), vs.pred$Category)
vs.pred <- vs.pred[-grep("#",vs.pred$vs.id),]

vs.pred$node <- gsub(pattern="VIRSorter_",replacement="",x=vs.pred$vs.id)
vs.pred$node <- gsub(pattern="-circular",replacement="",x=vs.pred$node)
vs.pred$node <- gsub(pattern="cov_(\\d+)_",replacement="cov_\\1.",x=vs.pred$node,perl=F)

input.meta$vs.phage1 <- rownames(input.meta) %in% vs.pred$node[vs.pred$Category=="Complete phage contigs 1"] + 0
input.meta$vs.phage2 <- rownames(input.meta) %in% vs.pred$node[vs.pred$Category=="Complete phage contigs 2"] + 0
input.meta$vs.phage3 <- rownames(input.meta) %in% vs.pred$node[vs.pred$Category=="Complete phage contigs 3"] + 0
input.meta$vs.prophage1 <- rownames(input.meta) %in% vs.pred$node[vs.pred$Category=="Prophages 1"] + 0
input.meta$vs.prophage2 <- rownames(input.meta) %in% vs.pred$node[vs.pred$Category=="Prophages 2"] + 0
input.meta$vs.prophage3 <- rownames(input.meta) %in% vs.pred$node[vs.pred$Category=="Prophages 3"] + 0

rm(vs.head,vs.cats,vs.num)

vs.circ <- read.table(virsortercircfile)
colnames(vs.circ) <- c("vs.id","length")
vs.circ$node <- gsub(pattern="VIRSorter_",replacement="",x=vs.circ$vs.id)
vs.circ$node <- gsub(pattern="-circular",replacement="",x=vs.circ$node)
vs.circ$node <- gsub(pattern="cov_(\\d+)_",replacement="cov_\\1.",x=vs.circ$node,perl=F)
vs.circ$circle <- 1

input.meta$vs.circ <- rownames(input.meta) %in% vs.circ$node + 0

virdef <- c("vf.phage1","vs.phage1","vs.prophage1")
virlike <- c("vf.phage2","vs.phage2","vs.phage3","vs.prophage2","vs.prophage3","vs.circ")
vircols <- c(virdef,virlike)

input.meta$virdef <- (rowSums(input.meta[,virdef])>0) + 0
input.meta$virlike <- (rowSums(input.meta[,virlike])>0) + 0


### FIND RNA CONTIGS

bacrib.in <- cbind(domain="Bacteria",read.table(bacrrnafile,header=F,sep="\t",quote="",comment="#",fill=T, stringsAsFactors=F))
arcrib.in <- cbind(domain="Archaea",read.table(arcrrnafile,header=F,sep="\t",quote="",comment="#",fill=T, stringsAsFactors=F))
eukrib.in <- cbind(domain="Eukarya",read.table(eukrrnafile,header=F,sep="\t",quote="",comment="#",fill=T, stringsAsFactors=F))
mitorib.in <- cbind(domain="Mitochondria",read.table(mitorrnafile,header=F,sep="\t",quote="",comment="#",fill=T, stringsAsFactors=F))
rib.in <- rbind(bacrib.in,arcrib.in,eukrib.in,mitorib.in)
colnames(rib.in) <- c("domain","orf","version","mol","bp.start","bp.end","evalue","strand","X","description")
rib.in$length <- rib.in$bp.end - rib.in$bp.start
rib.in$ribo16S <- grepl("16S",rib.in$description) + 0
rib.in$ribo23S <- grepl("23S",rib.in$description) + 0
rib.in$ribo5S <- grepl("5S",rib.in$description) + 0
rib.in <- rib.in[order(rib.in$evalue,rib.in$orf,-rib.in$length),]
rib.in <- rib.in[!duplicated(rib.in$orf),]

input.meta$rrna.mol <- NA
input.meta[input.meta$orf %in% rib.in$orf[rib.in$ribo16S > 0],"rrna.mol"] <- "16S"
input.meta[input.meta$orf %in% rib.in$orf[rib.in$ribo23S > 0],"rrna.mol"] <- "23S"
input.meta[input.meta$orf %in% rib.in$orf[rib.in$ribo5S > 0],"rrna.mol"] <- "5S"
input.meta[grep(patt="23S ribosomal RNA",input.meta$description),"rrna.mol"] <- "23S"
input.meta[grep(patt="16S ribosomal RNA",input.meta$description),"rrna.mol"] <- "16S"
input.meta[grep(patt="5S ribosomal RNA",input.meta$description),"rrna.mol"] <- "5S"

input.meta$ribo16S <- (input.meta$rrna.mol=="16S") + 0
input.meta$ribo23S <- (input.meta$rrna.mol=="23S") + 0
input.meta$ribo5S <- (input.meta$rrna.mol=="5S") + 0
input.meta$ribo16S[is.na(input.meta$ribo16S)] <- 0
input.meta$ribo23S[is.na(input.meta$ribo23S)] <- 0
input.meta$ribo5S[is.na(input.meta$ribo5S)] <- 0

input.meta$rrna <- !is.na(input.meta$rrna.mol) + 0

trna.orfs <- unique(input.meta$orf[grep(patt="tRNA-...\\(...\\)",input.meta$description)])
input.meta$trna <- input.meta$orf %in% trna.orfs + 0


### FIND HYPOTHETICAL CONTIGS

hypo.orfs <- unique(input.meta$orf[grep(patt="hypothetical",input.meta$description,invert=F)])
nonhypo.orfs <- unique(input.meta$orf[grep(patt="hypothetical",input.meta$description,invert=T)])

input.meta$hypo <- input.meta$orf %in% hypo.orfs + 0
input.meta$nonhypo <- input.meta$orf %in% nonhypo.orfs + 0

### ADD KAIJU TAXONOMY GROUPS
input.meta$tax.bacteria <- grepl("Bacteria",input.meta$taxonomy)+0
input.meta$tax.eukarya <- grepl("Eukaryota",input.meta$taxonomy)+0
input.meta$tax.archaea <- grepl("Archaea",input.meta$taxonomy)+0
input.meta$tax.viruses <- grepl("Viruses",input.meta$taxonomy)+0



### GC CONTENT

library(seqinr)

orf.seqs <- read.FASTA(file = scaffile)
orf.gc <- sapply(1:length(orf.seqs),function(x) GC.content(orf.seqs[x]))
names(orf.gc) <- names(orf.seqs)
input.meta$gc.content <- NA
input.meta[names(orf.gc),"gc.content"] <- orf.gc
rm(orf.seqs)
rm(orf.gc)

### AGGREGATE 
toagg <- c("marker","vf.phage1","vf.phage2","vs.phage1","vs.phage2","vs.phage3",
"vs.prophage1","vs.prophage2","vs.prophage3","virdef","virlike","hypo","nonhypo","trna",
"rrna","ribo5S","ribo16S","ribo23S","gc.content","contigcov","contigcount","contigsize","orfcount","tax.bacteria","tax.eukarya",
"tax.archaea","tax.viruses","vs.circ","tigr.mark","pfam.mark","bact.mark","arch.mark", "euk.mark")

agg.nodes <- aggregate(input.meta[,toagg], by=list(input.meta$node), FUN=sum, na.rm=TRUE)
rownames(agg.nodes) <- agg.nodes$Group.1
agg.nodes <- agg.nodes[,-1]
agg.nodes <- agg.nodes[rownames(cov.tsne),] #sort same as cov.tsne

#order by makers, then size
agg.sort <- agg.nodes[order(agg.nodes$marker,agg.nodes$contigsize,decreasing=T),]
marker.nodes.sort <- match(rownames(agg.sort), rownames(cov.tsne))

#nodes to exclude because viral or viral-like, but don't exclude contigs with markers
#row numbers in cov.tsne
vir.exclude <- unname(which(rowSums(agg.nodes[,vircols])>0))
vir.exclude <- setdiff(vir.exclude, which(agg.nodes$marker>0)) #don't exclude contigs with markers

#set taxonomy likelihoods
tax.mat <- as.matrix(agg.nodes[,c("tax.bacteria","tax.archaea","tax.eukarya","tax.viruses")])
tax.prop <- prop.table(tax.mat,1)
tax.prop[is.nan(tax.prop)] <- 0
tax.prop <- as.data.frame(tax.prop)

source("bin/ambwr-s40-binning.r",echo=T)
