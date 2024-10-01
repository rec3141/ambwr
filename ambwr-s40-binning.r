### START BINNING GENOMES

# idea for naive bayes classifier: for each genome in database, find agg.bin data and train to recognize bact, arch, euk, vir, megavir, plasmid, etc
# remove viral nodes until end
# calculate # of marker genes per node
# descending sort remaining nodes by # of marker genes
# for each node
#   for each ascending sorted cutoff less than 1unit
#     calculate completion and redundancy
#     choose cutoff where completion is high and redundancy is low
#	  pick max cutoff by finding highest cutoff where redundancy < 1.2
#	  find lowest cutoff at lowest cutoff at which redundancy

input.meta$bin <- NA
input.meta$binnumber <- NA
input.meta$bincutoff <- NA

#function to calculate N50 given a set of contig lengths
N50 <- function(x) {
	x <- rev(sort(x))
	x[which.max(cumsum(x) > sum(x)/2)]
}


    x <- cov.tsne      # bivariate normal n=100
	y1 <- rep.int(x[,1], round(agg.nodes$contigcov/min(agg.nodes$contigcov))) #weight by cov
	y2 <- rep.int(x[,2], round(agg.nodes$contigcov/min(agg.nodes$contigcov))) #weight by cov
# 	y1 <- rep.int(x[,1], round(agg.nodes$contigsize/min(agg.nodes$contigsize))) #weight by length
# 	y2 <- rep.int(x[,2], round(agg.nodes$contigsize/min(agg.nodes$contigsize))) #weight by length
	x <- cbind(y1,y2)
		
     ab <- matrix( c(-80,-80,80,80), 2, 2)       # interval [-5,5) x [-5,5)
     nbin <- c( 2000, 2000)                      # 400 bins
     bins <- bin2(x, ab, nbin)               # bin counts,ab,nskip


    m <- c(20,20)
     f <- ash2(bins,m)
     image(f$x,f$y,f$z)
#      contour(f$x,f$y,f$z,add=TRUE)
#      points(x,pch='.',cex=0.3)
     dev.off()
     
     
# this loops file in the input.meta$bin column with the appropriate bin
# exclude viruses in loop1
# only viruses for loop2
# rescue small bins in loop3
# break contaminated bins in loop4? determine if multiple genomes present using GC + coverage

maxredund <- 0.05 #max redundancy before starting new bin
redosize <- 100 #minimum orfs to avoid redo
max.cutoff <- 1 #for loop1
min.cutoff <- 0.3 #for bins with markers
binreport <- list()

for(loop in 1:4) {

	# main loop; maximize sizes of things with markers first, then from largest to smallest
	if(loop==1) {
		binnumber <- 0
		seen.contigs <- vir.exclude
	# give small bins another chance; counting orfs
	# break up mixed bins
	} else if(loop==2) {
		max.cutoff <- 0.5
		rle.bin <- as.data.frame(unclass(rle(sort(input.meta$binnumber))))
		bin.small <- rle.bin$values[which(rle.bin$lengths < redosize)]
		input.meta$bin[input.meta$binnumber %in% bin.small] <- NA
		input.meta$bincutoff[input.meta$binnumber %in% bin.small] <- NA
		input.meta$binnumber[input.meta$binnumber %in% bin.small] <- NA #do last		
		seen.contigs <- which(rownames(cov.tsne) %in% input.meta$node[!is.na(input.meta$bin)])
	# do viruses last with small cutoff
	} else if(loop==3) {
		seen.contigs <- setdiff(seen.contigs,vir.exclude)
		max.cutoff <- 0.5
	}
	
	dotcount <- 0
	#i is a row in the tsne contig matrix
	for(i in 1:length(marker.nodes.sort)) {		

		i.row <- marker.nodes.sort[i]

		#skip viruses and already binned contigs
		if(i.row %in% seen.contigs) next

		if(dotcount < 80) {	cat(".") } else { cat("\r"); dotcount <- 0 }
		dotcount <- dotcount + 1
		
		temp.contigs <- rep(0,nrow(input.meta))
		keep.contigs <- rep(0,nrow(input.meta))

		#get cophenetic distances
		dists <- NULL
		dists$xs <- tsne.tmat[i.row,]
		dists <- as.data.frame(dists)
		i.node <- rownames(dists)[i.row]

		cutoffs <- sort(unique(dists$xs)) # choose all cutoffs in hierarchical clustering
		cutoffs <- cutoffs[cutoffs <= max.cutoff] #enforce size limits

		for(j in 1:length(cutoffs)) {

			cutoff <- cutoffs[j]
		
			if(j==1) { dist.pick <- which(dists$xs <= cutoff) }
			if(j>1) { dist.pick <- which(dists$xs <= cutoff & dists$xs > cutoffs[j-1]) }

			dist.pick <- setdiff(dist.pick, c(seen.contigs,which(temp.contigs > 0))) #remove viral contigs from consideration in loop1; binned contigs in loop2; already in temp bin

			#temp.contigs is for each row of input.meta, cumulative over js
			temp.contigs[input.meta$node %in% rownames(dists)[dist.pick]] <- j

			#which domain are we working with
			i.domain <- as.data.frame(prop.table(t(colSums(tax.mat[input.meta$node[temp.contigs>0],,drop=F])),1))
			i.domain[is.nan(as.numeric(i.domain))] <- 0
			#has to be 70% one domain to use that marker set
			if(max(i.domain, na.rm=T) > 0.7) { i.max <- which.max(i.domain) } else { i.max <- 0 }

			if(length(i.max)==0) { markers <- input.meta$marker; refmarkcount <- nrow(markerinfo)
			} else if(i.max==1) { markers <- input.meta$bact.mark; refmarkcount <- bactmarkcount
			} else if(i.max==2) { markers <- input.meta$arch.mark; refmarkcount <- archmarkcount
			} else if(i.max==3) { markers <- input.meta$euk.mark; refmarkcount <- eukmarkcount
			} else { markers <- input.meta$marker; refmarkcount <- nrow(markerinfo) }

			# making this more conservative by changing 'markers' to 'input.meta$marker'
			# that way it won't include other domain's markers as well
			# ah but that won't work if there are legit hmm hits to other marker sets that just aren't domain-specific markers

			mark.hmm <- input.meta$mark.hmm[temp.contigs > 0 & markers > 0]
			completion <- length(unique(mark.hmm))/refmarkcount
			redund <- sum(rle(sort(mark.hmm))$lengths > 1)/refmarkcount
			
			#if bin has marker genes
			if(length(mark.hmm) > 0) {

				#add to kept if founder bin
				if(j==1) {
					cat("F")
					keep.contigs <- temp.contigs

				#add to kept if not too contaminated, or if still at small cutoff
				} else if(redund < maxredund | cutoff < min.cutoff) {

					cat("+")
					keep.contigs <- temp.contigs
					
				#if has marks this means we've exceeded the redundancy, in which case bin is done
				} else {
					cat("R")
					break
				}
				
			#if bin doesn't have marker genes
			} else {

				#add to kept if below cutoffs
				if(cutoff < max.cutoff) {
					cat("+")
					keep.contigs <- temp.contigs

				#if above cutoff, the bin is done
				} else {
					cat("C")
					break
				}
			}
		}

		bin.contigs <- which(rownames(input.meta) %in% unique(input.meta$node[keep.contigs > 0]))
		contigcount <- length(bin.contigs)
		orfcount <- sum(keep.contigs > 0)
		i.domain <- as.data.frame(prop.table(t(colSums(tax.mat[input.meta$node[bin.contigs],,drop=F])),1))
		i.domain[is.nan(as.numeric(i.domain))] <- 0
		names(i.domain)[which.max(i.domain)] <- toupper(names(i.domain)[which.max(i.domain)])
		i.taxa <- gsub("tax.","",ignore.case=T,paste(paste(names(i.domain),sprintf('%.2f',i.domain), sep=': '),collapse="  "))

		if(max(i.domain, na.rm=T) > 0.7) { i.max <- which.max(i.domain) } else { i.max <- 0 }

		if(length(i.max)==0) { markers <- input.meta$marker; refmarkcount <- nrow(markerinfo)
		} else if(i.max==1) { markers <- input.meta$bact.mark; refmarkcount <- bactmarkcount
		} else if(i.max==2) { markers <- input.meta$arch.mark; refmarkcount <- archmarkcount
		} else if(i.max==3) { markers <- input.meta$euk.mark; refmarkcount <- eukmarkcount
		} else { markers <- input.meta$marker; refmarkcount <- nrow(markerinfo) }

		mark.hmm <- input.meta$mark.hmm[keep.contigs > 0 & markers > 0]
		completion <- length(unique(mark.hmm))/refmarkcount
		redund <- sum(rle(sort(mark.hmm))$lengths > 1)/refmarkcount

		#if bin has marker genes but they aren't it's own kind, drop the whole bin in first loop
		all.hmm <- input.meta$mark.hmm[keep.contigs > 0 & input.meta$marker > 0]
# drop it for now (repeats for every contig in this bin, not efficient)
#		if(length(all.hmm) > 0 & completion == 0) { cat("0"); next }

		# otherwise process bin
		input.meta$bin[keep.contigs > 0] <- i.node
		input.meta$binnumber[keep.contigs > 0] <- binnumber
		input.meta$bincutoff[keep.contigs > 0] <- cutoff


		if(is.nan(completion) | is.nan(redund)) { binqual <- ""
		} else if(completion > 0.9 & redund < 0.05) { binqual <- "high"
		} else if (completion >= 0.5 & redund < 0.1) { binqual <- "medium"
		} else if (completion < 0.5 & completion > 0 & redund < 0.1) { binqual <- "low"
		} else if (redund > 0.1) { binqual <- "contam"
		} else { binqual <- "" }
		
		cat("\r")

		bincat <- paste(
		sprintf("%4.1f", round(sum(!is.na(input.meta$bin))/nrow(input.meta)*100, 3)),"%",
		"loop:", loop,
		"bin:", sprintf("%4d", binnumber), 
		"cutoff:", sprintf("%3.2f", round(cutoff,2)), 
		"contigs:",sprintf("%5d", contigcount),
		"orfs:",sprintf("%6d", orfcount),
		"|",i.taxa,"|",
		"N50:", sprintf("%6d",N50(input.meta$contiglen[bin.contigs])),
		"compl:", sprintf("%4.1f", 100*round(completion,3)),"%",
		"redun:",sprintf("%4.2f", 100*round(redund,3)),"%",
		binqual,
		"\n")
		
		cat(bincat)
		binreport[i] <- bincat
		
		#save seen contigs
		seen.contigs <- union(seen.contigs, which(rownames(cov.tsne) %in% input.meta$node[keep.contigs > 0]))
		binnumber <- binnumber+1

	}


}
	# BIN TABLE
	agg.nodes$bin <- input.meta[rownames(agg.nodes),"bin"]
	agg.nodes$binnumber <- input.meta[rownames(agg.nodes),"binnumber"]
	agg.nodes$binnumber[is.na(agg.nodes$bin)] <- 0
	agg.nodes$bin[is.na(agg.nodes$bin)] <- "unbinned"

	input.meta$shortbin <- paste0("bin",input.meta$binnumber,"_",do.call(rbind,strsplit(rownames(input.meta),split="_"))[,2])

	agg.nodes$shortbin <- paste0("bin",agg.nodes$binnumber,"_",do.call(rbind,strsplit(rownames(agg.nodes),split="_"))[,2])
	agg.bins$binname <- paste0("bin",agg.bins$binnumber)

	agg.nodes$node <- rownames(agg.nodes)
	rownames(agg.nodes) <- agg.nodes$shortbin
	agg.nodes$shortbin <- NULL

	input.meta$binnumber[is.na(input.meta$bin)] <- 0
	input.meta$bin[is.na(input.meta$bin)] <- "unbinned"

	agg.bins <- aggregate(input.meta[,toagg], by=list(input.meta$bin), FUN=sum, na.rm=TRUE)
	rownames(agg.bins) <- agg.bins$Group.1
	agg.bins <- agg.bins[,-1]
	agg.bins$meancov <- agg.bins$contigcov/agg.bins$contigcount
	agg.bins$contigcov <- NULL
	agg.bins$binnumber <- agg.nodes[rownames(agg.bins),"binnumber"]
	agg.bins <- agg.bins[order(agg.bins$binnumber),]

	# COMPLETENESS
	# definitions from MIMAG paper: https://www.nature.com/articles/nbt.3893/tables/1
	# Completion: ratio of observed single-copy marker genes to total single-copy marker genes in chosen marker gene set.
	# Contamination: ratio of observed single-copy marker genes in ≥2 copies to total single-copy marker genes in chosen marker gene set.
	# high: completion > 90%, contamination < 5%
	# medium: completion > 50%, contamination < 10%
	# low: completion < 50%, contamination < 10%

	# bacteria markers are also universal markers so other domains have lots of hits
	agg.bins$marktax.arch <- (agg.bins$arch.mark/agg.bins$bact.mark > 1) + 0
	agg.bins$marktax.euk <- (agg.bins$euk.mark/agg.bins$bact.mark > 1) + 0
	agg.bins$marktax.bact <- ((agg.bins$bact.mark > 0) & (!agg.bins$marktax.arch) & (!agg.bins$marktax.euk)) + 0
	agg.bins$marktax.arch[is.na(agg.bins$marktax.arch)] <- 0
	agg.bins$marktax.bact[is.na(agg.bins$marktax.bact)] <- 0
	agg.bins$marktax.euk[is.na(agg.bins$marktax.euk)] <- 0

	agg.bins$completion <- 0
	agg.bins$redundancy <- 0

	for(i in 1:nrow(agg.bins)) {
		markers <- NULL
		refmarkcount <- NULL
		if(agg.bins$marktax.arch[i] > 0) { markers <- input.meta$arch.mark; refmarkcount <- archmarkcount
		} else if(agg.bins$marktax.bact[i] > 0) { markers <- input.meta$bact.mark; refmarkcount <- bactmarkcount
		} else if(agg.bins$marktax.euk[i] > 0) { markers <- input.meta$euk.mark; refmarkcount <- eukmarkcount
		} else { next }

		#find all contigs in input.meta
		mark.hmm <- input.meta$mark.hmm[input.meta$bin %in% rownames(agg.bins)[i] & markers > 0]
		agg.bins$completion[i] <- length(unique(mark.hmm))/refmarkcount
		agg.bins$redundancy[i] <- sum(rle(sort(mark.hmm))$lengths>1)/refmarkcount

	}


	### BLACK CONTIG PLOT
	pdf(file=paste0("blackplots-",loop,".pdf"),width=24,height=24)

	par(bg = 'black')

	#plot hypotheticals
	plot(cov.tsne, cex=(sqrt(agg.nodes$hypo)+0)*0.2,pch=21,bg="grey",col=NULL)

	#add nonhypotheticals
	points(cov.tsne, cex=(sqrt(agg.nodes$nonhypo)+0)*0.2, pch=21, bg="white",col=NULL)

	#add markers
	points(cov.tsne, cex=sqrt(agg.nodes$marker)*0.4, pch=21, bg="blue",col=NULL)

	#add likely viruses
	vir.cex <- agg.nodes$virlike * 0.3
	points(cov.tsne, cex=vir.cex, pch=21, bg="yellow",col=NULL)

	#add definite viruses
	vir.cex <- agg.nodes$virdef * 0.4
	points(cov.tsne, cex=vir.cex, pch=21, bg="red",col=NULL)

	#add ribosomal RNAs
	points(cov.tsne[which(rownames(cov.tsne) %in% input.meta$node[grep(patt="23S ribosomal RNA",input.meta$description)]),], cex=1, pch=21, bg="green4",col=NULL)
	points(cov.tsne[which(rownames(cov.tsne) %in% input.meta$node[grep(patt="16S ribosomal RNA",input.meta$description)]),], cex=0.8, pch=21, bg="green3",col=NULL)
	points(cov.tsne[which(rownames(cov.tsne) %in% input.meta$node[grep(patt="5S ribosomal RNA",input.meta$description)]),], cex=0.5, pch=21, bg="green",col=NULL)
	points(cov.tsne[which(rownames(cov.tsne) %in% input.meta$node[grep(patt="tRNA-...\\(...\\)",input.meta$description)]),], cex=0.4, pch=21, bg="purple",col=NULL)


	#BLACK BIN plot
	par(bg="black")

	alpha <- rescale(log10(agg.nodes$contigcov/agg.nodes$contigcount),to=c(0,1))
	# Rgb <- sample(1:255,255,replace=T)/255 #this is why I was getting low diversity colors?
	Rgb <- sample(1:255,nrow(agg.bins),replace=T)/255
	rGb <- rescale(cov.tsne[,1],to=c(0,1))
	rgB <- rescale(cov.tsne[,2],to=c(0,1))
	rcols <- rgb(Rgb,rGb,rgB,alpha)
	names(rcols) <- rownames(cov.tsne)

	# plot(cov.tsne,pch=21,col=rcols[agg.nodes$bin],bg=NULL,cex=0.4)
	# text(cov.tsne, labels=agg.nodes$binnumber,cex=0.2,col="white")

	plot(cov.tsne,pch=21,bg=rcols[agg.nodes$bin],col=NULL,cex=sqrt(agg.nodes$orfcount)/10)
	text(cov.tsne[which(!is.na(match(rownames(cov.tsne),rownames(agg.bins)))),], labels=agg.bins$binnumber,cex=0.2,col="white")


	#BLACK ARCHAEA plot
	par(bg="black")

	alpha <- agg.bins$marktax.arch
	alpha <- with(agg.bins, tax.archaea/(tax.bacteria+tax.archaea+tax.eukarya+tax.viruses+1))
	Rgb <- sample(1:255,nrow(agg.bins),replace=T)/255
	rGb <- rescale(cov.tsne[rownames(agg.bins),1],to=c(0,1))
	rgB <- rescale(cov.tsne[rownames(agg.bins),2],to=c(0,1))
	rcols <- rgb(Rgb,rGb,rgB,alpha)
	names(rcols) <- rownames(agg.bins)

	# plot(cov.tsne,pch=21,col=rcols[agg.nodes$bin],bg=NULL,cex=0.4)
	# text(cov.tsne, labels=agg.nodes$binnumber,cex=0.2,col="white")

	plot(cov.tsne,pch=21,bg=rcols[agg.nodes$bin],col=NULL,cex=sqrt(agg.nodes$orfcount)/5, main="archaea", col.main="white", cex.main=5)
	text(cov.tsne[which(!is.na(match(rownames(cov.tsne),rownames(agg.bins)))),], labels=agg.bins$binnumber,cex=0.2,col="white")



	#BLACK BACTERIA plot
	par(bg="black")

	alpha <- agg.bins$marktax.bact
	alpha <- with(agg.bins, tax.bacteria/(tax.bacteria+tax.archaea+tax.eukarya+tax.viruses+1))
	Rgb <- sample(1:255,nrow(agg.bins),replace=T)/255
	rGb <- rescale(cov.tsne[rownames(agg.bins),1],to=c(0,1))
	rgB <- rescale(cov.tsne[rownames(agg.bins),2],to=c(0,1))
	rcols <- rgb(Rgb,rGb,rgB,alpha)
	names(rcols) <- rownames(agg.bins)

	# plot(cov.tsne,pch=21,col=rcols[agg.nodes$bin],bg=NULL,cex=0.4)
	# text(cov.tsne, labels=agg.nodes$binnumber,cex=0.2,col="white")

	plot(cov.tsne,pch=21,bg=rcols[agg.nodes$bin],col=NULL,cex=sqrt(agg.nodes$orfcount)/10, main="bacteria", col.main="white", cex.main=5)
	text(cov.tsne[which(!is.na(match(rownames(cov.tsne),rownames(agg.bins)))),], labels=agg.bins$binnumber,cex=0.2,col="white")


	#BLACK EUKARYOTA plot
	par(bg="black")

	alpha <- agg.bins$marktax.euk
	alpha <- with(agg.bins, tax.eukarya/(tax.bacteria+tax.archaea+tax.eukarya+tax.viruses+1))
	Rgb <- sample(1:255,nrow(agg.bins),replace=T)/255
	rGb <- rescale(cov.tsne[rownames(agg.bins),1],to=c(0,1))
	rgB <- rescale(cov.tsne[rownames(agg.bins),2],to=c(0,1))
	rcols <- rgb(Rgb,rGb,rgB,alpha)
	names(rcols) <- rownames(agg.bins)

	# plot(cov.tsne,pch=21,col=rcols[agg.nodes$bin],bg=NULL,cex=0.4)
	# text(cov.tsne, labels=agg.nodes$binnumber,cex=0.2,col="white")

	plot(cov.tsne,pch=21,bg=rcols[agg.nodes$bin],col=NULL,cex=sqrt(agg.nodes$orfcount)/5, main="eukarya", col.main="white", cex.main=5)
	text(cov.tsne[which(!is.na(match(rownames(cov.tsne),rownames(agg.bins)))),], labels=agg.bins$binnumber,cex=0.2,col="white")



	#BLACK VIRUS plot
	par(bg="black")

	alpha <- (agg.bins$virdef+agg.bins$virlike > 10)+0
	alpha <- with(agg.bins, tax.viruses/(tax.bacteria+tax.archaea+tax.eukarya+tax.viruses+1))
	Rgb <- sample(1:255,nrow(agg.bins),replace=T)/255
	rGb <- rescale(cov.tsne[rownames(agg.bins),1],to=c(0,1))
	rgB <- rescale(cov.tsne[rownames(agg.bins),2],to=c(0,1))
	rcols <- rgb(Rgb,rGb,rgB,alpha)
	names(rcols) <- rownames(agg.bins)

	# plot(cov.tsne,pch=21,col=rcols[agg.nodes$bin],bg=NULL,cex=0.4)
	# text(cov.tsne, labels=agg.nodes$binnumber,cex=0.2,col="white")

	plot(cov.tsne,pch=21,bg=rcols[agg.nodes$bin],col=NULL,cex=sqrt(agg.nodes$orfcount)/5, main="virus", col.main="white", cex.main=5)
	text(cov.tsne[which(!is.na(match(rownames(cov.tsne),rownames(agg.bins)))),], labels=agg.bins$binnumber,cex=0.2,col="white")



	#BLACK UNKNOWN plot
	par(bg="black")

	alpha <- (agg.bins$virdef+agg.bins$virlike > 10)+0
	alpha <- with(agg.bins, tax.viruses/(tax.bacteria+tax.archaea+tax.eukarya+tax.viruses+1))
	Rgb <- sample(1:255,nrow(agg.bins),replace=T)/255
	rGb <- rescale(cov.tsne[rownames(agg.bins),1],to=c(0,1))
	rgB <- rescale(cov.tsne[rownames(agg.bins),2],to=c(0,1))
	rcols <- rgb(Rgb,rGb,rgB,alpha)
	names(rcols) <- rownames(agg.bins)

	# plot(cov.tsne,pch=21,col=rcols[agg.nodes$bin],bg=NULL,cex=0.4)
	# text(cov.tsne, labels=agg.nodes$binnumber,cex=0.2,col="white")

	plot(cov.tsne,pch=21,bg=rcols[agg.nodes$bin],col=NULL,cex=sqrt(agg.nodes$orfcount)/5, main="virus", col.main="white", cex.main=5)
	text(cov.tsne[which(!is.na(match(rownames(cov.tsne),rownames(agg.bins)))),], labels=agg.bins$binnumber,cex=0.2,col="white")



	#BLACK BIN COMPLETENESS/REDUNDANCY PLOT
	#red=contam, green=good
	#alpha=completion

	alpha <- rescale(agg.bins$completion,to=c(0.1,1))
	Rgb <- rescale(agg.bins$redundancy,to=c(0,1))
	rGb <- rescale(agg.bins$redundancy,to=c(1,0))
	rgB <- rep(0,nrow(cov.tsne))
	rcols <- rgb(Rgb,rGb,rgB,alpha)
	names(rcols) <- rownames(cov.tsne)

	par(bg="black")
	plot(cov.tsne, pch=21, bg=rcols[agg.nodes$bin],col=NULL,cex=log10(agg.bins$contigsize)/5)

	#red=incomplete, green=good
	#alpha=redundancy

	alpha <- rescale(agg.bins$redundancy,to=c(1,0))
	Rgb <- rescale(agg.bins$completion,to=c(1,0))
	rGb <- rescale(agg.bins$completion,to=c(0,1))
	rgB <- rep(0,nrow(cov.tsne))
	rcols <- rgb(Rgb,rGb,rgB,alpha)
	names(rcols) <- rownames(cov.tsne)

	par(bg="black")
	plot(cov.tsne, pch=21, bg=rcols[agg.nodes$bin],col=NULL,cex=log10(agg.bins$contigsize)/5)

	#red=low cov, green=high cov
	#alpha=completion

	alpha <- rescale(agg.bins$completion,to=c(0.25,0.75))
	Rgb <- rescale(sqrt(sqrt(agg.bins$meancov)),to=c(1,0))
	rGb <- rescale(sqrt(sqrt(agg.bins$meancov)),to=c(0,1))
	rgB <- rep(0,nrow(cov.tsne))
	rcols <- rgb(Rgb,rGb,rgB,alpha)
	names(rcols) <- rownames(cov.tsne)

	par(bg="black")
	plot(cov.tsne, pch=21, bg=rcols[agg.nodes$bin],col=NULL,cex=log10(agg.bins$contigsize)/5)

	dev.off()



