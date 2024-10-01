# write out fasta files for each bin
# these will be used for mapping and taxonomy

library(seqinr)
bindir <- file.path("annotation/bins/",prefix)
scaf.fa <- read.fasta(file = scaffile, as.string = TRUE, forceDNAtolower=F)
for(bin in rownames(agg.bins)) {
	scaf.bin <- scaf.fa[ agg.nodes$node[ agg.nodes$bin == bin] ]
	write.fasta(scaf.bin, names(scaf.bin),file.path(bindir,paste0(agg.bins[bin,"binname"],".fasta")) , open = "w", nbchar = 1e7, as.string = TRUE)
	print(bin)
}
rm(scaf.fa)

# taxonomy classification from GTDB
# run bash scripts
gtdb.bacfile <- "annotation/bins/scaffolds-gt2000-2/gtdbtk/scaffolds-gt2000.bac120.summary.tsv"
gtdb.arcfile <- "annotation/bins/scaffolds-gt2000-2/gtdbtk/scaffolds-gt2000.ar122.summary.tsv"

gbac <- read.table(gtdb.bacfile,sep="\t",head=T,stringsAsFactors=F,row.names=1)
garc <- read.table(gtdb.arcfile,sep="\t",head=T,stringsAsFactors=F,row.names=1)
gtdb <- rbind(gbac,garc)
agg.bins[match(rownames(gtdb),agg.bins$binname),"gtdb"] <- gtdb$classification

# taxonomy classification from 
# run bash scripts
sk.file <- "annotation/bins/scaffolds-gt2000-2/sketchout.tsv"
sktab <- read.table(sk.file,sep="\t",head=F,stringsAsFactors=F,fill=T)
sktab <- sktab[,c(1,3,12)]
rownames(sktab) <- gsub(".fasta","",basename(sktab[,1]))

agg.bins[match(rownames(sktab),agg.bins$binname),"sketch1"] <- sktab[,2]
agg.bins[match(rownames(sktab),agg.bins$binname),"sketch2"] <- sktab[,3]

# concatenates them all together
# aggregate(input.meta[,"taxonomy"], by=list(input.meta$binnumber), FUN=toString, na.rm=TRUE)
kjtab <- aggregate(input.meta$taxonomy,by=list(input.meta$binnumber), 
	FUN=function(x) { 
          ux <- unique(x)
          ux[ux=="unclassified;"] <- NA
          ux[which.max(tabulate(match(x, ux)))] 
})
kjtab[is.na(kjtab)] <- "unclassified;"
rownames(kjtab) <- paste0("bin",kjtab[,1])
agg.bins[,"kaiju"] <- kjtab[match(agg.bins$binname,rownames(kjtab)),2]

// use gtdb + sketch as per-bin IDs, kaiju as per-orf (# per contig)

agg.taxa <- agg.bins[,c("binname","gtdb","kaiju","sketch1","sketch2")]


# GTDB
tmp.gtdb <- do.call(rbind,strsplit(agg.taxa$gtdb,split=";"))
rownames(tmp.gtdb) <- agg.bins$binname
uniq.gtdb <- unique(as.vector(as.matrix(tmp.gtdb)))
list.gtdb <- list()

for(i in 1:ncol(tmp.gtdb)) {
	uniq.gtdb <- unique(tmp.gtdb[,i])
	for(tax in uniq.gtdb) {
		tax2 <- gsub("^.__","",tax)
		print(tax2)
		list.gtdb[[tax2]] <- c(list.gtdb[[tax2]],rownames(tmp.gtdb)[which(tmp.gtdb[,i] %in% tax)])
	}
}

#JSON output
# {"taxon1": [bin1, bin2, bin3], "taxon2": [bin4, bin5]}
list.gtdb %>% toJSON() %>% write_lines("taxa_gtdb.json")



#KAIJU
tmp.kaiju <- do.call(rbind,strsplit(agg.taxa$kaiju,split=";"))
rownames(tmp.kaiju) <- agg.bins$binname
uniq.kaiju <- unique(as.vector(as.matrix(tmp.kaiju)))
list.kaiju <- list()

for(i in 1:ncol(tmp.kaiju)) {
	uniq.kaiju <- unique(tmp.kaiju[,i])
	for(tax in uniq.kaiju) {
		tax2 <- gsub("^ ","",tax)
		print(tax2)
		list.kaiju[[tax2]] <- c(list.kaiju[[tax2]],rownames(tmp.kaiju)[which(tmp.kaiju[,i] %in% tax)])
	}
}

#JSON output
# {"taxon1": [bin1, bin2, bin3], "taxon2": [bin4, bin5]}
list.kaiju %>% toJSON() %>% write_lines("taxa_kaiju.json")




#SKETCH
sketch.in <- read.table("annotation/bins/scaffolds-gt2000-2/sketchtax.tsv",sep="\t",fill=T)
agg.taxa$sketch <- sketch.in[match(agg.taxa$binname,sketch.in$V1),"V2"]
tmp.split <- strsplit(agg.taxa$sketch,split=";")
tmp.sketch <- data.frame(matrix(,nrow=nrow(agg.taxa),ncol=7))
colnames(tmp.sketch) <- c("sk","p","c","o","f","g","s")
for(i in 1:length(tmp.split)) {
	for(j in 1:length(tmp.split[[i]])) {
# 		print(tmp.split[[i]][j])[[1]];
		if(length(tmp.split[[i]])==0) next
		tm <- strsplit(x=tmp.split[[i]][j],split=":")[[1]]
		if(any(is.na(tm))) next
		if(tm[1] %in% colnames(tmp.sketch)) tmp.sketch[i,tm[1]] <- tm[2]
	}
}
rownames(tmp.sketch) <- agg.bins$binname
tmp.sketch[is.na(tmp.sketch)] <- ""
uniq.sketch <- unique(as.vector(as.matrix(tmp.sketch)))
list.sketch <- list()

for(i in 1:ncol(tmp.sketch)) {
	uniq.sketch <- unique(tmp.sketch[,i])
	for(tax in uniq.sketch) {
		tax2 <- gsub("^ ","",tax)
		print(tax2)
		list.sketch[[tax2]] <- c(list.sketch[[tax2]],rownames(tmp.sketch)[which(tmp.sketch[,i] %in% tax)])
	}
}

# {"taxon1": [bin1, bin2, bin3], "taxon2": [bin4, bin5]}
list.sketch %>% toJSON() %>% write_lines("taxa_sketch.json")


#COMBINED JSON output

list.combined <- list()
for(i in sort(unique(c(names(list.sketch),names(list.kaiju),names(list.gtdb))))) {
	list.combined[[i]] <- unique(c(list.sketch[[i]],list.kaiju[[i]],list.gtdb[[i]]))
}
list.combined %>% toJSON() %>% write_lines("taxa_combined.json")




### NOW FOR CONTIGS

# export metadata: descriptions

uniq.desc <- sort(unique(input.meta$desc))
list.desc <- list()

for(desc in uniq.desc) {
	list.desc[[desc]] <- sort(unique(input.meta$shortbin[input.meta$description == desc]))
}

#JSON output
# {"taxon1": [bin1, bin2, bin3], "taxon2": [bin4, bin5]}
list.desc %>% toJSON() %>% write_lines("desc_combined.json")

#combined
list.combined <- list()
for(i in sort(unique(c(names(list.sketch),names(list.kaiju),names(list.gtdb))))) {
	list.combined[[i]] <- unique(c(list.sketch[[i]],list.kaiju[[i]],list.gtdb[[i]]))
}
list.combined %>% toJSON() %>% write_lines("taxa_combined.json")




#### convert circle -> ellipse -> molleweide -> spherical -> cartesian

system("module load data/GDAL/2.2.3-pic-intel-2016b-Python-2.7.12")
system("LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/import/c1/MSL464/recollins/R/x86_64-redhat-linux-gnu-library/3.5/rgdal/libs/")
source("~/sph2car.R")
library(rgdal)
library(scales)

# get reverse molleweide mapping
# https://stackoverflow.com/questions/42174346/reverse-map-projection-how-to-get-lat-lon-coordinates-from-projected-coordinate
#The x-coordinate has a range of [−2R√2, 2R√2], and the y-coordinate has a range of [−R√2, R√2].
# agg.nodes$x <- rescale(cov.tsne[,"x"],to=c(7,346))
# agg.nodes$y <- rescale(cov.tsne[,"y"],to=c(-79,79))
# 
# theta = asin(agg.nodes$y/(rad.earth*sqrt(2)))
# latitude = asin((2*theta +sin(2*theta))/pi)
# meridian = 0
# longitude = meridian + (pi*agg.nodes$x)/(2*rad.earth*sqrt(2)*cos(theta))
# latitude = latitude/max(abs(latitude))*90
# longitude = longitude/max(abs(longitude))*180
# 
# latitude=180
# theta_0 = 2*asin((2*y)/pi)
# x	=	(2*sqrt(2))*(latitude)*cos(theta)/pi	
# y	=	sqrt(2)*sin(theta),
# write.table(agg.nodes,file="agg.nodes.tsv",sep="\t")
# agg.nodes <- read.table("agg.nodes.tsv",sep="\t")
# ref.points = data.frame(x=c(0,180),y=c(90,0))
# ref.points = SpatialPoints(ref.points, CRS('+proj=longlat'))
# ref.points = spTransform(ref.points, CRS('+proj=moll'))
# 
# sp.points = data.frame(x=agg.nodes$x * ref.points@coords[2,1]/180 , y=agg.nodes$y * ref.points@coords[1,2]/90)
# sp.points = SpatialPoints(sp.points, CRS('+proj=moll'))
# sp.points = as.data.frame(spTransform(sp.points, CRS('+proj=longlat')))
# 
# # write.table(agg.nodes,file="agg.nodes.tsv",sep="\t")
# # agg.nodes <- read.table("agg.nodes.tsv",sep="\t")

#direct inverse
#http://mathworld.wolfram.com/MollweideProjection.html
rad.earth <- 6378137
moll.r <- sqrt(2) #rad.earth*
agg.nodes$x <- rescale(cov.tsne[,"x"],to=c(-1.9*moll.r,1.9*moll.r))
agg.nodes$y <- rescale(cov.tsne[,"y"],to=c(-0.9*moll.r,0.9*moll.r))
# plot(agg.nodes$x,agg.nodes$y)

theta = asin(agg.nodes$y/sqrt(2))
latitude = asin((2*theta + sin(2*theta))/pi)
meridian = 0
longitude = meridian + (pi*agg.nodes$x)/(2*sqrt(2)*cos(theta))
# plot(latitude~longitude)
agg.nodes[,c("xpos","ypos","zpos")] <- sph2car(rescale(longitude,to=c(-179,179)),rescale(latitude,to=c(-89,89)),radius=rad.earth*1.01) # earth radius in meters


# #Hammer-Aitoff map projection
# #http://paulbourke.net/geometry/transformationprojection/
# #weirdness
# agg.nodes$x <- rescale(cov.tsne[,"x"],to=c(0.05,.95))
# agg.nodes$y <- rescale(cov.tsne[,"y"],to=c(-0.95,0.95))
# 
# z2 = 1 - agg.nodes$x^2/2 - agg.nodes$y^2/2
# longitude = 2 * atan(sqrt(2)* agg.nodes$x * sqrt(z2) / (2*z2 - 1))
# latitude = asin(sqrt(2)*agg.nodes$y*sqrt(z2))
# 
# agg.nodes[,c("xpos","ypos","zpos")] <- sph2car(rescale(longitude,to=c(-179,179)),rescale(latitude,to=c(-89,89)),radius=rad.earth*1.01) # earth radius in meters




# tsne.tmat is 54 Gb
# save.image(file=paste0(prefix,format(Sys.time(), "%Y%m%d"),".Rdata"))


#1 nodes as points
#     //[["series1", taxlevel, color, ["name", xpos, ypos, zpos, abundance, ... ],
#     // ["series2", taxlevel, color, ["name", xpos, ypos, zpos, abundance, ... ]]]

#decrease characters for faster loading in the browser
sedna.nodes <- agg.nodes
sedna.nodes$binname <- paste0("bin",sedna.nodes$binnumber)
for(i in 1:ncol(sedna.nodes)) { 
	if(is.numeric(sedna.nodes[1,i])) {
	if(all(sedna.nodes[,i]%%1==0)) { next }
	sedna.nodes[,i] <- sprintf('%.2f',sedna.nodes[,i])
	}
}




sedna.nodes$node <- NULL
sedna.nodes$binnumber <- NULL
sedna.nodes$xpos <- NULL
sedna.nodes$ypos <- NULL
sedna.nodes$zpos <- NULL

write.table(sedna.nodes,file="nodes.sedna.tsv",sep="\t",quote=F)

alpha <- rescale(log10(agg.nodes$contigcov/agg.nodes$contigcount),to=c(0,1))
Rgb <- sample(1:255, nrow(agg.bins), replace=T)/255
rGb <- rescale(cov.tsne[,1],to=c(0.2,1))
rgB <- rescale(cov.tsne[,2],to=c(0.2,1))
rcols <- rgb(Rgb,rGb,rgB,alpha)
names(rcols) <- rownames(cov.tsne)


sedna.bins <- agg.bins
rownames(sedna.bins) <- sedna.bins$binname
for(i in 1:ncol(sedna.bins)) { 
	if(is.numeric(sedna.bins[1,i])) {
	if(all(sedna.bins[,i]%%1==0)) { next }
	sedna.bins[,i] <- sprintf('%.2f',sedna.bins[,i])
	}
}
sedna.bins$binname <- NULL

sedna.bins[,c("gtdb","sketch1","sketch2","kaiju")] <- agg.bins[,c("gtdb","sketch1","sketch2","kaiju")]
write.table(sedna.bins,file="bins.sedna.tsv",sep="\t",quote=F)


# node.tab <- paste0('contig_points.add( {position : Cesium.Cartesian3(',agg.nodes$x,",",agg.nodes$y,",",agg.nodes$z,'), color : Cesium.Color.BLUE, pixelSize : 6, id : "',rownames(agg.nodes),'", translucencyByDistance : new Cesium.NearFarScalar(4, 1.0, 20, 0.0)});')
node_tab <- data.frame(binname=paste0("bin",agg.nodes$binnumber),agg.nodes$xpos,agg.nodes$ypos,agg.nodes$zpos,agg.nodes$contigcov/agg.nodes$contigcount,agg.nodes$contigsize)
rownames(node_tab) <- rownames(agg.nodes)
colnames(node_tab) <- c("binname","xpos","ypos","zpos","cov","size")

node_tab$cols <- rcols[agg.nodes$bin]




#for coverage
node_points <- NULL
for(mybin in 1:nrow(agg.bins)) {
# 	print(mybin)
	pick.rows <- which(node_tab$binname == agg.bins$binname[mybin])
	node_row <- paste0('"',rownames(node_tab)[pick.rows], '"', ",",sprintf('%.0f',node_tab$xpos[pick.rows]), ",",sprintf('%.0f',node_tab$ypos[pick.rows]), ",",sprintf('%.0f',node_tab$zpos[pick.rows]), ",",sprintf('%.2f',node_tab$cov[pick.rows]))
	node_row <- paste(node_row,collapse=",")
	node_row <- paste0('["',agg.bins$binname[mybin],'", "strain", "',node_tab$cols[pick.rows[1]],'", [',node_row,']]')
	node_points <- c(node_points,node_row)
}

node_points <- paste(node_points,collapse=",")

node_points <- paste0('[',node_points,']')
write.table(file="node_points_cov.js",node_points,sep="",eol="\n",quote=F,row.names=F,col.names=F)

#for coverage
node_points <- NULL
for(mybin in 1:nrow(agg.bins)) {
# 	print(mybin)
	pick.rows <- which(node_tab$binname == agg.bins$binname[mybin])
	node_row <- paste0('"',rownames(node_tab)[pick.rows], '"', ",",sprintf('%.0f',node_tab$xpos[pick.rows]), ",",sprintf('%.0f',node_tab$ypos[pick.rows]), ",",sprintf('%.0f',node_tab$zpos[pick.rows]), ",",sprintf('%.2f',node_tab$size[pick.rows]))
	node_row <- paste(node_row,collapse=",")
	node_row <- paste0('["',agg.bins$binname[mybin],'", "strain", "',node_tab$cols[pick.rows[1]],'", [',node_row,']]')
	node_points <- c(node_points,node_row)
}

node_points <- paste(node_points,collapse=",")

node_points <- paste0('[',node_points,']')
write.table(file="node_points_size.js",node_points,sep="",eol="\n",quote=F,row.names=F,col.names=F)


#2 bins as labels
#3 scaffolds as labels


# need to run GTDB on intermediary bins

# for contigs of re-assembled genomes shorter than 2000 bp, draw coverage from gaussian of good contigs rather than use rela data
# distance matrix between two bins as setdiff rather than union so that bins that complement each other cluster together
# then plot tsne single linkage distance vs complementation distance
