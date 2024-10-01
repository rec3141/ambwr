# wget http:#www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/tax_slv_ssu_123.txt
# cut -f3 tax_slv_ssu_123.txt | sort | uniq -c | sort -rn

# 	8275 genus
# 	  26 subfamily
# 	1412 family
# 	  10 superfamily
# 	  21 suborder
# 	1138 order
# 	  16 superorder
# 	   6 infraclass
# 	  58 subclass
# 	 541 class
# 	   1 superclass
# 	   1 infraphylum
# 	  33 subphylum
# 	 239 phylum
# 	   5 superphylum
# 	   1 infrakingdom
# 	   4 subkingdom
# 	  15 kingdom
# 	   2 superkingdom
# 	   4 major_clade
# 	   3 domain

# grep -E 'Eukaryota' tax_slv_ssu_123.txt | grep -E '\sdomain|\skingdom|\sphylum|\sclass|\sorder|\sfamily|\sgenus' > silva.nr_v123.euk.taxmap
# grep -E 'Eukaryota' tax_slv_ssu_123.txt > tax_slv_ssu_123_euk.txt

#!/usr/bin/R

map.in <- read.table("tax_slv_ssu_123_euk.txt",header=F,sep="\t",stringsAsFactors=F)
map.in <- map.in[,c(1,3)]
colnames(map.in) <- c("taxlabel","taxlevel")

tax.mat <- matrix(nrow=nrow(map.in),ncol=length(unique(map.in$taxlevel)))
colnames(tax.mat) <- c("domain","major_clade","superkingdom","kingdom","subkingdom","infrakingdom","superphylum","phylum","subphylum","infraphylum","superclass","class","subclass","infraclass","superorder","order","suborder","superfamily","family","subfamily","genus")
#colnames(tax.mat) <- c("domain","kingdom","phylum","class","order","family","genus")

for (i in 1:nrow(map.in)) {
	taxname <- unlist(strsplit(as.character(map.in[i,1]), split=';'))
	while ( length(taxname) > 0) {
		print(taxname)	
		g.hit <- grep(paste("^",paste(taxname,collapse=";"),";","$",sep=""),map.in$taxlabel)
		tax.mat[i,map.in[g.hit,2]] <- tail(taxname,1)		
		taxname <- head(taxname,-1)
	}

	for (j in 1:ncol(tax.mat)) {
		if(is.na(tax.mat[i,j])) { tax.mat[i,j] <- paste(tax.mat[i,j-1],".",sep="")}
	}
}

for (i in 1:nrow(tax.mat)) {
	map.in[i,"taxout"] <- paste(paste(tax.mat[i,c("domain","kingdom","phylum","class","order","family","genus")],collapse=";"),";",sep="")
}
write.table(map.in$output,file="output.csv",row.names=F,sep="\t")

tax.in <- read.table("silva.nr_v123.tax",header=F,stringsAsFactors=F,sep="\t")
colnames(tax.in) <- c("taxid","taxlabel")

tax.out <- merge(map.in,tax.in)

write.table(tax.out[,c("taxid","taxout")],file="silva.nr_v123.euk.taxonomy",sep="\t",row.names=F,quote=F,col.names=F)


  
  


