#make a kraken database of LSU and SSU sequences to screen for rRNA contamination in RNA-seq

wget http://www.arb-silva.de/fileadmin/silva_databases/release_123/Exports/taxonomy/taxmap_embl_lsu_ref_123.txt
wget http://www.arb-silva.de/fileadmin/silva_databases/release_123/Exports/taxonomy/taxmap_embl_ssu_ref_123.txt

wget http://www.arb-silva.de/fileadmin/silva_databases/release_123/Exports/SILVA_123_SSURef_tax_silva_trunc.fasta.gz
wget http://www.arb-silva.de/fileadmin/silva_databases/release_123/Exports/SILVA_123_LSURef_tax_silva_trunc.fasta.gz

cut -f1,6 taxmap_embl_ssu_ref_123.txt | tr -d '\t' > ssu.gi.acc
cut -f1,6 taxmap_embl_lsu_ref_123.txt | tr -d '\t' > lsu.gi.acc

#uggh weirdo lines
#CBMY010000039	447	3319	unclassified sequences;	Helicobacter pylori SA261C	 -1048507
#CBMY010000036	438	1920	unclassified sequences;	Helicobacter pylori SA261C	 -1048507

./replace-ids.pl ssu.gi.acc SILVA_123_SSURef_tax_silva_trunc.fasta > SILVA_123_SSURef_tax_silva_trunc.taxid.fasta
./replace-ids.pl lsu.gi.acc SILVA_123_LSURef_tax_silva_trunc.fasta > SILVA_123_LSURef_tax_silva_trunc.taxid.fasta

for file in *.taxid.fasta
do
   kraken-build --add-to-library $file --db rRNA
done

#spend three hours trying to build jellyfish in MacOSX
#give up and do it in the BiolinuxVM

kraken-build --build --db rRNA --jellyfish-hash-size 1000M --max-db-size 64

 
 
 
 
 
 
 
 
 
 
 
 
 
 