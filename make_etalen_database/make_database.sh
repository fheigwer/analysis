#!/usr/bin/bash

echo $1;
#gunzip $1*.gz;
#cat *.dat >> $1.all.gb;
#rm *.dat;
#perl gb_to_gff.pl $1.all.gb;
#cat *.fasta >> $1.all.dna.fa;
perl cpgi_for_all_files.pl /home/mount/fheigwer/make_etalen_database
rm *.fasta;
perl correct_cdna.pl $1.cdna.all.fa;
rm $1.cdna.all.fa;
mv $1.cdna.all.facorrected.fa $1.cdna.all.fa;
perl build_tree.pl $1;
rm *.gff;
rm *.gz.1
rm CHECKSUMS.*
rm README.*
perl include_cpg.pl /home/mount/fheigwer/make_etalen_database
#bowtie2-build $1.cdna.all.fa $1.cdna & bowtie2-build $1.all.dna.fa $1.dna;
echo "done";

