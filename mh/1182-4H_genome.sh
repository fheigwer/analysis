#!/bin/bash

### maternal haploid (mh) project
# analysis of 1182-4H cell line genomic sequence to identify mh mutation

## Final workflow

# create novoalign index for flies (here BDGP5.25)
if [ ! -f dmel_genome ]; then
    echo "Novoalign genome index not found. Creating new index."
	~/Genome/bin/novocraft/novoindex dmel_genome dmel_genome.fa
fi

# run novoalign
~/Genome/bin/novocraft/novoalign -d dmel_genome -f 1182-4H-1_S1_L001_R1_001.fastq 1182-4H-1_S1_L001_R2_001.fastq -o SAM > 1182_4H-1_out.sam
~/Genome/bin/novocraft/novoalign -d dmel_genome -f 1182-4H-3_S2_L001_R1_001.fastq 1182-4H-3_S2_L001_R2_001.fastq -o SAM > 1182-4H-3_out.sam

# since both 1182-4H samples are from the same cell line but run in different lanes/barcodes, we will merge the files and create a bam file
java -jar ~/Genome/bin/SortSam.jar INPUT=1182-4H-1_out.sam OUTPUT=1182-4H-1_out.bam SO=unsorted 
java -jar ~/Genome/bin/SortSam.jar INPUT=1182-4H-3_out.sam OUTPUT=1182-4H-3_out.bam SO=unsorted 
~/Genome/bin/samtools merge 1182-4H_out.bam 1182-4H-1_out.bam 1182-4H-3_out.bam 
~/Genome/bin/samtools view -h -o 1182-4H_out.sam 1182-4H_out.bam 

# filter for Q20
~/Genome/bin/samtools view -h -b -q 30 -o 1182-4H_out_filtered.bam 1182-4H_out.bam

# sort alignment file
~/Genome/bin/samtools sort 1182-4H_out_filtered.bam 1182-4H_out_filtered_sorted

# remove duplicates
~/Genome/bin/samtools rmdup 1182-4H_out_filtered_sorted.bam 1182-4H_out_filtered_sorted_dedup.bam 

# check whether .bam file is valid
java -jar ~/Genome/bin/ValidateSamFile.jar I=1182-4H_out_filtered_sorted_dedup.bam

# Fix .bam using Picard 
java -jar ~/Genome/bin/AddOrReplaceReadGroups.jar I=1182-4H_out_filtered_sorted_dedup.bam O=1182-4H_out_filtered_sorted_dedup_fixed.bam \
	LB=TruSeqLT PL=Illumina PU=1 SM=1182-4H

# Create index file for dmel_gemome.fa
java -jar ~/Genome/bin/CreateSequenceDictionary.jar R=dmel_genome.fa O=dmel_genome.dict
~/Genome/bin/samtools faidx dmel_genome.fa

# Fix order in .bam to reference genome
# java -Xmx2g -jar ~/Genome/bin/ReorderSam.jar I=1182-4H_out_filtered_sorted_dedup_fixed.bam \
#	O=1182-4H_out_filtered_sorted_dedup_fixed_ordered.bam \
#	REFERENCE=dmel_genome.fa

# Create .bai file
~/Genome/bin/samtools index 1182-4H_out_filtered_sorted_dedup_fixed.bam

# Calculate Depth of Coverage using GATK
java -Xmx2g -jar ~/Genome/bin/GenomeAnalysisTK.jar -R dmel_genome.fa \
	-T DepthOfCoverage \
	-o 1182-4H_DepthOfCoverage \
	-I 1182-4H_out_filtered_sorted_dedup_fixed.bam

# BAQ strategy - creating a .baq file (optional)
~/Genome/bin/samtools calmd -Abr 1182-4H_out_filtered_sorted_dedup_fixed.bam dmel_genome.fa > 1182-4H_out_filtered_sorted_dedup_fixed.baq.bam 

# Local re-align usign GATK
java -Xmx1g -jar ~/Genome/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator \
	-R dmel_genome.fa -I 1182-4H_out_filtered_sorted_dedup_fixed.bam \
	-o 1182-4H_out_filtered_sorted_dedup_fixed.bam.realign.intervals -L X:13121101-13621101

java -Xmx4g -jar ~/Genome/bin/GenomeAnalysisTK.jar -I 1182-4H_out_filtered_sorted_dedup_fixed.bam \
	-R dmel_genome.fa -T IndelRealigner \
	-targetIntervals 1182-4H_out_filtered_sorted_dedup_fixed.bam.realign.intervals -o 1182-4H_out_filtered_sorted_dedup_fixed.realigned.bam

# fix a formatting problem in the re-aligned BAM files 
java -jar ~/Genome/bin/AddOrReplaceReadGroups.jar \
	I=1182-4H_out_filtered_sorted_dedup_fixed.realigned.bam \
	O=1182-4H_out_filtered_sorted_dedup_fixed.realigned.fixed.bam \
	SORT_ORDER=coordinate RGID=1182-4H RGLB=TruSeqLT RGPL=Illumina RGPU=1 RGSM=1182-4H \
	CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT

# run GATK UnifiedGenotyper (normal ploidy)
java -jar ~/Genome/bin/GenomeAnalysisTK.jar -R dmel_genome.fa -T UnifiedGenotyper \
	-I 1182-4H_out_filtered_sorted_dedup_fixed.realigned.fixed.bam \
	-o 1182-4H_GATK.vcf \
	-stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 500 -L  X:13121101-13621101

# run GATK UnifiedGenotyper (with ploidy = 1 since 1182-H4 are believed to be haploid)
java -jar ~/Genome/bin/GenomeAnalysisTK.jar -R dmel_genome.fa -T UnifiedGenotyper \
	-I 1182-4H_out_filtered_sorted_dedup_fixed.realigned.fixed.bam \
	-o 1182-4H_GATK_ploidy.vcf \
	-stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 500 -ploidy 1 -L X:13121101-13621101

# Download database for SNPeff
java -Xmx4g -jar ~/snpEff/snpEff.jar download -c ~/snpEff/snpEff.config -v BDGP5.69

# annotate SNPs using snpEff
java -Xmx4g -jar ~/Genome/bin/snpEff-3.2/snpEff.jar eff -c ~/Genome/bin/snpEff-3.2/snpEff.config -v BDGP5.69 1182-4H_GATK_ploidy.vcf > 1182-4H_GATK_ploidy1.eff.vcf

open -a Safari snpEff_summary.html
open -a snpEff_genes.txt

# annotate previously known SNPs using SNPsif
java -Xmx4g -jar ~/Genome/bin/snpEff-3.2/SnpSift.jar annotate -v vcf_chr_X.vcf 1182-4H_GATK_ploidy.vcf > 1182-4H_GATK_ploidy.dbsnp.vcf
java -Xmx4g -jar ~/Genome/bin/snpEff-3.2/snpEff.jar eff -c ~/Genome/bin/snpEff-3.2/snpEff.config -v BDGP5.69 1182-4H_GATK_ploidy.dbsnp.vcf > 1182-4H_GATK_ploidy1.dbsnp.eff.vcf
java -Xmx4g -jar ~/Genome/bin/snpEff-3.2/SnpSift.jar filter -f 1182-4H_GATK_ploidy1.dbsnp.eff.vcf "! exists ID" > 1182-4H_GATK_ploidy1.dbsnp.eff.not_in_dbSnp.vcf


exit 0;

## worksheet

# Setup and file download
curl -O http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.10.1.zip

# Misc commands
open -a Safari fastqc_report.html 

# samtools index and stats (needs sorted .bam file)
~/Genome/bin/samtools idxstats 1182-4H_out.bam
~/Genome/bin/samtools index 1182-4H_out.bam

# create novoalign index for flies (here BDGP5.25)
~/Genome/bin/novocraft/novoindex dmel_genome dmel_genome.fa

# run novoalign
~/Genome/bin/novocraft/novoalign -d dmel_genome -f 1182-4H-1_S1_L001_R1_001.fastq 1182-4H-1_S1_L001_R2_001.fastq -o SAM > 1182_4H-1_out.sam
~/Genome/bin/novocraft/novoalign -d dmel_genome -f 1182-4H-3_S2_L001_R1_001.fastq 1182-4H-3_S2_L001_R2_001.fastq -o SAM > 1182-4H-3_out.sam

# since both 1182-4H samples are from the same cell line but run in different lanes/barcodes, we will merge the files and create a bam file
java -jar ~/Genome/bin/SortSam.jar INPUT=1182-4H-1_out.sam OUTPUT=1182-4H-1_out.bam SO=unsorted 
java -jar ~/Genome/bin/SortSam.jar INPUT=1182-4H-3_out.sam OUTPUT=1182-4H-3_out.bam SO=unsorted 
~/Genome/bin/samtools merge 1182-4H_out.bam 1182-4H-1_out.bam 1182-4H-3_out.bam 
~/Genome/bin/samtools view -h -o 1182-4H_out.sam 1182-4H_out.bam 
less 1182-4H_out.sam 

# filter for q20
~/Genome/bin/samtools view -h -b -q 20 -o 1182-4H_out_filtered.bam 1182-4H_out.bam

# sort alignment file
~/Genome/bin/samtools sort 1182-4H_out_filtered.bam -f 1182-4H_out_filtered_sorted.bam

# remove duplicates
~/Genome/bin/samtools rmdup 1182-4H_out_filtered_sorted.bam 1182-4H_out_filtered_sorted_dedup.bam 

# call SNP using samtools
# samtools mpileup -uD -r 2L:100,000-150,000 -f /data/drosophila/dmel-all-chromosome-r5.37.fasta ../RAL357_full_bwa.sorted.bam ../RAL391_full_bwa.sorted.bam | bcftools view -bvcg - > RAL_samtools.raw.bcf

~/Genome/bin/samtools mpileup -uD -f ~/Genome/db/dmel-all-chromosome-r5.51.fasta 1182-4H_out_filtered_sorted_dedup.bam | ~/Genome/bin/bcftools view -bvcg - > 1182-4H_out_filtered_sorted_dedup.bcf
~/Genome/bin/bcftools view 1182-4H_out_filtered_sorted_dedup.bcf | ~/Genome/bin/vcfutils.pl varFilter -D500 > 1182-4H_out_filtered_sorted_dedup.vcf

# call SNP using GATK
java -Xmx1g -jar ~/Genome/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator \
	-R ~/Genome/db/dmel-all-chromosome-r5.51.fasta -I ../RAL357_full_bwa.sorted.bam \
		-o RAL357.realign.intervals -L 2L:100000-150000



## Calculate Depth of coverage using GATK
# check whether .bam file is valid
java -jar ~/Genome/bin/ValidateSamFile.jar I=1182-4H_out_filtered_sorted_dedup.bam

# If not, add read group etc.
java -jar ~/Genome/bin/AddOrReplaceReadGroups.jar I=1182-4H_out_filtered_sorted_dedup.bam O=1182-4H_out_filtered_sorted_dedup_fixed.bam \
	LB=TruSeqLT PL=Illumina PU=1 SM=1182-4H

# Create sequence directory 
java -jar ~/Genome/bin/CreateSequenceDictionary.jar R=~/Genome/db/dmel-all-chromosome-r5.51.fasta O=~/Genome/db/dmel-all-chromosome-r5.51.dict
java -Xmx2g -jar ~/Genome/bin/GenomeAnalysisTK.jar -R ~/Genome/db/dmel-all-chromosome-r5.51.fasta -T DepthOfCoverage -o 1182-4H_Depth -I 1182-4H_out_filtered_sorted_dedup.bam


# template
java -Xmx2g -jar ~/Genome/bin/GenomeAnalysisTK.jar \
  -R ~/Genome/db/dmel-all-chromosome-r5.51.fasta \
  -T Coverage \
  -o 1182-4H_Depth \
  -I 1182-4H_out_filtered_sorted_dedup.bam
