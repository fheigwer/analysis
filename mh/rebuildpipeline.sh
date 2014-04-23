#!/usr/bin/bash

#quality trim the reads to ensure higher mapping quality, accuracy and spead
condetri -fastq1=1182_4H-1_S1_L001_R1_001.fastq -fastq2=1182_4H-1_S1_L001_R2_001.fastq & condetri -fastq1=1182_4H-3_S2_L001_R1_001.fastq -fastq2=1182_4H-3_S2_L001_R2_001.fastq

#align the reads to a reference genome with two different methods
novoalign -d dmel_genome -f 1182_4H-1_S1_L001_R1_001_trim1.fastq 1182_4H-1_S1_L001_R1_001_trim2.fastq -o SAM > 1182_4H-1_out.sam & novoalign -d dmel_genome -f 1182_4H-3_S2_L001_R1_001_trim1.fastq 1182_4H-3_S2_L001_R1_001_trim2.fastq -o SAM > 1182_4H-3_out.sam

#convert the alignment files from sam to bam to save space
samtools view -b -S 1182_4H-1_out.sam > 1182_4H-1_out.bam
samtools view -b -S 1182_4H-3_out.sam > 1182_4H-3_out.bam

# since both 1182_4H samples are from the same cell line but run in different lanes/barcodes, we will merge the files and create a bam file
samtools merge 1182_4H_out.bam 1182_4H-1_out.bam 1182_4H-3_out.bam

samtools view -h -o 1182_4H_out.sam 1182_4H_out.bam 

# filter for minimum mapping quality Q20
samtools view -h -b -q 20 -o 1182_4H_out_filtered.bam 1182_4H_out.bam

# sort the alignment file
samtools sort 1182_4H_out_filtered.bam 1182_4H_out_filtered_sorted

# remove duplicates
samtools rmdup 1182_4H_out_filtered_sorted.bam 1182_4H_out_filtered_sorted_dedup.bam 

# check whether .bam file is valid
java -jar /usr/bin/ValidateSamFile.jar I=1182_4H_out_filtered_sorted_dedup.bam

# Fix .bam using Picard 
java -jar /usr/bin/AddOrReplaceReadGroups.jar I=1182_4H_out_filtered_sorted_dedup.bam O=1182_4H_out_filtered_sorted_dedup_fixed.bam \
	LB=TruSeqLT PL=Illumina PU=1 SM=1182_4H
	
# Create index file for dmel_gemome.fa
java -jar /usr/bin/CreateSequenceDictionary.jar R=dmel_genome.69.fa O=dmel_genome.69.dict
samtools faidx dmel_genome.69.fa

# Create .bai file
samtools index 1182_4H_out_filtered_sorted_dedup_fixed.bam

# Calculate Depth of Coverage using GATK
java -Xmx2g -jar /usr/bin/GenomeAnalysisTK.jar -R dmel_genome.69.fa \
	-T DepthOfCoverage \
	-o 1182_4H_DepthOfCoverage \
	-I 1182_4H_out_filtered_sorted_dedup_fixed.bam

	
# BAQ strategy - creating a .baq file (optional)
samtools calmd -Abr 1182_4H_out_filtered_sorted_dedup_fixed.bam dmel_genome.69.fa > 1182_4H_out_filtered_sorted_dedup_fixed.baq.bam 

# Local re-align usign GATK
java -Xmx1g -jar /usr/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator \
	-R dmel_genome.69.fa -I 1182_4H_out_filtered_sorted_dedup_fixed.bam \
	-o 1182_4H_out_filtered_sorted_dedup_fixed.bam.realign.intervals -L X:10819664-13621237

java -Xmx4g -jar /usr/bin/GenomeAnalysisTK.jar -I 1182_4H_out_filtered_sorted_dedup_fixed.bam \
	-R dmel_genome.69.fa -T IndelRealigner \
	-targetIntervals 1182_4H_out_filtered_sorted_dedup_fixed.bam.realign.intervals -o 1182_4H_out_filtered_sorted_dedup_fixed.realigned.bam

# fix a formatting problem in the re-aligned BAM files 
java -jar /usr/bin/AddOrReplaceReadGroups.jar \
	I=1182_4H_out_filtered_sorted_dedup_fixed.realigned.bam \
	O=1182_4H_out_filtered_sorted_dedup_fixed.realigned.fixed.bam \
	SORT_ORDER=coordinate RGID=1182_4H RGLB=TruSeqLT RGPL=Illumina RGPU=1 RGSM=1182_4H \
	CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT

# run GATK UnifiedGenotyper (normal ploidy)
#java -jar /usr/bin/GenomeAnalysisTK.jar -R dmel_genome.69.fa -T UnifiedGenotyper \
#	-I 1182_4H_out_filtered_sorted_dedup_fixed.realigned.fixed.bam \
#	-o 1182_4H_GATK.vcf \
#	-stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 500 -L  X:10819664-13621237	

# run GATK UnifiedGenotyper (haploidy)
java -jar /usr/bin/GenomeAnalysisTK.jar -R dmel_genome.69.fa -T UnifiedGenotyper \
	-I 1182_4H_out_filtered_sorted_dedup_fixed.realigned.fixed.bam \
	-o 1182_4H_ploidy_GATK.vcf \
	-stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 500 -ploidy 1 -L  X:10819664-13621237

# Download database for SNPeff
java -Xmx4g -jar ~/snpEff/snpEff.jar download -c ~/snpEff/snpEff.config -v BDGP5.69

# annotate SNPs using snpEff
java -Xmx4g -jar ~/snpEff/snpEff.jar eff -c ~/snpEff/snpEff.config -v BDGP5.69 1182_4H_ploidy_GATK.vcf > 1182_4H_ploidy_GATK.eff.vcf
	
# annotate previously known SNPs using SNPsif
java -Xmx4g -jar ~/snpEff/SnpSift.jar annotate -v vcf_chr_X.vcf 1182_4H_ploidy_GATK.vcf > 1182_4H_ploidy_GATK.dbsnp.vcf
java -Xmx4g -jar ~/snpEff/SnpSift.jar filter -f 1182_4H_ploidy_GATK.dbsnp.vcf "! exists ID" > 1182_4H_ploidy_GATK.dbsnp.not_in_dbSnp.vcf
java -Xmx4g -jar ~/snpEff/snpEff.jar eff -c ~/snpEff-/snpEff.config -v 1182_4H_ploidy_GATK.dbsnp.not_in_dbSnp.vcf > 1182_4H_ploidy_GATK.dbsnp.not_in_dbSnp.eff.vcf

	
	
	