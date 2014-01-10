bowtie2 -p 2 --qc-filter --mm --very-sensitive-local -x /Users/b110-mm06/Desktop/Projects/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.73.genome -1 F2_S1_L001_R1_001.fastq -2 F2_S1_L001_R2_001.fastq -S F2_S1_L001.sam &
bowtie2 -p 2 --qc-filter --mm --very-sensitive-local -x /Users/b110-mm06/Desktop/Projects/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.73.genome -1 F3_S2_L001_R1_001.fastq -2 F3_S2_L001_R2_001.fastq -S F3_S2_L001.sam &
bowtie2 -p 2 --qc-filter --mm --very-sensitive-local -x /Users/b110-mm06/Desktop/Projects/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.73.genome -1 FM7_S3_L001_R1_001.fastq -2 FM7_S3_L001_R2_001.fastq -S FM7_S3_L001.sam &
bowtie2 -p 2 --qc-filter --mm --very-sensitive-local -x /Users/b110-mm06/Desktop/Projects/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.73.genome -1 OreR_S4_L001_R1_001.fastq -2 OreR_S4_L001_R2_001.fastq -S OreR_S4_L001.sam &

samtools view -u -q 10 -@ 8 -S F2_S1_L001.sam  | samtools sort -n -@ 8 -o - temp | samtools fixmate -r - - | samtools sort -o - tempo | samtools mpileup -uD -f /Users/b110-mm06/Desktop/Projects/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa - | bcftools call -O u -m -S snps -v - > F2_S1_L001.indels.bcf ;
samtools view -u -q 10 -@ 8 -S F3_S2_L001.sam  | samtools sort -n -@ 8 -o - temp | samtools fixmate -r - - | samtools sort -o - tempo | samtools mpileup -uD -f /Users/b110-mm06/Desktop/Projects/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa - | bcftools call -O u -m -S snps -v - > F3_S2_L001.indels.bcf ;
samtools view -u -q 10 -@ 8 -S FM7_S3_L001.sam  | samtools sort -n -@ 8 -o - temp | samtools fixmate -r - - | samtools sort -o - tempo | samtools mpileup -uD -f /Users/b110-mm06/Desktop/Projects/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa - | bcftools call -O u -m -S snps -v - > FM7_S3_L001.indels.bcf ;
samtools view -u -q 10 -@ 8 -S OreR_S4_L001.sam  | samtools sort -n -@ 8 -o - temp | samtools fixmate -r - - | samtools sort -o - tempo | samtools mpileup -uD -f /Users/b110-mm06/Desktop/Projects/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa - | bcftools call -O u -m -S snps -v - > OreR_S4_L001.indels.bcf &

bgzip OreR_S4_L001.indels.bcf
bgzip FM7_S3_L001.indels.bcf
bgzip F3_S2_L001.indels.bcf
bgzip F2_S1_L001.indels.bcf

bcftools index OreR_S4_L001.indels.bcf.gz
bcftools index FM7_S3_L001.indels.bcf.gz
bcftools index F3_S2_L001.indels.bcf.gz
bcftools index F2_S1_L001.indels.bcf.gz

bcftools merge -O b OreR_S4_L001.indels.bcf.gz FM7_S3_L001.indels.bcf.gz > wt.bcf.gz
bcftools merge -O b F3_S2_L001.indels.bcf.gz F2_S1_L001.indels.bcf.gz > mut.bcf.gz

bcftools index wt.bcf.gz
bcftools index mut.bcf.gz

bcftools isec  mut.bcf.gz wt.bcf.gz --complement 1 > mutant.only.variants.vcf




samtools view -u -q 10 -@ 8 -S F2_S1_L001.sam  | samtools sort -n -@ 8 -o - temp1 | samtools fixmate -r - - | samtools sort -o - temp5 > F2_S1_L001_clean.bam &
samtools view -u -q 10 -@ 8 -S F3_S2_L001.sam  | samtools sort -n -@ 8 -o - temp2 | samtools fixmate -r - - | samtools sort -o - temp6 > F3_S2_L001_clean.bam &
samtools view -u -q 10 -@ 8 -S FM7_S3_L001.sam  | samtools sort -n -@ 8 -o - temp3 | samtools fixmate -r - - | samtools sort -o - temp7 > FM7_S3_L001_clean.bam &
samtools view -u -q 10 -@ 8 -S OreR_S4_L001.sam  | samtools sort -n -@ 8 -o - temp4 | samtools fixmate -r - - | samtools sort -o - temp8 > OreR_S4_L001_clean.bam ;


java -jar /usr/bin/AddOrReplaceReadGroups.jar I=F2_S1_L001_clean.bam O=F2_S1_L001_clean_fixed.bam LB=TruSeqLT PL=Illumina PU=1 SM=F2 &
java -jar /usr/bin/AddOrReplaceReadGroups.jar I=F3_S2_clean.bam O=F3_S2_L001_clean_fixed.bam LB=TruSeqLT PL=Illumina PU=1 SM=F3 &
java -jar /usr/bin/AddOrReplaceReadGroups.jar I=FM7_S3_L001_clean.bam O=FM7_S3_L001_clean_fixed.bam LB=TruSeqLT PL=Illumina PU=1 SM=FM7 &
java -jar /usr/bin/AddOrReplaceReadGroups.jar I=OreR_S4_L001_clean.bam O=OreR_S4_L001_clean_fixed.bam LB=TruSeqLT PL=Illumina PU=1 SM=OreR ;

samtools index F2_S1_L001_sorted_fixed.bam
java -Xmx1g -jar /usr/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /Users/b110-mm06/Desktop/Projects/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa -I F2_S1_L001_sorted_fixed.bam -o F2_S1_L001.realign.intervals
java -Xmx4g -jar /usr/bin/GenomeAnalysisTK.jar -I F2_S1_L001_sorted_fixed.bam -R /Users/b110-mm06/Desktop/Projects/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa -T IndelRealigner -targetIntervals F2_S1_L001.realign.intervals -o F2_S1_L001_sorted_fixed_realigned.bam

