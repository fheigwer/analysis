bowtie2 -p 2 --qc-filter --mm --very-sensitive-local -x /Users/b110-mm06/Desktop/Projects/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.73.genome -1 F2_S1_L001_R1_001.fastq -2 F2_S1_L001_R2_001.fastq -S F2_S1_L001.sam &
bowtie2 -p 2 --qc-filter --mm --very-sensitive-local -x /Users/b110-mm06/Desktop/Projects/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.73.genome -1 F3_S2_L001_R1_001.fastq -2 F3_S2_L001_R2_001.fastq -S F3_S2_L001.sam &
bowtie2 -p 2 --qc-filter --mm --very-sensitive-local -x /Users/b110-mm06/Desktop/Projects/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.73.genome -1 FM7_S3_L001_R1_001.fastq -2 FM7_S3_L001_R2_001.fastq -S FM7_S3_L001.sam &
bowtie2 -p 2 --qc-filter --mm --very-sensitive-local -x /Users/b110-mm06/Desktop/Projects/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.73.genome -1 OreR_S4_L001_R1_001.fastq -2 OreR_S4_L001_R2_001.fastq -S OreR_S4_L001.sam &

vcftools --vcf merged_variants.raw.vcf --keep-only-indels --min-meanDP 5 --minQ --recode-to-stream > merged_variants.raw.indel.minqual.vcf

java -jar /usr/bin/AddOrReplaceReadGroups.jar I=F2_S1_L001_sorted.bam O=F2_S1_L001_sorted_fixed.bam LB=TruSeqLT PL=Illumina PU=1 SM=F2
samtools index F2_S1_L001_sorted_fixed.bam
java -Xmx1g -jar /usr/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /Users/b110-mm06/Desktop/Projects/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa -I F2_S1_L001_sorted_fixed.bam -o F2_S1_L001.realign.intervals
java -Xmx4g -jar /usr/bin/GenomeAnalysisTK.jar -I F2_S1_L001_sorted_fixed.bam -R /Users/b110-mm06/Desktop/Projects/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa -T IndelRealigner -targetIntervals F2_S1_L001.realign.intervals -o F2_S1_L001_sorted_fixed_realigned.bam



samtools view -u -q 10 -@ 8 -S F2_S1_L001.sam  | samtools sort -n -@ 8 -o - temp | samtools fixmate -r - - | samtools sort -o - tempo | samtools mpileup -uD -f /Users/b110-mm06/Desktop/Projects/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa - | bcftools call -O u -m -S snps -v - > F2_S1_L001.indels.bcf ;
samtools view -u -q 10 -@ 8 -S F3_S2_L001.sam  | samtools sort -n -@ 8 -o - temp | samtools fixmate -r - - | samtools sort -o - tempo | samtools mpileup -uD -f /Users/b110-mm06/Desktop/Projects/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa - | bcftools call -O u -m -S snps -v - > F3_S2_L001.indels.bcf ;
samtools view -u -q 10 -@ 8 -S FM7_S3_L001.sam  | samtools sort -n -@ 8 -o - temp | samtools fixmate -r - - | samtools sort -o - tempo | samtools mpileup -uD -f /Users/b110-mm06/Desktop/Projects/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa - | bcftools call -O u -m -S snps -v - > FM7_S3_L001.indels.bcf ;
samtools view -u -q 10 -@ 8 -S OreR_S4_L001.sam  | samtools sort -n -@ 8 -o - temp | samtools fixmate -r - - | samtools sort -o - tempo | samtools mpileup -uD -f /Users/b110-mm06/Desktop/Projects/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa - | bcftools call -O u -m -S snps -v - > OreR_S4_L001.indels.bcf &
