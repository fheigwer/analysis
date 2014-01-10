#!/usr/bin/bash

condetri -fastq1=$1R1_001.fastq -fastq2=$1R2_001.fastq -sc=33 -lq=5 -hq=15 -rmN -minlen=130;

#fastq- $1R1_001_trim1.fastq $1R1_001_trim2.fastq -o $1.fastq

#./parsefortrans $1.fastq ;
perl /Users/b110-mm06/Desktop/analysis/transposon/parsefortrans_paired.pl $1R1_001_trim1.fastq $1R1_001_trim2.fastq $1R1_001_trim_unpaired.fastq fastq $1 

bowtie2 -p 8 -x /Users/b110-mm06/Desktop/Projects/transposons/database/Homo_sapiens.GRCh37.73.genome -1 $1_5_R1.fastq -2 $1_5_R2.fastq -U $1_5_U.fastq -S $1_5prime.sam 
bowtie2 -p 8 -x /Users/b110-mm06/Desktop/Projects/transposons/database/Homo_sapiens.GRCh37.73.genome -1 $1_3_R1.fastq -2 $1_3_R2.fastq -U $1_3_U.fastq  -S $1_3prime.sam

samtools view -b -q 10 -S $1_5prime.sam > $1_5prime.bam 
samtools view -b -q 10 -S $1_3prime.sam > $1_3prime.bam 

samtools sort  $1_5prime.bam $1_5prime_sorted
samtools sort  $1_3prime.bam $1_3prime_sorted

samtools index $1_5prime_sorted.bam
samtools index $1_3prime_sorted.bam

macs14 -t $1_5prime_sorted.bam -f BAM -g hs -n $1_5prime_sorted
macs14 -t $1_3prime_sorted.bam -f BAM -g hs -n $1_3prime_sorted