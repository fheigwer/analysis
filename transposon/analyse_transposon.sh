#!/usr/bin/bash

condetri -fastq1=$1.R1_001.fastq -fastq2=$1.R2_001.fastq -q33 -rmN -hq=15 -lq=5 -minlen=50;

fastq-join $1.R1_001_trim1.fastq $1.R1_001_trim2.fastq -o $1.join.fastq

./parsefortrans $1.join.fastq ;

bowtie2 -p 8 -x dmel_genome -U $1.join_35prime.fastq -S $1.join_35prime.sam --qc-filter --time
bowtie2 -p 8 -x dmel_genome -U $1.join_53prime.fastq -S $1.join_53prime.sam --qc-filter --time
bowtie2 -p 8 -x dmel_genome -U $1.join_other.fastq -S $1.join_other.sam --qc-filter --time

samtools view -b -S $1.join_35prime.sam > $1.join_35prime.bam 
samtools view -b -S $1.join_53prime.sam > $1.join_53prime.bam 
samtools view -b -S $1.join_other.sam > $1.join_other.bam 

samtools sort  $1.join_35prime.bam $1.join_35prime_sorted
samtools sort  $1.join_53prime.bam $1.join_53prime_sorted
samtools sort  $1.join_other.bam $1.join_other_sorted

samtools index $1.join_35prime_sorted.bam
samtools index $1.join_53prime_sorted.bam
samtools index $1.join_other_sorted.bam

macs14 -t $1.join_35prime_sorted.bam -f BAM -g dm -n $1.join_35prime_sorted
macs14 -t $1.join_53prime_sorted.bam -f BAM -g dm -n $1.join_53prime_sorted
macs14 -t $1.join_other_sorted.bam -f BAM -g dm -n $1.join_other_sorted