#!/usr/bin/bash

condetri -fastq1=TraDIS-Drosophila-3-5_S1_L001_R1_001.fastq -fastq2=TraDIS-Drosophila-3-5_S1_L001_R2_001.fastq -q33 -rmN -hq=15 -lq=5 -minlen=50;

fastq-join TraDIS-Drosophila-3-5_S1_L001_R1_001_trim1.fastq TraDIS-Drosophila-3-5_S1_L001_R1_001_trim2.fastq -o TraDIS-Drosophila-3-5_S1_L001_join.fastq

./parsefortrans TraDIS-Drosophila-3-5_S1_L001_join.fastq ;

bowtie2 -p 8 -x dmel_genome -U TraDIS-Drosophila-3-5_S1_L001_join_35prime.fastq -S TraDIS-Drosophila-3-5_S1_L001_join_35prime.sam --qc-filter --time
bowtie2 -p 8 -x dmel_genome -U TraDIS-Drosophila-3-5_S1_L001_join_53prime.fastq -S TraDIS-Drosophila-3-5_S1_L001_join_53prime.sam --qc-filter --time
bowtie2 -p 8 -x dmel_genome -U TraDIS-Drosophila-3-5_S1_L001_join_other.fastq -S TraDIS-Drosophila-3-5_S1_L001_join_other.sam --qc-filter --time

samtools view -b -S TraDIS-Drosophila-3-5_S1_L001_join_35prime.sam > TraDIS-Drosophila-3-5_S1_L001_join_35prime.bam 
samtools view -b -S TraDIS-Drosophila-3-5_S1_L001_join_53prime.sam > TraDIS-Drosophila-3-5_S1_L001_join_53prime.bam 
samtools view -b -S TraDIS-Drosophila-3-5_S1_L001_join_other.sam > TraDIS-Drosophila-3-5_S1_L001_join_other.bam 

samtools sort  TraDIS-Drosophila-3-5_S1_L001_join_35prime.bam TraDIS-Drosophila-3-5_S1_L001_join_35prime_sorted
samtools sort  TraDIS-Drosophila-3-5_S1_L001_join_53prime.bam TraDIS-Drosophila-3-5_S1_L001_join_53prime_sorted
samtools sort  TraDIS-Drosophila-3-5_S1_L001_join_other.bam TraDIS-Drosophila-3-5_S1_L001_join_other_sorted

samtools index TraDIS-Drosophila-3-5_S1_L001_join_35prime_sorted.bam
samtools index TraDIS-Drosophila-3-5_S1_L001_join_53prime_sorted.bam
samtools index TraDIS-Drosophila-3-5_S1_L001_join_other_sorted.bam

macs14 -t TraDIS-Drosophila-3-5_S1_L001_join_35prime_sorted.bam -f BAM -g dm -n TraDIS-Drosophila-3-5_S1_L001_join_35prime_sorted
macs14 -t TraDIS-Drosophila-3-5_S1_L001_join_53prime_sorted.bam -f BAM -g dm -n TraDIS-Drosophila-3-5_S1_L001_join_53prime_sorted
macs14 -t TraDIS-Drosophila-3-5_S1_L001_join_other_sorted.bam -f BAM -g dm -n TraDIS-Drosophila-3-5_S1_L001_join_other_sorted