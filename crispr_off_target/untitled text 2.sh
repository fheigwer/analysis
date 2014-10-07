bowtie2 -p 2 --qc-filter --local --very-sensitive-local -x /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.genome -1 MUT1_ACAGTG_L005_R1_complete_filtered.fastq -2 MUT1_ACAGTG_L005_R2_complete_filtered.fastq | samtools view -u -q 10 -@ 1 - | samtools sort -n -@ 1 -o - temp1 | samtools fixmate -r - - | samtools sort -o - temp5 | samtools mpileup -puD -f /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa - | bcftools call -O b -m -S snps - > MUT1_S1_L001.indels.bcf.gz &
bowtie2 -p 2 --qc-filter --local --very-sensitive-local -x /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.genome -1 MUT2_GCCAAT_L005_R1_complete_filtered.fastq -2 MUT2_GCCAAT_L005_R2_complete_filtered.fastq | samtools view -u -q 10 -@ 1 - | samtools sort -n -@ 1 -o - temp2 | samtools fixmate -r - - | samtools sort -o - temp6 | samtools mpileup -puD -f /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa - | bcftools call -O b -m -S snps - > MUT2_S2_L001.indels.bcf.gz &
bowtie2 -p 2 --qc-filter --local --very-sensitive-local -x /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.genome -1 FM7_CTTGTA_L002_R1_complete_filtered.fastq -2 FM7_CTTGTA_L002_R2_complete_filtered.fastq | samtools view -u -q 10 -@ 1 - | samtools sort -n -@ 1 -o - temp3 | samtools fixmate -r - - | samtools sort -o - temp7 | samtools mpileup -puD -f /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa - | bcftools call -O b -m -S snps - > FM7_S3_L001.indels.bcf.gz &
bowtie2 -p 2 --qc-filter --local --very-sensitive-local -x /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.genome -1 OreR_GTGAAA_L005_R1_complete_filtered.fastq -2 OreR_GTGAAA_L005_R2_complete_filtered.fastq | samtools view -u -q 10 -@ 1 - | samtools sort -n -@ 1 -o - temp4 | samtools fixmate -r - - | samtools sort -o - temp8 | samtools mpileup -puD -f /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa - | bcftools call -O b -m -S snps - > OreR_S4_L001.indels.bcf.gz ;
bowtie2 -p 8 --qc-filter --local --very-sensitive-local -x /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.genome -1 MUT3_GCCAAT_L002_R1_complete_filtered.fastq -2 MUT3_GCCAAT_L002_R2_complete_filtered.fastq | samtools view -u -q 10 -@ 1 - | samtools sort -n -@ 1 -o - temp10 | samtools fixmate -r - - | samtools sort -o - temp9 | samtools mpileup -puD -f /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa - | bcftools call -O b -m -S snps - > MUT3_S1_L001.indels.bcf.gz ;

bcftools index OreR_S4_L001.indels.bcf.gz
bcftools index FM7_S3_L001.indels.bcf.gz
bcftools index MUT2_S2_L001.indels.bcf.gz
bcftools index MUT1_S1_L001.indels.bcf.gz
bcftools index MUT3_S1_L001.indels.bcf.gz

bcftools view -O v OreR_S4_L001.indels.bcf.gz > OreR_S4_L001.indels.vcf &
bcftools view -O v FM7_S3_L001.indels.bcf.gz > FM7_S3_L001.indels.vcf &
bcftools view -O v MUT2_S2_L001.indels.bcf.gz > MUT2_S2_L001.indels.vcf &
bcftools view -O v MUT1_S1_L001.indels.bcf.gz > MUT1_S1_L001.indels.vcf &
bcftools view -O v MUT3_S1_L001.indels.bcf.gz > MUT3_S1_L001.indels.vcf ;

perl ~/Desktop/analysis/crispr_off_target/intersect.pl MUT1_S1_L001.indels.vcf FM7_S3_L001.indels.vcf OreR_S4_L001.indels.vcf > MUT1_woFM7_woOreR.vcf &
perl ~/Desktop/analysis/crispr_off_target/intersect.pl MUT2_S2_L001.indels.vcf FM7_S3_L001.indels.vcf OreR_S4_L001.indels.vcf > MUT2_woFM7_woOreR.vcf &
perl ~/Desktop/analysis/crispr_off_target/intersect.pl MUT3_S1_L001.indels.vcf FM7_S3_L001.indels.vcf OreR_S4_L001.indels.vcf > MUT3_woFM7_woOreR.vcf ;

perl ~/Desktop/analysis/crispr_off_target/final_filter_vcf.pl MUT1_woFM7_woOreR.vcf > MUT1_woFM7_woOreR.filtered.vcf &
perl ~/Desktop/analysis/crispr_off_target/final_filter_vcf.pl MUT2_woFM7_woOreR.vcf > MUT2_woFM7_woOreR.filtered.vcf &
perl ~/Desktop/analysis/crispr_off_target/final_filter_vcf.pl MUT3_woFM7_woOreR.vcf > MUT3_woFM7_woOreR.filtered.vcf ;

perl ~/Desktop/analysis/crispr_off_target/intersect.pl MUT1_S1_L001.indels.vcf MUT2_S2_L001.indels.vcf MUT3_S1_L001.indels.vcf > MUT1_woFM7_woOreR_others.vcf &
perl ~/Desktop/analysis/crispr_off_target/intersect.pl MUT2_S2_L001.indels.vcf MUT1_S1_L001.indels.vcf MUT3_S1_L001.indels.vcf > MUT2_woFM7_woOreR_others.vcf &
perl ~/Desktop/analysis/crispr_off_target/intersect.pl MUT3_S1_L001.indels.vcf MUT1_S1_L001.indels.vcf MUT2_S2_L001.indels.vcf > MUT3_woFM7_woOreR_others.vcf ;