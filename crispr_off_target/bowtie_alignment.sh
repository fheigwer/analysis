bowtie2 -p 2 --qc-filter --local --very-sensitive-local -x /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.genome -1 F2_ACAGTG_L005_R1_complete_filtered.fastq -2 F2_ACAGTG_L005_R2_complete_filtered.fastq | samtools view -u -q 10 -@ 1 - | samtools sort -n -@ 1 -o - temp1 | samtools fixmate -r - - | samtools sort -o - temp5 | samtools mpileup -puD -f /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa - | bcftools call -O b -m -S snps - > F2_S1_L001.indels.bcf.gz &
bowtie2 -p 2 --qc-filter --local --very-sensitive-local -x /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.genome -1 F3_GCCAAT_L005_R1_complete_filtered.fastq -2 F3_GCCAAT_L005_R2_complete_filtered.fastq | samtools view -u -q 10 -@ 1 - | samtools sort -n -@ 1 -o - temp2 | samtools fixmate -r - - | samtools sort -o - temp6 | samtools mpileup -puD -f /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa - | bcftools call -O b -m -S snps - > F3_S2_L001.indels.bcf.gz &
bowtie2 -p 2 --qc-filter --local --very-sensitive-local -x /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.genome -1 FM7_CTTGTA_L002_R1_complete_filtered.fastq -2 FM7_CTTGTA_L002_R2_complete_filtered.fastq | samtools view -u -q 10 -@ 1 - | samtools sort -n -@ 1 -o - temp3 | samtools fixmate -r - - | samtools sort -o - temp7 | samtools mpileup -puD -f /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa - | bcftools call -O b -m -S snps - > FM7_S3_L001.indels.bcf.gz &
bowtie2 -p 2 --qc-filter --local --very-sensitive-local -x /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.genome -1 OreR_GTGAAA_L005_R1_complete_filtered.fastq -2 OreR_GTGAAA_L005_R2_complete_filtered.fastq | samtools view -u -q 10 -@ 1 - | samtools sort -n -@ 1 -o - temp4 | samtools fixmate -r - - | samtools sort -o - temp8 | samtools mpileup -puD -f /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa - | bcftools call -O b -m -S snps - > OreR_S4_L001.indels.bcf.gz ;
bowtie2 -p 8 --qc-filter --local --very-sensitive-local -x /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.genome -1 F5_GCCAAT_L002_R1_complete_filtered.fastq -2 F5_GCCAAT_L002_R2_complete_filtered.fastq | samtools view -u -q 10 -@ 1 - | samtools sort -n -@ 1 -o - temp10 | samtools fixmate -r - - | samtools sort -o - temp9 | samtools mpileup -puD -f /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa - | bcftools call -O b -m -S snps - > F5_S1_L001.indels.bcf.gz ;

bowtie2 -p 2 --qc-filter --local --very-sensitive-local -x /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.genome -1 F2_ACAGTG_L005_R1_complete_filtered.fastq -2 F2_ACAGTG_L005_R2_complete_filtered.fastq | samtools view -u -q 10 -@ 1 - > F2_S1_L001.sam ;
bowtie2 -p 2 --qc-filter --local --very-sensitive-local -x /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.genome -1 F3_GCCAAT_L005_R1_complete_filtered.fastq -2 F3_GCCAAT_L005_R2_complete_filtered.fastq | samtools view -u -q 10 -@ 1 - > F3_S2_L001.sam ;
bowtie2 -p 2 --qc-filter --local --very-sensitive-local -x /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.genome -1 FM7_CTTGTA_L002_R1_complete_filtered.fastq -2 FM7_CTTGTA_L002_R2_complete_filtered.fastq | samtools view -u -q 10 -@ 1 - > FM7_S3_L001.sam ; 
bowtie2 -p 2 --qc-filter --local --very-sensitive-local -x /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.genome -1 OreR_GTGAAA_L005_R1_complete_filtered.fastq -2 OreR_GTGAAA_L005_R2_complete_filtered.fastq | samtools view -u -q 10 -@ 1 - > OreR_S4_L001.sam;
bowtie2 -p 8 --qc-filter --local --very-sensitive-local -x /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.genome -1 F5_GCCAAT_L002_R1_complete_filtered.fastq -2 F5_GCCAAT_L002_R2_complete_filtered.fastq | samtools view -u -q 10 -@ 1 - > F5_S1_L001.sam ;

bowtie2 -p 8 --qc-filter --local --very-sensitive-local -x mh/Drosophila_melanogaster.BDGP5.74.genome -1 F5_S1_L001_R1_001.fastq -2 F5_S1_L001_R2_001.fastq | samtools view -u -q 10 -@ 1 - | samtools sort -n -@ 1 -o - temp10 | samtools fixmate -r - - | samtools sort -o - temp9 | samtools mpileup -puD -f mh/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa - | bcftools call -O b -m -S snps - > F5_S1_L001.indels.bcf.gz ;



#bcftools merge -O b OreR_S4_L001.indels.bcf.gz FM7_S3_L001.indels.bcf.gz > wt.bcf.gz
#bcftools merge -O b F3_S2_L001.indels.bcf.gz F2_S1_L001.indels.bcf.gz > mut.bcf.gz

#bcftools index wt.bcf.gz
#bcftools index mut.bcf.gz

#bcftools isec  mut.bcf.gz wt.bcf.gz --complement 1 -p . 
#sed '/LowQual/d' 0000.vcf > mutants_only_highqual_samtools.vcf

#samtools view -u -q 10 -@ 8 -S F2_S1_L001.sam  | samtools sort -o - temp5 > F2_S1_L001_clean.bam &
#samtools view -u -q 10 -@ 8 -S F3_S2_L001.sam  | samtools sort -o - temp6 > F3_S2_L001_clean.bam &
#samtools view -u -q 10 -@ 8 -S FM7_S3_L001.sam  | samtools sort -o - temp7 > FM7_S3_L001_clean.bam &
#samtools view -u -q 10 -@ 8 -S OreR_S4_L001.sam  | samtools sort -o - temp8 > OreR_S4_L001_clean.bam ;

#java -jar /Users/b110-mm06/Desktop/gatk-protected/dist/AddOrReplaceReadGroups.jar I=F2_S1_L001_clean.bam O=F2_S1_L001_clean_fixed.bam LB=TruSeqLT PL=Illumina PU=1 SM=F2 &
#java -jar /Users/b110-mm06/Desktop/gatk-protected/dist/AddOrReplaceReadGroups.jar I=F3_S2_L001_clean.bam O=F3_S2_L001_clean_fixed.bam LB=TruSeqLT PL=Illumina PU=1 SM=F3 &
#java -jar /Users/b110-mm06/Desktop/gatk-protected/dist/AddOrReplaceReadGroups.jar I=FM7_S3_L001_clean.bam O=FM7_S3_L001_clean_fixed.bam LB=TruSeqLT PL=Illumina PU=1 SM=FM7 &
#java -jar /Users/b110-mm06/Desktop/gatk-protected/dist/AddOrReplaceReadGroups.jar I=OreR_S4_L001_clean.bam O=OreR_S4_L001_clean_fixed.bam LB=TruSeqLT PL=Illumina PU=1 SM=OreR ;

#samtools index F2_S1_L001_clean_fixed.bam & 
#samtools index F3_S2_L001_clean_fixed.bam &
#samtools index FM7_S3_L001_clean_fixed.bam &
#samtools index OreR_S4_L001_clean_fixed.bam ;

#java -Xmx1g -jar /Users/b110-mm06/Desktop/gatk-protected/dist/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa -I F2_S1_L001_clean_fixed.bam -o F2_S1_L001.realign.intervals &
#java -Xmx1g -jar /Users/b110-mm06/Desktop/gatk-protected/dist/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa -I F3_S2_L001_clean_fixed.bam -o F3_S2_L001.realign.intervals &
#java -Xmx1g -jar /Users/b110-mm06/Desktop/gatk-protected/dist/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa -I FM7_S3_L001_clean_fixed.bam -o FM7_S3_L001.realign.intervals &
#java -Xmx1g -jar /Users/b110-mm06/Desktop/gatk-protected/dist/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa -I OreR_S4_L001_clean_fixed.bam -o OreR_S4_L001.realign.intervals ;

#java -Xmx4g -jar /Users/b110-mm06/Desktop/gatk-protected/dist/GenomeAnalysisTK.jar -I F2_S1_L001_clean_fixed.bam -R /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa -T IndelRealigner -targetIntervals F2_S1_L001.realign.intervals -o F2_S1_L001_clean_fixed_realigned.bam &
#java -Xmx4g -jar /Users/b110-mm06/Desktop/gatk-protected/dist/GenomeAnalysisTK.jar -I F3_S2_L001_clean_fixed.bam -R /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa -T IndelRealigner -targetIntervals  F3_S2_L001.realign.intervals -o F3_S2_L001_clean_fixed_realigned.bam;
#java -Xmx4g -jar /Users/b110-mm06/Desktop/gatk-protected/dist/GenomeAnalysisTK.jar -I FM7_S3_L001_clean_fixed.bam -R /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa -T IndelRealigner -targetIntervals FM7_S3_L001.realign.intervals -o FM7_S3_L001_clean_fixed_realigned.bam &
#java -Xmx4g -jar /Users/b110-mm06/Desktop/gatk-protected/dist/GenomeAnalysisTK.jar -I OreR_S4_L001_clean_fixed.bam -R /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa -T IndelRealigner -targetIntervals OreR_S4_L001.realign.intervals -o OreR_S4_L001_clean_fixed_realigned.bam;

#java -jar /Users/b110-mm06/Desktop/gatk-protected/dist/AddOrReplaceReadGroups.jar I=F2_S1_L001_clean_fixed_realigned.bam O=F2_S1_L001_clean_fixed_realigned_fixed.bam SORT_ORDER=coordinate RGID=F2_S1 RGLB=TruSeqLT RGPL=Illumina RGPU=1 RGSM=F2_S1 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT &	
#java -jar /Users/b110-mm06/Desktop/gatk-protected/dist/AddOrReplaceReadGroups.jar I=F3_S2_L001_clean_fixed_realigned.bam O=F3_S2_L001_clean_fixed_realigned_fixed.bam SORT_ORDER=coordinate RGID=F3_S2 RGLB=TruSeqLT RGPL=Illumina RGPU=1 RGSM=F3_S2 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT &
#java -jar /Users/b110-mm06/Desktop/gatk-protected/dist/AddOrReplaceReadGroups.jar I=FM7_S3_L001_clean_fixed_realigned.bam O=FM7_S3_L001_clean_fixed_realigned_fixed.bam SORT_ORDER=coordinate RGID=FM7_S3 RGLB=TruSeqLT RGPL=Illumina RGPU=1 RGSM=FM7_S3 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT &
#java -jar /Users/b110-mm06/Desktop/gatk-protected/dist/AddOrReplaceReadGroups.jar I=OreR_S4_L001_clean_fixed_realigned.bam O=OreR_S4_L001_clean_fixed_realigned_fixed.bam SORT_ORDER=coordinate RGID=OreR_S4 RGLB=TruSeqLT RGPL=Illumina RGPU=1 RGSM=OreR_S4 CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT ;

#java -jar /Users/b110-mm06/Desktop/gatk-protected/dist/GenomeAnalysisTK.jar -R /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa -T UnifiedGenotyper -I F2_S1_L001_clean_fixed_realigned_fixed.bam -o F2_S1_bwt2_gatk.bcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 500 & 
#java -jar /Users/b110-mm06/Desktop/gatk-protected/dist/GenomeAnalysisTK.jar -R /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa -T UnifiedGenotyper -I F3_S2_L001_clean_fixed_realigned_fixed.bam -o F3_S2_bwt2_gatk.bcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 500 ;
#java -jar /Users/b110-mm06/Desktop/gatk-protected/dist/GenomeAnalysisTK.jar -R /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa -T UnifiedGenotyper -I FM7_S3_L001_clean_fixed_realigned_fixed.bam -o FM7_S3_bwt2_gatk.bcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 500 &
#java -jar /Users/b110-mm06/Desktop/gatk-protected/dist/GenomeAnalysisTK.jar -R /Volumes/IMAGES/databasefiles/Drosophila_melanogaster.BDGP5.73/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa -T UnifiedGenotyper -I OreR_S4_L001_clean_fixed_realigned_fixed.bam -o OreR_S4_bwt2_gatk.bcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 500 ;

#vcftools --bcf F2_S1_bwt2_gatk.bcf --recode-to-stream | bgzip -c  >  F2_S1_bwt2_gatk.vcf.gz &
#vcftools --bcf F3_S2_bwt2_gatk.bcf --recode-to-stream | bgzip -c  >  F3_S2_bwt2_gatk.vcf.gz &
#vcftools --bcf FM7_S3_bwt2_gatk.bcf --recode-to-stream | bgzip -c  >  FM7_S3_bwt2_gatk.vcf.gz &
#vcftools --bcf OreR_S4_bwt2_gatk.bcf --recode-to-stream | bgzip -c  >  OreR_S4_bwt2_gatk.vcf.gz ;

#tabix -p vcf F2_S1_bwt2_gatk.vcf.gz
#tabix -p vcf F3_S2_bwt2_gatk.vcf.gz
#tabix -p vcf FM7_S3_bwt2_gatk.vcf.gz
#tabix -p vcf OreR_S4_bwt2_gatk.vcf.gz

#vcf-merge -d -t F2_S1_bwt2_gatk.vcf.gz F3_S2_bwt2_gatk.vcf.gz | bgzip -c  > mutants_gatk.vcf.gz &
#vcf-merge -d -t FM7_S3_bwt2_gatk.vcf.gz OreR_S4_bwt2_gatk.vcf.gz | bgzip -c  > wt_gatk.vcf.gz

#tabix -p vcf mutants_gatk.vcf.gz &
#tabix -p vcf wt_gatk.vcf.gz ;

#vcf-isec -c mutants_gatk.vcf.gz wt_gatk.vcf.gz -o -f > mutants_only_gatk.vcf

#sed '/LowQual/d' mutants_only_gatk.vcf > mutants_only_highqual_gatk.vcf

#bgzip  mutants_only_highqual_gatk.vcf 
#bgzip  mutants_only_highqual_samtools.vcf 

#tabix -p vcf mutants_only_highqual_gatk.vcf.gz
#tabix -p vcf mutants_only_highqual_samtools.vcf.gz

#vcf-merge mutants_only_highqual_samtools.vcf.gz mutants_only_highqual_gatk.vcf.gz > mutants_only_highqual_merge_samtools_gatk.vcf
#vcftools --vcf mutants_only_highqual_merge_samtools_gatk.vcf --keep-only-indels --recode-to-stream --minQ 80 > mutants_only_merged_qual80_indelonly.vcf

#java -Xmx4g -jar ~/snpEff/snpEff.jar eff -c ~/snpEff/snpEff.config -v BDGP5.69  mutants_only_merged_qual80_indelonly.vcf >  mutants_only_merged_qual80_indelonly_snpEff.vcf








