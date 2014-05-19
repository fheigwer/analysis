#!/usr/bin/perl
use strict;
use warnings;
use Bio::Tools::GFF;
use Set::IntervalTree;

my $file=$ARGV[0]; 
use Cwd;
my $dir = getcwd;
opendir INDIR, $dir;

foreach my $infile (readdir(INDIR)){
    if($infile=~m/(\S+)\.gff/){
	my $outfile=$1."_indexed.mygff";
	my %chroms=();
	my %trees=();
	my $file=$infile;
	open INFILE ,$file;
	open OUTFILE ,">".$outfile;
	my $count="";
	while(<INFILE>) {
		my $line=$_;
		chomp $line;
		my @line= split("\t", $line);
		my $chrom=$line[0];
		$chroms{$chrom}++;
		if($chroms{$chrom}==1){
			$trees{$chrom}=Set::IntervalTree->new;
		}else{
			my $gene_name="";						
			my $start=$line[3];
			my $end=$line[4];
		    	if($line[2]=~/gene/ ){
				if($line=~/.*\s(\S+)\slocus_tag= (\S+);.*/){
					my $temp="gene_".$1."::".$2;
					print OUTFILE $temp,"\t",$start,"\t",$end,"\n";
				}
			}
		}
	}
	close(INFILE);
	close(OUTFILE);
	}
}
close INDIR;

open INFILE, $file.".gtf";
my %chrom=();
my $count=0;
foreach my $line (<INFILE>){
    chomp $line;
    my @line=split("\t",$line);
    $chrom{$line[0]}++;
    
    if($chrom{$line[0]}==1){
        if ($count==0) {
            open OUTFILE, ">>".$line[0]."_indexed.mygff";
        }else{
            close OUTFILE;
            open OUTFILE, ">>".$line[0]."_indexed.mygff";
        }
    }else{
        if($line=~m/start_codon/ || $line=~m/stop_codon/ || $line=~m/CDS/ || $line=~m/exon/){ 
            if($line=~m/transcript_id \"(\S+)\"\; exon_number \"(\S+)\"\;/){            
                print OUTFILE $line[2]."::".$1."::".$2."\t".$line[3]."\t".$line[4]."\n";
            }
        }  
    }
    $count++;
}
close OUTFILE;
