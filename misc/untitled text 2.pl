#!/usr/bin/perl
use strict;
use warnings;
#use Bio::Restriction::Analysis;
#use Bio::PrimarySeq;
#use Parallel::ForkManager; #important package to enable mutlithreading of the script
open (my $outfiletab, "<", $temp_dir . "/all_results_together.tab");
open (my $libraryfasta, ">", $temp_dir . "/library.fasta");
foreach my $line (<$outfiletab>){
	my @line=split("\t",$line);
	my $seq=$line[5];
	$seq=~s/\w{3}$//;
	my $cuts = 0;
	$cuts = $cuts+(length($seq=~/GAAGAC/));
	$cuts = $cuts+(length($seq=~/GAATTC/));
	$cuts = $cuts+(length($seq=~/CTTAAG/));
	$cuts = $cuts+(length($seq=~/CAATTG/));
	$cuts = $cuts+(length($seq=~/GTTAAC/));
	$cuts = $cuts+(length($seq=~/CTCGAG/));
	$cuts = $cuts+(length($seq=~/GAGCTC/));
	if (!($cuts > 0)){
		print $libraryfasta ">".$line[0]."\n"."CTGAGCTCATAGAAGACCTCACC".$seq."GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTGGGTCTTCGTTCG\n";
		};
}
close $outfiletab;
close $libraryfasta;
