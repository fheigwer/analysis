#!/usr/bin/perl
use strict;
use warnings;

my $file=shift;
open INFILE, $file;
foreach my $line (<INFILE>){
	my @line=split("\t",$line);
	my $seq=$line[5];
	$seq=~s/.GG$//;
	print ">".$line[0]."\n"."CTGAGCTCATAGAAGACCTCACC".$seq."GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTGGGTCTTCGTTCG\n";
	#print ">".$line[0]."\n".$seq."\n";
}

close INFILE;
	
