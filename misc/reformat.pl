#!/usr/bin/perl

use strict;
use warnings;

#my $infile="/Volumes/UNTITLED/flo_IC50.TXT";
my $last="A01";
my %hash;
open my $infile, "<", "/Volumes/UNTITLED/flo_IC50.TXT";

while(<$infile>){
	my $line = $_;
	chop $line;
	my @line=split("\t",$line);
	if($line[1]=~m/^(\w)(\w+)/){
		${$hash{$1}}{$2}=$line[2];
		${$hash{$1}}{$2}=~s/\s+//;
	}
}
close($infile);

foreach my $key (sort keys %hash){
	foreach my $subkey (sort keys $hash{$key}){		
		if($subkey<24){
			print ${$hash{$key}}{$subkey}."\t";
		}else{
			print ${$hash{$key}}{$subkey}."\n";
		}
	}
}

