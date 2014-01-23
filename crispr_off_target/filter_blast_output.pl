#!/usr/bin/perl -s -w

open $infile, shift;
while(<$infile>){
	if($_=~m/^\w/){
		my @line=split("\t",$_);
		if($line[7]>22){
			if($line[8]>900 && $line[8]<1100){
				print $_;							
			}
		}
	}
}
close $infile;
		