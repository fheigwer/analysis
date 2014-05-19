#!/usr/bin/perl
use strict;
use warnings;
my $file=shift;
open INFILE, $file;
foreach my $line (<INFILE>){
	my @line=split("\t",$line);
	print "\>".$line[0]."\n".$line[5]."\n";	
}
close INFILE;
