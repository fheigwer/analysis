#!/usr/bin/perl

use strict;
use warnings;

my $name="";
my $length="";
my $coverage="";

opendir INDIR, shift;
foreach my $file (readdir(INDIR)){
	if ($file=~m/(.*)_stats.tab/){
		open OUTFILE, ">$1.txt";
		open INFILE, $file;
		foreach my $line (<INFILE>){
			if($line=~m/(\S+)   Query length: (\d+)/){
				$name=$1;
				$length=$2;
			}elsif($line=~m/Number of successful designs = (\d+)/){
				$coverage=$1;
				print OUTFILE "$name\t$length\t$coverage\n";
			}
		}
		close OUTFILE;
		close INFILE;
	}
}
closedir INDIR;
 
