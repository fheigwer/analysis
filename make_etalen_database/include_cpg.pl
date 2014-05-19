#!/usr/bin/perl
use strict;
use warnings;

opendir INDIR, shift;
chdir INDIR;
foreach my $file (readdir(INDIR)){
#    print $file."\n";
    if ($file=~m/(.*)\.csv/) {
	open OUTFILE, ">>".$1."_indexed.mygff";
	open INFILE, $file;
	foreach my $line (<INFILE>){
	    my @line=split(",",$line);
	    print OUTFILE "CpG_".$line[1]."_".$line[2]."\t".$line[1]."\t".$line[2]."\n";
	}
	close INFILE;
	close OUTFILE; 
    }
}
closedir INDIR;
	
