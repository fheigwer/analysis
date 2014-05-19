#!/usr/bin/perl
use warnings;
use strict;
my $filename=shift;
open INFILE, $filename;
open OUTFILE, ">".$filename."corrected.fa";

foreach my $line (<INFILE>){
	$line=~s/^>(\S+).*gene:(\S+).*$/>$2 transcript:$1/;
	print OUTFILE $line;
}
close INFILE;
close OUTFILE;
#unlink $filename;
exit;
