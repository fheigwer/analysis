#!/usr/bin/perl

open $infile, shift;
open($outfile, ">", shift ) or die $!;
%hash=();
while(<$infile>){
    if($_=~m/.*gene_id \"(.+?)\"\;.*?gene_name \"(.+?)\"\;/){
        $1=~s/\'//ig;
        $2=~s/\'//ig;
        $hash{$1}=$2;
    }
}

print $outfile 'library=['."\n";
foreach $name (keys %hash){
    print $outfile '{ENS_ID: \''.$name.'\',GENE_SYMBOL: \''.$hash{$name}.'\'},'."\n";
}
print $outfile '{}];'."\n";
close $infile;
close $outfile;