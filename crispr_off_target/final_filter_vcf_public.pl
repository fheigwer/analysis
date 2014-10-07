#!/usr/bin/perl -s -w

open $infile, shift or die $!;
while(<$infile>){
	if($_=~m/^X\s/){
		my @line=split("\t",$_);
		if($line[5]>40){#filter for QUAL greater 40
			if($line[7]=~m/DP=(\d+)/){	#filter for depth greater 40
				if($1>40){
					if($line[7]=~m/IDV=(\d+)/){ #filter for number of supportive reads
						if($1>5){
							if($line[9]=~m/([\d\.]+)[\/\|]([\d\.]+)/){  #filter to get heterozygous mutations only
								$S11=$1;$S12=$2;
								if($S11 ne $S12){
									print $_;
								}
							}
						}
					}
				}
			}
		}
	}elsif($_=~m/^#/){
		print $_;
	}
}
close $infile;
