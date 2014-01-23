#!/usr/bin/perl -s -w

open $infile, shift;
while(<$infile>){
	if($_=~m/^X/){
		my @line=split("\t",$_);
		if($line[5]>40){
			if($line[7]=~m/DP=(\d+)/){
				if($1>4){
					if($line[9]=~m/([\d\.]+)[\/\|]([\d\.]+)/){
						$S11=$1;$S12=$2;
						if($line[10]=~m/([\d\.]+)[\/\|]([\d\.]+)/){
						$S21=$1;$S22=$2;
							#if($S11 ne $S12 || $S21 ne $S22 ){
								print $_;
							#}
						}
					}
				}
			}
		}
	}
}
close $infile;
		