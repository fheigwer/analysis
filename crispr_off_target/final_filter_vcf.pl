#!/usr/bin/perl -s -w

open $infile, shift or die $!;
while(<$infile>){
	if($_=~m/^X\s/){
		my @line=split("\t",$_);
		if($line[5]>40){#filter for QUAL greater 40
			if($line[7]=~m/DP=(\d+)/){	#filter for depth greater 4
				if($1>4){
					if($line[9]=~m/([\d\.]+)[\/\|]([\d\.]+)/){  #filter to get heterozygous mutations only
						$S11=$1;$S12=$2;
						if($S11 ne $S12){
							if($line[3]!~m/AAAAA+/gi && $line[4]!~m/AAAAA+/gi){ #filter out sequences with homopolymers
								if($line[3]!~m/CCCCC+/gi && $line[4]!~m/CCCCC+/gi){
									if($line[3]!~m/GGGGG+/gi && $line[4]!~m/GGGGG+/gi){
										if($line[3]!~m/TTTTT+/gi && $line[4]!~m/TTTTT+/gi){
											if(!is_repetitive($line[3]) && !is_repetitive($line[4])){ #filter out low complex repetitive sequences
												print $_;												
											}
										}
									}
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

sub is_repetitive {
	my $test=$_[0];
	my $count=0;
	foreach my $first ("A","G","T","C"){
		foreach my $second ("A","G","T","C"){
			my @matches = $test=~ m/($first$second)/g;
			if(scalar(@matches)>$count){
				$count=scalar(@matches);
			}			
		}
	}
	if($count>=4){
		return 1;
	}else{
		return 0;
	}
}
