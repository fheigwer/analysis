#!/usr/bin/perl
use strict;
use warnings;
my $file=shift;
open INFILE, $file;
my $name=0;
my $length=0;
my $loc=0;
my $poss=0;
my $succ=0;
my $invar=0;
my $unspec=0;
my $spec=0;
print "NAME\tLENGTH\tLocation\tPossible\tSuccessful\tSpecific\tInvariable\tUnspecific\n";
foreach my $line (<INFILE>){
	chomp $line;
	if($line=~m/Query name: (\S+)/){
		print "$name\t$length\t$loc\t$poss\t$succ\t$spec\t$invar\t$unspec\n";
		$line=~m/Query name: (\S+)/;		
		$name=$1;
		$line=~m/Query length: (\d+)/;
		$length=$1;
		$line=~m/Query location: (\S+)/;
		$loc=$1;
		$poss=0;
		$succ=0;
		$invar=0;
		$unspec=0;
		$spec=0;
		
	}else{
		if($line=~m/possible designs = (\d+)/){
			$poss=$1;
		}
		if($line=~m/successful designs = (\d+)/){
			$succ=$1;
		}
		if($line=~m/specific target = (\d+)/){
			$spec=$1;
		}
		if($line=~m/to invariable = (\d+)/){
			$invar=$1;
		}
		if($line=~m/or none = (\d+)/){
			$unspec=$1;
		}
	}
	
}
close INFILE;
	
