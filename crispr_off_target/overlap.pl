#!/usr/bin/perl
use strict;
use warnings;
use Set::IntervalTree;
 my $tree = Set::IntervalTree->new;
  

open my $file1, "<".shift;
open my $file2, "<".shift;
open my $file3, "<".shift;

my %positions=();

while(<$file2>){
	if(!($_=~m/#/)){
		my @temp=split("\t",$_);
		#$positions{$temp[1]}++;
		if(length($temp[3])>length($temp[4])){
			$tree->insert($temp[3],$temp[1],($temp[1]+length($temp[3])));
		}else{
			$tree->insert($temp[3],$temp[1],($temp[1]+length($temp[4])));
		}
	}
}
while(<$file3>){
	if(!($_=~m/#/)){
		my @temp=split("\t",$_);
		#$positions{$temp[1]}++;
		if(length($temp[3])>length($temp[4])){
			$tree->insert($temp[3],$temp[1],($temp[1]+length($temp[3])));
		}else{
			$tree->insert($temp[3],$temp[1],($temp[1]+length($temp[4])));
		}
	}
}

while(<$file1>){
	if(!($_=~m/#/)){
		my @temp=split("\t",$_);
		my $results ="";
		if(length($temp[3])>length($temp[4])){
			$results = $tree->fetch($temp[1],($temp[1]+length($temp[3])));
		}else{
			$results = $tree->fetch($temp[1],($temp[1]+length($temp[4])));
		}
  		#print." intervals found.\n";
		#foreach my $tol (-20..20){
			#if(exists($positions{$temp[1]+$tol})){
			#	$count++;
		#	}
		#}
		if( scalar(@$results)>0){
			print $_;
		}
	}else{
		print $_;
	}
}

close $file1;
close $file2;
close $file3;
