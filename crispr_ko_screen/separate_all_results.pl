#!/usr/bin/perl
use strict;
use warnings;
my %result=();
my $lastexon="";
my $lastscore="";
my $gene="";
my $lastcpg="";
my %controls=();
open CONTROLS, "controls.tab";

foreach my $line (<CONTROLS>){
	chomp $line;
	$controls{$line}++;
}
close CONTROLS;

open INFILE, shift;
open OUTFILE, ">out.txt";
open OUTFILECONTROLS, ">out_controls.txt";
my @file=(<INFILE>);
	foreach my $line (@file){
	   my @line=split("\t",$line);
			if(!($line[6]=~m/_/)){
			my $curgene=$line[6];
			if(!exists($result{$curgene})){$result{$curgene}=0;}
			if($line[0]=~m/([a-zA-Z0-9\.\-]+)_/){
				 $gene=$1;
			}else{
				 $gene="";
			}
			my $curexon=$line[8];
			my $curcpg=$line[9];
			my $curscore=$line[11];
			if(exists $controls{$gene}){
				if($curexon ne $lastexon && $result{$curgene}< 16){
					if($curcpg ne $lastcpg && $result{$curgene}< 16){
						if($curscore ne $lastscore && $result{$curgene}< 16){
							print OUTFILE $line;
							$result{$curgene}++;
						}
					}
				}
			}else{
				if($result{$curgene}< 16 && $curexon ne "NA"){
					if($result{$curgene}< 16 && $curcpg eq "NA"){
						if($curscore ne $lastscore && $result{$curgene}< 16){
							print OUTFILE $line;
							$result{$curgene}++;
						}
					}
				}
			}
			$lastexon=$curexon;
			$lastscore=$curscore;
			$lastcpg=$curcpg;
			}
		}
		foreach my $line (@file){
		   my @line=split("\t",$line);
			if(!($line[6]=~m/_/)){
			my $curgene=$line[6];
			if(!exists($result{$curgene})){$result{$curgene}=0;}			
			if($line[0]=~m/([a-zA-Z0-9\.\-]+)_/){
				 $gene=$1;
			}else{
				 $gene="";
			}
			my $curexon=$line[8];
			my $curcpg=$line[9];
			my $curscore=$line[11];
			if(exists $controls{$gene}){
				if($curexon ne $lastexon && $result{$curgene}< 16){
					if($curcpg ne $lastcpg && $result{$curgene}< 16){
						if($curscore ne $lastscore && $result{$curgene}< 16){
							print OUTFILE $line;
							$result{$curgene}++;
						}
					}
				}
			}else{
				if($result{$curgene}< 16 && $curexon ne "NA"){
					if($result{$curgene}< 16){
						if($curscore ne $lastscore && $result{$curgene}< 16){
							print OUTFILE $line;
							$result{$curgene}++;
						}
					}
				}
			}
			$lastscore=$curscore;
			$lastcpg=$curcpg;
			}
		}
		foreach my $line (@file){
		   my @line=split("\t",$line);
			if(!($line[6]=~m/_/)){
			my $curgene=$line[6];
			if(!exists($result{$curgene})){$result{$curgene}=0;}
			if($line[0]=~m/([a-zA-Z0-9\.\-]+)_/){
				 $gene=$1;
			}else{
				$gene="";
			}
			my $curexon=$line[8];
			my $curcpg=$line[9];
			my $curscore=$line[11];
			if(exists $controls{$gene}){
				if($curexon ne $lastexon && $result{$curgene}< 16){
					if($curcpg ne $lastcpg && $result{$curgene}< 16){
						if($curscore ne $lastscore && $result{$curgene}< 16){
							print OUTFILE $line;
							$result{$curgene}++;
						}
					}
				}
			}else{
				if($result{$curgene}< 16 && $curexon ne "NA"){
					if($result{$curgene}< 16){
						print OUTFILE $line;
						$result{$curgene}++;						
					}
				}
			}
			$lastexon=$curexon;
			$lastscore=$curscore;
			$lastcpg=$curcpg;
			}
		}
		foreach my $line (@file){
		   my @line=split("\t",$line);
			if(!($line[6]=~m/_/)){
			my $curgene=$line[6];
			if(!exists($result{$curgene})){$result{$curgene}=0;}			
			my $curexon=$line[8];
			my $curcpg=$line[9];
			my $curscore=$line[11];
			if($result{$curgene}< 16){
				print OUTFILE $line;
				$result{$curgene}++;					
			}			
			$lastexon=$curexon;
			$lastscore=$curscore;
			$lastcpg=$curcpg;
			}
		}



		close INFILE;	
		close OUTFILE;
		close OUTFILECONTROLS;


foreach my $key (sort(keys(%result))){
	if($result{$key}>0){
		print $key."\t".$result{$key}."\n";	
	}
}
