#!/usr/bin/perl
use strict;
use warnings;
#use Bio::Restriction::Analysis;
#use Bio::PrimarySeq;
#use Parallel::ForkManager; #important package to enable mutlithreading of the script
my $file=shift;
open INFILE, $file;
#my $pm = new Parallel::ForkManager(8);
foreach my $line (<INFILE>){
	#$pm->start and next;
	my @line=split("\t",$line);
	my $seq=$line[5];
	$seq=~s/.GG$//;
	#my $seqo = Bio::PrimarySeq->new
      	#(-seq =>"CACC".$seq."GTT",
      	##-primary_id => 'synopsis',
       #	-molecule => 'dna');
	#my $ra = Bio::Restriction::Analysis->new(-seq=>$seqo);
	my $cuts = 0;
	$cuts = $cuts+(length($seq=~/GAAGAC/));
	$cuts = $cuts+(length($seq=~/GAATTC/));
	$cuts = $cuts+(length($seq=~/CTTAAG/));
	$cuts = $cuts+(length($seq=~/CAATTG/));
	$cuts = $cuts+(length($seq=~/GTTAAC/));
	$cuts = $cuts+(length($seq=~/CTCGAG/));
	$cuts = $cuts+(length($seq=~/GAGCTC/));
	if (!($cuts > 0)){print $line};
	#$pm->finish();
}
#$pm->wait_all_children();
#print "DONE!\n";
close INFILE;
