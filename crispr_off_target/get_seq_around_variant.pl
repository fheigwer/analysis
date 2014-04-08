#!/usr/bin/perl
use strict;
use warnings;
use Bio::DB::Fasta;
my $infilename=shift;
my $db = Bio::DB::Fasta->new( '/Volumes/IMAGES/databasefiles/ensemble_fasta_files/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa', -makeid => \&make_my_id );
open my $infile , "$infilename";
	while (<$infile>){
		if($_=~m/^\w/){
			my @line=split("\t",$_);
			my $obj        = $db->get_Seq_by_id($line[0]);
			my $target = substr $obj->seq, ($line[1]-1000) , (2000);
   		   	print ">$line[0]:($line[1]-1000)..($line[1]+1000):$line[3]:$line[4]\n$target\n";
  		}
      }
close $infile;
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
#########################################################################################
#name:      make_my_id
#function:  return searchable ids for generation of the FASTA index
#input:     (FASTA header)
#output:    array of searchable ids (strings)
#########################################################################################
sub make_my_id {
      $_[0] =~m/^>(\S+) locus_tag= (\S+);/;
      return ( $1, $2);
}