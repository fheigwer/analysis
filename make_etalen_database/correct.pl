#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my $file =shift;
my $seqio_obj=Bio::SeqIO->new(-file=>$file ,-format=>"fasta");
while(my $seq_obj=$seqio_obj->next_seq){
	open OUTFILE, ">".$seq_obj->display_id().".fa";
	my $description=$seq_obj->description();
	my $seq=$seq_obj->seq();
	my $chrom=$seq_obj->display_id();
	my $length=$seq_obj->length();
	my $number=int($length/20);
	my $cut=0;
	while( $cut < $length ){
		print ">".$chrom." ;chrom:".$chrom.":".$cut."..".($cut+$number)."\n";
		print OUTFILE ">".$chrom." ;chrom:".$chrom.":".$cut."..".($cut+$number)."\n".substr($seq , $cut , $number )."\n";
		$cut=$cut+$number;
	}	
}
close OUTFILE;
