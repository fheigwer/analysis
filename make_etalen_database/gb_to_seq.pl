#!/usr/bin/perl
use strict;

use Bio::Tools::GFF;
use Bio::SeqIO;

my ($seqfile) = @ARGV;
die("must define a valid seqfile to read") unless ( defined $seqfile && -r $seqfile);

my $seqio = new Bio::SeqIO(-format => 'genbank', -file => $seqfile);
while( my $seq = $seqio->next_seq ) {
#   $count++;
    # defined a default name
    my $fname = $seq->display_id;
	open SEQOUT, ">$fname.fasta";
	print SEQOUT ($seq->seq())."\n";					
	close SEQOUT;
}
#unlink $seqfile;
exit;
