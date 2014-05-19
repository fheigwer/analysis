#!/usr/bin/perl
use strict;

use Bio::Tools::GFF;
use Bio::SeqIO;
use Parallel::ForkManager;
my ($seqfile) = @ARGV;
die("must define a valid seqfile to read") unless ( defined $seqfile && -r $seqfile);
$pm = Parallel::ForkManager->new(8);
my $seqio = new Bio::SeqIO(-format => 'genbank', -file => $seqfile);
my $count = 0;
while( my $seq = $seqio->next_seq ) {
    my $pid = $pm->start and next;
    $count++;
    # defined a default name
    my $fname = $seq->display_id;
	open SEQOUT, ">$fname.fasta";
	open GFFOUT, ">$fname.gff";
	    foreach my $feature ( $seq->top_SeqFeatures() ) {		
		if($feature->primary_tag()=~m/gene/){			
			my $chrom=$fname.":".$feature->start()."..".$feature->end();
			my $name="";
			for my $tag ($feature->get_all_tags) {         
	     			 for my $value ($feature->get_tag_values($tag)) {
					if($tag eq "gene"){                
	      		   			$name.="$value ";
					}else{  
						$name.="$tag= $value;";
					}           
	  			    }		          
	  		}
			print GFFOUT "$fname\tENSEMBL\tgene\t".$feature->start()."\t".$feature->end()."\t.\t+\t.\t$name\n";
			$name.="chrom:$chrom";
			$name=~s/>gene= (.+);/>$1 /;  		
			print SEQOUT ">$name\n".substr($seq->seq(),$feature->start(),($feature->end()-$feature->start()))."\n";					
		}
	    }
	close SEQOUT;
	close GFFOUT;
    $pm->finish;
}
$pm->wait_all_children;
#unlink $seqfile;
exit;
