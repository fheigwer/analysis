use Bio::SeqIO;
use strict;
use warnings;
# get command-line arguments, or die with a usage statement
my $usage         = "parsefortrans_paired.pl R1 R2 U infileformat outfilebase outfileformat\n";
my $infile_R1        = shift or die $usage;
my $infile_R2        = shift or die $usage;
#my $unpairedfile       = shift or die $usage;
my $infileformat  = shift or die $usage;
my $outfilebase       = shift or die $usage;
my %valid_3_ids=();
my %valid_5_ids=();
my %invalid_ids=();
my $phredstring='!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~';
my @phredarray=split("",$phredstring);
my %phredhash=();
my %hash=();
my $match="";
my $start=0;
foreach my $n (0..104){
	$phredhash{$n}=$phredarray[$n];
}

# create one SeqIO object to read in,and another to write out
my $seq_forward = Bio::SeqIO->new(
                             -file   => "<$infile_R1",
                             -format => $infileformat
                             );
my $seq_reverse = Bio::SeqIO->new(
                             -file   => "<$infile_R2",
                             -format => $infileformat,
                             );
#my $unpaired = Bio::SeqIO->new(
#                             -file   => "<$unpairedfile",
 #                            -format => $infileformat,
 #                            );
open seq_forward_3_out , ">$outfilebase"."_3_R1.fastq"; 
open seq_reverse_3_out , ">$outfilebase"."_3_R2.fastq";                       
open seq_forward_5_out , ">$outfilebase"."_5_R1.fastq"; 
open seq_reverse_5_out , ">$outfilebase"."_5_R2.fastq";

# write each entry in the input file to the output file
my $done=0;
FORWARD: while (my $data = $seq_forward->next_dataset) {
	%hash=%$data;
	if($hash{"-seq"}=~m/.*CGTCAATTTT.*TTAA(.+)$/){
		$match=$1;$start=$-[1];
		if(length($match)>2){
			print seq_forward_3_out "\@".$hash{"-id"}.":3 ".$hash{"-desc"}."\n".$match."\n\+\n".substr($hash{"-raw_quality"},$start)."\n"; #
			$valid_3_ids{$hash{"-id"}}++;
			$done++;
		}
    }elsif($hash{"-seq"}=~m/.*GTACGTCACAAT.*TTAA(.+)$/){
    	$match=$1;$start=$-[1];
    	if(length($match)>2){
    		print seq_forward_5_out "\@".$hash{"-id"}.":5 ".$hash{"-desc"}."\n".$match."\n\+\n".substr($hash{"-raw_quality"},$start)."\n"; #
    		$valid_5_ids{$hash{"-id"}}++;
    		$done++;
    	}
    		
    }else{
    	$invalid_ids{$hash{"-id"}}++;    
    }	
   #if($done>10000){
   	#last FORWARD;
   #	}
}
$done=0;
REVERSE:while (my $data = $seq_reverse->next_dataset) {
	my %hash=%$data;
	if(exists $valid_3_ids{$hash{"-id"}}){
		print seq_reverse_3_out  "\@".$hash{"-id"}.":3 ".$hash{"-desc"}."\n".substr($hash{"-seq"},0,100)."\n\+\n".substr($hash{"-raw_quality"},0,100)."\n";
		$done++;
	}elsif(exists $valid_5_ids{$hash{"-id"}}){
		print seq_reverse_5_out  "\@".$hash{"-id"}.":5 ".$hash{"-desc"}."\n".substr($hash{"-seq"},0,100)."\n\+\n".substr($hash{"-raw_quality"},0,100)."\n";
		$done++;
	}
	#if($done>100000){
   	#	last REVERSE;
   	#}
}	
close (seq_forward_3_out);
close (seq_reverse_3_out);
close (seq_forward_5_out);
close (seq_reverse_5_out);


#open seq_U_3_out , ">$outfilebase"."_3_U.fastq"; 
#open seq_U_5_out , ">$outfilebase"."_5_U.fastq"; 

#while (my $data = $unpaired->next_dataset) {
#		my %hash=%$data;
 #   if($hash{"-seq"}=~m/.*CGTCAATTTT.+TTAA(.*)/){
  #  	print seq_U_3_out  "\@".$hash{"-id"}.":3 ".$hash{"-desc"}."\n$1\n\+\n".substr($hash{"-raw_quality"},$-[1],$+[1])."\n";
   # 	$valid_3_ids{$hash{"-id"}}++;
	#}elsif($hash{"-seq"}=~m/.*GTACGTCACAAT.+TTAA(.*)/){
    #	print seq_U_5_out "\@".$hash{"-id"}.":3 ".$hash{"-desc"}."\n$1\n\+\n".substr($hash{"-raw_quality"},$-[1],$+[1])."\n";
    #	$valid_5_ids{$hash{"-id"}}++;
    #}else{
    #	$invalid_ids{$hash{"-id"}}++;    
    #}    	
#}
#close(seq_U_3_out);
#close(seq_U_5_out);

print "3_prime count =".keys(%valid_3_ids)." 5_prime count =".keys(%valid_5_ids)." invalid id count =".keys(%invalid_ids)."\n";
