#!/usr/bin/perl
use strict;
use warnings;

#my $string="ATGCTAGCAGCTGATGCTGCTGCTAGCTGCTGATGCTAGTGGGGATGCGTAAGATGCTAGATGCTAGTGGTAGTGCTAGCTAGCTGCATGATCGTACTACGACTG";

#my $string="ATGCTAATGGATCGAATGCTAATGGATCGATCGTAGCTATGCGACGATGCTAATGGATCGAATGCTAATGGATCGATCGTAGCTATGCGACGACATCGACGATCGTAGCTATGCGACGACATCGACGAATGCTAATGGATCGAATGCTAATGGATCGATCGTAGCTATGCGACGACATCGACGATCGTAGCTATGCGACGACATCGACGAACATCGACGATCGTAGCTATGCGACGACATCGACGA";
	
foreach my $stu (1..100000){
	my $string="";
	my %hash=('0'=>"A",'1'=>"G",'2'=>"C",'3'=>"T");
	foreach my $number (1..1000){
		$string.=$hash{int(rand(3))};
	}
	#print "string done\n";
	my $position=300;
	my $threshold=20;
	my $length=5;
	my %postitions;
	%{$postitions{"A"}}=make_pos_index(\$string,"A");
	%{$postitions{"G"}}=make_pos_index(\$string,"G");
	%{$postitions{"C"}}=make_pos_index(\$string,"C");
	%{$postitions{"T"}}=make_pos_index(\$string,"T");
	print score_micro_homolgy(\%postitions,$threshold,$position,$length,\$string)."\n";
}





sub score_micro_homolgy {
	my $score=0;
	my $count=0;
	my $lengthhom=1;
	my $seq;
	#my $position=$_[2];
	#my $threshold=$_[1];
	#my $length=$_[3];
	#my $postitions=%{$_[0]};
	my $stuff=0;
	my $outframe=1;
	my $inframe=1;
	RIGTHSTART:foreach my $rightarmstart ($_[2]..($_[2]+$_[1])){		
		LENGTH:foreach my $length (2..$_[3]){
			if($lengthhom>0){
				$lengthhom=0;
				my @right_seq=split("",substr(${$_[4]},$rightarmstart,$length));
				$stuff=$length+1;
				LEFTARM: while($stuff<$_[1]){
					$stuff++;
					$count=0;
					$seq="";
					foreach my $letter (@right_seq){
						if(exists(${${$_[0]}{$letter}}{($_[2]-$stuff+$count)})){
							$seq.=$letter;
							$count++;												
						}else{
							next LEFTARM;					
						}
					}
					if($count==scalar(@right_seq)){
						$lengthhom=1;
						my $gaplength=($stuff+($rightarmstart-$_[2]));
						if($gaplength%3 == 0){
							$inframe=$inframe+($count*exp(0.1*(-$gaplength)));
						}else{
							$outframe=$outframe+($count*exp(0.1*(-$gaplength)));
						}
						#print "(@right_seq): $seq ".scalar(@right_seq)." ".$count." $stuff $position $rightarmstart hit"."\n";
						#print substr($string,($position-20),40)."\n";
						#print substr($string,($position-20),(20-$stuff+$count));
						#print '-' x ($stuff+($rightarmstart-$position-$count)); 
						#print substr($string,($rightarmstart),20-($rightarmstart-$position))."\n";
					}else{
						next LEFTARM;
					}
				}	
			}else{
				$lengthhom=1;
				next RIGTHSTART;		
			}	
		}
	}
	return log($outframe/$inframe);
}

sub make_pos_index {
      my %pos     = ();
      my $result  = index( ${$_[0]}, $_[1], 0);
      while ( $result != -1 ) {
            $pos{$result}++;
            $result = index( ${$_[0]}, $_[1], ($result + 1) );
      }
      return %pos;
}