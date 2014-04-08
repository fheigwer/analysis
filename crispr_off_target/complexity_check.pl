#!/usr/bin/perl -s -w
my @test_strings=("AACAAACAATAAG","ATATATATATATA","GTCGTCGTCGTCG","AAAAAAAAAAAAA");
foreach my $test (@test_strings){
	if(!is_repetitive($test)){
		print $test."\n";
	}
}

sub is_repetitive {
	my $test=$_[0];
	my %randletter=("0"=>"A","1"=>"G","2"=>"C","3"=>"T");
	my $stringo="";
	foreach (1..length($test)){
		$stringo.=$randletter{int(rand(3))};
	}
	my $stringi="";
	foreach (1..length($test)){
		$stringi.=$randletter{"1"};
	}
	my $randomscore=(complexity_check($stringo,2)+complexity_check($stringo,3)+complexity_check($stringo,4))/3;
	my $monoscore=(complexity_check($stringi,2)+complexity_check($stringi,3)+complexity_check($stringi,4))/3;
	my $result=(complexity_check($test,2)+complexity_check($test,3)+complexity_check($test,4))/3;
	if($result<=(($randomscore+$monoscore)/2)){
		return 1;
	}else{
		return 0;
	}
}
sub complexity_check {
	my $string=$_[0];
	my $k=$_[1];
	my %result=();
	my $temp=0;
	foreach(0..(length($string)-$k)){
		$offset=$_;
		while(($offset+$k)<=(length($string))){			
			$result{substr($string,$offset,$k)}++;
			$offset=$offset+$k;
		}
	}	
	foreach my $key (sort keys %result){
		$temp=$temp+$result{$key};
	}
	my $number= (keys %result);
	return $temp*$number;
}

