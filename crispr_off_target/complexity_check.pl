#!/usr/bin/perl -s -w
my @test_strings=("AACAAACAATAAGAACAAGAA","ATATATATATATA","GTCGTCGTCGTCGTC");

foreach my $test (@test_strings){
	complexity_check($test);
}



sub complexity_check {
	my $string=$_[0];
	if($string=~m/(\w{1,4})(\w{1,4})(\w{1,4})(\w{1,4})/){
		print "$1\n";
	}
}

