#!/usr/bin/perl

$word="AGAGAGAGA";
$count=0;
foreach my $first ("A","G","T","C"){
	foreach my $second ("A","G","T","C"){
		my @matches = $word =~ m/($first$second)/g;
		my $count = scalar(@matches);
		if($count>=4){
			print "is not\n";	
		}
	}
}

