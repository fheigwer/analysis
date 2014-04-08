#!/usr/bin/perl -s -w
 use Image::Info qw(image_info image_type);

 if(-e "B8fld1wvCy5Cy5.tif"){
 	if(!(-z "B8fld1wvCy5Cy5.tif")){
 			print "blub\n";
 	}
 }
#my %Q=();
#foreach my $letter ("A","B","C"){
#	foreach my $number (1..4){
#		my @temp=($letter,$number);
#		$Q{$letter.$number} = \@temp;
#	}
#}
#while(%Q){
#	foreach my $key (sort keys(%Q)){
#		print $key."\n";
#		
#		delete $Q{$key};
#	}
	#sleep 1;
#}
#my $letter="A";
#my $number="1";
#my $count=0;
#my @names=();
#my @result=();
#open $result_tab , "A3.tab";
#	while(<$result_tab>){
#		chomp $_;
#		if($count<=0){
#			@names=split("\t",$_);
#			$count++;
#		}elsif($_=~m/^\d/){
#			@result=split("\t",$_);
#		}
		
#	}
#close $result_tab;

#print join("_",@names)."\n".join("_",@result)."\n";