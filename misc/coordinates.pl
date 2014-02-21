$test="";
$x=50;
$y=50;
foreach my $number (1..24){
	foreach my $letter ("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"){	
			$test.=$letter.$number."_".$x."_".$y."::";
			$y+=55;
			
		
	}
	$x+=55;
	$y=50;
}

print $test;
#my %coords=();

#foreach $element (split("::",$test)){
#	my @temp=split("_",$element);
#	my $coord=shift(@temp);
#	$coords{$coord}=\@temp;
#}

#foreach $key (keys(%coords)){
#	print $key." ".join("_",@{$coords{$key}})."\n";
#}