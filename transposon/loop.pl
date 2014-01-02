my $count=0;
foreach my $G  (0..100){
	foreach my $C  (0..100){
		foreach my $T  (0..100){
			foreach my $A  (0..100){
				if(($A+$T+$C+$G)==100){
					$count++;
				}
			}
		}
	}
}
print $count."\n";