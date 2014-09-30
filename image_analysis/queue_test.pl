my %Q=();
foreach my $letter ("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"){		
        foreach my $number (1..24){
                foreach my $field (1..4){
                        my @temp=($letter,$number,$field);
                        $Q{$letter.$number."_".$field} = \@temp;
                }
        }
}
        
        
while(%Q){
    foreach my $key (sort keys(%Q)){
        print "$key\n";
        delete $Q{$key};
    }
}
print "done\n";