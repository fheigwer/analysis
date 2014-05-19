open (my $outfiletab_filtered, ">", $temp_dir . "/libraryfile.tab");
open (my $outfiletab, "<", $temp_dir . "/all_results_together.tab");

my @file=(<$outfiletab>);
while(<$outfiletab>){
	push @file, $_;
}
close $outfiletab;

print $outfiletab_filtered $file[0];
my $current_size=0;
my $last_gene="";
my @new_file=();
my $limit=0;
if(($something{"library_size"}-$something{"random_seqs"})<= scalar(@file)){
	$limit=$something{"library_size"}-$something{"random_seqs"};
}else{
	$limit=(scalar(@file)-1);
}
while($current_size<($something{"library_size"}-$something{"random_seqs"})){
	my $count=0;
	foreach my $element (@file){
		$element=~m/^(.+)_\d+_\d+/;
		my  $current_gene=$1;
		if($current_gene ne $last_gene){
			push @new_file, $element;
			$current_size++;
			$last_gene=$current_gene;
			splice(@file,$count,1);			
		}else{
			next;
		}
		$count++;
	}
	$last_gene="";
}
foreach my $element (@new_file){
		print $outfiletab_filtered $element;
	}
close $outfiletab_filtered;	