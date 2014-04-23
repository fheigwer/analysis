$infile=shift;
%hash=();

open INFILE , "<".$infile;

while(<INFILE>){
	my @line=split("\t",$_);
	my @gene=split("_",$line[0]);
	$hash{$gene[0]}++; #."\n";
}

foreach my $key (sort keys %hash){
	print $key."\t".$hash{$key}."\n";
}
close INFILE;
