my %hash=('0'=>"A",'1'=>"G",'2'=>"C",'3'=>"T");
open(my $rand_crispr_file ,">", $temp_dir . '/rand_crispr.fasta' );
foreach my $stu (1..10000){
	my $string="";
	
	foreach my $number (1..21){
		$string.=$hash{int(rand(3))};
	}
	print $rand_crispr_file ">random_$stu\n$string"."GG\n";
}
close($rand_crispr_file);

system( '/usr/bin/bowtie2 -p 4 -f -k 30 --'.$something{"bowtie_mode"}.' --end-to-end --no-hd -x ' . $databasepath.".genome" . ' -U ' . $temp_dir . '/rand_crispr.fasta > ' . $temp_dir . '/rand_crispr.bwt' );

my %non_hit;
open(my $rand_crispr_bwt ,"<", $temp_dir . '/rand_crispr.bwt' );
while(<$rand_crispr_bwt>){
	my @line=split("\t",$_);
	if($line[2] eq "*"){
		$non_hit{$line[0]}++;
	}	
}
close($rand_crispr_bwt);

open(my $rand_crispr_file ,"<", $temp_dir . '/rand_crispr.fasta' );
open(my $rand_crispr_filtered_file ,">", $temp_dir . '/rand_crispr_filtered.fasta' );
while(<$rand_crispr_file>){
	if($_=~m/^>(.+)$/){
		if(exists($non_hit{$1})){
			$allowed=1;
		}
	}elsif($allowed==1){
		chomp $_;
		$_=~s/.GG$//;
		if($crisprcount<$something{"random_seqs"}){
			print $rand_crispr_filtered_file ">$1\n$_\n";
		}
	}
}
close($rand_crispr_file);
close($rand_crispr_filtered_file);

unlink $temp_dir . '/rand_crispr.fasta';
unlink $temp_dir . '/rand_crispr.bwt';
