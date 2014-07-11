my %hash=('0'=>"A",'1'=>"G",'2'=>"C",'3'=>"T");
open(my $rand_crispr_file ,">", 'rand_crispr.fasta' );
foreach my $stu (1..100000){
	my $string="";	
	foreach my $number (1..20){
		$string.=$hash{int(rand(3))};
	}
	print $rand_crispr_file ">random_$stu\n$string"."AGG\n";
	print $rand_crispr_file ">random_$stu\n$string"."TGG\n";
	print $rand_crispr_file ">random_$stu\n$string"."CGG\n";
	print $rand_crispr_file ">random_$stu\n$string"."GGG\n";
}
close($rand_crispr_file);

system( 'bowtie /Volumes/IMAGES/DATABASEFILES/Homo_sapiens.GRCh37.73/Homo_sapiens.GRCh37.73.genome -S -5 3 -v 2 -f -p 8 rand_crispr.fasta  > rand_crispr.bwt' );

my %non_hit;
open(my $rand_crispr_bwt ,"<", 'rand_crispr.bwt' );
while(<$rand_crispr_bwt>){
	my @line=split("\t",$_);
	if($line[2] eq "*"){
		$non_hit{$line[0]}++;
	}	
}
close($rand_crispr_bwt);

open(my $rand_crispr_file ,"<", 'rand_crispr.fasta' );
my $allowed=0;
my $crisprcount=0;
my $name="";
my %fin=();
open $libraryfasta, ">random_seq.fasta";
foreach my $line (<$rand_crispr_file>){
	if($line=~m/^>(.+)$/){
		$name=$1;
		if(exists($non_hit{$name}) && $non_hit{$name}==4){		  
			$allowed=1;
		}
	}elsif($allowed==1){
		chomp $line;
		$line=~s/.GG$//;
		$line=~s/^./G/;
		if(($crisprcount<2000) && !exists($fin{$name})){
			print $libraryfasta ">".$name."\n"."CTGAGCTCATAGAAGACCTCACC".$line."GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTGGGTCTTCGTTCG\n";
			$fin{$name}++;
			$crisprcount++;
		}
		$allowed=0;
	}
}
close($rand_crispr_file);
close($libraryfasta);