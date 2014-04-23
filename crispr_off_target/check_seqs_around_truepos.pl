open $infile, shift or die $!;
@crispr=split("","GATCAGGAGCTATTAATTCGCGG");
@crispr_rev=split("","CCGCGAATTAATAGCTCCTGATC");
while(<$infile>){
	$line=$_;
	if($_=~m/^>(.*)\n/){
		$name=$1;
	}else{	
		$length=0;
			while($line=~m/(.{21}GG)/){
				@cur=split("",$1);
				$count=0;
				
				foreach $x (0..22){
						if(@crispr[$x] eq @cur[$x]){
							$count++;
						}
				}
				print $name."\t".$1."\t".abs((length($`)+$length)-1000)."\t".($count-2)."\n";
				$line=substr($line,length($`)+length($1));
				$length=length($`)+$length;
			
			}
			while($line=~m/(CC.{21})/){
				@cur=split("",$1);
				$count=0;
				
				foreach $x (0..22){
						if(@crispr_rev[$x] eq @cur[$x]){
							$count++;
						}
				}
				print $name."\t".$1."\t".abs((length($`)+$length)-1000)."\t".($count-2)."\n";
				$line=substr($line,length($`)+length($1));
				$length=length($`)+$length;
			
			}
		
	}
}
close $infile;