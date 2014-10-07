#!/usr/bin/perl -s -w
my $indir=shift;
opendir $outerdir, $indir;
while(my $subdir = readdir($outerdir)){
	if(-d $indir."/".$subdir && !($subdir=~m/\./)){
		my %ctrl=();
		my %sample=();
		open $wellnames, $indir."/".$subdir."/"."wellnames.txt" or die "cannot open $wellnames\n" ;
			while(<$wellnames>){
				my @line=split("\t",$_);				
				if($line[2]=~m/pos/ || $line[2]=~m/neg/ ){
					$ctrl{$line[0]}++;
					
				}elsif($line[2]=~m/sample/){
					$sample{$line[0]}++;
				}
			}
		close $wellnames;
		foreach my $key (keys(%ctrl)){
			print $key."\n";
		}
		opendir $innerdir, $indir."/".$subdir;
		while(my $file = readdir($innerdir)){
			if($file=~m/(DefaultOUT_Nuclei)\.csv/){
				my $name=$1;
				open $ctrlfile, ">".$indir."/".$subdir."/".$name."_ctrl.tab";
				open $allctrlfile, ">>".$indir."/All_ctrl_non_rand.tab";
				open $allctrlfilerand, ">>".$indir."/All_ctrl.tab";
				open $infile, $indir."/".$subdir."/".$file;
				while(<$infile>){
					my @line=split("\t",$_);
					if(exists($ctrl{$line[0]}) && $line[2]!=1  ){
						print $ctrlfile $_;
						print $allctrlfile $_;
					}elsif(exists($sample{$line[0]})){
					}
					if((exists($ctrl{$line[0]}) && $line[2]!=1) || int(rand(50))==25){
						print $allctrlfilerand $_;
					}elsif(exists($sample{$line[0]})){
					}
				}
				close $ctrlfile;
				close $allctrlfile;
				close $allctrlfilerand;
				close $infile;
			}
		}
		closedir $innerdir;
	}	
}
closedir $outerdir;