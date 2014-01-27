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
			if($file=~m/(.+)\.csv/){
				my $name=$1;
				open $samplefile, ">".$indir."/".$subdir."/".$name."_sample.tab";
				open $ctrlfile, ">".$indir."/".$subdir."/".$name."_ctrl.tab";
				open $infile, $indir."/".$subdir."/".$file;
				while(<$infile>){
					my @line=split("\t",$_);
					#print $line[0]."\n";
					if(exists($ctrl{$line[0]}) && $line[2]!=1){
						print $ctrlfile $_;
					}elsif(exists($sample{$line[0]})){
						print $samplefile $_;
					}
				}
				close $samplefile;
				close $ctrlfile;
				close $infile;
			}
		}
		closedir $innerdir;
	}	
}
closedir $outerdir;