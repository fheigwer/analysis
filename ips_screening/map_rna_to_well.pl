#!/usr/bin/perl -s -w
my $indir=shift;
my %RNA=();
my $done=0;
my $count=1;
opendir $outerdir, $indir;
while(my $subdir = readdir($outerdir)){
	#komisch if($subdir=~m/Pilotkinases_384W_2012_corrected_04-26-12_4cellHTS2/){
	if($subdir=~m/PilotkinaseB_plate 3 adapted_03-27-12_4cellHTS2.txt/){
		open $plateconf, "<".$indir."/".$subdir;
		while(<$plateconf>){
			my @line=split("\t",$_);
			${$RNA{$line[0]}}{$line[1]}=$line[6];
		}
		close $plateconf;
		$done=1;
	}
}
closedir  $outerdir;

opendir $outerdir, $indir;
while(my $subdir = readdir($outerdir)){
	print $done;
	if(-d $indir."/".$subdir && !($subdir=~m/\./ ) && $done==1){
		$subdir=~m/\d{4}(\d)/;
		my $folder=$1;
		print $folder."\n";		
		open $wellnames_rna,">".$indir."/".$subdir."/wellnames_rna.txt";
		$count=1;
		foreach my $letter ("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"){		
			foreach my $number ("03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24"){
				print  $wellnames_rna $count."\t".$letter.$number."\t".${$RNA{$folder}}{$letter.$number}."\n";
				$count++;
			}
		}		
		close $wellnames_rna;
	}
}
			