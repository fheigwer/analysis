#!/usr/bin/perl
use Parallel::ForkManager;
$pm = Parallel::ForkManager->new(8);
opendir INDIR, shift;

foreach $file (readdir(INDIR)){
    if($file=~m/\.fasta/){
	my $pid = $pm->start and next;
	 system("perl cpgi130.pl $file");
	$pm->finish; 
    }
}
$pm->wait_all_children;
closedir INDIR;
