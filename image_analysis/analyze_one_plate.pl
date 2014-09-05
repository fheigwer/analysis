#!/usr/bin/perl
use strict;
use warnings;

my $dir=shift;
my $newdir=shift;
my %name_hash=();
opendir (my $indir, $dir );
foreach my $file (sort(readdir($indir))){
	if($file=~m/(.+)\.tif/){
		my $name= my $newname=$1;
		$name=~s/\s/\\ /g;
		$name=~s/([\(\)])/\\$1/g;
		$newname=~s/\W//g;	
		$name_hash{$name}++;		
		system("cp $dir/$name.tif $newdir/$newname.tif");
	}
}
print "done\n";

#use Parallel::ForkManager;
#my $pm = new Parallel::ForkManager(8);
#my $string="ABCDEFGHIJKLMNOP";
#my @letters=split("",$string);
#
#
#
#foreach $letter (@letters){
#	foreach $number (1..24){
#		$pm->start and next;
#			system('/usr/bin/R -f /Users/b110-mm06/Desktop/Projects/image_analysis/KristinasStuff/R_Files/EBImage_pipeline.R --slave --args '.$letter.$number.'fld1wvDAPIDAPI.tif '.$letter.$number.'fld1wvCy3Cy3.tif '.$letter.$number.'fld1wvCy5Cy5.tif'); 
#		$pm->finish();
#	}
#}
#$pm->wait_all_children();
