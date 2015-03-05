#!/usr/bin/perl
use strict;
use warnings;
use File::Find;

my $dir=shift;
sub wanted {
	if($_=~m/(.+)\.tif/){
		my $newname=$1;
		$newname=~s/\W//g;
		$newname=~s/^_//g;
		$newname=~s/fld/_/g;
		$newname=~s/wv/_/g;
		$newname=~s/DAPIDAPI/DAPI/g;
		$newname=~s/Cy3Cy3/Cy3/g;
		$newname=~s/Cy5Cy5/Cy5/g;
		$newname=~s/FITCFITC/FITC/g;
		$newname=~s/0(\d)/$1/g;
		$File::Find::dir=~s/.+\/(\S+)_.+$/$1/g;
		if (!(-e $File::Find::dir."_$newname.tif") && !($_=~m/$File::Find::dir/)) {
			rename($_,$File::Find::dir."_$newname.tif");
		}
	}
}
find(\&wanted,$dir);


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
