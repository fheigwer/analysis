#!/usr/bin/perl

use strict;
use warnings;

my $dir=shift;
use Parallel::ForkManager;
my $pm = new Parallel::ForkManager(8);
my $string="ABCDEFGHIJKLMNOP";
my @letters=split("",$string);


foreach my $letter (@letters){
    foreach my $number (1..24){
	$pm->start and next;
	if(-e $dir.'/'.$letter.$number.'fld1wvCy3Cy3_segmented.png'){}else{
            print $dir.'/'.$letter.$number."\n";            
	    system('/home/mount/fheigwer/R-3.0.2/bin/R -f /home/mount/fheigwer/image_analysis/segmentedRData_to_png.R --slave --args '.$dir.'/'.$letter.$number.'fld1wvCy3Cy3.tif');
	}
	$pm->finish();
    }
}
$pm->wait_all_children();
