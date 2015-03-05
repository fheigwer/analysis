#!/usr/bin/perl -s -w
use Parallel::ForkManager;
use Sys::Info;
use Sys::Info::Constants qw( :device_cpu );
use Capture::Tiny ':all';
use File::Find;


my $infolder=$ARGV[0];
opendir $bigindir, $infolder;
open(my $bigout, ">", $infolder."/all_results.tab") or die $!;
    my $bighead=0;
    foreach my $subfolder (readdir($bigindir)){
        print $infolder."/".$subfolder."\n";
        if (    -d $infolder."/".$subfolder
                && !($subfolder=~m/^\./)
            ) 
        {
            my $barcode= $subfolder;
            $barcode=~s/.+\/(\S+)_.+$/$1/g;
            opendir my $indir, $infolder."/".$subfolder;
                open(my $outfile,">", $infolder."/".$subfolder."/".$barcode."_results.tab" ) or die $!;
                    my $heading=0;
                    while(readdir ($indir)) {
                        if($_=~/.*\.tab/){
                            open(my $infile, "<", $infolder."/".$subfolder."/".$_) or die $!;
                                $_=~m/_([A-Z0-9]+?)_([1,2,3,4]).tab/;
                                my $well=$1;
                                 my $field=$2;
                                while (<$infile>) {
                                    if ($_=~m/^\d/) {
                                        print $outfile $well."\t".$field."\t".$barcode."\t".$_;
                                         print $bigout $well."\t".$field."\t".$barcode."\t".$_;
                                    }elsif($heading==0){
                                        print $outfile "well\tfield\tbarcode\t".$_;
                                        if($bighead==0){
                                            print $bigout "well\tfield\tbarcode\t".$_;
                                            $bighead=1;
                                        }
                                        $heading=1;
                                    }
                                }
                            close $infile;
                        }
                    }
                close $outfile;
            closedir $indir;
        }
    }
close $bigout;
closedir $bigindir;