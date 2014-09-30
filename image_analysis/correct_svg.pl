#!/usr/bin/perl -s -w
use Parallel::ForkManager;
use Sys::Info;
use Sys::Info::Constants qw( :device_cpu );
use Capture::Tiny ':all';
use File::Find;

open ($input, "<".$output_folder."/".$name.".svg") or die $!;		
		while(<$input>){
			if($_=~m/ShowTooltip\(evt\,\'(\-*\d+[\.\d+]*).*/){
				my $tempval=$1;				
				if($tempval>$temp[0]){					
					$temp[0]=$tempval;
				}elsif($tempval<$temp[1]){
					$temp[1]=$tempval;
				}
			}
		}