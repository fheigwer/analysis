#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;
use Sys::Info;
use Sys::Info::Constants qw( :device_cpu );
use Capture::Tiny ':all';
use File::Find;

my $infolder=$ARGV[0];

opendir(my $bigindir, $infolder);

foreach my $subfolder (sort readdir($bigindir)){
	if (-d $infolder."/".$subfolder
		&& !($subfolder=~m/^\./)
	    ) 
	{
		my $outfolder="/data/results/Dmel_HP_20x_4tiles_ERC_HD3A_S13_S14/".$subfolder;
		if (!(-d $outfolder)) {
			system('mkdir '.$outfolder);
		}
		start_off_image_analysis ("/data/scripts/genome_wide_image_analysis/image_analysis_ERC.v2.R",$infolder."/".$subfolder,$outfolder,384,"DAPI.tif","Cy3.tif","FITC.tif",40);
		print $subfolder."\n";
	}	
}

#########################################################################################
#name:      start image analysis
#function:  closes an browser readable svg in a given file if closed already than unclose
#input:     (path/to/filename string) 
#output:    a closed or open svg
#########################################################################################
sub start_off_image_analysis {
	my $path_to_script=$_[0];
	my $input_folder=my $barcode= $_[1];
	$barcode=~s/.+\/(\S+)_.+$/$1/g;
	my $output_folder=$_[2];
	my $format=$_[3];
	my $nuclei_scheme=$_[4];	#fld1wvDAPIDAPI.tif
	my $body_scheme=$_[5];		#fld1wvCy3Cy3.tif
	my $extra_scheme=$_[6];	
	my $tilesize=$_[7];	
	my %Q=();
	foreach my $letter ("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"){		
		foreach my $number (1..24){
			foreach my $field (1..4){
				my @temp=($letter,$number,$field);
				$Q{$letter.$number."_".$field} = \@temp;
			}
		}
	}
	
	#system('echo \'perl /data/scripts/genome_wide_image_analysis/waiter.quick.pl '.$barcode.' '.$output_folder.' '.$tilesize.' '.$format.'\' | qsub -l walltime=12:00:00 ');
	#my $sub_pm = new Parallel::ForkManager(6);#((Sys::Info->new)->device( CPU => %options ))->count);
	my %smQ;
	my $test=0;
	my $counter=0;
	my $test_injob=0;
		QITER: while(%Q){				
				$counter=0;
				CHECKITER: foreach my $key (keys(%Q)){
					my $letter=@{$Q{$key}}[0];
					my $number=@{$Q{$key}}[1];
					my $field=@{$Q{$key}}[2];
					
					#print $output_folder."/".$barcode."_".$letter.$number."_".$field.".tab"."\n";
					#print $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme."\n";
					my $queue_length=capture_stdout {system("qstat -B | grep b ")};					
					$queue_length=~s/\S+\s+\S+\s+\S+\s+(\d+).*$/$1/;
					if(	(-e $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme) 
						&& (-e $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$body_scheme) 
						&& (-e $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$extra_scheme)
						&& !(-z $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme) 
						&& !(-z $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$body_scheme) 
						&& !(-z $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$extra_scheme)
						&& !(-e $output_folder."/".$barcode."_".$letter.$number."_".$field.".tab")
						&& $queue_length<20
					){
						
						$smQ{$key}=$Q{$key};						
						delete $Q{$key};
						$test++;
						$counter++;
						if ($counter>=16) {
									#do the raw analysis either with R ,CP, or Hcell
									my $cmd_start='/opt/software/R-3.1.2/bin/R -f '.$path_to_script.' --slave --args';
									my $filename="";
									foreach my $smkey (keys(%smQ)){
										$filename=$filename." ".$input_folder.'/'.$barcode."_".@{$smQ{$smkey}}[0].@{$smQ{$smkey}}[1]."_".@{$smQ{$smkey}}[2]."_".$nuclei_scheme."::".$input_folder.'/'.$barcode."_".@{$smQ{$smkey}}[0].@{$smQ{$smkey}}[1]."_".@{$smQ{$smkey}}[2]."_".$body_scheme."::".$input_folder.'/'.$barcode."_".@{$smQ{$smkey}}[0].@{$smQ{$smkey}}[1]."_".@{$smQ{$smkey}}[2]."_".$extra_scheme;
									$test_injob++;
									}
									my $cmd_tail=$output_folder;
									system('echo "'.$cmd_start.$filename.' '.$cmd_tail.'" | qsub ');
									sleep 1;
									%smQ=();
									next QITER;
						}					
					}elsif(-e $output_folder."/".$barcode."_".$letter.$number."_".$field.".tab"){						
						delete $Q{$key};
					}
				}				
			}
		if ($counter!=0) {
			my $cmd_start='/opt/software/R-3.1.2/bin/R -f '.$path_to_script.' --slave --args';
			my $filename="";
			foreach my $smkey (keys(%smQ)){
				$filename=$filename." ".$input_folder.'/'.$barcode."_".@{$smQ{$smkey}}[0].@{$smQ{$smkey}}[1]."_".@{$smQ{$smkey}}[2]."_".$nuclei_scheme."::".$input_folder.'/'.$barcode."_".@{$smQ{$smkey}}[0].@{$smQ{$smkey}}[1]."_".@{$smQ{$smkey}}[2]."_".$body_scheme."::".$input_folder.'/'.$barcode."_".@{$smQ{$smkey}}[0].@{$smQ{$smkey}}[1]."_".@{$smQ{$smkey}}[2]."_".$extra_scheme;
			$test_injob++;
			}
			my $cmd_tail=$output_folder;
			system('echo "'.$cmd_start.$filename.' '.$cmd_tail.'" | qsub ');
		}
					
	#open(my $LOGFILE, ">>", "log.txt") or die $!;	
		print "$input_folder\t$test\t$test_injob\n";
		#die;
	#close($LOGFILE);
}
		

