#!/usr/bin/perl -s -w
use Parallel::ForkManager;
use Sys::Info;
use Sys::Info::Constants qw( :device_cpu );
use Capture::Tiny ':all';
use File::Find;

my $infolder=$ARGV[0];
$pid=0;
opendir $bigindir, $infolder;

foreach my $subfolder (readdir($bigindir)){
	if (-d $infolder."/".$subfolder
	    #&& !(-d "/Users/b110-mm06/Sites/test_images_2/".$subfolder)
	#    && !(kill 0, $pid)
		&& !($subfolder=~m/^\./)
	    ) 
	{
		my $outfolder="/data/results/set13_14/".$subfolder;
		if (!(-d $outfolder)) {
			system('mkdir '.$outfolder);
		}
		start_off_image_analysis ("/data/scripts/image_analysis/image_analysis_ERC.v2.short.R",$infolder."/".$subfolder,$outfolder,384,"DAPI.tif",40);
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
	my $tilesize=$_[5];
	my %Q=();
	foreach my $letter ("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"){		
		foreach my $number (1..24){
			foreach my $field (1..4){
				my @temp=($letter,$number,$field);
				$Q{$letter.$number."_".$field} = \@temp;
			}
		}
	}
	
		system('echo \'perl /data/scripts/image_analysis/waiter.quick.pl '.$barcode.' '.$output_folder.' '.$tilesize.' '.$format.'\' | qsub -l walltime=12:00:00 ');
	#my $sub_pm = new Parallel::ForkManager(3);#((Sys::Info->new)->device( CPU => %options ))->count);
		while(%Q){
			foreach my $key (sort keys(%Q)){
				my $letter=@{$Q{$key}}[0];
				my $number=@{$Q{$key}}[1];
				my $field=@{$Q{$key}}[2];
				my $queue_length=capture_stdout {system("qstat -B | grep b ")};
				#print $output_folder."/".$barcode."_".$letter.$number."_".$field.".tab"."\n";
				$queue_length=~s/\S+\s+\S+\s+\S+\s+(\d+).*$/$1/;
				if(	(-e $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme) 
					&& !(-z $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme) 
					&& !(-e $output_folder."/".$barcode."_".$letter.$number."_".$field.".qtab")
					&& $queue_length<200
					
				){	
				  delete $Q{$key};
					#$sub_pm->start and next;
						#do the raw analysis either with R ,CP, or Hcell
							 system('echo "/opt/software/R-3.1.2/bin/R -f '.$path_to_script.' --slave --args '.$input_folder.'/'.$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme.' '.$output_folder.' '.$barcode."_".$letter.$number."_".$field.'" | qsub ' );
							#system(" echo '/usr/local/bin/R -f $path_to_script --slave --args ".$input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme.' '.
							#							$input_folder."/".$barcode."_".$letter.$number."_".$field."_".$body_scheme.' '.
							#							$input_folder."/".$barcode."_".$letter.$number."_".$field."_".$extra_scheme.' '.
							#							$output_folder.' '.
							#							$barcode."_".$letter.$number."_".$field."' | qsub "); 
							#system('/usr/bin/hcell -f /Users/b110-mm06/Desktop/Projects/image_analysis/KristinasStuff/R_Files/EBImage_pipeline.R --slave --args '.$letter.$number.'fld1wvDAPIDAPI.tif '.$letter.$number.'fld1wvCy3Cy3.tif '.$letter.$number.'fld1wvCy5Cy5.tif'); 
							#system('/usr/bin/CellProfiler -f /Users/b110-mm06/Desktop/Projects/image_analysis/KristinasStuff/R_Files/EBImage_pipeline.R --slave --args '.$letter.$number.'fld1wvDAPIDAPI.tif '.$letter.$number.'fld1wvCy3Cy3.tif '.$letter.$number.'fld1wvCy5Cy5.tif'); 
					#$sub_pm->finish();
				}elsif(-e $output_folder."/".$barcode."_".$letter.$number."_".$field.".qtab"){
				
					delete $Q{$key};
				}
				#sleep 3;
			}
			sleep 1;			
			
		}
		#$sub_pm->wait_all_children();
}
