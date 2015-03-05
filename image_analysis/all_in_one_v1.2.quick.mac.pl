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
		my $outfolder="~/Desktop/test/results/".$subfolder;
		if (!(-d $outfolder)) {
			system('mkdir '.$outfolder);
		}
		start_off_image_analysis ("~/Desktop/analysis/image_analysis/image_analysis_ERC.v2.short.R",$infolder."/".$subfolder,$outfolder,384,"DAPI.tif",40);
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
	
	#system('perl ~/Desktop/analysis/image_analysis/waiter.quick.pl '.$barcode.' '.$output_folder.' '.$tilesize.' '.$format);
	my $sub_pm = new Parallel::ForkManager(6);#((Sys::Info->new)->device( CPU => %options ))->count);
		my %smQ;
		QITER: while(%Q){
			if (scalar(keys(%Q))<16) {
				%smQ=();
				CHECKITER: foreach my $key (sort keys(%Q)){
					my $letter=@{$Q{$key}}[0];
					my $number=@{$Q{$key}}[1];
					my $field=@{$Q{$key}}[2];
					if(	(-e $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme) 
						&& !(-z $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme) 
						&& !(-e $output_folder."/".$barcode."_".$letter.$number."_".$field.".qtab")					
					){
						$smQ{$key}=$Q{$key};					
						
						if (scalar(keys(%smQ))==scalar(keys(%Q))) {
									#do the raw analysis either with R ,CP, or Hcell
									my $cmd_start='/usr/bin/R -f '.$path_to_script.' --slave --args';
									my $filename="";
									foreach my $smkey (sort keys(%smQ)){
										$filename=$filename." ".$input_folder.'/'.$barcode."_".@{$smQ{$smkey}}[0].@{$smQ{$smkey}}[1]."_".@{$smQ{$smkey}}[2]."_".$nuclei_scheme;
									}
									my $cmd_tail=$output_folder;
									system($cmd_start.' '.$filename.' '.$cmd_tail);
						}else{
							next CHECKITER
						}
					}elsif(-e $output_folder."/".$barcode."_".$letter.$number."_".$field.".qtab"){								
						delete $Q{$key};
					}else{
						next CHECKITER;
					}
				}
			}else{
			
				%smQ=();
				CHECKITER: foreach my $key (sort keys(%Q)){
					my $letter=@{$Q{$key}}[0];
					my $number=@{$Q{$key}}[1];
					my $field=@{$Q{$key}}[2];
					if(	(-e $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme) 
						&& !(-z $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme) 
						&& !(-e $output_folder."/".$barcode."_".$letter.$number."_".$field.".qtab")					
					){
						$smQ{$key}=$Q{$key};						
						delete $Q{$key};
						if (scalar(keys(%smQ))==16) {
							$sub_pm->start and next QITER;
									#do the raw analysis either with R ,CP, or Hcell
									my $cmd_start='/usr/bin/R -f '.$path_to_script.' --slave --args';
									my $filename="";
									foreach my $smkey (keys(%smQ)){
										$filename=$filename." ".$input_folder.'/'.$barcode."_".@{$smQ{$smkey}}[0].@{$smQ{$smkey}}[1]."_".@{$smQ{$smkey}}[2]."_".$nuclei_scheme;
									}
									my $cmd_tail=$output_folder;
									system($cmd_start.' '.$filename.' '.$cmd_tail);
							$sub_pm->finish();
							
						}						
					}elsif(-e $output_folder."/".$barcode."_".$letter.$number."_".$field.".qtab"){
						
						delete $Q{$key};
					}else{
						next CHECKITER;
					}
					#sleep 3;
				}
			}
			#sleep 1;
		}
		$sub_pm->wait_all_children();
	
}
		

