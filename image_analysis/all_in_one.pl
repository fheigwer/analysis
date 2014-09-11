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
	    && !(-d "/var/www/live_images/test_images/".$subfolder)
	#    && !(kill 0, $pid)
	    ) 
	{
		my $outfolder="/var/www/live_images/test_images/".$subfolder;
		if (!(-d $outfolder)) {
			system('mkdir '.$outfolder);
		}
		start_off_image_analysis ("/home/mount/fheigwer/analysis/image_analysis/image_analysis_ERC.R",$infolder."/".$subfolder,$outfolder,384,"DAPI.tif","Cy3.tif","FITC.tif",40);
	}	
}

#########################################################################################
#name:      remove all image files in folder
#function:  renames any tif file in a given folder in a way that each non word character 
#			or white space gets removed
#input:     (path/to/folder string) 
#output:    renamed files if open was possible
#########################################################################################
sub rename_tifs_in_folder {
	my $indir=$_[0];
	opendir $dir, $indir or die print "folder ".$indir." not found";
	foreach $file (sort(readdir($dir))){
		if($file=~m/(.+)\.tif/ && -e $indir.'/'.$file && !(-z $indir.'/'.$file)){
			$name=$newname=$1;
			$name=~s/\s/\\ /g;
			$name=~s/([\(\)])/\\$1/g;
			$newname=~s/\W//g;		
			if($name ne $newname){
				system('mv '.$indir.'/'.$name.'.tif '.$indir.'/'.$newname.'.tif');
			}
		}
	}
	closedir $dir;
	return(1);      
}
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
		rename($_,$File::Find::dir."_$newname.tif");
		
	}
}
#########################################################################################
#name:      make a new svg and write svg header
#function:  writes an browser readable svg header in a given file
#input:     (path/to/filename format tilesize [string] [6|12|24|96|384] [integer] ) 
#output:    a new svg file with header and platedesign
#########################################################################################
sub make_new_plate_svg {
	my $filename=$_[0];		#complete path to file including filename prefix "*.svg"
	my $format=$_[1];		#may be 6, 12, 24, 96, 384 well plate
	my $tilesize=$_[2];		#may be whatever
	my $id=$_[3];
	open my $output ,'>'.$filename.'.svg' or die $!;
	print $output '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">';
	print $output '<svg version="1.1" onload="init(evt)" baseProfile="full" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:ev="http://www.w3.org/2001/xml-events">
	<script type="text/ecmascript">
		<![CDATA[
			function init(evt)
			{
				if ( window.svgDocument == null )
				{
					svgDocument = evt.target.ownerDocument;
				}
				tooltip = svgDocument.getElementById(\'tooltip\');
				tooltip_bg = svgDocument.getElementById(\'tooltip_bg\');
			}
			function ShowTooltip(evt, mouseovertext)
			{
				tooltip.setAttributeNS(null,"x",evt.clientX+11);
				tooltip.setAttributeNS(null,"y",evt.clientY+35);
 				tooltip.firstChild.data = mouseovertext;
 				
 				tooltip_bg.setAttributeNS(null,"x",evt.clientX+8);
				tooltip_bg.setAttributeNS(null,"y",evt.clientY+16);
				length = tooltip.getComputedTextLength();
				tooltip_bg.setAttributeNS(null,"width",length+8);			
				
 				tooltip_bg.setAttributeNS(null,"visibility","visible");
 				tooltip.setAttributeNS(null,"visibility","visible");
			}
			function HideTooltip()
			{
				tooltip.setAttributeNS(null,"visibility","hidden");
				tooltip_bg.setAttributeNS(null,"visibility","hidden");
			}
			function putToTop(evt, whom) {
				//get node reference
				var element = svgDocument.getElementById(\'whom\');
				//appendChild after the last child
				element.parentNode.appendChild(element);
			}
		]]>
	</script>
	<style>
		.tooltip{
			font-size: 20px;
		}
		.tooltip_bg{
			fill: white;
			stroke: black;
			stroke-width: 1;
			opacity: 0.85;
		}
	</style>
	<g id="'.$id.'">'."\n";
	my $x=int($tilesize/2);
	my $y=$x+$tilesize;
	foreach my $letter ("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"){
		print $output	'<text x="'.$x.'" y="'.$y.'" font-family="Verdana" font-size="12" fill="black" >'.$letter.'</text>'."\n";
		$y+=$tilesize;		
	}
	$y=int($tilesize/2);
	$x=$y+$tilesize;		
	foreach my $number (1..24){
		print $output	'<text x="'.$x.'" y="'.$y.'" font-family="Verdana" font-size="12" fill="black" >'.$number.'</text>'."\n";
		$x+=$tilesize;
	}
	close $output;
	return(1);      
}
#########################################################################################
#name:      add image to svg
#function:  adds an base64 encoded image to the prepared svg
#input:     (path/to/filename base64image /path/to/folder [string] [string] [string]) 
#output:    svg
#########################################################################################
sub add_image_to_svg {
	my $filename=$_[0];		#complete path to file including filename 
	my $b64=$_[1];
	my $output_folder=$_[2];
	my $tilesize=$_[3];
	my $format=$_[4];
	my $barcode=$_[5];
	my %coords=();
	my $x=int($tilesize/2);
	my $y=int($tilesize/2);
	foreach my $letter ("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"){
		$y+=int($tilesize/2);
		$x=int($tilesize/2);			
		foreach my $number (1..24){
			foreach my $field (1..2){
				$x+=int($tilesize/2);
				my @temp=($x,$y);
				$coords{$letter.$number."_".$field}=\@temp;
			}
		}
		$y+=int($tilesize/2);
		$x=int($tilesize/2);
		foreach my $number (1..24){
			foreach my $field (3..4){
				$x+=int($tilesize/2);
				my @temp=($x,$y);
				$coords{$letter.$number."_".$field}=\@temp;
			}
		}
	}	
	open ($output, ">>".$output_folder."/plate_image.svg") or die $!;
		my $file=$barcode."_".$filename."_segmented.tif";
		$x=@{$coords{$filename}}[0];
		$y=@{$coords{$filename}}[1];
		print $output '<image id="'.$filename.'" x="'.$x.'" y="'.$y.'" onmousemove="ShowTooltip(evt,\''.$filename.'\')" onmouseout="HideTooltip()" onclick="window.open(\''.$file.'\',\'_blank\');" onmouseover="this.style.cursor=\'pointer\'" width="'.(int($tilesize/2)-4).'px" height="'.(int($tilesize/2)-4).'px" xlink:href="data:image/png;base64,'.$b64.'"></image>'."\n";
	close $output;
	return(1);      
}
#########################################################################################
#name:      add rect to svg
#function:  adds an base64 encoded image to the prepared svg
#input:     (path/to/filename digital_value /path/to/folder [string] [double] [string]) 
#output:    svg
#########################################################################################
sub add_rect_to_svg {
	my $filename=$_[0];		#complete path to file including filename 
	my $value=$_[1];
	my $output_folder=$_[2];
	my $name=$_[3];
	my $tilesize=$_[4];
	my $format=$_[5];
	my $barcode=$_[6];
	my %coords=();
	my $x=int($tilesize/2);
	my $y=int($tilesize/2);
	foreach my $letter ("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"){
		$y+=int($tilesize/2);
		$x=int($tilesize/2);			
		foreach my $number (1..24){
			foreach my $field (1..2){
				$x+=int($tilesize/2);
				my @temp=($x,$y);
				$coords{$letter.$number."_".$field}=\@temp;
			}
		}
		$y+=int($tilesize/2);
		$x=int($tilesize/2);
		foreach my $number (1..24){
			foreach my $field (3..4){
				$x+=int($tilesize/2);
				my @temp=($x,$y);
				$coords{$letter.$number."_".$field}=\@temp;
			}
		}
	}
	open ($output, ">>".$output_folder."/".$name.".svg") or die $!;
		my $file=$barcode."_".$filename."_segmented.tif";
		$x=@{$coords{$filename}}[0];
		$y=@{$coords{$filename}}[1];
		print $output '<rect x="'.$x.'" y="'.$y.'" width="'.(int($tilesize/2)-4).'px" height="'.(int($tilesize/2)-4).'px" style="fill:rgb(0,0,0);" onmousemove="ShowTooltip(evt,\''.$value.'\')" onmouseout="HideTooltip()" onclick="window.open(\''.$file.'\',\'_blank\');" onmouseover="this.style.cursor=\'pointer\'" />'."\n";
	close $output;
	return(1);      
}
#########################################################################################
#name:      close svg
#function:  closes an browser readable svg in a given file if closed already than unclose
#input:     (path/to/filename string) 
#output:    a closed or open svg
#########################################################################################
sub close_open_svg {
	my $filename=$_[0];		#complete path to file including filename 
	open ($FH, "+< $filename") or die "can't update $file: $!";
	while ( <$FH> ) {
		$addr = tell($FH) unless eof($FH);
		if(eof($FH) && $_=~m/svg/){
			truncate($FH, $addr) or die "can't truncate $file: $!";
			return (1);
		} 
	}
	close $FH;	              
	open my $output ,'>>'.$filename or die $!;
		print $output '</g><rect class="tooltip_bg" id="tooltip_bg" x="0" y="0" rx="4" ry="4" width="30" height="30" visibility="hidden"/><text class="tooltip" id="tooltip" x="0" y="0" visibility="hidden">Tooltip</text></svg>';
	close $output;
	return(1);      
}

#########################################################################################
#name:      create total svg 
#function:  collects all other svgs as grouped elements and combines them into one single svg
#input:     (/path/to/folder [string]) 
#output:    svg
#########################################################################################
sub combine_to_svg {
	my $output_folder=$_[0];	
	open ($output, ">".$output_folder."/total.svg") or die $!;
	print $output '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">';
	print $output '<svg version="1.1" onload="init(evt)" baseProfile="full" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:ev="http://www.w3.org/2001/xml-events">
	<script type="text/ecmascript">
		<![CDATA[
			function init(evt)
			{
				if ( window.svgDocument == null )
				{
					svgDocument = evt.target.ownerDocument;
				}
				tooltip = document.getElementById(\'tooltip\');
				tooltip_bg = document.getElementById(\'tooltip_bg\');
			}
			function ShowTooltip(evt, mouseovertext)
			{
				tooltip.setAttributeNS(null,"x",evt.clientX+11);
				tooltip.setAttributeNS(null,"y",evt.clientY+35);
 				tooltip.firstChild.data = mouseovertext;
 				
 				tooltip_bg.setAttributeNS(null,"x",evt.clientX+8);
				tooltip_bg.setAttributeNS(null,"y",evt.clientY+16);
				length = tooltip.getComputedTextLength();
				tooltip_bg.setAttributeNS(null,"width",length+8);			
				
 				tooltip_bg.setAttributeNS(null,"visibility","visible");
 				tooltip.setAttributeNS(null,"visibility","visible");
			}
			function HideTooltip()
			{
				tooltip.setAttributeNS(null,"visibility","hidden");
				tooltip_bg.setAttributeNS(null,"visibility","hidden");
			}
			function putToTop(evt, whom) {
				//get node reference
				var element = svgDocument.getElementById(whom);
				//appendChild after the last child
				element.parentNode.appendChild(element);
				element.parentNode.appendChild(tooltip_bg);
				element.parentNode.appendChild(tooltip);
								
				others = svgDocument.getElementsByClassName("clickable");	
					for(var i=0 ; i<others.length; i++){
						others[i].setAttributeNS(null,"fill","black");	
					}
				var that = evt.target;
				that.setAttributeNS(null,"fill","blue");			
			}
		]]>
	</script>
	<style>
		.tooltip{
			font-size: 20px;
		}
		.tooltip_bg{
			fill: white;
			stroke: black;
			stroke-width: 1;
			opacity: 0.85;
		}
		.clickable { cursor: pointer;  }
		.clickable:hover { color: blue;  }
		.clickable:visited {
			color: red;
			font-family: Verdana;
			text-decoration:none;
		}
	</style>'."\n";
	opendir $outdir, $output_folder;
	my $y=25;
	my $firstfile=0;
	foreach my $file (readdir($outdir)){		
		my $bool=0;
		if($file=~m/(.*)\.svg/ && !($file=~m/total.svg/ )){
			my $name=$1;
			print $output '<text x="1380" y="'.$y.'"  class="clickable" font-family="Verdana" font-size="12" fill="black" onclick="putToTop(evt, \''.$name.'\')" >'.$name.'</text>'."\n";
			$y+=14;
			open $FH, "<".$output_folder."/".$file;
				while (<$FH>){
					my $line =$_;
					if($line=~m/\<g id.*/){
						print $output $line;
						$bool=1;
					}elsif($bool==1){
						if($line=~m/^(.*)<\/svg>/){
							print $output $1;
						}else{
							if($firstfile==0 && $line=~m/\<text .*/){
								print $output $line;
							}elsif(!($line=~m/\<text .*/)){
								print $output $line;
							}
						}
					}	
				}
			close $FH;
			$firstfile++;
		}
	}
	closedir $outdir;
	print $output '<rect class="tooltip_bg" id="tooltip_bg" x="0" y="0" rx="4" ry="4" width="30" height="30" visibility="hidden"/><text class="tooltip" id="tooltip" x="0" y="0" visibility="hidden">Tooltip</text></svg>'."\n";
	close $output;
	return(1);      
}
#########################################################################################
#name:      create total svg 
#function:  collects all other svgs as grouped elements and combines them into one single svg
#input:     (/path/to/folder [string]) 
#output:    svg
#########################################################################################
sub combine_to_html {
	my $output_folder=$_[0];	
	open ($output, ">".$output_folder."/total.html") or die $!;
	print $output '<!DOCTYPE html>
			<html lang="enc">
			<head>
			<meta http-equiv="Content-Type" content="text/html;charset=utf-8">
			<meta name="Author" content="Florian Heigwer"><title>Image Surveillance</title>
			<script type="text/javascript">
			    function init(evt)
			    {
				    if ( window.svgDocument == null )
				    {
					    svgDocument = evt.target.ownerDocument;
				    }
				    tooltip = svgDocument.getElementById(\'tooltip\');
				    tooltip_bg = svgDocument.getElementById(\'tooltip_bg\');
			    }
			    function ShowTooltip(evt, mouseovertext)
			    {
				    tooltip.setAttributeNS(null,"x",evt.clientX+11);
				    tooltip.setAttributeNS(null,"y",evt.clientY+35);
				    tooltip.firstChild.data = mouseovertext;
				    
				    tooltip_bg.setAttributeNS(null,"x",evt.clientX+8);
				    tooltip_bg.setAttributeNS(null,"y",evt.clientY+16);
				    length = tooltip.getComputedTextLength();
				    tooltip_bg.setAttributeNS(null,"width",length+8);			
				    
				    tooltip_bg.setAttributeNS(null,"visibility","visible");
				    tooltip.setAttributeNS(null,"visibility","visible");
			    }
			    function HideTooltip()
			    {
				    tooltip.setAttributeNS(null,"visibility","hidden");
				    tooltip_bg.setAttributeNS(null,"visibility","hidden");
			    }
			    function putToTop(evt, whom) {
				    //get node reference
				    var element = svgDocument.getElementById(whom);
				    //appendChild after the last child
				    element.parentNode.appendChild(element);
				    element.parentNode.appendChild(tooltip_bg);
				    element.parentNode.appendChild(tooltip);
								    
				    others = svgDocument.getElementsByClassName("clickable");	
					    for(var i=0 ; i<others.length; i++){
						    others[i].setAttributeNS(null,"fill","black");	
					    }
				    var that = evt.target;
				    that.setAttributeNS(null,"fill","blue");			
			    }
			    function change_plot(feature){
					var el = document.getElementById(\'plate_container\');
					el.innerHTML = \'<object id="plate_plot" type="image/svg+xml" width="100%" height="1000" data="\'+feature+\'.svg"></object>\'
			    
			    }
			</script>
			<style>
				.tooltip{
					font-size: 20px;
				}
				.tooltip_bg{
					fill: white;
					stroke: black;
					stroke-width: 1;
					opacity: 0.85;
				}
				.clickable { cursor: pointer;  }
				.clickable:hover { color: blue;  }
				.clickable:visited {
					color: red;
					font-family: Verdana;
					text-decoration:none;
				}
			</style>
			</head>
			<body>
			<select name="feature" id="feature" onchange="change_plot(document.getElementById(\'feature\').value)">'."\n";
	opendir $outdir, $output_folder;
	foreach my $file (sort readdir($outdir)){
		if($file=~m/(.*)\.svg/ && !($file=~m/total.svg/ )){
			my $name=$1;
			print $output '<option value=\''.$name.'\' >'.$name.'</option>'."\n";			
		}
	}
	closedir $outdir;
	print $output '</select>
			<br>
			<div id="plate_container">
			<object id="plate_plot" type="image/svg+xml" width="100%" height="1000" data="plate_image.svg"></object>
			</div>
			</body>
			</html>'."\n";
	close $output;
	return(1);      
}

#########################################################################################
#name:      update tile colors on plate
#function:  parses an svg value plate plot for the values and adjusts the colors 
#			according to the values
#input:     (/path/to/folder name [string] [string]) 
#output:    svg
#########################################################################################
sub adjust_colors_in_svg {
	my $output_folder=$_[0];
	my $name=$_[1];
	my @temp=(0,"inf");
	my @color_vector=('rgb(103,0,31)','rgb(178,24,43)','rgb(214,96,77)','rgb(244,165,130)','rgb(253,219,199)','rgb(247,247,247)','rgb(209,229,240)','rgb(146,197,222)','rgb(67,147,195)','rgb(33,102,172)','rgb(5,48,97)');
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
	close $input;
	open ($input, "<".$output_folder."/".$name.".svg") or die $!;		
	open $tempf ,">".$output_folder."/".$name."temp";
		while(<$input>){
			my $line=$_;
			my $color_val=0;
			my $color="rgb(0,0,0)";
			if($line=~m/ShowTooltip\(evt\,\'(\-*\d+[\.\d+]*).*/){
				my $val=$1;	
				if($val>$temp[0]){					
					$temp[0]=$val;
				}elsif($val<$temp[1]){
					$temp[1]=$val;
				}
				if(($temp[0]-$temp[1])!=0){			
					$color_val=10*(($val-$temp[1])/($temp[0]-$temp[1]));
				}else{
					$color_val=0
				}			
				$color=$color_vector[int($color_val)];
			}
			$line=~s/(style\=\"fill\:)rgb\(\d+\,\d+\,\d+\)/$1$color/;
			print $tempf $line;
		}		
		close $tempf;
	close $input;
	rename $output_folder."/".$name."temp",$output_folder."/".$name.".svg";
	return(1);      
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
	my $extra_scheme=$_[6];		#fld1wvCy5Cy5.tif
	my $tilesize=$_[7];
	my $pm = new Parallel::ForkManager((((Sys::Info->new)->device( CPU => %options ))->count));
	my %Q=();
	foreach my $letter ("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"){		
		foreach my $number (1..24){
			foreach my $field (1..4){
				my @temp=($letter,$number,$field);
				$Q{$letter.$number."_".$field} = \@temp;
			}
		}
	}
		
	unless($pid=fork){
		#my $sub_pm = new Parallel::ForkManager((((Sys::Info->new)->device( CPU => %options ))->count));
		while(%Q){
			#find(\&wanted,$input_folder);
			#print join("_",keys(%Q))."\n";
			foreach my $key (sort keys(%Q)){
				my $letter=@{$Q{$key}}[0];
				my $number=@{$Q{$key}}[1];
				my $field=@{$Q{$key}}[2];
				my $queue_length=capture_stdout {system("qstat -B | grep b ")};
				$queue_length=~s/\S+\s+\S+\s+\S+\s+(\d+).*$/$1/;
				if(	(-e $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme) 
					&& (-e $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$body_scheme) 
					&& (-e $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$extra_scheme)
					&& !(-z $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme) 
					&& !(-z $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$body_scheme) 
					&& !(-z $input_folder."/".$barcode."_".$letter.$number."_".$field."_".$extra_scheme)
					&& !(-e $output_folder."/".$barcode."_".$letter.$number.".tab")
					&& $queue_length<20
					
				){
					delete $Q{$key};
					#$sub_pm->start and next;
						#do the raw analysis either with R ,CP, or Hcell
							#system("/usr/bin/R -f $path_to_script --slave --args ".	$input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme.' '.
							#							$input_folder."/".$barcode."_".$letter.$number."_".$field."_".$body_scheme.' '.
							#							$input_folder."/".$barcode."_".$letter.$number."_".$field."_".$extra_scheme.' '.
							#							$output_folder.' '.
							#							$barcode."_".$letter.$number."_".$field);
							system(" echo '/usr/local/bin/R -f $path_to_script --slave --args ".$input_folder."/".$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme.' '.
														$input_folder."/".$barcode."_".$letter.$number."_".$field."_".$body_scheme.' '.
														$input_folder."/".$barcode."_".$letter.$number."_".$field."_".$extra_scheme.' '.
														$output_folder.' '.
														$barcode."_".$letter.$number."_".$field."' | qsub "); 
							
							
							#system('/usr/bin/hcell -f /Users/b110-mm06/Desktop/Projects/image_analysis/KristinasStuff/R_Files/EBImage_pipeline.R --slave --args '.$letter.$number.'fld1wvDAPIDAPI.tif '.$letter.$number.'fld1wvCy3Cy3.tif '.$letter.$number.'fld1wvCy5Cy5.tif'); 
							#system('/usr/bin/CellProfiler -f /Users/b110-mm06/Desktop/Projects/image_analysis/KristinasStuff/R_Files/EBImage_pipeline.R --slave --args '.$letter.$number.'fld1wvDAPIDAPI.tif '.$letter.$number.'fld1wvCy3Cy3.tif '.$letter.$number.'fld1wvCy5Cy5.tif'); 
					#$sub_pm->finish();		
					
				}elsif(-e $output_folder."/".$barcode."_".$letter.$number."_".$field.".tab"){
					delete $Q{$key};
				}
				sleep 3;
			}
			sleep 10;
		}
		#$sub_pm->wait_all_children();
	}else{
		
		while(%Q){
			foreach my $key (sort keys(%Q)){				
				my $letter=@{$Q{$key}}[0];
				my $number=@{$Q{$key}}[1];
				my $field=@{$Q{$key}}[2];
				if(	(-e $output_folder."/".$barcode."_".$letter.$number."_".$field.".tab")){
					#convert the R result
						$base64_string = capture_stdout {system('convert '.$output_folder."/".$barcode."_".$letter.$number."_".$field.'_segmented.tif -thumbnail 100x100 png:- | openssl enc -base64  ')}; #-out '.$output_folder.'/'.$letter.$number.'.b64'
					#read in the tab-delim line for this image
						my $count=0;
						my @names=();
						my @result=();
						open $result_tab , $output_folder."/".$barcode."_".$letter.$number."_".$field.".tab";
							while(<$result_tab>){
								chomp $_;
								if($count<=0){
									@names=split("\t",$_);
									$count++;
								}elsif($_=~m/^\d/){
									@result=split("\t",$_);
								}		
							}
						close $result_tab;								
						if(!(-e $output_folder."/plate_image.svg")){
							make_new_plate_svg($output_folder."/plate_image",$format,$tilesize,"plate_image");	#start new plate svg file
							add_image_to_svg($letter.$number."_".$field,$base64_string,$output_folder,$tilesize,$format,$barcode);
							close_open_svg($output_folder."/plate_image.svg");
						}else{
							close_open_svg($output_folder."/plate_image.svg");
							add_image_to_svg($letter.$number."_".$field,$base64_string,$output_folder,$tilesize,$format,$barcode);
							close_open_svg($output_folder."/plate_image.svg");
						}
						my $element=0;
						foreach my $key (@names){
							my $value=$result[$element];							
							$element++;
							if(!(-e $output_folder."/".$key.".svg")){
								make_new_plate_svg($output_folder."/".$key,$format,$tilesize,$key);	#start new plate svg file
								add_rect_to_svg($letter.$number."_".$field,$value,$output_folder,$key,$tilesize,$format,$barcode);
								adjust_colors_in_svg($output_folder,$key);
								close_open_svg($output_folder."/".$key.".svg");
							}else{
								close_open_svg($output_folder."/".$key.".svg");
								add_rect_to_svg($letter.$number."_".$field,$value,$output_folder,$key,$tilesize,$format,$barcode);
								adjust_colors_in_svg($output_folder,$key);
								close_open_svg($output_folder."/".$key.".svg");
							}
						}		
					delete $Q{$key};
				}
			}
			sleep 10;
			#combine_to_svg($output_folder);
			combine_to_html($output_folder);
			system("chmod -R a+r ".$output_folder."")
		}
	}
}












































