#!/usr/bin/perl -s -w

my $indir=shift;

opendir dir, $indir or die $!;
foreach my $subdir (readdir(dir)){
	if (!($subdir=~m/^\.+$/) && (opendir subdir, $indir.'/'.$subdir)  ){
		foreach my $subsubdir (readdir(subdir)){
			if (!($subsubdir=~m/^\.+$/) && (opendir subsubdir, $indir.'/'.$subdir.'/'.$subsubdir)  ){
				my $header=0;
				open outfile, '>'.$indir.'/'.$subdir.'/'.$subsubdir.'/'.$subdir.'_'.$subsubdir.'.tab' or die $!;
				foreach my $file (readdir(subsubdir)){
					if($file=~m/(.*)\.tab/){
						my $name=$1;
						open infile,  $indir.'/'.$subdir.'/'.$subsubdir.'/'.$file;
							while(<infile>){
								if($_=~/^\d/){
									print outfile $name."\t".$_;
								}elsif($header==0){
									print outfile "well\t".$_;
									$header++;
								}
							}
						close infile;
					}
				}
				close outfile;
			}
			close subsubdir;
		}
		closedir subdir;
	}
	
}
closedir dir;