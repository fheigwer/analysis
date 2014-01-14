opendir dir, shift ;
foreach $file (sort(readdir(dir))){
	if($file=~m/(.+)\.tif/){
		$name=$newname=$1;
		$name=~s/\s/\\ /g;
		$name=~s/([\(\)])/\\$1/g;
		$newname=~s/\W//g;	
		$name_hash{$name}++;
		print "$name\n$newname\n";
		system("mv $name.tif $newname.tif");
	}
}

use Parallel::ForkManager;
my $pm = new Parallel::ForkManager(8);
my $string="ABCDEFGHIJKLMNOP";
my @letters=split("",$string);



foreach $letter (@letters){
	foreach $number (1..24){
		$pm->start and next;
			system('/usr/bin/R -f /Users/b110-mm06/Desktop/Projects/image_analysis/KristinasStuff/R_Files/EBImage_pipeline.R --slave --args '.$letter.$number.'fld1wvDAPIDAPI.tif '.$letter.$number.'fld1wvCy3Cy3.tif '.$letter.$number.'fld1wvCy5Cy5.tif'); 
		$pm->finish();
	}
}
$pm->wait_all_children();
