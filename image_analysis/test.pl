use File::Find;

$input_folder=$barcode=shift;
$letter="A";
$number="2";
$field="3";
$nuclei_scheme="DAPI.tif";
$barcode=~s/.+\/(\S+)_.+$/$1/g;
if(-e $input_folder.$barcode."_".$letter.$number."_".$field."_".$nuclei_scheme){print "blub\n";}



sub wanted {
    $exp=$letter.$number."_".$field."_".$nuclei_scheme;
    if($_=~/$exp/){$shit=1}
}
#############
					