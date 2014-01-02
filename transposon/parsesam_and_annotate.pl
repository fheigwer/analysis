#!/usr/bin/perl
use strict;
use warnings;
use Set::IntervalTree;
my $file         = shift; # get the file name, somehow
open INFILE , $file;
$file=~/(.*)\.+/;
my $filebase=$1;
my %trees=();
my @line=();
my %results=();
foreach my $line (<INFILE>){
	if(!($line=~m/^@/)){
	    chomp $line;
	@line = split( "\t", $line );
	if ( exists $trees{$line[2]} ) {
	my $annotations = $trees{$line[2]}->fetch( int($line[3]), int($line[3]+1) );
	$results{$line[2]."_".substr($line[3],0,4)}++;
	}else{
        print $line[2]."\n";
	$trees{$line[2]} = build_tree( "/home/mount/fheigwer/transposon_sequencing/dmel_1/DATABASE/DMEL/" . $line[2] . "_indexed.mygff" );
	my $annotations = $trees{$line[2]}->fetch( int($line[3]), int($line[3]+1) );
	print $line[0]."\t".$line[2]."\t".$line[3].
	}
	}
}
close INFILE;

foreach my $key (sort { $results{$a} <=> $results{$b}} keys %results){
    print $key."\t".$results{$key}."\n";
}






sub build_tree {
	my  $file  = $_[0];
	my  $trees = Set::IntervalTree->new;
	open FILE, $file;
	while (<FILE>) {
		  my  $line = $_;
		  chomp $line;
		  my  @line = split( "\t", $line );
		  my  $object = $line[0] . "_" . $line[1] . "_" . $line[2];
		  if($line[1]!=$line[2]){
		      $trees->insert( $object, $line[1], $line[2] );
		      
		  }
	}
	close(FILE);
	return $trees
}
