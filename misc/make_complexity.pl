#!/usr/bin/perl
use strict;
use warnings;

use List::Util qw(sum);
use Parallel::ForkManager;
use JSON::XS qw(encode_json decode_json); #library to encode perl objects in a file
use File::Slurp qw(read_file write_file); #library to write perl objects to a file
my $pm = new Parallel::ForkManager(8);
open INFILE, shift;
my $name="";
my $seq="";
my $var=0;
my @count=();
my @names=();
my $i=0;
my %seq_vars=();
foreach my $line (<INFILE>){
	if($line=~/^>(.*)/){
		$name=$1;
		push @names, $name;
		
	}else{	
		$pm->start and next;
		$i=1;
		$var=0;
		while (($i+21)<length($line)){
			$seq=substr($line,$i,21);
			@count=find_base_count(\$seq);
			if(variance(\@count)!=0){
				$var=$var+21/variance(\@count);
			}else{
				$var=$var+21/0.00001;
			}
			$i++;
		}
		$seq_vars{$name}=$var/$i;
		{
			my $json = JSON::XS::encode_json(\%seq_vars);
			write_file( $name . '.json', { binmode => ':raw' }, $json );
		}
		$pm->finish();			
		print $name."\n";
	}
}

$pm->wait_all_children();
foreach  my $name (@names) {
				 my $json = read_file( $name . '.json', { binmode => ':raw' } );
				%seq_vars = ( %seq_vars, %{ decode_json $json } );
				unlink $name.".json";				
			}
close INFILE;

foreach my $key (sort(keys(%seq_vars))){
	print $key."\t".$seq_vars{$key}."\n";
}

sub find_base_count {
	  my  $seqref  = $_[0];
	  my  $seq     = $$seqref;
	  my  $A_count = $seq =~ tr/A//;
	  my  $C_count = $seq =~ tr/C//;
	  my  $T_count = $seq =~ tr/T//;
	  my  $G_count = $seq =~ tr/G//;
	  my  $N_count = $seq =~ tr/N//;
	  my  $count   = length($seq);
	  if ($count==0) {
		$count=1;
	  }
	$A_count = int( $A_count * 100 / $count );
	$C_count = int( $C_count * 100 / $count );
	$T_count = int( $T_count * 100 / $count );
	$G_count = int( $G_count * 100 / $count );
	$N_count = int( $N_count * 100 / $count );
	return $A_count, $C_count, $T_count, $G_count, $N_count;
}
sub mean {
	return sum(@_) / @_;
}
sub variance {
	return ( sum( map { ( $_ - mean( @{ $_[0] } ) )**2 } @{ $_[0] } ) / @{ $_[0] } );
}
