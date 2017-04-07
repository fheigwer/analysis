
BEGIN { $ENV{CLUSTALDIR} = '/usr/bin/' }

use warnings; #throw warnings for debugging the perl syntax mistakes
use strict; #require strict definition of every single variable, obey the accidential reuse of global variables

use CGI qw(:standard VARS);
use CGI::Carp qw ( fatalsToBrowser );

use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::SeqFeature::Generic; #important package to handle sequence formats and objects
use Bio::Tools::Run::Alignment::Clustalw; #library used to align DNA sequences
use Bio::SimpleAlign;
use Bio::Graphics;
use List::Util;
use Parallel::ForkManager; #important package to enable mutlithreading of the script
use Bio::Location::Split; #library to make splitted location objects
use Set::IntervalTree; #library providing methods for using interval trees

######################################################################################
my $query         = new CGI;
my @params		= $query->param;

my $process		= param("PROCESS");
my $temp_dir	= param("ID");
my $waiting_time	= param("WAITINGTIME");
my $process_id    = param("PID");
my $progress      = param("PRG");

my @seq_array           = ();
my %something           = ();
my $seqsrc              = "";
my %trees               = ();
my $seqio_obj           = "";
my $starttime           = time();

chdir("/var/www/E-CRISP/workdir");

if ( $temp_dir eq "") {
      $temp_dir = localtime.time;
      $temp_dir =~ s/\s+/\_/ig;
      mkdir( "$temp_dir") or die $!;
}

system('chmod -R o+rwx /var/www/E-CRISP/workdir/'.$temp_dir.';');

foreach my $element (@params) {
      $something{$element} = $query->param($element);
}

      
if (!($something{"peakdet"} eq "on")) {
      $something{"peakdet_level"}="strict";   
}
      

my $result = "/var/www/E-CRISP/workdir/$temp_dir/fertig.txt";
my $error = "/var/www/E-CRISP/workdir/$temp_dir/error.html";

if ( $waiting_time eq "") {
	$waiting_time = 0;
} elsif($waiting_time<8) {
	$waiting_time += 2;
}
#
print "Content-type:text/html\r\n\r\n";
open(HEADER,"header_CRISPR_nogoog.txt");
while(<HEADER>){
	print $_;
}
close HEADER;

my $filename    = $something{"input_file"};
my $seq_in_file = "";
if ($filename) {
	$CGI::POST_MAX = 1024 * 20000;
	my $safe_filename_characters = "a-zA-Z0-9_.-";
	my $upload_dir = $temp_dir;
	my ( $name, $path, $extension ) = fileparse( $filename, '\..*' );
	$filename = $name . $extension;
	$filename =~ tr/ /_/;
	$filename =~ s/[^$safe_filename_characters]//g;
	if ( $filename =~ /^([$safe_filename_characters]+)$/ ) {
		$filename = $1;
	} else {
		print_error_html( $temp_dir, "Filename contains invalid characters");
		die;
	}
	my $upload_filehandle = $query->upload("input_file");
	open( UPLOADFILE, ">$upload_dir/$filename" ) or die print_error_html( $temp_dir, "Could not write the uploaded file");  #Not needed
		binmode UPLOADFILE;
		while (<$upload_filehandle>) {
			print UPLOADFILE;
			}
	close UPLOADFILE;
	$seq_in_file = $upload_dir . "/" . $filename;
} else {
	$seq_in_file = "";
}

my $databasepath = "/data/DATABASEFILES/" . $something{"ref_organism"} . "/" . $something{"ref_organism"};

if ( $process ne "") {
	my $result = "/var/www/E-CRISP/workdir/$temp_dir/fertig.txt";
	my $error = "/var/www/E-CRISP/workdir/$temp_dir/error.html";
     
	if ( -e $error ) {
		print "<meta http-equiv=\"refresh\" content=\"0; URL=workdir/$temp_dir/error.html\">";
		print end_html;
	}elsif ( -e $result ) {
		print "<meta http-equiv=\"refresh\" content=\"0; URL=workdir/$temp_dir/results.html\">";
		print end_html;
	}elsif(!(kill 0, $process_id)){
                  print_error_html( $temp_dir, "SORRY! <br><br>Some serious error occured during the process.<br>In case this happens again do not hesitate to contact us at crispr\@dkfz.de .\n" );
                  print "<meta http-equiv=\"refresh\" content=\"0; URL=workdir/$temp_dir/error.html\">";
                  print end_html;
      } else {
		
            open FACTFILE,"/var/www/E-CRISP/workdir/funfacts.txt" or die $!;	
		my @array=<FACTFILE>;
		my $randomline=$array[rand @array];
		close FACTFILE;
            
            open (my $statsfile, "<", "/var/www/E-CRISP/workdir/$temp_dir/statsfile.txt") or die $!;
                  my @stats=();
                  while(<$statsfile>){
                        @stats=split("\t",$_);
                  }
            close $statsfile;
			
            my @circle = (-710+$progress*7.1,-710+$stats[0]*7.1);
		
            print '
                <tr>
                  <br><br>
                  <td style="width:800px" align=center>
                  <svg version="1.1" id="Layer_0" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0px" y="0px" width="240px" height="240px" viewBox="0 0 240 240" enable-background="new 0 0 240 240" xml:space="preserve">
							<image x="20" y="20" width="200" height="200" xlink:href="/E-CRISP/workdir/Align.svg"/>
							
							<circle cx="120" cy="120" r="114" fill="none" stroke="#2662C3" stroke-width="10" stroke-dasharray="800" stroke-dashoffset="'.$circle[0].'" transform="rotate(-90, 120, 120)">
                                          	<animate
                                                      attributeName="stroke-dashoffset"
                                                      from="'.$circle[0].'"
                                                      to="'.$circle[1].'"
                                                      begin="0s"
                                                      dur="1s"
                                                      values="'.$circle[0].';'.$circle[1].'"
                                                      keySplines="0.1 0.8 0.2 1; " 
                                                      keyTimes="0;1" 
                                                      calcMode="spline"
                                                      fill="freeze"
                                                />                      
                                          </circle>
						</svg>
					</td>
				</tr><tr>
					<td align=center>
						<font size=3><h3><br>'.$stats[1].'</h3>Search for '.$stats[0].'% done</font>
					</td>
				</tr><tr style="height:200px">
					<td align=center style=" vertical-align=top">      
						<br>
						<font size=3>
							While you are waiting, why don\'t you visit E-RNAi or E-TALEN?
						<br><br>Did you know that:<br>
						'.$randomline.'

				
			</font>
			</td>

                </tr>
		<tr>
		<td>
			<br>This page will be automatically updated every '.$waiting_time.' seconds until search is done
			<br>
		</td>
                </tr>
		
		';
		open(FOOTER,"/var/www/E-CRISP/workdir/footer.txt");
				while(<FOOTER>){
					print $_;
				}
				close FOOTER;
		print "<meta http-equiv=\"refresh\" content=\"" . $waiting_time . "; URL=SeqAlign.pl?PROCESS=1&ID=$temp_dir&&WAITINGTIME=$waiting_time&PID=$process_id&PRG=".$stats[0]."\">";

	}
} else {
      
      unless( $process_id = fork){
            
               close(STDOUT);
               close(STDERR);
               open STDERR, ">/var/www/E-CRISP/workdir/$temp_dir/error.log";
               open STDOUT, ">/var/www/E-CRISP/workdir/$temp_dir/out.log";
               
               if ($something{string_length}>5000) {
                        print_error_html( $temp_dir, "Please provide a peak length shorter than 5000bp.\n" );
                        die;
               }   
                  
               my $ctfac = Bio::Tools::Run::Alignment::Clustalw->new();    #set up the clustalw factory
               my $pm = Parallel::ForkManager->new(15); #set up forkmanager for a maximum of 15 sequences
               
               
                                    
            if ( exists $something{"GENE.SYMB"}) {
            
            #################################################################################################################################################################################
            # For ENSEMBLE: converts the ENSEMBLE data to a temporary FASTA file
            #################################################################################################################################################################################
             
                  my $db = Bio::DB::Fasta->new( $databasepath . ".all.dna.fa", -makeid => \&make_my_id );
                  
                  if ( !$filename ) {
                        my $temp = $something{"pasted_seq"};
                        $temp =~ s/\s+//ig;
                        my @ids = split( ";", $temp );
                              
                        make_temp_fasta_file(\@ids, \%trees, \%something, $db, $temp_dir, 0);
                        $seqio_obj = Bio::SeqIO->new( -file => $temp_dir . "/tempfile.fasta", -format => "fasta" ); #read the temporary fasta file

                  } else {
                        my @ids = ();
                        open (my $infile, "<", $seq_in_file);
                  
                              while (<$infile>) {
                                    my $line = $_;
                                    chomp $line;
                                    push @ids, $line;
                              }
                        
                        close $infile;
                              
                        make_temp_fasta_file(\@ids, \%trees, \%something, $db, $temp_dir, 1);
                        $seqio_obj = Bio::SeqIO->new( -file => $temp_dir . "/tempfile.fasta", -format => "fasta" ) or die print_error_html($temp_dir, "No gene symbols have been entered.<br>"); #read the temporary fasta file
                              
                  }
                       
                  $seqsrc = $temp_dir."/tempfile.fasta";       # Newly created temporary FASTA file will be used for alignment process
                       
            } else {
                  
                  open (my $statsfile, ">", "/var/www/E-CRISP/workdir/$temp_dir/statsfile.txt") or die $!;
                  print $statsfile "0\tChecking sequences...\n";
			close $statsfile;
                        
            #################################################################################################################################################################################
            # For FASTA: checks if the fasta sequence/file is in the right format
            #################################################################################################################################################################################
                        
                  if ( !$filename ) { #check if the pasted sequence should be used
              
                        my $pasted = $something{"pasted_seq"}; #get the sequence(s) from the textentry window
              
                        if (length($pasted)>201000) {
                              print_error_html($temp_dir, "Please enter a sequence that is shorter than 200 kilo base pairs");
                              die;
                        }
              
                        if ($pasted=~m/^(>[^\n]+)$/) {
                              print_error_html( $temp_dir, "\\\"$1\\\" is not a FASTA format sequence.<br>FASTA sequences need a header starting with &quot;&gt;&quot;, which is new line separated from the sequence.<br><br> e.g.: &gt;some Sequence<br>AGCTGATCGATCTAGCTAGCTGCTAGCTAGTCGATCGATGCTAGCTAGCTAGCTGCTAG<br>");
                              die;
                        }
                              
                        my @pasted=split("\n",$pasted);
                        my $count=0;
              
                        foreach my $line (@pasted){
                              if ($line=~m/^(>.+)/) {
                                    $count++;
                              }elsif ($line=~m/([^ACGTUN\s]+)/){
                                    print_error_html( $temp_dir, "\\\"$line\\\" is not a FASTA format file because it contains &quot;$1&quot; as bases (only ATCGUN are allowed).<br>");
                                    die;
                              }
                        }
              
                        if ($count==0){
                              print_error_html($temp_dir, "No FASTA sequences have been entered.<br>FASTA sequences need a header starting with &quot;&gt;&quot;, which is new line separated from the sequence.<br><br> e.g.: &gt;some Sequence<br>AGCTGATCGATCTAGCTAGCTGCTAGCTAGTCGATCGATGCTAGCTAGCTAGCTGCTAG<br>");
                                    die;
                              }   
              
                        if ($count>15) {
                              print_error_html( $temp_dir, "Your input is more than 15 Sequences\n" );
                              die;
                        }
                              
                        open (my $tempfile, ">", $temp_dir . "/tempfile.fasta"); #write it to an temporary fasta format file
                              print $tempfile $pasted;
                        close $tempfile;
                              
                        $seqsrc = $temp_dir."/tempfile.fasta";
                        $seqio_obj = Bio::SeqIO->new( -file => $temp_dir . "/tempfile.fasta", -format => "fasta" ) or die print_error_html($temp_dir, "No FASTA sequences have been entered.<br>FASTA sequences need a header starting with &quot;&gt;&quot;, which is new line separated from the sequence.<br><br> e.g.: &gt;some Sequence<br>AGCTGATCGATCTAGCTAGCTGCTAGCTAGTCGATCGATGCTAGCTAGCTAGCTGCTAG<br>");

                        
                  } else {
                        
                        my $count=0;
                        my $temp="";
                        
                        open(my $infile, "<", $seq_in_file);
                              while (my $line = <$infile>){
                                    if ($line=~m/^(>.+)/) {
                                          $count++;
                                    }elsif ($line=~m/([^ACGTUN\s]+)/){
                                          print_error_html( $temp_dir, "\\\"$seq_in_file\\\" is not a FASTA format file because it contains &quot;$1&quot; as bases\n" );
                                          die;
                                    }
                                    $temp=$temp.$line;
                              }
                        close $infile;
                              
                        if (length($temp)>201000) {
                              print_error_html($temp_dir, "Please enter a sequence that is shorter than 200 kilo base pairs");
                              die;
                        }
                        
                        if ($temp=~m/^(>[^\n]+)$/) {
                              print_error_html( $temp_dir, "\\\"$1\\\" is not a FASTA format sequence.<br>FASTA sequences need a header starting with &quot;&gt;&quot;, which is new line separated from the sequence.<br><br> e.g.: &gt;some Sequence<br>AGCTGATCGATCTAGCTAGCTGCTAGCTAGTCGATCGATGCTAGCTAGCTAGCTGCTAG<br>");
                              die;
                        }
                        
                        if($count==0){
                              print_error_html( $temp_dir, "A FASTA format sequence needs to have a header starting with &quot;&gt;&quot;<br><br> e.g.: &gt;some Sequence<br>AGCTGATCGATCTAGCTAGCTGCTAGCTAGTCGATCGATGCTAGCTAGCTAGCTGCTAG<br>" );
                              die;
                        }   
                        
                        if ($count>50) {
                              print_error_html( $temp_dir, "Your input is more than 50 Sequences\n" );
                              die;
                        }
                                                      
                        $seqsrc = $temp_dir.$filename;   #uploaded FASTA file will be used for the alignment process
                        $seqio_obj = Bio::SeqIO->new( -file => $seq_in_file, -format => "fasta" ) or die print_error_html( $temp_dir, "Your input wasn't FASTA format <br><br> e.g.: &gt;some Sequence<br>AGCTGATCGATCTAGCTAGCTGCTAGCTAGTCGATCGATGCTAGCTAGCTAGCTGCTAG<br>" ); #if neither online or pasted sequences are used it will use the input file as sequence input
                        
                  }
            }
            
            open (my $statsfile, ">", "/var/www/E-CRISP/workdir/$temp_dir/statsfile.txt") or die $!;
            print $statsfile "10\tPerforming alignment...\n";
		close $statsfile;
      
      
            #################################################################################################################################################################################
            # Align sequences using clustalw            
            #################################################################################################################################################################################
      
            my $aln = $ctfac->align($seqsrc) or die print_error_html("An error happened during the alignment process.<br>Please make sure your sequences are properly formatted.");
            $aln -> sort_alphabetically();
            
            #################################################################################################################################################################################
            # Detect homology peaks
            #################################################################################################################################################################################

            my @conservation = $aln->consensus_conservation();

            my %peakobj;

            my @peakfeatures;
            
            open (my $statsfile, ">", "/var/www/E-CRISP/workdir/$temp_dir/statsfile.txt") or die $!;
            print $statsfile "50\tSearching for homology peaks...\n";
		close $statsfile;
                  
            if ($something{"peakdet_level"} eq "strict") {
                  
                  %peakobj = detect_peaks_strict(\@conservation,$something{"string_length"},$something{merge_peaks}); 
            
            }  else  {
                  
                  %peakobj = detect_peaks_dynamic(\@conservation,$something{"string_length"},$something{merge_peaks},$something{"peakdet_level"});
            
            }
            

            
            foreach my $peak (0..($peakobj{"peak_count"}-1)){
               
               my $peakfeature = Bio::SeqFeature::Generic->new(-start => $peakobj{"peaks"}[$peak]{pos},
                                                               -end   => ($peakobj{"peaks"}[$peak]{pos}+$peakobj{"peaks"}[$peak]{"size"}),
                                                               -display_name => "Peak".($peak+1)."::consensus:".$peakobj{"peaks"}[$peak]{"pos"}."..".($peakobj{"peaks"}[$peak]{"pos"}+$peakobj{"peaks"}[$peak]{"size"}-1),
                                                               -primary_tag =>  "Peak".($peak+1)."::consensus:".$peakobj{"peaks"}[$peak]{"pos"}."..".($peakobj{"peaks"}[$peak]{"pos"}+$peakobj{"peaks"}[$peak]{"size"}-1),
                                                           );               
            
               push @peakfeatures, $peakfeature;
               
            }   
                  
            #################################################################################################################################################################################
            # Fetch annotations, locate gaps for each sequence
            #################################################################################################################################################################################              
                  
            my @seq_array;
            my %seq_coords;
      
            while( my $seq = $seqio_obj->next_seq() ) {     #fetch all seq_obj sequences to get their annotations
                  
                  push(@seq_array,$seq);
                  
            }
            
            
                  
          foreach my $seq_obj ( @seq_array ) {       #for each provided sequence a child process is started
               
               $pm->start and next;
               
               my $fname; 
               my $chrom               = ""; #create an empty name chromosome
               my $location_offset     = 0;
               my $location_end        = 0; #set the location offset value to 0
               if ( ( $seq_obj->description ) =~ m/chrom:([\w\.]+):(\d+)..(\d+)/ig) { #if the descrition of the sequence contains information of the form loc=chr1:0..13000000
                     #print the result of the pattern matching
                     $chrom = $1; #save the location information in the certain objects
                     $location_offset = $2;
                     $location_end = $3; #the location offset is the start of the sequence on that certain chromosome, for the whole chromosome this is of cause 0
               }
                    
               $fname = ( $seq_obj->display_id ); #current name is the sequence' id
               if ( $chrom eq "" ) { $chrom = $fname; } #if $chrom is still empty fill it with the sequence' id
               my $whole_seq = $seq_obj->seq; #deduce the complete nucleotide sequence as a alphanumerical string from the SeqIo object
               
                           
               open (my $statsfile, ">", "/var/www/E-CRISP/workdir/$temp_dir/statsfile.txt") or die $!;
               print $statsfile "70\tMapping $fname...\n";
               close $statsfile;
                    
               my %transcripts_hash=();
               my %CDS_hash        =();
                    
               print "<p>\n\nfetching annotations for  ".$fname."...\n";
                                          
               if ( exists $trees{$chrom} ) {
                                    
                        my $annotations   = $trees{$chrom}->fetch( int($location_offset), int($location_end) );     

                        foreach my $anno ( sort( @{$annotations} ) ) {
                              if ( $anno =~ m/exon::(\S+)::(\d)::(.+)_(\d+)_(\d+)$/ig ) {
                                                
                                    my @pair = ( $4-$location_offset, $5-$location_offset );
                                    ${$transcripts_hash{$1}}{$2}=\@pair;
                              
                              } elsif ( $anno =~ m/CDS::(\S+)::(\d)::(.+)_(\d+)_(\d+)$/ig ) {
                                    
                                    my @pair = ( $4, $5 );
                                    push @{ $CDS_hash{$1} }, \@pair;
                                    
                              }
                                    
                        }
					
			}
			     
			
                    # Sort the exons by their start position and put them into an array
			
			my %sorter;
			my %transcripts_hash_sorted;
               my %transcriptfeatures;
			
			foreach my $transcript ( sort( keys (%transcripts_hash) ) ){
    
				 $sorter{$transcript} = [sort { $transcripts_hash{$transcript}{$a}[0] <=> $transcripts_hash{$transcript}{$b}[0] } keys $transcripts_hash{$transcript} ];
    
				 foreach my $exon ( @{$sorter{$transcript}} ){
	  
					 my @pair=($transcripts_hash{$transcript}{$exon}[0],$transcripts_hash{$transcript}{$exon}[1]);
					 push @{$transcripts_hash_sorted{$transcript}}, \@pair;
 
				 }

			}
               
               my %strand=();
               foreach my $key ( sort( keys(%transcripts_hash) ) ) {
                    my $curstart = 0;
                    my $diff = 0;
                    $strand{$key} = 1;
                    
                    foreach my $pair_ref (sort(keys %{$transcripts_hash{$key}})) {
                         $diff = ${$transcripts_hash{$key}}{$pair_ref}->[0]-$curstart;
                         $curstart = ${$transcripts_hash{$key}}{$pair_ref}->[0];
                    }
                    
                    if ($diff<0) {
                         $strand{$key} = -1;
                    }
               }
						
			# locate the gaps in the alignment sequences      
			
			my $sequence=$aln->get_seq_by_id($fname)->seq();
			my @gaps;
    
			while ($sequence =~ /\.+/g) {
                    
				 my $pos=$-[0];
				 my $size=$+[0]-$-[0];

				 push @gaps, [$pos,$size];
	  
			}
               
               my $seqoffset = 0;         
               my $seqend = length($sequence);
                
               unless ($gaps[0][0]!=0){ $seqoffset=$gaps[0][1]; }
               
               (reverse $sequence) =~ /\.+/;
               
               my $pos=$-[0];
               my $size=$+[0]-$-[0];
               
               unless ($pos!=0){ $seqend=(length($sequence)-$size); }
               
               print "\n".$fname." coords: ".$seqoffset."..".$seqend;
               
               $seq_coords{$fname}=[$seqoffset,$seqend];

                # for each transcript, map exons with located gaps
				
			foreach my $transcript ( sort( keys (%transcripts_hash_sorted) ) ){
                    
                    print "\n\nProcessing transcript ".$transcript."... ";
                    
                    my $exon = 0;
                    my $last_exon = 0;
                    my $transcript_end=$transcripts_hash_sorted{$transcript}[scalar @{$transcripts_hash_sorted{$transcript}}-1][1];

                    foreach my $gap (0..(scalar (@gaps)-1)) { 
                         
                         my $gap_pos=$gaps[$gap][0];
                         my $gap_size=$gaps[$gap][1];
                         
                         $exon=$last_exon;
                         my $loc=0;
                         
                         if ($gap_pos > $transcript_end) {
						last;
					}
					
                         
                         if($gap_pos <= $transcripts_hash_sorted{$transcript}[0][0]){                              
                              $loc=1; #if gap is located before first exon
                         }
                         
                         while ($transcripts_hash_sorted{$transcript}[$exon]&&$loc!=1) {   #loop for locating the gap
						
                              if ($transcripts_hash_sorted{$transcript}[$exon][0]<$gap_pos&&$gap_pos<$transcripts_hash_sorted{$transcript}[$exon][1]) {  #if gap is located in an exon, cut it at the gaps position
							
                                  	my @new_exon = ($gap_pos, $transcripts_hash_sorted{$transcript}[$exon][1]);
                                         
                                   $transcripts_hash_sorted{$transcript}[$exon][2]=$gap_size; 								
							$transcripts_hash_sorted{$transcript}[$exon][1]=$gap_pos; 
							
							splice(@{$transcripts_hash_sorted{$transcript}},$exon+1,0,\@new_exon);
                                   $loc=1;
                                   $last_exon=$exon;
                                   
						} elsif ($transcripts_hash_sorted{$transcript}[$exon+1]) {
                              
                                   if($transcripts_hash_sorted{$transcript}[$exon][1]<$gap_pos&&$gap_pos<$transcripts_hash_sorted{$transcript}[$exon+1][0]) {
                                        
                                        $loc=1;
                                        $last_exon=$exon;

                                   }                                   
                                   
                              }     
                              
                              $exon++;
                              
					}
                         
                         if ($loc) { #if gap was located, shift the following exons
                              
                              $transcript_end+=$gap_size;
						
                              foreach my $e (($exon)..(scalar (@{$transcripts_hash_sorted{$transcript}})-1)){
                              
                                   $transcripts_hash_sorted{$transcript}[$e][0]+=$gap_size;
                                   $transcripts_hash_sorted{$transcript}[$e][1]+=$gap_size;
                              
                              }
					
                         }
                         
                    }                                       				
                    
                    my @features;
				my $splitlocation = Bio::Location::Split->new();
                    my $firstblock=1;

				foreach my $exon (0..@{$transcripts_hash_sorted{$transcript}}-1){     #generate features for the graph
                         
                         my $subloc= 	Bio::Location::Simple->new(
														-start => int( $transcripts_hash_sorted{$transcript}[$exon][0] ),
														-end => int( $transcripts_hash_sorted{$transcript}[$exon][1] ),
														-strand => 0,
                                                                       );																
					$splitlocation->add_sub_Location($subloc);                              				
					
					if (defined $transcripts_hash_sorted{$transcript}[$exon][2] ) {
                              
                              my $feature = Bio::SeqFeature::Generic->new (
															-location => $splitlocation,
															-score => 0,
                                                                           );
                              
                              if($firstblock){
                                   
                                   if($strand{$transcript}<0)  {
                                        $feature->set_attributes(-display_name => "< Transcript::" . $transcript,
                                                                 -primary_tag => "< Transcript::" . $transcript,
                                                                 -strand => -1);
                                        
                                   } else {
                                        $feature->set_attributes(-display_name => "Transcript::" . $transcript." >",
                                                                 -primary_tag => "Transcript::" . $transcript,);
                                   }
                                        
                                   $firstblock=0;
                                   
                              } else { $feature->set_attributes(-display_name => " ",
												-primary_tag => " ",
                                                           );   }
                              
                              push @features, $feature;
						
						if (defined $transcripts_hash_sorted{$transcript}[$exon+1][0] ){	#create a spacer feature for the gap
							
							my $gaplocation = Bio::Location::Split->new();
							$gaplocation->add_sub_Location(Bio::Location::Simple->new   (
																			-start => int( $transcripts_hash_sorted{$transcript}[$exon][1] ),
																			-end => int( $transcripts_hash_sorted{$transcript}[$exon][1] ),
																			));
                                   
							$gaplocation->add_sub_Location(Bio::Location::Simple->new   (
																			-start => int( $transcripts_hash_sorted{$transcript}[$exon+1][0]),
																			-end => int( $transcripts_hash_sorted{$transcript}[$exon+1][0]),
																			));
							
							my $gapfeature = Bio::SeqFeature::Generic->new(  
                                                                                     -location => $gaplocation,
                                                                                     -score => 1, #score indicates this is a gap
                                                                                     -display_name => " ",
                                                                                     -primary_tag => " ",
                                                                                  );
                                   
							push @features, $gapfeature;
                                   						
						}
						
                              $splitlocation = undef;
						$splitlocation = Bio::Location::Split->new();
											
					}
                         
				}
                    
                    my $feature = Bio::SeqFeature::Generic->new (
															-location => $splitlocation,
															-score => 0,
                                                                           -display_name => " ",
                                                                           -primary_tag => " ",
                                                                           );
                    
                    if($strand{$transcript}==1)  {    $feature->set_attributes( -strand => 1 );  }
                    
                    if($firstblock){
                                                                     
                                   if($strand{$transcript}<0)  {
                                        $feature->set_attributes(-display_name => "< Transcript::" . $transcript,
                                                                 -primary_tag => "Transcript::" . $transcript,
                                                                 -strand => -1);
                                        
                                   } else {
                                        $feature->set_attributes(-display_name => "Transcript::" . $transcript." >",
                                                                 -primary_tag => "Transcript::" . $transcript,);
                                   }
                                   
                    }
                                                                 
                    
                    push @features, $feature;
                    
				$transcriptfeatures{$transcript} = \@features;                  
                    
                    
               }
                
               my $gaprawlocation = Bio::Location::Split->new();
               foreach my $gap(@gaps){
                         
               $gaprawlocation->add_sub_Location(Bio::Location::Simple->new   (
																-start => $gap->[0]+1,
																-end => $gap->[0]+$gap->[1],
																));
                         
               }
                    
               my $gapraw = Bio::SeqFeature::Generic->new(
                                                               
                    -location => $gaprawlocation
                                                               
               ); 
						 
               #####################################################################################################################################################################
               # Draw the Image of the sequence
               #####################################################################################################################################################################
               
               open (my $statsfile, ">", "/var/www/E-CRISP/workdir/$temp_dir/statsfile.txt") or die $!;
               print $statsfile "80\tGenerating alignment image...\n";
               close $statsfile;
               
               my $consensus_string = $aln->consensus_string();
               
               my $consensus_feature = Bio::SeqFeature::Generic->new( -start => int(1),
                           -end => length ($aln->consensus_string()),
                           -display_name => ( $seq_obj->display_name )
               );
               
               my $wholeseq = Bio::SeqFeature::Generic->new( -start => int(1),
                           -start => $seqoffset,
                           -end => $seqend,
                           -display_name => $fname,
               );
               
               my $panel = Bio::Graphics::Panel->new(
                           -length => (length ($aln->consensus_string())),
                           -key_style => 'bottom',
                           -width     => 770,
                           -pad_left  => 10,
                           -pad_right => 10,
                           -grid      => 1,
                           -spacing   => 10,
                           -key_spacing   => 10,
                           -image_class=> 'GD::SVG',
                           -bump_limit => 0,
                              );
               
               $panel->add_track( $consensus_feature,
                           -glyph  => 'arrow',
                           -bump   => 0,
                           -double => 1,
                           -tick   => 2
                              );
               
               $panel->add_track( $consensus_feature,
                           -glyph   => 'generic',
                           -bgcolor => '#88CCFF',
                           -fgcolor => 'transparent',
                           -label   => 'Consensus',
                           -key => 'Consensus'
                    
                              );
               
               $panel->add_track( $wholeseq,
                           -glyph   => 'generic',
                           -bgcolor => '#2662C3',
                           -fgcolor => 'transparent',
                           -label   => 1,
                           -description => $seq_obj->description(),
                           -key     => 'Gene',
                              );
               
               $panel->add_track( $gapraw,
                           -glyph   => 'generic',
                           -bgcolor => 'green',
                           -fgcolor => 'transparent',
                           -key     => 'Alignment gap',
                              );
               
               my $key=0;
               
               foreach my $transcript ( sort( keys (%transcripts_hash_sorted) ) ){

                    unless($key){$panel->add_track(\@{$transcriptfeatures{$transcript}},
                              -bgcolor     => '#FF9900',
                              -fgcolor => 'transparent',
                              -glyph => 'generic',
                              -description   => 1,
                              -connector=> sub    {
                                                       my ($feature,$p) = @_;
                                                       if ($feature->can("score")) {
                                                            if ($feature->score==1){return 'solid';} else {return 'hat';}
                                                       }     
                                                  },
                              -connector_color=>'black',
                              -strand_arrow=>'ends',
                              -bump => 0,
                              -bump_limit => 0,
                              -key  => 'Transcript',
                              -label => 1,

                         );
                                 
                                 $key++;
                                 
                    }    else {
                         
                         $panel->add_track(\@{$transcriptfeatures{$transcript}},
                              -bgcolor     => '#FF9900',
                              -fgcolor => 'transparent',
                              -glyph => 'generic',
                              -description   => 1,
                              -connector=> sub    {
                                                       my ($feature,$p) = @_;
                                                       if ($feature->can("score")) {
                                                            if ($feature->score==1){return 'solid';} else {return 'hat';}
                                                       }     
                                                  },
                              -connector_color=>'black',
                              -strand_arrow=>'ends',
                              -bump => 0,
                              -bump_limit => 0,
                              -label => 1,

                         );
                         
                    }
                    
               }
               
               $panel->add_track( \@peakfeatures,
                           -glyph   => 'generic',
                           -bgcolor => '#EE3300',
                           -fgcolor => 'transparent',
                           -key     => 'Homology peak',
                           -label   => 1
                              );
                     
               open (my $svgfile, ">", $temp_dir . "/SVG_".$fname.".svg");
                    print $svgfile $panel->svg;
               close $svgfile;
               
               $pm->finish; # do the exit in the child process
                                    	
          }
          
          #####################################################################################################################################################################
          # Get the coordinates of the Sequences in the Align
          #####################################################################################################################################################################
           
          foreach my $seq_obj (@seq_array){
               
               my $fname = ( $seq_obj->display_id );
               
               my $sequence=$aln->get_seq_by_id($fname)->seq();
    
			$sequence =~ /\.+/g;
               
               my $pos=$-[0];
               my $gap=$+[0]-$-[0];
               
               my $seqoffset = 0;         
               my $seqend = length($sequence);
                
               unless ($-[0]!=0){ $seqoffset=$gap; }
               
               (reverse $sequence) =~ /\.+/;
               
               my $pos=$-[0];
               my $gap=$+[0]-$-[0];
               
               unless ($pos!=0){ $seqend=(length($sequence)-$gap); }
               
               $seq_coords{$fname}=[$seqoffset,$seqend];
               
          }
          
          $pm->wait_all_children; #wait for all mapping processes to finish
          
          my $peak_chrom = "N/A";   # get data for a fasta header for e-crisp
          my $magnification = 1;
          my $seq_start = 0;
          my $seq_mapped_start = 0;
          
          if ( ( $seq_array[0]->description ) =~ m/chrom:([\w\.]+):(\d+)..(\d+)/ig) { #if the descrition of the sequence contains information of the form loc=chr1:0..13000000
                    
                     $peak_chrom = $1; #save the location information in the certain objects
                     my $location_offset = $2;
                     my $location_end = $3; 
                     
                     my $fname = ( $seq_array[0]->display_id ); #current name is the sequence' id
                     $magnification = ($location_end - $location_offset)/($seq_coords{$fname}[1]-$seq_coords{$fname}[0]);
                     
                     $seq_start = $location_offset;
                     $seq_mapped_start = $seq_coords{$fname}[0];
                     
          }
          
          ##########################################################################################################################################
          # draw image of the alignment
          ##########################################################################################################################################
          
          my $consensus_string = $aln->consensus_string();
               
          my $consensus_feature = Bio::SeqFeature::Generic->new( -start => int(1),
                           -end => length ($aln->consensus_string()),
          );
          
          my $panel = Bio::Graphics::Panel->new(
                           -length => (length ($aln->consensus_string())),
                           -key_style => 'bottom',
                           -width     => 770,
                           -pad_left  => 10,
                           -pad_right => 10,
                           -grid      => 1,
                           -spacing   => 10,
                           -key_spacing   => 10,
                           -image_class=> 'GD::SVG',
                           -bump_limit => 0,
                              );
               
          $panel->add_track( $consensus_feature,
                           -glyph  => 'arrow',
                           -bump   => 0,
                           -double => 1,
                           -tick   => 2
                              );
               
          $panel->add_track( $consensus_feature,
                           -glyph   => 'generic',
                           -bgcolor => '#88CCFF',
                           -fgcolor => 'transparent',
                           -label   => 'Consensus',
                              );
                    
          foreach my $seq_obj ( @seq_array ){
               
               my $fname = ( $seq_obj->display_id );
                             
               my $seq = Bio::SeqFeature::Generic->new( 
                           -start => $seq_coords{$fname}[0],
                           -end => $seq_coords{$fname}[1],
                           -display_name => $fname,
               );
               
               $panel->add_track( $seq,
                           -glyph   => 'generic',
                           -bgcolor => '#2662C3',
                           -fgcolor => 'transparent',
                           -label   => $fname,
                              );
          }
          
          $panel->add_track( \@peakfeatures,
                           -glyph   => 'generic',
                           -bgcolor => '#EE3300',
                           -fgcolor => 'transparent',
                           -label   => 1
                              );
          
          open (my $tempfile, ">", $temp_dir . "/SVG_Alignment.svg");
               print $tempfile $panel->svg;
          close $tempfile; 

                       
            
            
            #################################################################################################################################################################################
            # Print the Information as HTML document
            #################################################################################################################################################################################
            
            my $seqID         = "";                        
            my $seqstring     = "";
            my $consensus_string = $aln->consensus_string();
            my @sequences; #will start from sequences[1]!
            my $num           = $aln->num_sequences();
            
            open (my $statsfile, ">", "/var/www/E-CRISP/workdir/$temp_dir/statsfile.txt") or die $!;
            print $statsfile "90\tCreating report file...\n";
            close $statsfile;
            
            for(1..($num)){
                             
                  $sequences[$_] = $aln->get_seq_by_pos($_);
                  
            }

            
            open (my $file, ">", $temp_dir . "/results.html");
            
                  open(my $header, "<", "/var/www/E-CRISP/workdir/header_Align.txt");
                        while(<$header>){
                              print $file $_;
                        }                                                                                                                                                                          
                  close $header;                                                                                                                                                                  
            
            print $file "<style type='text/css'>\n html { overflow-y: scroll; }    div.container      {width: 795px; max-width:795px; margin-right:auto; margin-left:auto; clear:left; clear:right; overflow-x:hidden;}\n div.sequences     {width: 634px; overflow-x: scroll; }\n table.head{table-layout:fixed; font-family:arial; border:1px; margin-right:auto; margin-left:auto;}\n table.main{table-layout:fixed; font-family:arial; border:1px; float:left; text-align:center;}\n tr    {height:24px; font-size:15;}\n tr.cons    {height:132px;} th    {border-bottom: 1px solid #000000; color: #FFFFFF; padding: 2px; background-color: #2662C3;font-size: 14px;}\n td.b    {max-width:13px; min-width:13px; width:13px; vertical-align:center; text-align:center;}\n td:hover    {color:#000000;} th.tab:hover    {background-color: #56a2f3;}\n</style>"; 
            
            print $file "<table class='main' style='width:100%'><tr><td colspan=2 style='text-align:left;'><h3>Sequence alignment results</h3>(after ".(time()-$starttime)."s)</td></tr><tr><td colspan=2><embed style='width: 790px; height: 10px;' alt='there should appear a  line' src='http://e-crisp-test.dkfz.de/E-CRISP/workdir/gradientline.svg' type='image/svg+xml'></td></tr><tr     style='text-align:left'><td>Alignment score: ".$aln->score()."</td><td>Average homology: ".round_digits((List::Util::sum(@conservation)/@conservation),2)."%</td></tr>      <tr style='text-align:left'><td>Homology peaks detected: ".$peakobj{"peak_count"}."</td><td>Length of consensus string: ".length($consensus_string)."bp</tr><tr></tr></table>\n<img src='SVG_Alignment.svg'>\n<p align='center'>
               <div style='display:block; height:32px; min-height:32px'><table class='main' style='width:800px; height:24px; border-spacing: 0px; table-layout: fixed'>
               <tr style='height:8px'><th colspan=4 style='border-bottom-width: 0px; width: auto'></th></tr>
               <tr><th> </th>
               <th class='tab' id='tab_peaks' style='background-color: #ffffff; color: #000000; border-bottom-color: #ffffff' onclick='switchToPeaks()'><a onclick='switchToPeaks()'>
               Peaks</a>
               </th>
               <th class='tab' id='tab_sequences' onclick='switchToSeqs()'><a onclick='switchToSeqs()'>
               Alignment Map</a>
               </th><th></th></tr><tr><td></td></tr>
               </table></div>
               </p><br><div id='page_peaks' style='display:block; width: 800px; max-width: 800px'>";
                  
            #Print the table for every peak:
             
            foreach my $o (1..$peakobj{"peak_count"}){
                  
                  my $peak_num=$o-1;
                  my $peak_pos=$peakobj{"peaks"}->[$peak_num]->{"pos"};
                  my $peak_size=$peakobj{"peaks"}->[$peak_num]->{"size"};
                  my $avg = (List::Util::sum(@conservation[$peak_pos..($peak_pos+$peak_size-1)]))/$peak_size;
                  
                  my $start = cutoffdigits($peak_pos*$magnification+$seq_start-$seq_mapped_start);
                  my $end = cutoffdigits(($peak_pos+$peak_size)*$magnification+$seq_start-$seq_mapped_start);
                  
                  my $peak_tag = "consensus chrom:".$peak_chrom.":".$start."..".$end;
                  
                  print $file "<div class='container'>\n";      
                  print $file "<table style='width:100%; height:24px; horizontal-align:left; font-family:arial;'><tr><th> Peak $o - pos ".$peak_pos."-".($peak_pos+$peak_size-1)." on alignment consensus   string</th></tr>\n";
                  
                  if ($peak_size>$something{"string_length"}) {
                  
                         print $file "<tr><td style='max-width:100%; min-width:100%; width:100%; background-color:#FD8E51;'> This peak was extended due to multiple peaks of the same score overlapping</td></tr>\n";
                  
                  }
            
            
                  print $file "</table><table class='main' cellpadding=2>\n";
            
                  for(1..$num){
                                  
                        $seqID      = $sequences[$_]->display_id();

                        print $file "<tr>\n";           
                        print $file "<th style='minimum-width:150; width:150px;'>".$seqID."</th>\n";            
                        print $file "</tr>\n";
            
                  }
      
                  print $file "<tr style='height:16px; '></tr>\n";
                  print $file "<tr><th style='minimum-width:150 width: 150px;'> Consensus seq.: </th></tr>\n";      
                  print $file "<tr class='cons' style='height:132; minimum-height:132; font-size: 10px; vertical-align:center;'><th> Conservation <br> (hover for %)<br><br> Average<br> conservation: <br> ".round_digits($avg)."%</th></table>\n";      
      
                  print $file "<div class='sequences'><table class='main'>\n";
      
                  for(1..$num){
                             
                        my $seq = $aln->get_seq_by_pos($_);
                        my $seqstring= $seq->seq();
      
                        print $file "<tr>\n";
              
                        print $file colorcode_DNA(substr ($seqstring, $peak_pos, $peak_size));
                      
                        print $file "\n</tr>\n";
            
                  }
      
                  print $file "<tr style='height:16px; maximum-height:16px;font-size: 10px;'>\n";

      
                  for(($peak_pos)..($peak_pos+$peak_size-1)) {
    
                        if (($_%5)==0) {
        
                              print $file "<td class='b'>$_</td>";
    
                        }     else  {
                  
                        print $file "<td class='b'></td>";
      
                        }
      
                  }                                  
                                       
                  print $file "\n</tr><tr >\n";
                  print $file colorcode_DNA(substr($consensus_string,$peak_pos,$peak_size));
                  print $file "\n</tr><tr style='height:116; font-size: 10px; vertical-align:bottom; color:#FFFFFF;'>\n";      

                  foreach my $i (($peak_pos)..($peak_pos+$peak_size-1)){
                  print $file "<td class='b' style='border-bottom:".$conservation[$i]."px solid #8FBEDD;'>".cutoffdigits($conservation[$i])."</td>";
                  }
            
                  print $file "</table></div>
                  <table style='width:100%; height:24px; horizontal-align:left; font-family:arial;'><tr><th class='tab'><a target='_blank' style='color: white;' href='../../designcrispr.html?tag=>".$peak_tag." &pasted_seq=".substr($consensus_string,$peak_pos,$peak_size)."&ref_org=".$something{"ref_organism"}."'>>> Import '$peak_tag' into E-CRISP << </a></th></tr></table></div><br><p align=center>\n";
      
            }
      
            
            print $file "</div><div id='page_sequences' style='display:none; width: 7795px; max-width: 795px; margin-left:auto; margin-right:auto'><table style='width:100%; horizontal-align:left; height:24px;'>";
               
            for(1..($num)){
                                  
               $seqID      = $sequences[$_]->display_id();
                        
               print $file "<tr><th>".$seqID."</th></tr>";
               print $file "<tr><td><img src='SVG_".$seqID.".svg'><br></td></tr><tr style'height:24px'></tr>"
            
            }
               
            print $file "</th></tr></table></div>";
            
            open(my $footer, "<", "/var/www/E-CRISP/workdir/footer.txt");
                  while(<$footer>){
                        print $file $_;
                  }
            close $footer;
      
            close ($file);
      
            chmod 0755, $temp_dir . "/results.html";
            
            open(LOG,">>", "//var/log/talecrisp.log");
            print LOG $temp_dir."\tE-CRISP\tALIGN\n";
            close(LOG);
            system( 'touch ' . $temp_dir . '/fertig.txt' );
            print end_html;
      
      } else {
		
            print "
			<tr>
				<td align=\"left\"><font size=\"3\">WAITING for results</font></td>
			</tr>
            	";
		open(FOOTER,"footer.txt");
				while(<FOOTER>){
					print $_;
				}
		close FOOTER;
		print "<meta http-equiv=\"refresh\" content=\"" . $waiting_time . "; URL=SeqAlign.pl?PROCESS=1&ID=$temp_dir&WAITINGTIME=$waiting_time&PID=$process_id\">"; 
		print end_html;
                
      }
      
}

                  
#########################################################################################
#name:      make_my_id
#function:  return searchable ids for generation of the FASTA index
#input:     (FASTA header)
#output:    array of searchable ids (strings)
#########################################################################################
sub make_my_id {
      $_[0] =~m/^>(\S+) locus_tag= (\S+);/;
      return ( $1, $2);
}

#########################################################################################
#name:      make_temp_fasta_file
#function:  creates a temporary fasta file
#input:     (given id-Array, tree reference, something-Hashreference,
#           database, temp_dir, 1/0 if file or not)
#output:    N/A
#########################################################################################
sub make_temp_fasta_file {
      if ($_[5] == 0 && scalar(@{$_[0]})>15) {
            print_error_html( $_[4], "Your input is more than 15 Sequences.\n" );
             die;
      }
      if ($_[5] == 0 && scalar(@{$_[0]})<2) {
            print_error_html( $_[4], "Please provide at least 2 different sequences.\n" );
             die;
      }
      if ($_[5] == 1 && scalar(@{$_[0]})>15) { #@Flo die 500 sind hier Absicht?
            print_error_html( $_[4], "Your input is more than 15 lines with IDs.<br> Please shorten the list, or maybe change to option to FASTA.<br>" );
            die;
      }
      open (my $tempfile, ">", $_[4] . "/tempfile.fasta");
            foreach my $id (@{$_[0]}) { 
                  $id =~ s/\s//ig;
                  my $seq_obj = $_[3]->get_Seq_by_id($id); # get a PrimarySeq obj
                  if ($seq_obj) {
                        my $header = $_[3]->header($id); # get the header, or description line
				print "the header should go here ".$header;
                        $header =~ m/chrom:([\w\.]+):(\d+)..(\d+)/ig;
                        my $chrom = $1;
                        my $location_offset = $2;
                        my $location_end = $3;
                        if ( !exists $_[1]->{$chrom} ) {
                              $_[1]->{$chrom} = build_tree( "/data/DATABASEFILES/" . $_[2]->{"ref_organism"} . "/" . $chrom . "_indexed" );
                        }
				print "\nSequence $id: length: ".length($seq_obj->seq())."<br>";
                        print $tempfile ">", $seq_obj->display_id(), " ", $header, "\n", $seq_obj->seq(), "\n";
                  } else {
                  open (my $failfile, ">>", $_[4] . "/failfile.tab");
                         print $failfile substr($id,0)."\n";
                  close($failfile);
                  }
            }
      close $tempfile;
      if ( open (my $failfile, "<", $_[4] . "/failfile.tab") && !( $_[2]->{"ignore_missing_id"} eq "ignore")) {
            my $error="";
            while (<$failfile>) {
                  $error.=$_."<br>";
            }
            print_error_html( $_[4], "No Database entry found for <br> \"". $error."\" in the \" ".$_[2]->{"ref_organism"}."\" genome.<br> Please enter a valid ensembl ID or gene symbol (case sensitive) or change the input option above to FASTA sequence.<br>" );
            close($failfile);
            die;
      }
      
}



#########################################################################################
#name:      print_error_html
#function:  print any error messages as a html file for convenience
#input:     (/path/to/workdir as string, message as string) 
#output:    readable html file displaying the error message
#########################################################################################
sub print_error_html {
      my $temp_dir      = $_[0];
      my $message       = $_[1];
      open (my $errorfile, ">", $temp_dir . "/error.html");
            open(my $header, "<", "/var/www/E-CRISP/workdir/header_CRISPR.txt");
                  while(<$header>){
                        print $errorfile $_;
                  }
            close $header;
            print $errorfile '<tr><td align="left" ><span class="errormessage"><br><br><br><big>'. $message . '</big></span></td></tr>';		
		
            open(my $footer, "<", "/var/www/E-CRISP/workdir/footer.txt");
                  while(<$footer>){
                        print $errorfile $_;
                  }
            close $footer;
      close $errorfile;
      open(my $log,">>", "/var/log/talecrisp.log");
            print $log $temp_dir."\tE-CRISP\tERROR\t$message\n";
      close($log);
      chmod 0755, $temp_dir . "/error.html";
}

#########################################################################################
#name:      round_digits
#function:  cutoff after certain number of digits
#input:     (float number as string) 
#output:    shortened float
#########################################################################################
sub round_digits{
      if($_[0]=~m/(\d+\.\d{1,2})/){
            return $1;
      }
      return $_[0];
}

#########################################################################################
#name:      colorcode_DNA
#function:  color-codes a DNA sequence
#input:     Sequence as string
#output:    string Rigged with html tags for colors
#########################################################################################

sub colorcode_DNA{
      
      my %colors;      
      my $seq = $_[0];
      

            %colors = (
            
                  "A" => "<td class='b' bgcolor=#90D055>A</td>", #A => green      aaff88
                  "C" => "<td class='b' bgcolor=#8FBEDD>C</td>", #C => blue       88aaff
                  "G" => "<td class='b' bgcolor=#F3ED63>G</td>", #G => yellow     ffee88
                  "T" => "<td class='b' bgcolor=#FD8E51>T</td>", #T => red        ffaa88
                  "." => "<td class='b' bgcolor=#FFFFFF>-</td>",            
            
            );
            

      $seq=~s/([ACGT\.])/$colors{$1}/g;
      
      return $seq;
      
}

#########################################################################################
#name:      cutoffdigits
#function:  cuts off all digits behind the comma
#input:     float
#output:    float with only a single digit
#########################################################################################

sub cutoffdigits{
      
      if($_[0]=~/(\d+)\.\d/){return $1;} else {return $_[0];}
      
}

#########################################################################################
#name:      detect_peaks_strict
#function:  non-dynamic (but accurate) search for local score peaks
#input:     array containing conservation percentages, peak width to look for
#output:    array containing positions of local peaks
#########################################################################################

sub detect_peaks_strict{
      
      my $maxscore = 0;
      my %peakobj = ();
      my $localsum = 0;
      my $width = $_[1];
      my @scores = @{$_[0]};
      my $end = @scores-$width;
      my $peakcount = 0;
      my $lastpos = undef;
      my $merge = $_[2];

      $width--;

      my $i = 0;
      while ($i <= $end){      
      
            $localsum=List::Util::sum(@scores[$i..($i+$width)]);
               
            if($localsum>=$maxscore){

                  if ($localsum>$maxscore) {
      
                        $maxscore=$localsum;
                        %peakobj=();
                        $peakcount=0;
                        $lastpos = undef;
                        
                        $peakobj{"peaks"}->[$peakcount]->{"pos"}=$i;
                        $peakobj{"peaks"}->[$peakcount]->{"size"}=$something{"string_length"};                  
                        $peakcount++;
                        $lastpos=$i;
                                                  
                  } elsif($merge&&$i<($lastpos+$width)&&$peakcount>0) {
                              
                        $peakobj{"peaks"}->[$peakcount-1]->{"size"}=(($i+$width)-$lastpos);
                        
                  } elsif($i==$lastpos+1&&$peakcount>0) {
                        
                        $peakobj{"peaks"}->[$peakcount-1]->{"size"}++;
                        $lastpos++
                  
                  }  else  {                  
                  
                        $peakobj{"peaks"}->[$peakcount]->{"pos"}=$i;
                        $peakobj{"peaks"}->[$peakcount]->{"size"}=$something{"string_length"};                  
                        $peakcount++;
                        $lastpos=$i;

                  }
                                  
            }
            
            $i++;
    
      }
      
      $peakobj{"peak_count"}=$peakcount;
         
      return %peakobj;
      
}

#########################################################################################
#name:      detect_peaks_dynamic
#function:  Less accurate but much faster search for local score peaks
#input:     array containing conservation percentages, peak width to look for, integer speed level (>=1, the higher the slower and more accurate)
#output:    peak object containing positions and sizes of local peaks
#########################################################################################

sub detect_peaks_dynamic{

      my $maxscore = 0;
      my %peakobj = ();
      my $localsum = 0;
      my $width = $_[1];
      my @scores = @{$_[0]};
      my $end = @scores-$width;
      my @speedlevel = 6-$_[3];
      my $peakcount = 0;
      my $speed = 1;
      my $lastpos = -$width;
      my $lastsum = 0;
      my $merge = $_[2];


      $width--;

      my $i = 0;
      while ($i <= $end){      
      
            $localsum=List::Util::sum(@scores[$i..($i+$width)]);
               
            if($localsum>=$maxscore){

                  if ($localsum>$maxscore) {
      
                        $maxscore=$localsum;
                        %peakobj=();
                        $peakcount=0;
                        $lastpos = undef;
                        
                        $peakobj{"peaks"}->[$peakcount]->{"pos"}=$i;
                        $peakobj{"peaks"}->[$peakcount]->{"size"}=$something{"string_length"};                  
                        $peakcount++;
                        $lastpos=$i;
                                                  
                  } elsif($merge&&$i<($lastpos+$width)&&$peakcount>0) {
                              
                        $peakobj{"peaks"}->[$peakcount-1]->{"size"}=(($i+$width)-$lastpos);
                        
                  } elsif($i==$lastpos+1&&$peakcount>0) {
                        
                        $peakobj{"peaks"}->[$peakcount-1]->{"size"}++;
                        $lastpos++
                  
                  }  else  {                  
                  
                        $peakobj{"peaks"}->[$peakcount]->{"pos"}=$i;
                        $peakobj{"peaks"}->[$peakcount]->{"size"}=$something{"string_length"};                  
                        $peakcount++;
                        $lastpos=$i;

                  }
                                  
            }
    
            if ($localsum>=$lastsum) {
        
                  $speed-=@speedlevel;
    
            }   else    {
        
                  $speed++;
        
            }
    
            if ($speed < 1) {
                  
                  $speed=1;
                  
            }
    
            $lastsum=$localsum;

            $i+=$speed;
    
      }
      
      $peakobj{"peak_count"}=$peakcount;
      
      return %peakobj;
      
}

#########################################################################################
#name:      build_tree
#function:  if a pre.build tree exists, load this tree into memory if not build it 
#           from source and save to it to file and return the tree anyways
#input:     (path to data without file ending)
#output:    augmented black/red tree (Set::IntervalTree->new)
#########################################################################################
sub build_tree {
      my  $trees = Set::IntervalTree->new;
      if (-e $_[0].".otree") {
            $trees->LoadTree($_[0].".otree");
      }elsif(-e $_[0].".mygff")  {
            open (my $infile, "<", $_[0].".mygff");
                  foreach my $line (<$infile>) {
                        chomp $line;
                        my @line=split("\t",$line);
                        my $object= $line[0] . "_" . $line[1] . "_" . $line[2];
                        $trees->insert(  $object, $line[1], $line[2] );
                  }
            close($infile);
            $trees->SaveTree($_[0].".otree");
      }
      return $trees
}

