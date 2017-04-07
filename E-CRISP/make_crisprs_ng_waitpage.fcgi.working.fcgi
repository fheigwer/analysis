#!/usr/bin/perl

###################################################################################################################################################################################################
# setup the software's infrastructure - with FCGI this will be executed exactly 1 time
###################################################################################################################################################################################################

use warnings; #throw warnings for debugging the perl syntax mistakes
use strict; #require strict definition of every single variable, obey the accidential reuse of global variables
use FCGI; #Imports the library; required line
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::SeqFeature::Generic; #important package to handle sequence formats and objects
use Bio::Location::Split; #library to make splitted location objects
use Set::IntervalTree; #library providing methods for using interval trees
use JSON::XS qw(encode_json decode_json); #library to encode perl objects in a file
use File::Slurp qw(read_file write_file); #library to write perl objects to a file
use File::Copy::Recursive qw(dircopy);
use File::Basename;
use List::MoreUtils qw{
      any all none notall true false
      firstidx first_index lastidx last_index
      insert_after insert_after_string
      apply indexes
      after after_incl before before_incl
      firstval first_value lastval last_value
      each_array each_arrayref
      pairwise natatime
      mesh zip uniq distinct minmax part}; #some math and list utilities, later needed to partition a list of entries
use List::Util qw(sum);
use Archive::Zip;
use Bio::Graphics;
use Parallel::ForkManager; #important package to enable mutlithreading of the script
use CGI qw(:standard VARS);
use CGI::Carp qw ( fatalsToBrowser );

###################################################################################################################################################################################################

my %enzymelibrary       = (   "sevencutters"    => "enzymecollections/sevencutters.enz",
                              "sixcutters"      => "enzymecollections/sixcutters.enz",
                              "fivecutters"     => "enzymecollections/fivecutters.enz",
                              "fourcutters"     => "enzymecollections/fourcutters.enz",
                              "threecutters"    => "enzymecollections/threecutters.enz",
                              "twocutters"      => "enzymecollections/twocutters.enz",
                              "ecocutters"      => "enzymecollections/ecocutters.enz",
                              "hindcutters"     => "enzymecollections/hindcutters.enz",
                              "allcutters"      => "enzymecollections/enzymecollection.enz"
);
my $query               = "";
my $process             = "";
my $process_id          = "";
my $temp_dir            = "";
my $waiting_time        = "";
my $req                 = FCGI::Request(); 

###############################################################################################################################################################################################
# Response loop - will be executed everytime a request is send to the server
###################################################################################################################################################################################################

while ($req->Accept() >= 0) {
      CGI::initialize_globals();
      my %trees               = ();
      my $seqio_obj           = "";
      my %something           = ();
      my $parallel_number     = 2;
      my @enzarray            = ();
      my $query     = new CGI;
      $ENV{QUERY_STRING}=~m/PROCESS=(.*)&ID=(.*)&WAITINGTIME=(.*)&PID=(.*)/;
      $process       = $1;
      $temp_dir      = $2;
      $waiting_time  = $3;
      $process_id    = $4;
	$progress	   = $5;
      chdir("/var/www/E-CRISP/workdir");
      if ( $temp_dir eq "stuff") {
            $temp_dir = localtime.time;
            $temp_dir =~ s/\s+/\_/ig;
            mkdir( "$temp_dir") or die $!;
      } 
      system('chmod -R o+rwx /var/www/E-CRISP/workdir/'.$temp_dir.';');
      my $result = "/var/www/E-CRISP/workdir/$temp_dir/fertig.txt";
      my $error = "/var/www/E-CRISP/workdir/$temp_dir/error.html";
      if ( $waiting_time eq "stuff") {
            $waiting_time = 2;
      } else {
            $waiting_time = $waiting_time+1;
      }
    
      print "Content-type:text/html\r\n\r\n";
      open(my $header, "<", "header_CRISPR_nogoog.txt");
            while(<$header>){
                  print $_;
            }
      close $header;
      
      my $id_tag="";
      %something = ();
      foreach my $element (sort $query->param) {
            $something{$element} = $query->param($element);
            $id_tag.= $something{$element}."=".$query->param($element).";;";
            $id_tag=~s/\W/0/ig;
      }
      
      my $filename    = $something{"input_file"};
      my $seq_in_file = "";
      if ($filename) {
            $CGI::POST_MAX = 1024 * 20000;
            my $safe_filename_characters = "a-zA-Z0-9_.-";
            my ( $name, $path, $extension ) = fileparse( $filename, '\..*' );
            $filename = $name . $extension;
            $filename =~ tr/ /_/;
            $filename =~ s/[^$safe_filename_characters]//g;
            if ( $filename =~ /^([$safe_filename_characters]+)$/ ) {
                  $filename = $1;
            } else {
                  die "Filename contains invalid characters";
            }
            my $upload_filehandle = $query->upload("input_file");
            open( my $uploadfile, ">", "$temp_dir/$filename" ) or die "$!";
                  binmode $uploadfile;
                  while (<$upload_filehandle>) {
                        print $uploadfile $_;
                  }
            close $uploadfile;
            $seq_in_file = $temp_dir . "/" . $filename;
      } else {
            $seq_in_file = "";
      }
      
      #############################################################################################################################################################################################
      # If it is not the first call of the script just refresh the site and load a new funfact
      #############################################################################################################################################################################################
      
      if ( $process ne "stuff") {
            if ( -e $result ) {
                  print "<meta http-equiv=\"refresh\" content=\"0; URL=workdir/$temp_dir/index.html\">";
                  print end_html;
            } elsif ( -e $error ) {
                  print "<meta http-equiv=\"refresh\" content=\"0; URL=workdir/$temp_dir/error.html\">";
                  print end_html;
            }elsif(!(kill 0, $process_id)){
                  print_error_html( $temp_dir, "SORRY! <br><br>Some serious error occured during the process.<br>In case this happens again do not hesitate to contact us at crispr\@dkfz.de .\n" );
                  print "<meta http-equiv=\"refresh\" content=\"0; URL=workdir/$temp_dir/error.html\">";
                  print end_html;
            } else {
                  open (my $factfile, "<", "/var/www/E-CRISP/workdir/funfacts.txt") or die $!;
                        my @array=<$factfile>;
                        my $randomline=$array[rand @array];
                  close $factfile;
                  open (my $statsfile, "<", "/var/www/E-CRISP/workdir/$temp_dir/statsfile.txt") or die $!;
                        my @stats=(0,"");
                        while(<$statsfile>){
                              @stats=split("\t",$_);
                        }
                  close $statsfile;
			
			my $circle = (-710+$stats[0]*7.1);
			
                  print '<tr>
					<td style="width:800px" align=center>
						<br><br>
						<svg version="1.1" id="Layer_0" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" x="0px" y="0px" width="240px" height="240px" viewBox="0 0 240 240" enable-background="new 0 0 240 240" xml:space="preserve">
							<image x="20" y="20" width="200" height="200" xlink:href="/E-CRISP/workdir/E-CRISP.svg"/>
							
							<circle cx="120" cy="120" r="114" fill="none" stroke="#2662C3" stroke-width="10" stroke-dasharray="800" stroke-dashoffset="'.$circle.'" transform="rotate(-90, 120, 120)">
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
                                    <br>This page will be automatically updated every '.$waiting_time.' seconds until search is done.<br>
                              </td>
                        </tr>';
                  open(my $footer, "<", "/var/www/E-CRISP/workdir/footer.txt");
                        while(<$footer>){
                              print $_;
                        }
                  close $footer;
                  print "<meta http-equiv=\"refresh\" content=\"" . $waiting_time . "; URL=make_crisprs_ng_waitpage.fcgi?PROCESS=1&ID=$temp_dir&WAITINGTIME=$waiting_time&PID=$process_id\">";
                  print end_html;
            }
      } else {
            
            #######################################################################################################################################################################################
            # Start the logic only for the first call -> fork new childs for the calculations (each will get a part of the sequence)
            #######################################################################################################################################################################################
            
            my $databasepath = "/var/www/E-CRISP/DATABASEFILESNEW/" . $something{"ref_organism"} . "/" . $something{"ref_organism"};
            
            $req->Detach();
            unless( $process_id = fork){
                  $SIG{CHLD} = 'DEFAULT'; #reactivation of this default signal handling needed to avoid zombies and deadend processes (will be set to "IGNORE" in the else branche by the parent process)
                  #close needed, or content will be written into the wrong file
                  close(STDOUT);
                  close(STDERR);
                  open STDERR, ">", "/var/www/E-CRISP/workdir/$temp_dir/error.log"; #TODO eventuell an zentraler Stelle speichern und nicht im temp-dir (Idee: http://perl-howto.de/2008/06/sicheres-offnen-von-dateien-mi.html)
                  open STDOUT, ">", "/var/www/E-CRISP/workdir/$temp_dir/out.log";
                  
                  #################################################################################################################################################################################
                  # check if the calculation already have been done in the past and is stored in the buffer - if yes, just copy the existing folder and end the process
                  #################################################################################################################################################################################
                  
                  open(my $index, "<", "/var/www/E-CRISP/workdir/buffer/buffer.idx") or die $!;
                        my $tag_count  = 0;
                        my $line_count = 0;
                        my $tag_folder = "";
                        my $line_found = "";
                        my $access     = 0;
                        my $new_line   = "";
                        while (my $line = <$index>){
                              chomp $line;
                              $line_count++;
                              my @line=split("\t",$line);
                              if($line[0] eq $id_tag) {
                                    $tag_count++;
                                    $tag_folder = $line[1];
                                    $line_found = $line_count."d";
                                    $access = $line[2] + 1;
                                    $new_line = $line[0]."\t".$line[1]."\t".$access."\n";
                                    last; # end loop if found
                              }
                        }
                  close($index);
                  if( $tag_count>0) {
                        system("sed -i $line_found /var/www/E-CRISP/workdir/buffer/buffer.idx;");
                        open(my $indexfile, ">>", "/var/www/E-CRISP/workdir/buffer/buffer.idx") or die $!;
                              print $indexfile $new_line;
                        close($indexfile) or die $!;
                        
                        dircopy("/var/www/E-CRISP/workdir/buffer/$tag_folder",$temp_dir);
                        open(my $log, ">>", "/var/log/talecrisp.log");
                            print $log $temp_dir."\tE-CRISP\tDESIGN\tcopied\n";
                        close($log);
                        last;
                  }
                  my $runtime=time;
                  
                  #################################################################################################################################################################################
                  # upload a file and save it in the temp directory for bowtie
                  #################################################################################################################################################################################
                  open (my $statsfile, ">", "/var/www/E-CRISP/workdir/$temp_dir/statsfile.txt") or die $!;
                  print $statsfile "10\tsearch started\n";
		  close $statsfile;
                  if ( exists $something{"do_restriction_analysis"} ) { #check wether restriction enzyme anlysis is wanted at all
                        open (my $enzymefile, "<", "/var/www/E-CRISP/".$enzymelibrary{ $something{"enzymes"} }); #open the library file
                              while (<$enzymefile>) { #parse the file
                                    my @array = split( "\t", $_ ); #split the line after tab stops
                                    push @enzarray, \@array; #save a references to this line array in the enzyme library array
                              }
                        close $enzymefile;
                  }
                  if (exists $something{"sec_off_target"}) {
                        my $temp_count=0;
                        #built index for secondary off targets and store it for bowtie
                        my @sec_off_target_ids=$something{"sec_off_targets"};
                        my $sec_db = Bio::DB::Fasta->new( '/var/www/E-CRISP/DATABASEFILESNEW/secondary_off_targets.fasta');
                        open(my $tempsec, ">", $temp_dir."/temp_sec.fasta") or die $!;
                              foreach my $sec_off_target_id (@sec_off_target_ids){
                                    my $sec_seq_obj = $sec_db->get_Seq_by_id($sec_off_target_id);
                                    print $tempsec ">".$sec_seq_obj->display_id()."\n".$sec_seq_obj->seq()."\n";
                                    $temp_count++;
                              }
                              my $pasted_sec_seq = $something{"secondary_targets_seq"}."\n";
                              if($pasted_sec_seq=~m/\S/ig){
                                    if ($pasted_sec_seq=~m/^(>.+)\n([ACGT\s\n]+)\n$/) {
                                          my $header=$1;
                                          my $sec_seq=$2;
                                          $sec_seq=~s/[^ATGC]//ig;
                                          print $tempsec  "$header\n$sec_seq\n";
                                          $temp_count++;
                                    }else{
                                          print_error_html( $temp_dir, "\\\"$pasted_sec_seq\\\" is not a FASTA format sequence\n" );
                                          die;
                                    }
                              }
                        close $tempsec;
                        if ($temp_count==0) {
                              print_error_html( $temp_dir, "No secondary off-target sequence has been specified.\n Please choose one from the presets or paste one in the text area.\n" );
                              die;
                        }
                        system('bowtie2-build '.$temp_dir.'/temp_sec.fasta '.$temp_dir.'/temp_sec ;');
                  }
                  
                  #################################################################################################################################################################################
                  # For ENSEMBLE: define the path to the bowtie index and do checks if an database entry is found for the ensemble accesion number - if all checks passed, create the $seqio_obj
                  #################################################################################################################################################################################
                  #print( "shittistiff".$something{"GENE.SYMB"}."stuff\n");
                  if ( exists $something{"GENE.SYMB"}) {
                        my $db = Bio::DB::Fasta->new( $databasepath . '.all.dna.fa', -makeid => \&make_my_id );
                        if ( !$filename ) {
                              my $temp = $something{"pasted_seq"};
                              $temp =~ s/\s+//ig;
                              my @ids = split( ";", $temp );
                              
                              make_temp_fasta_file(\@ids, \%trees, \%something, $db, $temp_dir, 0);
                              $seqio_obj = Bio::SeqIO->new( -file => $temp_dir . "/tempfile.fasta", -format => "fasta" ); #read the temporary fasta file
                        } 
                  } else { 
                        
                        #################################################################################################################################################################################
                        # For FASTA: define the path to the bowtie index and do checks if the fasta sequence/file is in the right format - if all checks passed, create the $seqio_obj
                        #################################################################################################################################################################################
                        
                        if ( !$filename ) { #check if the pasted sequence should be used
                              my $pasted = $something{"pasted_seq"}; #get the sequence(s) from the textentry window
                              $pasted=~tr/acgt/ACGT/;
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
                                          print_error_html( $temp_dir, "$line is not a FASTA format file because it contains &quot;$1&quot; as bases (only ATCGUN are allowed).<br>");
                                          die;
                                    }
                              }
                              if ($count==0){
                                    print_error_html($temp_dir, "No FASTA sequence has been entered.<br>FASTA sequences need a header starting with &quot;&gt;&quot;, which is new line separated from the sequence.<br><br> e.g.: &gt;some Sequence<br>AGCTGATCGATCTAGCTAGCTGCTAGCTAGTCGATCGATGCTAGCTAGCTAGCTGCTAG<br>");
                                    die;
                              }   
                              if ($count>50) {
                                   print_error_html( $temp_dir, "Your input is more than 50 Sequences\n" );
                                    die;
                              }
                              if ( !( $something{"specific_transcript"} eq "any") && $count >1) {
                                    print_error_html( $temp_dir, "Transcript specificity is only defined for single gene analyses.\n" );
                                    die;
                              }
                              
                              open (my $tempfile, ">", $temp_dir . "/tempfile.fasta"); #write it to an temporary fasta format file
                                    print $tempfile $pasted;
                              close $tempfile;
                              $seqio_obj = Bio::SeqIO->new( -file => $temp_dir . "/tempfile.fasta", -format => "fasta" ) or die print_error_html($temp_dir, "No FASTA sequence has been entered.<br>FASTA sequences need a header starting with &quot;&gt;&quot;, which is new line separated from the sequence.<br><br> e.g.: &gt;some Sequence<br>AGCTGATCGATCTAGCTAGCTGCTAGCTAGTCGATCGATGCTAGCTAGCTAGCTGCTAG<br>");
                        } 
                  }
                  
                  #################################################################################################################################################################################
                  # Start the creation of the report (index.html)
                  #################################################################################################################################################################################
                  
                  #define empty object to be filled in the process
                  my $fname               = "";
                  my $dont_asses_context  = 0;
                  open( my $report, ">", "$temp_dir/index.html" ) or die "can not open report file";
                        open(my $header, "<", "/var/www/E-CRISP/workdir/header_CRISPR.txt");
                              while(<$header>){
                                    print $report $_;
                              }
                        close $header;
                        print $report '   <tr>
                                                <td>
                                                      <br><br>
                                                      <a href="all_results_together.tab" target="_blank" download="all_results_together.tab"><input type="button" align="center" value="Download a tabular report for all query sequences together"></a>
                                                      <br>
                                                      <a href="all_results_together.gff" target="_blank" download="all_results_together.gff"><input type="button" align="center" value="Download a GFF-File for all query sequences together"></a>
                                                      <br>';
		  if(-e "$temp_dir/failfile.tab"){ print $report '<a href="failfile.tab" target="_blank" download="failfile.tab"><input type="button" align="center" value="Get the table of genes that could not be found"></a>
                                                       ';}
			print $report '<br><br></td>
                                          </tr>';
                        print $report '   <tr>
                                                <td colspan="2" rowspan="1" style="vertical-align: top;">
                                                      <embed style="width: 100%; height: 10px;" alt="there should appear a  line" src="/E-CRISP/gradientline.svg" type="image/svg+xml"><br>
                                                </td>
                                          </tr>
                                          ';
                  close $report;
                  chmod 0755, $temp_dir . "/index.html";
                  
                  #################################################################################################################################################################################
                  #start the main loop which is looping through the sequence found in the SeqIO object, may be one or many more
                  #################################################################################################################################################################################
                  
                  my %statistics    = ();
                  my %CRISPR_hash   = ();
                  my %CRISPR_cnt    = ();
                  my @seq_array     = ();
                  
                  while( my $seq = $seqio_obj->next_seq() ) {
                        push(@seq_array,$seq);
                  }
                  open ($statsfile, ">", "/var/www/E-CRISP/workdir/$temp_dir/statsfile.txt") or die $!;
                  print $statsfile "20\tsequences were read in\n";
		  close $statsfile;
                  foreach my $seq_obj ( @seq_array ) {
                        #find the genomic coordinates of the sequence either directly from the fasta sequence or from other smyces
                        my $chrom               = ""; #create an empty name chromosome
                        my $location_offset     = 0;
                        my $location_end        = 0; #set the location offset value to 0
                        if ( ( $seq_obj->description ) =~ m/chrom:([\w\.]+):(\d+)..(\d+)/ig) { #if the descrition of the sequence contains information of the form loc=chr1:0..13000000
                              #print the result of the pattern matching
                              $chrom = $1; #save the location information in the certain objects
                              $location_offset = $2;
                              $location_end = $3; #the location offset is the start of the sequence on that certain chromosome, for the whole chromosome this is of cause 0
                              if ( !exists $trees{$chrom} ) {
                                    $trees{$chrom} = build_tree( "/var/www/E-CRISP/DATABASEFILESNEW/" . $something{"ref_organism"} . "/" . $chrom . "_indexed" );
                              }
                        } else {
                              $dont_asses_context = 1;
                        }
                        #define the name for the current sequence as its display id
                        $fname                                                = ( $seq_obj->display_id );#current name is the sequence' id
                        $statistics{$fname}{"seq_name"}                       = $fname; 
                        $statistics{$fname}{"seq_length"}                     = $seq_obj->length;
                        $statistics{$fname}{"seq_location"}                   = $chrom."::".$location_offset."::".$location_end;
                        $statistics{$fname}{"Number of successful designs"}   = 0;
                        if ( $chrom eq "" ) { $chrom = $fname; } #if $chrom is still empty fill it with the sequence' id
                        my $whole_seq = $seq_obj->seq; #deduce the complete nucleotide sequence as a alphanumerical string from the SeqIo object
                        if (int(log($statistics{$fname}{"seq_length"})) <= 8) {
                              if (int(log($statistics{$fname}{"seq_length"}))>=2) {
                                     $parallel_number=int(log($statistics{$fname}{"seq_length"}));
                              }else{
                                     $parallel_number=2;
                              }
                        }else{
                               $parallel_number=8;
                        }
                        
                        ###########################################################################################################################################################################
                        #create the hashes for the CRISPR and the statistics
                        ###########################################################################################################################################################################
                        open (my $statsfile, ">", "/var/www/E-CRISP/workdir/$temp_dir/statsfile.txt") or die $!;
                  print $statsfile "20\t$fname will now be scanned for target sites\n";
		  close $statsfile;
                        if (!exists $CRISPR_hash{$fname}) { #only execute following calucaltions if it is not already done for the sequence
                              my @findings = find_and_print_CRISPRS(    \$seq_obj,
                                                                        $chrom,
                                                                        $location_offset,
                                                                        \%trees,
                                                                        $dont_asses_context,
                                                                        $temp_dir,
                                                                        $parallel_number,
                                                                        \%something);
                              %{$CRISPR_hash{$fname}} = %{$findings[0]};
                              %{$statistics{$fname}} = (%{$statistics{$fname}},%{$findings[1]});
			      open (my $statsfile, ">", "/var/www/E-CRISP/workdir/$temp_dir/statsfile.txt") or die $!;
                  print $statsfile "50\ttarget search finished\n";
		  close $statsfile;
                        }
                        
                        $CRISPR_cnt{$fname} = 0;
                  } # end the loop after Hash and statistic were created, so the Bowtie Magic will be executed only once
                  
                  #################################################################################################################################################################################
                  #Bowtie Magic for the found CRISPR
                  #################################################################################################################################################################################
                  
                  if ( exists $something{"exclude_unspecific"} ) {
                        open (my $statsfile, ">", "/var/www/E-CRISP/workdir/$temp_dir/statsfile.txt") or die $!;
                  print $statsfile "50\tdesigns are now evaluated for off-target effects\n";
		  close $statsfile;
                        ###########################################################################################################################################################################
                        #Bowtie for single sequence
                        ###########################################################################################################################################################################
                        
                        if ($something{"kind"} eq "single") {
                              open (my $crisprs, ">", $temp_dir . "/temp_CRISPRS.fasta");
                                    foreach my $seq ( sort( keys(%CRISPR_hash) ) ) {
                                          foreach my $key ( sort( keys(%{$CRISPR_hash{$seq}}) ) ) {
                                                foreach my $letter ("A","C","G","T"){
                                                           print $crisprs "\>" . $key . "\n";
                                                            if(${ ${ $CRISPR_hash{$seq} } {$key} }{"strand"} eq "minus"){
                                                                  print $crisprs substr(reverse_comp(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"}),$something{"unspecific_leading_bases"},(length(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"})-3-$something{"unspecific_leading_bases"})).$letter.substr(reverse_comp(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"}),(length(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"})-2)). "\n"; #$whole_CRISPR_seq
                                                            }else{
                                                                 print $crisprs substr(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"},$something{"unspecific_leading_bases"},(length(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"})-3-$something{"unspecific_leading_bases"})).$letter.substr(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"},(length(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"})-2)). "\n"; #$whole_CRISPR_seq
                                                            }
                                                }
                                          }
                                    }
                              close $crisprs;
                              
                              #####################################################################################################################################################################
                              #teemp_sec
                              if ($something{"bowtie_version"} eq "bowtie2") {
                                    if ($something{"offtargetdb"} eq "gDNA") {
                                          system( '/usr/bin/bowtie2 -p 4 -f -k 30 --'.$something{"bowtie_mode"}.' --end-to-end --no-hd -x ' . $databasepath .".dna". ' -U ' . $temp_dir . '/temp_CRISPRS.fasta > ' . $temp_dir . '/temp_out.bwt' );
                                    }elsif($something{"offtargetdb"} eq "cDNA"){
                                          system( '/usr/bin/bowtie2 -p 4 -f -k 30 --'.$something{"bowtie_mode"}.' --end-to-end --no-hd -x ' . $databasepath.".cdna" . ' -U ' . $temp_dir . '/temp_CRISPRS.fasta > ' . $temp_dir . '/temp_out.bwt' );
                                    }else{
                                          system( '/usr/bin/bowtie2 -p 4 -f -k 30 --'.$something{"bowtie_mode"}.' --end-to-end --no-hd -x ' . $databasepath.".genome" . ' -U ' . $temp_dir . '/temp_CRISPRS.fasta > ' . $temp_dir . '/temp_out.bwt' );
                                    }
                              }else{
                                    if ($something{"offtargetdb"} eq "gDNA") {
                                          system( '/usr/bin/bowtie ' . $databasepath .".dna". ' ' . $temp_dir . "/" . 'temp_CRISPRS.fasta -f -v 3 -y -k 30 -S --sam-nohead --sam-nosq -p 4  > ' . $temp_dir . '/temp_out.bwt' );
                                    }elsif($something{"offtargetdb"} eq "cDNA"){
                                          system( '/usr/bin/bowtie ' . $databasepath .".cdna". ' ' . $temp_dir . "/" . 'temp_CRISPRS.fasta -f -v 3 -y -k 30 -S --sam-nohead --sam-nosq -p 4  > ' . $temp_dir . '/temp_out.bwt' );
                                    }else{
                                          if (-e $databasepath.".genome.1.ebwtl") {
                                                system( '/usr/bin/bowtie ' . $databasepath.".genome" . ' ' . $temp_dir . "/" .'temp_CRISPRS.fasta -f -v 3 -y -k 30 -S --sam-nohead --large-index  --sam-nosq -p 4  > ' . $temp_dir .'/temp_out.bwt' );
                                          }else{
                                                system( '/usr/bin/bowtie ' . $databasepath.".genome" . ' ' . $temp_dir . "/" .'temp_CRISPRS.fasta -f -v 3 -y -k 30 -S --sam-nohead --sam-nosq -p 4  > ' . $temp_dir .'/temp_out.bwt' );
                                          }
                                          
                                          
                                    }
                              }
                              open (my $bowtie, "<", $temp_dir . "/temp_out.bwt");
                                    while (my $line = <$bowtie>) {
                                          chomp $line;
                                          my @line = split( "\t", $line );
                                          if ($line[2] eq "*") {
                                                $line[0] =~m/(\S+)_(\S+)_/;
                                                my $seq = $1;
                                                ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=";;".$line[2]."い0い0いNAい0いNA";
                                          }else{
                                                $line=~m/NM:i:(\d+)/;
                                                my $edit_distance=$1;
                                                if (($edit_distance <= $something{"edit_distance_allowed"})) {
                                                      my @line = split( "\t", $line );
                                                      #decide if it is forward (fw) or backward (bw) query-sequence
                                                      my $direction = "rc";
                                                      if ($line[1] == 0 || $line[1] == 256) {
                                                            $direction = "fw";
                                                      }
                                                      if ( $line[0] =~ m/([^_]+_[^_]+_[^_]+)/ig ) {
                                                            $line[0] =~m/(\S+)_(\S+)_/;
                                                            my $seq = $1;
                                                            my @matchstringo=make_mismatch_string (\$line,$something{"unspecific_leading_bases"}, $direction);
                                                            if ( ($direction eq "fw" && $matchstringo[scalar(@matchstringo)-1] ne "X" && $matchstringo[scalar(@matchstringo)-2] ne "X" && $matchstringo[scalar(@matchstringo)-3] ne "X")
                                                                || ($direction eq "rc" && $matchstringo[0] ne "X" && $matchstringo[1] ne "X" && $matchstringo[2] ne "X")
                                                                ) {
                                                                  my $startcoordinate=0;
                                                                  if ($something{"offtargetdb"} eq "genomicDNA") {
                                                                        my $namestuff="";
                                                                        if ( !exists $trees{$line[2]} ) { #TODO if und else fast identisch, aber relativ kurzer Part und daher Funktion performancetechnisch nachteilig
                                                                              $trees{$line[2]} = build_tree( "/var/www/E-CRISP/DATABASEFILESNEW/" . $something{"ref_organism"} . "/" . $line[2] . "_indexed" );
                                                                        }
                                                                        my $annotations = $trees{$line[2]}->fetch( int($line[3]), int(($line[3])) );
                                                                        foreach  my $anno ( @{$annotations} ) {
                                                                              if ( $anno =~ m/gene_(\S+)_([0-9]+)_([0-9]+)/ ) {
                                                                                    $namestuff=$1;
                                                                                    $startcoordinate=$2;
                                                                              }
                                                                        }
                                                                        if ($namestuff eq "" && exists($something{"ignore_intergenic"})) {                                                                              
                                                                        }elsif($namestuff ne ""){
                                                                              ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=";;".$namestuff."い".($line[3]-$startcoordinate)."い".($line[3]+@matchstringo-$startcoordinate)."い".join("",@matchstringo)."い".$edit_distance."い".$direction;
                                                                        }elsif($namestuff eq "" && !exists($something{"ignore_intergenic"})){
                                                                              ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=";;".$line[2]."い".($line[3]-$startcoordinate)."い".($line[3]+@matchstringo-$startcoordinate)."い".join("",@matchstringo)."い".$edit_distance."い".$direction;
                                                                        }
                                                                  }else{
                                                                        ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=";;".$line[2]."い".($line[3])."い".($line[3]+@matchstringo)."い".join("",@matchstringo)."い".$edit_distance."い".$direction;
                                                                  }
                                                            }
                                                      }    
                                                }
                                          }
                                    }
                              close $bowtie;
                              unlink $temp_dir . "/temp_CRISPRS.fasta";
                              unlink $temp_dir . "/temp_out.bwt";
                              open (my $statsfile, ">", "/var/www/E-CRISP/workdir/$temp_dir/statsfile.txt") or die $!;
                                    print $statsfile "80\tsequences were aligned\n";
                              close $statsfile;
                              if ( exists $something{"sec_off_target"} ) { #ckeck if this is wanted
                                    open (my $crisprs, ">", $temp_dir . "/temp_CRISPRS.fasta");
                                          foreach my $seq ( sort( keys(%CRISPR_hash) ) ) {
                                                foreach my $key ( sort( keys(%{$CRISPR_hash{$seq}}) ) ) {
                                                      print $crisprs "\>" . $key . "\n";
                                                      print $crisprs ${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"} . "\n"; #$whole_CRISPR_seq
                                                }
                                          }
                                    close $crisprs;
                                    
                                    ###############################################################################################################################################################
                                    
                                    #do send a bowtie2 job fot the two temporary written fasta files as if they were paired seqencing reads and save the result in a ~out.bwt file
                                    system( '/usr/bin/bowtie2 -p 4  -f -k 30 --'.$something{"bowtie_mode"}.' --end-to-end --no-hd -x /var/www/E-CRISP/workdir/' . $temp_dir .'/temp_sec -U ' . $temp_dir . '/temp_CRISPRS.fasta > ' . $temp_dir . '/temp_out.bwt' );
                                    open (my $bowtie, "<", $temp_dir . "/temp_out.bwt");
                                          my $edit_distance=0;
                                          while (my $line = <$bowtie>) {
                                                chomp $line;
                                                my @line = split( "\t", $line );
                                                if($line=~m/NM:i:(\d+)/){
                                                      $edit_distance=$1;
                                                }
                                                if (($edit_distance <= $something{"edit_distance_allowed"}) && ($line[2] ne "*")) {
                                                      if ( $line[0] =~ m/([^_]+_[^_]+_[^_]+)/ig ) {
                                                            my $key = $1;
                                                            $line[0] =~m/(\S+)_(\S+)_/;
                                                            my $seq = $1;
                                                            push @{ ${ ${ $CRISPR_hash{$seq} } {$key} }{"sec_hits"} }, $line[2];
                                                      }
                                                }
                                          }
                                    close $bowtie;
                                    unlink $temp_dir . "/temp_CRISPRS.fasta";
                                    unlink $temp_dir . "/temp_out.bwt";
                              }
                        }else{
                              
                              #####################################################################################################################################################################
                              #Bowtie for double sequence
                              #####################################################################################################################################################################
                              
                              open (my $leftcrisprs, ">", $temp_dir . "/temp_LEFTCRISPRS.fasta");
                                    open (my $rightcrisprs, ">", $temp_dir . "/temp_RIGHTCRISPRS.fasta");
                                          foreach my $seq ( sort( keys(%CRISPR_hash) ) ) {
                                                foreach my $key ( sort( keys(%{$CRISPR_hash{$seq}}) ) ) {
                                                      my $stringleft = @{${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"}}[0];
                                                      my $stringright = @{${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"}}[1];
                                                      foreach my $left_PAM ("CCA","CCC","CCT","CCG","CTA","CTC","CTT","CTG"){
                                                            foreach my $right_PAM ("CCA","CCC","CCT","CCG","CTA","CTC","CTT","CTG"){
                                                                  print $leftcrisprs "\>" . $key . "\n";
                                                                  print $leftcrisprs $left_PAM.substr(@{${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"}}[0], 3, (length($stringleft) -3 - $something{"unspecific_leading_bases"})) . "\n"; #TODO am Ende entfernen (...,0,leading)
                                                                  print $rightcrisprs "\>" . $key . "\n";
                                                                  print $rightcrisprs $right_PAM.substr(reverse_comp(@{${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"}}[1]), 3, (length($stringright) -3 - $something{"unspecific_leading_bases"})) . "\n"; #TODO am Ende entfernen (...,0,leading)
                                                            }
                                                      }
                                                }
                                          }
                                    close ($leftcrisprs);
                              close ($rightcrisprs);
                              
                              #####################################################################################################################################################################
                              #teemp_sec
                              if ($something{"offtargetdb"} eq "gDNA") {
                                    system( '/usr/bin/bowtie2 -p 4 -f -k 30 --'.$something{"bowtie_mode"}.' --end-to-end --no-mixed --no-discordant --no-hd -x ' . $databasepath .".dna". ' -I '.$something{"minspacerlength"}.' -X 100 -1 ' . $temp_dir . '/temp_LEFTCRISPRS.fasta'. ' -2 ' . $temp_dir . '/temp_RIGHTCRISPRS.fasta'.' > ' . $temp_dir . '/temp_out.bwt' );
                              }elsif($something{"offtargetdb"} eq "cDNA"){
                                    system( '/usr/bin/bowtie2 -p 4 -f -k 30 --'.$something{"bowtie_mode"}.' --end-to-end --no-mixed --no-discordant --no-hd -x ' . $databasepath.".cdna". ' -I '.$something{"minspacerlength"}.' -X 100 -1 ' . $temp_dir . '/temp_LEFTCRISPRS.fasta'. ' -2 ' . $temp_dir . '/temp_RIGHTCRISPRS.fasta'.' > ' . $temp_dir . '/temp_out.bwt' );
                              }else{
                                    system( '/usr/bin/bowtie2 -p 4 -f -k 30 --'.$something{"bowtie_mode"}.' --end-to-end --no-mixed --no-discordant --no-hd -x ' . $databasepath.".genome".'  -I '.$something{"minspacerlength"}.' -X 100 -1 ' . $temp_dir . '/temp_LEFTCRISPRS.fasta'. ' -2 ' . $temp_dir . '/temp_RIGHTCRISPRS.fasta'.' > ' . $temp_dir . '/temp_out.bwt' );
                              }
                              open (my $bowtie, "<", $temp_dir . "/temp_out.bwt");
                              my $edit_distance=0;
                              my $was_hit=0;
                                    while (my $line = <$bowtie>) {
                                          chomp $line;
                                          my @line = split( "\t", $line );
                                          if ($line[2] ne "*") {
                                                #decide if it is forward (fw) or backward (rc) query-sequence - only for the first pair, the 2nd will be the oposite - Workaround: oposite interpretation of the flag-numbers
                                                my $direction           = "fw";
                                                if ($line[1]==97 || $line[1]==99 || $line[1]==355 || $line[1]==161 || $line[1]==163 || $line[1]==419)  {
                                                      $direction        = "rc";
                                                }
                                                if ( $line[0] =~ m/([^_]+_[^_]+_[^_]+)/ig ) {
                                                      $line[0] =~m/(\S+)_(\S+)_/;
                                                      my $seq = $1;
                                                      my @matchstringo=make_mismatch_string (\$line,$something{"unspecific_leading_bases"}, $direction);
                                                      my $startcoordinate=0;
                                                      my $spacer = abs($line[8]) - ((abs($line[8]) - abs($line[3] - $line[7])) * 2);
                                                      if ($something{"offtargetdb"} eq "genomicDNA") {                                                      
                                                            if ($line[1]==81 || $line[1]==83 || $line[1]==97 || $line[1]==99 || $line[1]==339 || $line[1]==355) {
                                                                  $line=~m/NM:i:(\d+)/;
                                                                  $edit_distance=$1;
                                                                  if (($edit_distance <= $something{"edit_distance_allowed"})) {
                                                                        if ( (($direction eq "fw" && $matchstringo[scalar(@matchstringo)-1] ne "X" && $matchstringo[scalar(@matchstringo)-2] ne "X" && $matchstringo[scalar(@matchstringo)-3] ne "X")
                                                                              || ($direction eq "rc" && $matchstringo[0] ne "X" && $matchstringo[1] ne "X" && $matchstringo[2] ne "X"))
                                                                           ) {
                                                                              my $namestuff="";
                                                                              if ( !exists $trees{$line[2]} ) { #TODO if und else fast identisch, aber relativ kurzer Part und daher Funktion performancetechnisch nachteilig
                                                                                    $trees{$line[2]} = build_tree( "/var/www/E-CRISP/DATABASEFILESNEW/" . $something{"ref_organism"} . "/" . $line[2] . "_indexed" );
                                                                              }
                                                                              my $annotations = $trees{$line[2]}->fetch( int($line[3]), int(($line[3])) );
                                                                              foreach  my $anno ( @{$annotations} ) {
                                                                                    if ( $anno =~ m/gene_(\S+)_([0-9]+)_([0-9]+)/ ) {
                                                                                          $namestuff=$1;
                                                                                          $startcoordinate=$2;
                                                                                    } 
                                                                              }
                                                                              if ($namestuff eq "") {
                                                                                    ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=";;".$line[2]."い".($line[3]-$startcoordinate)."い".($line[3]+@matchstringo-$startcoordinate)."い".join("",@matchstringo)."い".$edit_distance."い".$direction."い".$spacer;
                                                                              }else{
                                                                                    ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=";;".$namestuff."い".($line[3]-$startcoordinate)."い".($line[3]+@matchstringo-$startcoordinate)."い".join("",@matchstringo)."い".$edit_distance."い".$direction."い".$spacer;
                                                                              }
                                                                              $was_hit=1;
                                                                        }
                                                                  }
                                                            }elsif($line[1]==161 || $line[1]==163 || $line[1]==145 || $line[1]==147 || $line[1]==403 || $line[1]==419){
                                                                  if ($was_hit==1) {                                                                  
                                                                        $line=~m/NM:i:(\d+)/;
                                                                        $edit_distance=$1;
                                                                        if ( (($direction eq "fw" && $matchstringo[scalar(@matchstringo)-1] ne "X" && $matchstringo[scalar(@matchstringo)-2] ne "X" && $matchstringo[scalar(@matchstringo)-3] ne "X")
                                                                              || ($direction eq "rc" && $matchstringo[0] ne "X" && $matchstringo[1] ne "X" && $matchstringo[2] ne "X"))
                                                                              && ($edit_distance <= $something{"edit_distance_allowed"})
                                                                        ) {
                                                                                    my @hitarray=split(";;",${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"});
                                                                                    my @lasthitarray=split("い",$hitarray[-1]);
                                                                                    $lasthitarray[3].="-".join("",@matchstringo);
                                                                                    $hitarray[-1]=join("い",@lasthitarray);
                                                                                    ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}=join(";;",@hitarray);
                                                                        }else{
                                                                                    my @hitarray=split(";;",${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"});
                                                                                    pop @hitarray;
                                                                                    ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}=join(";;",@hitarray);
                                                                        }
                                                                        $was_hit=0;
                                                                  }
                                                            }                                                  
                                                      }else{
                                                            if ($line[1]<147) {
                                                                  ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=";;".$line[2]."い".($line[3]-$startcoordinate)."い".($line[3]+@matchstringo-$startcoordinate)."い".join("",@matchstringo)."い".$edit_distance."い".$direction."い".$spacer;
                                                            }elsif($line[1]==147){
                                                                  my @hitarray=split(";;",${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"});
                                                                  my @lasthitarray=split("い",$hitarray[-1]);
                                                                  $lasthitarray[3].="-".join("",@matchstringo);
                                                                  $hitarray[-1]=join("い",@lasthitarray);
                                                                  ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}=join(";;",@hitarray);
                                                            }else{
                                                                  ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=";;".$line[2]."い".($line[3]-$startcoordinate)."い".($line[3]+@matchstringo-$startcoordinate)."い".join("",@matchstringo)."い".$edit_distance."い".$direction."い".$spacer;
                                                            }
                                                            $was_hit=1;
                                                      }
                                                }                                                     
                                          }else{
                                                #$line[0] =~m/(\S+)_(\S+)_/;
                                                #my $seq = $1;
                                                #${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=";;".$line[2]."い0い0いNAい0いNA";
                                          }
                                    }
                              close $bowtie;
                              unlink $temp_dir . "/temp_LEFTCRISPRS.fasta";
                              unlink $temp_dir . "/temp_RIGHTCRISPRS.fasta";
                              unlink $temp_dir . "/temp_out.bwt";
                              open (my $statsfile, ">", "/var/www/E-CRISP/workdir/$temp_dir/statsfile.txt") or die $!;
                              print $statsfile "80\tsequences were aligned\n";
                              close $statsfile;
                              #####################################################################################################################################################################
                              #  if ( exists $something{"sec_off_target"} ) needed for double! @Florian
                              #####################################################################################################################################################################
                        }
                  }
                  
                  ###########################################################################################################################################################################
                  #Evaluate the infos of the CRISPR-Hash (additional information, changes)
                  ###########################################################################################################################################################################
                  
                  my $number_of_hits = 1;
                  SEQLOOP: foreach my $fname ( keys(%CRISPR_hash) ) {
                        CRISPRHASHLOOP: foreach my $key ( keys(%{$CRISPR_hash{$fname}}) ) {
                              if ( exists $something{"exclude_unspecific"} ){
                                    if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"hits"} ) {
                                          $number_of_hits = scalar(split(";;",${ ${ $CRISPR_hash{$fname} } {$key} }{"hits"}));
                                          ${ ${ $CRISPR_hash{$fname} } {$key} }{"number_of_hits"} = $number_of_hits-1;
                                          if (( $number_of_hits > $something{"off-targets-allowed"} + 2 || $number_of_hits < 1 ) ) {
                                                delete $CRISPR_hash{$fname}{$key};
                                                $statistics{$fname}{"Number of designs excluded because they hit multiple targets or none"}++;
                                                next CRISPRHASHLOOP;
                                          }else{
                                                $statistics{$fname}{"Number of designs that hit a specific target"}++;
                                          }
                                    }else{
                                          delete $CRISPR_hash{$fname}{$key};
                                          next CRISPRHASHLOOP;
                                    }
                              }
                              
                              if ( exists $something{"sec_off_target"} ){
                                    if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"sec_hits"} ) {
                                          $number_of_hits = @{ ${ ${ $CRISPR_hash{$fname} } {$key} }{"sec_hits"} };
                                          if (( $number_of_hits > 0 ) ) {
                                                delete $CRISPR_hash{$fname}{$key};
                                                $statistics{$fname}{"Number of designs excluded because they hit a secondary offtarget"}++;
                                                next CRISPRHASHLOOP;
                                          }else{
                                                $statistics{$fname}{"Number of designs that do not hit a secondary target"}++;
                                          }
                                    }else{
                                          $statistics{$fname}{"Number of designs that do not hit a secondary target"}++;
                                    }
                              }
                              
                              my $whole_crisp_seq = "";
                              if($something{"kind"} ne "single"){
                                    $whole_crisp_seq = join( "", @{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}});
                              }else{
                                    $whole_crisp_seq = ${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"};
                              }
                              ${ ${ $CRISPR_hash{$fname} } {$key} }{"Nuc_comp"} = join( " ", my @basecomp = find_base_count( $whole_crisp_seq ) ); #store the nucleotide composition as a string object in the hash
                              ${ ${ $CRISPR_hash{$fname} } {$key} }{"score"} = 120;
                              ${ ${ $CRISPR_hash{$fname} } {$key} }{"exon"} = join( "_", keys( %{ ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} } {"exon"} } ) );
                        }
                        
                        my %exon_count = ();
                        EXONLOOP: foreach my $key ( sort { $CRISPR_hash{$fname}{$a}->{"exon"} cmp $CRISPR_hash{$fname}{$b}->{"exon"} } sort { $CRISPR_hash{$fname}{$b}->{"score"} <=> $CRISPR_hash{$fname}{$a}->{"score"} } keys(%{$CRISPR_hash{$fname}}) ) {
                              my $exon = ${ ${ $CRISPR_hash{$fname} } {$key} }{"exon"};
                              $exon_count{$exon}++;
                              if ( $exon_count{$exon} > $something{"max_per_exon"} ) {
                                    delete $CRISPR_hash{$fname}{$key};
                                    $statistics{$fname}{"Number of designs excluded because the maximum of designs per exon was exceeded"}++;
                                    next EXONLOOP;
                              }
                        }
                  }
                  
                  #################################################################################################################################################################################
                  #reopen the main loop which is looping through the sequence found in the SeqIO object, may be one or many more
                  #################################################################################################################################################################################
                  
                  foreach my $seq_obj ( @seq_array ) {
                        #reinitiate values...
                        my $chrom               = ""; #create an empty name chromosome
                        my $location_offset     = 0;
                        my $location_end        = 0; #set the location offset value to 0
                        if ( ( $seq_obj->description ) =~ m/chrom:([\w\.]+):(\d+)..(\d+)/ig) { #if the descrition of the sequence contains information of the form loc=chr1:0..13000000
                              #print the result of the pattern matching
                              $chrom = $1; #save the location information in the certain objects
                              $location_offset = $2;
                              $location_end = $3; #the location offset is the start of the sequence on that certain chromosome, for the whole chromosome this is of cause 0
                        } else {
                              $dont_asses_context = 1;
                        }
                        #define the name for the current sequence as its display id
                        $fname = ( $seq_obj->display_id ); #current name is the sequence' id
                        if ( $chrom eq "" ) { $chrom = $fname; } #if $chrom is still empty fill it with the sequence' id
                        my $whole_seq = $seq_obj->seq; #deduce the complete nucleotide sequence as a alphanumerical string from the SeqIo object
                        
                        ###########################################################################################################################################################################
                        #Build te image of the gene
                        ###########################################################################################################################################################################
                        
                        if ($CRISPR_cnt{$fname} == 0) { #only execute following calucaltions if it is not already done for the sequence
                              my %transcripts_hash=();
                              my %CDS_hash      =();
                              if ( exists $trees{$chrom} ) {
                                    my $annotations   = $trees{$chrom}->fetch( int($location_offset), int($location_end) );
                                    my %cpg_hash      =();                                   
                                    foreach my $anno ( sort( @{$annotations} ) ) {
                                          if ( $anno =~ m/gene_(.+)_(\d+)_(\d+)$/ig ) {
                                                my $feat = Bio::SeqFeature::Generic->new(
                                                            -display_name => "gene",
                                                            -primary_tag => "gene::" . $1,
                                                            -start => $2 - $location_offset,
                                                            -strand => 1,
                                                            -end => $3 - $location_offset
                                                );
                                                $seq_obj->add_SeqFeature($feat);
                                          } elsif ( $anno =~ m/exon::(\S+)::(.+)_(\d+)_(\d+)$/ig ) {
                                                my @pair = ( $3, $4 );
                                                ${$transcripts_hash{$1}}{$2}=\@pair;
                                          } elsif ( $anno =~ m/CDS::(\S+)::(.+)_(\d+)_(\d+)$/ig ) {
                                                my @pair = ( $3, $4 );
                                                push @{ $CDS_hash{$1} }, \@pair;
                                          } elsif ( $anno =~ m/start_codon::(\S+)::(.+)_(\d+)_(\d+)$/ig ) {
                                                my $feat = Bio::SeqFeature::Generic->new(
                                                            -display_name => "start_codon",
                                                            -primary_tag => "start_codon::" . $1,
                                                            -strand => 1,
                                                            -start => $3 - $location_offset,
                                                            -end => $4 - $location_offset
                                                );
                                                $seq_obj->add_SeqFeature($feat);
                                          } elsif ( $anno =~ m/stop_codon::(\S+)::(.+)_(\d+)_(\d+)$/ig ) {
                                                my $feat = Bio::SeqFeature::Generic->new(
                                                            -display_name => "stop_codon",
                                                            -primary_tag => "stop_codon::" . $1,
                                                            -strand => 1,
                                                            -start => $3 - $location_offset,
                                                            -end => $4 - $location_offset
                                                );
                                                $seq_obj->add_SeqFeature($feat);
                                          } elsif ( $anno =~ m/CpG_(\d+)_(\d+).*/ig ) {
                                                my @pair=($1,$2);
                                                push @{$cpg_hash{"CpG"} }, \@pair;
                                          }
                                    }
                                    if ($something{"specific_transcript"} ne "any" && !exists $transcripts_hash{$something{"specific_transcript"} }) {
                                          print_error_html( $temp_dir, "\\\"".$something{"specific_transcript"}."\\\" is not a transcript in the specified gene\n" );
                                          die;
                                    }
                                    foreach my $key ( sort( keys(%cpg_hash) ) ) {
                                          my $splitlocation = Bio::Location::Split->new();
                                          foreach my $pair_ref ( @{ $cpg_hash{$key} } ) {
                                                my $t1=int( @{$pair_ref}[0] - $location_offset);
                                                my $t2=int( @{$pair_ref}[1] - $location_offset);
                                                $splitlocation->add_sub_Location( Bio::Location::Simple->new(
                                                                  -start => int($t1), #leftsplit start is the CRISPR's start
                                                                  -end => int($t2), #left splits end is the CRISPR's start + the left tale length
                                                                  -strand => 1 ) );
                                          }
                                          
                                          
                                          my $feat = Bio::SeqFeature::Generic->new(
                                                      -display_name => $key,
                                                      -primary_tag => "CpGIsland::" . $key,
                                                      -strand => 1,
                                                      -location => $splitlocation
                                          );
                                          $seq_obj->add_SeqFeature($feat);
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
                                    foreach my $key ( sort( keys(%transcripts_hash) ) ) {
                                                my $splitlocation = Bio::Location::Split->new();
                                                foreach my $pair_ref (sort(keys %{$transcripts_hash{$key}})) {
                                                      $splitlocation->add_sub_Location( Bio::Location::Simple->new(
                                                                        -start => int( ${$transcripts_hash{$key}}{$pair_ref}->[0] ) - $location_offset, #leftsplit start is the CRISPR's start
                                                                        -end => int( ${$transcripts_hash{$key}}{$pair_ref}->[1] ) - $location_offset, #left splits end is the CRISPR's start + the left tale length
                                                                        -strand => $strand{$key} ) );
                                                }
                                                my $feat = Bio::SeqFeature::Generic->new(
                                                      -display_name => $key,
                                                      -primary_tag => "Transcript::" . $key,
                                                      -location => $splitlocation
                                                );
                                                $seq_obj->add_SeqFeature($feat);
                                    }
                                    foreach my $key ( sort( keys(%CDS_hash) ) ) {
                                          my $splitlocation = Bio::Location::Split->new();
                                          foreach my $pair_ref ( @{ $CDS_hash{$key} } ) {
                                                $splitlocation->add_sub_Location( Bio::Location::Simple->new(
                                                                  -start => int( $pair_ref->[0] ) - $location_offset, #leftsplit start is the CRISPR's start
                                                                  -end => int( $pair_ref->[1] ) - $location_offset, #left splits end is the CRISPR's start + the left tale length
                                                                  -strand => $strand{$key} ) );
                                          }
                                          my $feat = Bio::SeqFeature::Generic->new(
                                                      -display_name => $key,
                                                      -primary_tag => "CDS::" . $key,
                                                      -location => $splitlocation
                                          );
                                          $seq_obj->add_SeqFeature($feat);
                                    }
                              }
                              if ( exists $something{"do_restriction_analysis"} ) {
                                    ENZYMES: foreach my $line (@enzarray) { #do something for every elemet in the enzarray, thus for every restriction enzyyme motif
                                          MATCHES: while ( $whole_seq =~ m/($line->[1])/ig ) { #for every occurance of the motif in the region do something
                                                my $cutsite = ( length($`) + $line->[2] );
                                                my %tag = ();
                                                $tag{"Restriction_Enzyme"} = $line->[0];
                                                my $feat = Bio::SeqFeature::Generic->new(
                                                            -display_name => "Restriction_Site",
                                                            -primary_tag => "Restriction Site::" . $line->[0],
                                                            -tag => \%tag,
                                                            -strand => 1,
                                                            -start => $cutsite,
                                                            -end => $cutsite
                                                );
                                                $seq_obj->add_SeqFeature($feat);
                                          }
                                    }
                              }
                              foreach my $key ( keys(%{$CRISPR_hash{$fname}}) ) {
                                    my @targets = split(";;", ${ ${ $CRISPR_hash{$fname} } {$key} }{"hits"} );
                                    my $number_off_targets=1;
                                    foreach my $hit (@targets){
                                          if ($hit ne "") {
                                                my @splithit=split("い",$hit);
                                                if (${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}-((20-(100/${ ${ $CRISPR_hash{$fname} } {$key} }{"length"}*$splithit[4]))/$number_off_targets)>0) {
                                                      ${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}=${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}-((20-(100/${ ${ $CRISPR_hash{$fname} } {$key} }{"length"}*$splithit[4]))/$number_off_targets);
                                                      $number_off_targets++;
                                                }else{
                                                      ${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}=0;
                                                }
                                          }
                                    }
                                    if(exists(${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"new_score"})){
                                          @{${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"new_score"}}[0]=${ ${ $CRISPR_hash{$fname} } {$key} }{"score"};
                                          ${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}=${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"new_score"};
                                          @{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[1]=@{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[1]*100/((5*(scalar(keys(%CDS_hash))))+(scalar(keys(%CDS_hash)))+(5*(scalar(keys(%transcripts_hash))))+1);
                                          @{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[2]=100*(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[2]-(-21))/(40-(-21));
                                    }
                                    foreach my $i (0..2){
                                          @{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[$i]=round_digits(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[$i],4);
                                    }
                                    ${ ${ $CRISPR_hash{$fname} } {$key} }{"spec_score"}= @{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[0];
                                    ${ ${ $CRISPR_hash{$fname} } {$key} }{"anno_score"}= @{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[1];
                                    ${ ${ $CRISPR_hash{$fname} } {$key} }{"eff_score"}= @{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[2];
                              }
                              
                              #####################################################################################################################################################################
                              # Create the Table for the findings
                              #####################################################################################################################################################################
                              
                              open (my $outfiletab, ">", $temp_dir . "/" . $fname . "_" . "table.tab");
                                    if ($something{"kind"} eq "single") {
                                          print $outfiletab "Name\tLength\tStart\tEnd\tStrand\tNucleotide sequence\tGene Name\tTranscripts\tTranscript:: Exon\tNumber of Cpg Islands hit\tSequence around the cutside\t%A %C %T %G\tS-Score\tA-Score\tE-Score\tpercent of total transcripts hit\tTarget\tMatch-start\tMatch-end\tMatchstring\tEditdistance\tNumber of Hits\tDirection\n";
                                    }
                                    else {
                                          print $outfiletab "Name\tLength\tStart\tEnd\tStrand\tNucleotide sequence\tGene Name\tTranscripts\tTranscript:: Exon\tNumber of Cpg Islands hit\tSequence around the cutside\t%A %C %T %G\tS-Score\tA-Score\tE-Score\tpercent of total transcripts hit\tTarget\tMatch-start\tMatch-end\tMatchstring\tEditdistance\tNumber of Hits\tDirection\tSpacer\n";
                                    }
                                    PRINTLOOP: foreach my $key ( sort { $CRISPR_hash{$fname}{$b}->{"spec_score"} <=> $CRISPR_hash{$fname}{$a}->{"spec_score"} } sort { $CRISPR_hash{$fname}{$b}->{"anno_score"} <=> $CRISPR_hash{$fname}{$a}->{"anno_score"} } sort { $CRISPR_hash{$fname}{$b}->{"eff_score"} <=> $CRISPR_hash{$fname}{$a}->{"eff_score"} } sort { $CRISPR_hash{$fname}{$b}->{"exon"} cmp $CRISPR_hash{$fname}{$a}->{"exon"} } keys(%{$CRISPR_hash{$fname}}) ) {
                                          $statistics{$fname}{"Number of successful designs"}++;
                                          my @targets=split(";;",${ ${ $CRISPR_hash{$fname} } {$key} }{"hits"} );
                                          #write the tab-delimited file
                                          HITLOOP: foreach my $hit (@targets){
                                                if ($hit ne "") {
                                                      #print the candidates name
                                                      print $outfiletab "$key\t";
                                                      my @splithit = split("い",$hit);
                                                      #print its length on this whole sequence these are not genomic coordinates
                                                      if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"length"} ) {
                                                            print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"length"} . "\t";
                                                      }
                                                      #print its start on this whole sequence these are not genomic coordinates
                                                      if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"} ) {
                                                            print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"} . "\t";
                                                      }
                                                      #print its end on this whole sequence these are not genomic coordinates
                                                      if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"} ) {
                                                            print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"} . "\t";
                                                      }
                                                      if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} ) {
                                                            print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} . "\t";
                                                      }
                                                      if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"} ) {
                                                            if ($something{"kind"} eq "single") {
                                                                  if (${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} eq "minus") {
                                                                        print $outfiletab reverse_comp(substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},0,2)."N ".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},3)) . "\t";
                                                                  }else{
                                                                        print $outfiletab substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},0,length(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"})-3)." N".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},length(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"})-2) . "\t";
                                                                  }
                                                            }else{
                                                                  print $outfiletab reverse_comp(
                                                                                                 substr(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[0],0,2)."N ".substr(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[0],3))
                                                                  ."_" .
                                                                                               substr(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[1],0,length(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[1])-3)." N".substr(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[1],length(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[1])-2). "\t";
                                                            }
                                                      }
                                                      #print the gene name it overlaped with if there was any
                                                      if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"gene"} ) {
                                                            print $outfiletab join( "_", keys( %{ ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"gene"} } ) ) . "\t";
                                                      } else {
                                                            print $outfiletab "NA\t";
                                                      }
                                                      my $percent_of_transcripts_hit = "NA";
                                                      my $number_of_transcripts = scalar(keys(%transcripts_hash));
                                                      my $number_of_target_transcripts = scalar(keys(%{ ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transexons"}}));
                                                      if ($number_of_transcripts != 0) {
                                                            $percent_of_transcripts_hit = $number_of_target_transcripts*100/$number_of_transcripts;
                                                      }
                                                      if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transcripts"} ) {
                                                            print $outfiletab ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transcripts"}."\t";
                                                      } else {
                                                            print $outfiletab "NA\t";
                                                      }
                                                      #print the exons it overlaped with if there was any
                                                      if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transexons"} ) {
                                                            print $outfiletab join( " ", keys(%{ ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transexons"}}))  . "\t";
                                                      } else {
                                                            print $outfiletab "NA\t";
                                                      }
                                                      #print the number of CpG islands it overlapped with if there was any
                                                      if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"CpG"} ) {
                                                            print $outfiletab ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"CpG"} . "\t";
                                                      } else {
                                                            print $outfiletab "NA\t";
                                                      }
                                                      #write the CRISPR as annotation to the original sequence
                                                      my $whole_crisp_seq = substr( $whole_seq, ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"}, ( ${ ${ $CRISPR_hash{$fname} } {$key} }{"length"} + 2 ) );
                                                      #print homology matrix also sub divided into left and right sequence
                                                      if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} }{"left"} ) {
                                                            print $outfiletab "left::" . ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} }{"left"} . "::right::" . ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} }{"right"} . "\t";
                                                      } else {
                                                            print $outfiletab "NA\t";
                                                      }
                                                      ${ ${ $CRISPR_hash{$fname} } {$key} }{"Nuc_comp"} = join( " ", my @basecomp = find_base_count( $whole_crisp_seq ) ); #store the nucleotide composition as a string object in the hash
                                                      
                                                      print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"Nuc_comp"}, "\t", join("\t",@{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}});
                                                      print $outfiletab "\t$percent_of_transcripts_hit\t";
                                                      print $outfiletab $splithit[0]."\t".$splithit[1]."\t".$splithit[2]."\t".$splithit[3]."\t".$splithit[4]."\t";
                                                      print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"number_of_hits"}, "\t";
                                                      print $outfiletab $splithit[5]."\t";
                                                      if (!($something{"kind"} eq "single")) {
                                                            print $outfiletab $splithit[6]."\t";
                                                      }
                                                      
                                                      
                                                      #print the end of the line as a new line
                                                      print $outfiletab "\n";
                                                      #make a featureanntotaion for that CRISPR
                                                }
                                          }
                                          #print the candidates name
                                          if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} }{"context"} }{"gene"} ) {
                                                ${ ${ $CRISPR_hash{$fname} } {$key} }{"gene"} = join( "_", keys( %{ ${ ${ ${ $CRISPR_hash{$fname} } {$key} }{"context"} } {"gene"} } ) );
                                          } else {
                                                ${ ${ $CRISPR_hash{$fname} } {$key} }{"gene"} = "NA";
                                          }
                                          my $percent_of_transcripts_hit="NA";
                                          my $number_of_transcripts = scalar(keys(%transcripts_hash));
                                          my $number_of_target_transcripts = scalar(keys(%{ ${ ${ ${ $CRISPR_hash{$fname} } {$key} }{"context"} }{"transexons"}}));
                                          if ($number_of_transcripts != 0) {
                                                $percent_of_transcripts_hit=$number_of_target_transcripts*100/$number_of_transcripts;
                                          }
                                          if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transcripts"} ) {
                                                ${ ${ $CRISPR_hash{$fname} } {$key} }{"transcripts"} = ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transcripts"} ;
                                          } else {
                                                ${ ${ $CRISPR_hash{$fname} } {$key} }{"transcripts"} = "NA";
                                          }
                                          #print the exons it overlaped with if there was any
                                          if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transexons"} ) {
                                                ${ ${ $CRISPR_hash{$fname} } {$key} }{"transexons"} = join( " ", keys(%{ ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} } {"transexons"} }));
                                          } else {
                                                ${ ${ $CRISPR_hash{$fname} } {$key} }{"transexon"} = "NA";
                                          }
                                          #print the number of CpG islands it overlapped with if there was any
                                          if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"CpG"} ) {
                                                ${ ${ $CRISPR_hash{$fname} } {$key} }{"CpG"} = ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"CpG"};
                                          } else {
                                                ${ ${ $CRISPR_hash{$fname} } {$key} }{"CpG"} = "NA";
                                          }
                                          #write the CRISPR as annotation to the original sequence
                                          my $whole_crisp_seq = substr( $whole_seq, ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"}, ( ${ ${ $CRISPR_hash{$fname} } {$key} }{"length"} + 2 ) );
                                          #print homology matrix also sub divided into left and right sequence
                                          if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} }{"left"} ) {
                                                ${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} = "left::" . ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} }{"left"} . "::right::" . ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} }{"right"};
                                          } else {
                                                ${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} = "NA";
                                          }
                                          undef ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"}; #delete the "context" hash of the crispr hash it is not longer needed
                                          ${ ${ $CRISPR_hash{$fname} } {$key} }{"Nuc_comp"} = join( " ", my @basecomp = find_base_count( $whole_crisp_seq ) ); #store the nucleotide composition as a string object in the hash
                                    } #end of printloop
                              close $outfiletab;
                              
                              #####################################################################################################################################################################
                              #Add some other features to the seq_obj - adjust the image
                              #####################################################################################################################################################################
                              
                              PRINTLOOP: foreach my $key ( sort { $CRISPR_hash{$fname}{$b}->{"spec_score"} <=> $CRISPR_hash{$fname}{$a}->{"spec_score"} } sort { $CRISPR_hash{$fname}{$b}->{"anno_score"} <=> $CRISPR_hash{$fname}{$a}->{"anno_score"} } sort { $CRISPR_hash{$fname}{$b}->{"eff_score"} <=> $CRISPR_hash{$fname}{$a}->{"eff_score"} } sort { $CRISPR_hash{$fname}{$b}->{"exon"} cmp $CRISPR_hash{$fname}{$a}->{"exon"} } keys(%{$CRISPR_hash{$fname}}) ) {
                                    my $feat="";
                                    if ($something{"kind"} eq "single") {
                                          if (${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} eq "plus") {
                                                $feat = Bio::SeqFeature::Generic->new(
                                                      -start => ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"},
                                                      -end => ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"},
                                                      -strand => 1,
                                                      -display_name => $key,
                                                      -primary_tag => "CRISPR::" . join( "_", ${ ${ $CRISPR_hash{$fname} } {$key} }{"gene"} ) . "::" . $key,
                                                      -tag => $CRISPR_hash{$fname}{$key}
                                                );
                                          }else{
                                                $feat = Bio::SeqFeature::Generic->new(
                                                      -start => ${ ${ $CRISPR_hash{$fname}} {$key} }{"start"},
                                                      -end => ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"},
                                                      -strand => -1,
                                                      -display_name => $key,
                                                      -primary_tag => "CRISPR::" . join( "_", ${ ${ $CRISPR_hash{$fname} } {$key} }{"gene"} ) . "::" . $key,
                                                      -tag => $CRISPR_hash{$fname}{$key}
                                                );
                                          }
                                          $seq_obj->add_SeqFeature($feat);
                                    }else{
                                          my $splitlocation = Bio::Location::Split->new(); #initialize a new splited location object
                                          $splitlocation->add_sub_Location(
                                                            Bio::Location::Simple->new(
                                                                  -start => int( ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"} ), #leftsplit start is the CRISPR's start
                                                                  -end => int( ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"} ) + @{${ ${ $CRISPR_hash{$fname} } {$key} }{"lengthcombo"}}[0], #left splits end is the CRISPR's start + the left tale length
                                                                  -strand => 1 )
                                                            );
                                          $splitlocation->add_sub_Location(
                                                            Bio::Location::Simple->new(
                                                               -start => int( ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"} ) - @{${ ${ $CRISPR_hash{$fname} } {$key} }{"lengthcombo"}}[0], #right splits start is the CRISPR's start + the left tale length +spacer length
                                                               -end => int( ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"}), #right split end is the CRISPR's end
                                                               -strand => 1 )
                                                            );
                                          my $feat = Bio::SeqFeature::Generic->new(
                                                      -display_name => $key,
                                                      -primary_tag => "CRISPR::" . join( "_", ${ ${ $CRISPR_hash{$fname} } {$key} }{"gene"} ) . "::" . $key,
                                                      -tag => $CRISPR_hash{$fname}{$key},
                                                      -location => $splitlocation
                                          );
                                          $seq_obj->add_SeqFeature($feat);
                                    }
                              }
                               
                              #####################################################################################################################################################################
                              #If asked, write all the stuff to a gff format file
                              #####################################################################################################################################################################
                              
                              if ( exists $something{"out_gff"}) { 
                                    open(my $gfffile, ">",$temp_dir . "/" . $fname . ".gff" ) or die $!;
                                    print $gfffile "##gff-version 3\n";
                                    PRINTLOOP: foreach my $key ( sort { $CRISPR_hash{$fname}{$b}->{"spec_score"} <=> $CRISPR_hash{$fname}{$a}->{"spec_score"} } sort { $CRISPR_hash{$fname}{$b}->{"anno_score"} <=> $CRISPR_hash{$fname}{$a}->{"anno_score"} } sort { $CRISPR_hash{$fname}{$b}->{"eff_score"} <=> $CRISPR_hash{$fname}{$a}->{"eff_score"} } sort { $CRISPR_hash{$fname}{$b}->{"exon"} cmp $CRISPR_hash{$fname}{$a}->{"exon"} } keys(%{$CRISPR_hash{$fname}}) ) {
                                          my @locus=split("::",$statistics{$fname}{"seq_location"});
                                          print $gfffile $locus[0]."\tE-CRISP\tCRISPRtarget\t".(${ ${ $CRISPR_hash{$fname} } {$key} }{"start"}+$locus[1])."\t";
                                          print $gfffile (${ ${ $CRISPR_hash{$fname} } {$key} }{"end"}+$locus[1])."\t".sum(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}})."\t";
                                          if (${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} eq "minus") {
                                                print $gfffile "-"."\t."."\t";
                                          }else{
                                                print $gfffile "+"."\t."."\t";
                                          }
                                          print $gfffile "id=".$key."; ";
                                          print $gfffile "spec_score=".${ ${ $CRISPR_hash{$fname} } {$key} }{"spec_score"}."; ";
                                          print $gfffile "anno_score=".${ ${ $CRISPR_hash{$fname} } {$key} }{"anno_score"}."; ";
                                          print $gfffile "eff_score=".${ ${ $CRISPR_hash{$fname} } {$key} }{"eff_score"}."; ";
                                          if ($something{"kind"} eq "single") {
                                                if (${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} eq "minus") {
                                                      print $gfffile "seq=".reverse_comp(substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},0,2)."N_".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},3)) . "; ";
                                                }else{
                                                      print $gfffile "seq=".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},0,length(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"})-3)."_N".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},length(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"})-2) . "; ";
                                                }
                                          }else{
                                                print $gfffile "seq=".reverse_comp(
                                                      substr(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[0],0,2)."N_".substr(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[0],3))
                                                      ."_" .
                                                      substr(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[1],0,length(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[1])-3)."_N".substr(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[1],length(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[1])-2). "; ";
                                          }
                                          print $gfffile "offtargetcount=".${ ${ $CRISPR_hash{$fname} } {$key} }{"number_of_hits"}.";\n";
                                    }
                                    close $gfffile;
                              }
                              open (my $statsfile, ">", "/var/www/E-CRISP/workdir/$temp_dir/statsfile.txt") or die $!;
                                    print $statsfile "90\tresults were written\n";
                              close $statsfile;
                              #####################################################################################################################################################################
                              #Not implemented anymore at current state (might come back in the future)
                              #####################################################################################################################################################################
                              
                              if ( exists $something{"out_seq"} ) { #should all the stuff be written to a genbank format file
                                    #my $seq_out = Bio::SeqIO->new( -file => ">" . $temp_dir . "/" . $fname . "_anno.seq", -format => "genbank" ); #open a new genbank format file
                                    #$seq_out->write_seq($seq_obj); #write the annotated sequence into that seqformat file
                              }
                              
                              #####################################################################################################################################################################
                              # Draw the Image of the sequence
                              #####################################################################################################################################################################
                              
                              my $seq = $seq_obj;
                              if ( exists $something{"draw_image"} &&  $statistics{$fname}{"Number of successful designs"}<200 ){
                                    open (my $statsfile, ">", "/var/www/E-CRISP/workdir/$temp_dir/statsfile.txt") or die $!;
                                          print $statsfile "90\tresult image is beeing drawn\n";
                                    close $statsfile;
                                    my $wholeseq = Bio::SeqFeature::Generic->new( -start => int(1),
                                                -end => int( $seq->length ),
                                                -display_name => ( $seq->display_name )
                                    );
                                    my $panel = Bio::Graphics::Panel->new(
                                                -length => int( $seq->length ),
                                                -key_style => 'bottom',
                                                -width     => 750,
                                                -pad_left  => 10,
                                                -pad_right => 10,
                                                -grid      => 1,
                                                -spacing   => 10,
                                                -key_spacing   => 10,
                                                -image_class=> 'GD::SVG'
                                    );
                                    $panel->add_track( $wholeseq,
                                                -glyph  => 'arrow',
                                                -bump   => 0,
                                                -double => 1,
                                                -tick   => 2 );
                                    $panel->add_track( $wholeseq,
                                                -glyph   => 'generic',
                                                -bgcolor => 'blue',
                                                -label   => 1,
                                                -description => $seq_obj->description()
                                    );
                                    # sort features by their primary tags
                                    my %sorted_features;
                                    for my $f ($seq->all_SeqFeatures) {
                                          my $tag = $f->primary_tag;
                                          $tag =~ m/^(\w+)::.+/;
                                          push @{ $sorted_features{$1} }, $f;
                                    }
                                    my $exon    = 0;
                                    my $CpG     = 0;
                                    my $CDS     = 0;
                                    my $res     = 0;
                                    my $trans   = 0;
                                    my $gene    = 0;
                                    my $design  = 0;
                                    for my $f ($seq->all_SeqFeatures) {
                                          my $tag = $f->primary_tag;
                                          $tag =~ m/^(\w+)::.+/;
                                          $tag = $1;
                                          if ( $tag =~ m/Transcript/ && $trans == 0 ) {
                                                $panel->add_track( $sorted_features{$tag},
                                                            -glyph => 'transcript2',
                                                            -bgcolor => 'orange',
                                                            -fgcolor => 'black',
                                                            -font2color => 'red',
                                                            -key => 'Transcript', 
                                                            -truetype=> '1',
                                                            -bump => +1,
                                                            -height => 8,
                                                            -label => \&gene_label
                                                );
                                                $trans=1;
                                                delete $sorted_features{$tag};
                                          }
                                          if ( $tag =~ m/CpG/ && $CpG == 0 ) {
                                                $panel->add_track( $sorted_features{$tag},
                                                            -glyph => 'segments',
                                                            -bgcolor => 'lightyellow',
                                                            -fgcolor => 'black',
                                                            -font2color => 'red',
                                                            -truetype=> '1',
                                                            -key => 'CpG Island',
                                                            -bump => +1,
                                                            -height => 8,
                                                            -label => \&gene_label
                                                );
                                                delete $sorted_features{$tag};
                                                $CpG = 1;
                                          }
                                          if ( $tag =~ m/CDS/ && $CDS==0) {
                                                $panel->add_track( $sorted_features{$tag},
                                                            -glyph => 'transcript',
                                                            -bgcolor => 'lightgreen',
                                                            -fgcolor => 'black',
                                                            -font2color => 'red',
                                                            -truetype=> '1',
                                                            -key => 'CDS',
                                                            -bump => '+1',
                                                            -height => 8,
                                                            -label => \&gene_label
                                                );
                                                $CDS=1;
                                                delete $sorted_features{$tag};
                                          }
                                          if ( $tag =~ m/gene/ && $gene==0) {
                                                $panel->add_track( $sorted_features{$tag},
                                                            -glyph => 'generic',
                                                            -bgcolor => 'lightblue',
                                                            -fgcolor => 'black',
                                                            -font2color => 'red',
                                                            -truetype=> '1',
                                                            -key => 'gene',
                                                            -bump => +1,
                                                            -height => 8,
                                                            -label => \&gene_label
                                                );
                                                $gene=1;
                                                delete $sorted_features{$tag};
                                          }
                                          if ( $tag =~ m/CRISPR/ && $design==0 ) {
                                                $panel->add_track( $sorted_features{$tag},
                                                            -glyph => 'generic',
                                                            -bgcolor => 'chartreuse',
                                                            -fgcolor => 'black',
                                                            -font2color => 'red',
                                                            -truetype=> '1',
                                                            -key => 'CRISPR',
                                                            -bump => +1,
                                                            -height => 12,
                                                            -label => \&gene_label
                                                );
                                                $design=1;
                                                delete $sorted_features{$tag};
                                          }
                                          if ( $tag =~ m/Restriction Site/ig && !( $res eq $tag ) ) {
                                                $panel->add_track( $sorted_features{$tag},
                                                            -glyph => 'generic',
                                                            -bgcolor => 'green',
                                                            -fgcolor => 'black',
                                                            -font2color => 'red',
                                                            -truetype=> '1',
                                                            -key => 'Restriction Site',
                                                            -bump => +1,
                                                            -height => 8,
                                                            -label => \&gene_label
                                                );
                                                delete $sorted_features{$tag};
                                                $res = $tag;
                                          }
                                    }
                                    open (my $tempfile, ">", $temp_dir . "/temp.svg");
                                          print $tempfile $panel->svg;
                                    close $tempfile;
                                    open ($tempfile, "<", $temp_dir . "/temp.svg");
                                          open (my $myimage, ">", $temp_dir . "/" . $fname . ".svg");
                                                while (my $line = <$tempfile>){
                                                      if($line =~ m/^\<svg .*$/){
                                                            $line = $line.'<script xlink:href="/E-CRISP/SVGPan.js"/>'."\n".'<g id="viewport" transform="translate(0,0)">'."\n";
                                                      }elsif($line =~ m/\<\/svg\>/){
                                                            $line = "\n".'</g>'."\n".$line;
                                                      }
                                                      if($line =~ s/(\>($fname\_\d+_\d+)\<\/text)/ onmouseover="this.style.cursor='pointer'" onclick="window.parent.svgElementClicked(&quot;$2&quot;)" $1/){
                                                            my $temp = $2;
                                                            $line =~ s/id="\S+"/id="$temp"/;
                                                      }
                                                      print $myimage $line;
                                                }
                                          close $myimage;
                                    close $tempfile;
                                    unlink $temp_dir . "/temp.svg";
                              }else{
                                    open( my $report, ">>", "$temp_dir/index.html" ) or die "can not open report file";
                                          print $report '<tr><td>The image could not be drawn, because to many designs were found.</td></tr>';
                                    close $report;
                              }
                              
                              #####################################################################################################################################################################
                              # ZIP the stuff
                              #####################################################################################################################################################################
                              
                              unlink $temp_dir . "/tempfile.fasta";
                              my $zip    = Archive::Zip->new();
                              my $member = "";
                              if ( -e $temp_dir . "/" . $fname . ".gff") { $member = $zip->addFile( $temp_dir . "/" . $fname . ".gff", $fname . "_CRISPR.gff" ); }
                              if ( -e $temp_dir . "/" . $fname . "_anno.seq") { $member = $zip->addFile( $temp_dir . "/" . $fname . "_anno.seq", $fname . "_CRISPR.gb" ); }
                              $member = $zip->addFile( $temp_dir . "/" . $fname . "_" . "table.tab", $fname . "_CRISPR.tab" );
                              if( -e $temp_dir . "/" . $fname . ".svg"){$member = $zip->addFile( $temp_dir . "/" . $fname . ".svg", $fname . "_CRISPR.svg" );} 
                              $zip->writeToFileNamed( $temp_dir . '/' . $fname . '.zip' );
                              
                              #####################################################################################################################################################################
                              # Print the report site with the results (header was created earlier in the loop)
                              #####################################################################################################################################################################
                              if ( exists $something{"draw_image"} && scalar(keys %{ $CRISPR_hash{$fname} })<200) {
                                    open( my $report, ">>", "$temp_dir/index.html" ) or die "can not open report file";
                                    #print the statistics
                                    print $report '   <tr><td>
                                                            <span class="main3"><B>Query name: '.$statistics{$fname}{"seq_name"}.'   Query length: '.$statistics{$fname}{"seq_length"}.'   Query location: '.$statistics{$fname}{"seq_location"}.'</B><br><br>';
                                    foreach my $statistic(sort {$b cmp $a} keys(%{$statistics{$fname}})){
                                          if ($statistic =~ m/Number/ig) {
                                                print $report $statistic.' = '.$statistics{$fname}{$statistic}.'<br>';
                                          }
                                    }
                                    print $report '         <br>S: Specificity score A: Annotation score E: Efficiency score<br>for more information please see the <a href="../../aboutpage.html" style="color: blue;">Help</a> pages
                                    </span><br><br>
                                                      </td></tr>';
                                    #print the table
                                    print $report '   <tr><td>
                                                            <table class="talenHits" style="text-align: center; width: 800px; table-layout: fixed; border-bottom: 2px solid black;border-top: 1px solid black; padding: 2px; background-color: white;">';
                                    open (my $tabfile, "<", $temp_dir . "/" . $fname . "_" . "table.tab");
                                         my $count = 0;
                                                while (my $line = <$tabfile>){
                                                      chomp $line;
                                                      my @line=split("\t",$line);
                                                      print $report '<tr class="" id="'.$line[0].'">';
                                                      foreach my $element (0,5,12,16,19,21){
                                                            if ($count == 0){
                                                                  if ($element==12) {
                                                                       print $report '<th style="text-align: center; border-bottom: 1px solid black; color: white; padding: 2px; background-color: #2662C3;"> SAE-Score </th>';
                                                                  }else{
                                                                        print $report '<th style="text-align: center; border-bottom: 1px solid black; color: white; padding: 2px; background-color: #2662C3;"> '.$line[$element].' </th>';
                                                                  }
                                                            }else{
                                                                  my $hittype = "bad";
                                                                  if ($line[21] == 1) {
                                                                        $hittype = "good";
                                                                  }
                                                                  if ($element == 19) {
                                                                        #create popup for the matchstring-info if the option is chosen, otherwise only the matchstring itself
                                                                        if (exists $something{"match_info"} ) {
                                                                              my $popup = "";
                                                                              if ($line[16] ne "*") {                                                                              
                                                                                    if ($something{"kind"} eq "single") {
                                                                                          $popup = create_popup(\@line, $databasepath, $something{"unspecific_leading_bases"}, 10);
                                                                                    }else {
                                                                                          $popup = create_popupD(\@line, $databasepath, $something{"unspecific_leading_bases"}, 5);
                                                                                    }
                                                                                    $line[0] =~ s/\W/_/ig; #for html
                                                                                    $line[16] =~ s/\W/_/ig; #for html
                                                                                    print $report '<td class="'.$hittype.' main-button" style="text-align: center; word-wrap: break-word; font-family: arial, monospace; color: black; font-size:12px; background-color: white;">'
                                                                                                      .'<div id="dialog_'.$line[0].'__'.$line[16].'__'.$line[19].'" class="dialog" title="Matchstring Info for '.$line[0].' on Target '.$line[16].'">'.$popup.'</div>';
                                                                                    print $report     '<button id="'.$line[0].'__'.$line[16].'__'.$line[19].'" class="opener">Matchstring Info</button>'
                                                                                                .'</td>';
                                                                              }else{
                                                                                    print $report '<td class="'.$hittype.' main-button" style="text-align: center; word-wrap: break-word; font-family: arial, monospace; color: black; font-size:12px; background-color: white;">NO-TARGET</td>';
                                                                              }                  
                                                                        }else {
                                                                              print $report '   <td class="'.$hittype.' main" style="text-align: center; word-wrap: break-word; font-family: arial, monospace; color: black; font-size:12px; background-color: white;">'
                                                                                                      . print_offtarget_string($line[$element]) .
                                                                                                '</td>';
                                                                        }
                                                                  }elsif( $element == 12 ){
                                                                        print $report '<td class="'.$hittype.'" >';
                                                                        print $report hor_bar_chart($line[12],$line[13],$line[14]);
                                                                        print $report '</td>';
                                                                  }else{
                                                                        print $report '   <td class="'.$hittype.' main" style="text-align: center; word-wrap: break-word; font-family: arial, monospace; color: black; font-size:12px; background-color: white;"> '
                                                                                                .$line[$element].
                                                                                          '</td>';
                                                                  }
                                                            }
                                                      }
                                                      $count++;
                                                      print $report '</tr>';
                                                }                                          
                                    close $tabfile;
                                    print $report '         </table>
                                                      </td></tr>
                                                      <tr><td>
                                                            <br><br>
                                                      </td></tr>';
                                    #print the Image, the ZIP-Download and the bottom-line
                                    print $report '   <table colspan="2" rowspan="2" style="vertical-align: top;border: 1px dashed black; width: 800px;">
                                                            <tr><td><br><br>';
                                    if ( exists $something{"draw_image"} &&  $statistics{$fname}{"Number of successful designs"}<200) {
                                          print $report '         <object data="' . $fname . '.svg" type="image/svg+xml" /><br>';
                                    }else{
                                          print $report 'The image could not be provided. Please see the downloadable results file for your designs.';
                                    }
                                    print $report '               <br><br>
                                                            </td></tr>
                                                            <tr><td>
                                                                  <a href="' . $fname . '.zip' . '" target="_blank" download="'.$fname.'_CRISPRS.zip" class="main"><input type="button" value="Download  results as zipped folder"></a><br><br><br>
                                                            </td></tr>
                                                      </table>';
                                    print $report '   <tr><td>
                                                            <br><br>
                                                      </td></tr>
                                                      <tr><td colspan="2" rowspan="1" style="vertical-align: top;">
                                                            <embed style="width: 100%; height: 10px;" alt="there should appear a  line" src="/E-CRISP/gradientline.svg" type="image/svg+xml"><br>
                                                      </td></tr>
                                                      <tr><td>
                                                            <br><br>
                                                      </td></tr>';
                              }else{
				  print $report '<tr><td>The table could not be drawn, because to many designs were found. Please appreciate the downloadable table to find your results.</td></tr>';
			      }
                              close $report;
                              
                              $runtime=time-$runtime;
                              open(my $log,">>", "/var/log/talecrisp.log");
                                    print $log $temp_dir."\tE-CRISP\tDESIGN\t$fname\t$runtime\n";
                              close($log);
                              
                              $CRISPR_cnt{$fname}++;
                        }
                  } #end Sequence loop
                  open ($statsfile, ">", "/var/www/E-CRISP/workdir/$temp_dir/statsfile.txt") or die $!;
                  print $statsfile "100\timage was drawn. You will be redirected to your results soon.\n";
		  close $statsfile;
                  open( $report, ">>", "$temp_dir/index.html" ) or die "can not open report file";
                        open( my $footer, "<", "/var/www/E-CRISP/workdir/footer.txt");
                              while(<$footer>){
                                    print $report $_;
                              }
                        close $footer;
                  close $report;
                  print end_html;
                  
                  #################################################################################################################################################################################
                  # End process with creation of "fertig.txt" telling the waiting-process to terminate
                  # Also write an entry into the buffer-index and copy the findings into the buffer for further requests
                  #################################################################################################################################################################################
                  my $head=1;
                  opendir(TEMPDIR,$temp_dir);
                        open (my $outfile, ">>", $temp_dir."/all_results_together.tab");
                        open (my $outfile2, ">>", $temp_dir."/all_results_together.gff");
                              foreach my $filename (readdir(TEMPDIR)){
                                    if($filename=~m/table.tab/){
                                          open (my $file, "<", $temp_dir."/".$filename);
                                                while (<$file>){
                                                      if ($_=~m/Name/ && $head==1) {
                                                            print $outfile $_;
                                                            $head=0;
                                                      }elsif($_=~m/Name/ && $head!=1){
                                                            next;
                                                      }else{
                                                             print $outfile $_;
                                                      }      
                                                }
                                          close $file;
                                    }elsif($filename=~m/.*\.gff/){
                                          open (my $file, "<", $temp_dir."/".$filename);
                                                while (<$file>){
                                                      print $outfile2 $_;
                                                }
                                          close $file;
                                    }
                              }
                        close $outfile;
                        close $outfile2;                        
                  closedir TEMPDIR;
                  system( 'sort -r ' . $temp_dir . '/all_results_together.gff | uniq > ' . $temp_dir . '/temp; mv ' . $temp_dir . '/temp ' . $temp_dir . '/all_results_together.gff;');
                  system( 'touch ' . $temp_dir . '/fertig.txt' );                  
                  open(my $indexfile, "<", "/var/www/E-CRISP/workdir/buffer/buffer.idx") or die $!;
                        $tag_count         = 0;
                        $line_count     = 0;
                        $line_found     = "";
                        my $min_access     = 2000;
                        my $delete_line    = "";
                        while (my $line = <$indexfile>){
                              chomp $line;
                              $line_count++;
                              my @line=split("\t",$line);
                              if($line[0] eq $id_tag) {
                                    $tag_count++;
                                    last; # end loop if found
                              }
                              if ($min_access > $line[2]) { # store the line with the lowest access
                                    $min_access = $line[2];
                                    $line_found = $line_count."d";
                                    $delete_line = $line[1];
                              }
                        }
                  close($indexfile);
                  if ($line_count>500) {
                        system("sed -i $line_found /var/www/E-CRISP/workdir/buffer/buffer.idx;");
                        system("rm -rf /var/www/E-CRISP/workdir/buffer/$delete_line");
                  }
                  open($indexfile, ">>", "/var/www/E-CRISP/workdir/buffer/buffer.idx") or die $!;
                        if( $tag_count==0) { # should be allways executed, if the program works
                              dircopy($temp_dir,"/var/www/E-CRISP/workdir/buffer/$temp_dir");
                              print $indexfile "$id_tag\t$temp_dir\t1\n";
                        }
                  close($indexfile) or die $!;
                  
                  last; # end forked process
            } else {
                  $req->Attach();
                  $SIG{CHLD} = 'IGNORE'; #do not wait for child processes, so that they will not end as zombie BUT requires reactivation in child processes (in the unless branch of this else)
                  print "<tr>
                              <td align=\"left\"><font size=\"3\">WAITING for results</font></td>
                        </tr>
                        ";
                  open(my $footer, "<", "footer.txt");
                        while(<$footer>){
                              print $_;
                        }
                  close $footer;
                  print "<meta http-equiv=\"refresh\" content=\"" . $waiting_time . "; URL=make_crisprs_ng_waitpage.fcgi?PROCESS=1&ID=$temp_dir&WAITINGTIME=$waiting_time&PID=$process_id\">";
                  print end_html;
            }
      }
}

#########################################################################################
#name:      create_popup
#function:  creates popup with. the target-, match- and querry-string adjusted to each
#input:     (line of the outfile.tab as string, databasepath as string,
#           $something{"unspecific_leading_bases"}, int overflow before and after Sequence)
#output:    HTML formated string
#########################################################################################
sub create_popup {
      my @line          = @{$_[0]};
      my $report        = "";
      my $target        = "";
      my $match         = $line[19];
      my $querry        = "";
      my $adjust        = 0;
      my $adjustleft    = 0;
      my $start         = 0;
      my $end           = 0;
      my $insert        = 0;
      
      if ($line[22] eq "rc") { #turn and translate the match and querry string if direction is reverse complementary
            $querry = reverse_comp($line[5]);
            $querry =~s/\s+//ig;
            $querry =~s/^(\w{2})\w(\w+)/$1N$2/ig;
      }
      else {
            $querry = $line[5];
            $querry =~s/\s+//ig;
             $querry =~s/(\w+)\w(\w{2})$/$1N$2/ig;
            $adjust = $_[2];
      }
      
      #search the target-string in database - substr workaround used, because the Bioperl function does not work proper with large numbers
      my $db            = Bio::DB::Fasta->new( $_[1] . '.all.dna.fa', -makeid => \&make_my_id );
      my @target_name   = split("::", $line[16]);
      my $obj           = $db->get_Seq_by_id($target_name[0]);
      if (!defined $obj) { #have to search in whole chromosom
            $db         = Bio::DB::Fasta->new( $_[1] . '.dna.toplevel.fa', -makeid => \&make_my_id );
            $obj        = $db->get_Seq_by_id($line[16]);
      }
      $start   = $line[17]-$adjust-$_[3]-1;
      if ($start < 0 ) {
            $adjustleft = abs($start);
            $start = 0;
      }
      
      $end     = ($line[18]-$adjust+$_[3]);
      $target  = substr $obj->seq, $start, $end - $start - 1;
      #count insertions and adjust the end property to display
      $insert  = $line[17] =~ tr/I/x/;
      
      $report .= '<table>';
      $report .= '<tr><td>Target:</td><td>|'.$start.'*|</td>'.create_popup_string($target, $_[3], 0, 0, $match, 0, 0, $adjustleft).'<td>|'.($end-$insert).'*|</td></tr>';
      $report .= '<tr><td>Matchstring:</td><td></td>'.create_popup_string($match, $_[3], 1, 1, $match, 0, 0, $adjustleft).'</tr>';
      $report .= '<tr><td>Query:</td><td></td>'.create_popup_string($querry, $_[3], 1, 0, $match, 0, 0, $adjustleft).'</tr>';
      $report .= '</table>';
      $report .= '<hr /><p style="font-size:50%">*Start and End-point in the target</p>';
      
      return $report;
}

#########################################################################################
#name:      create_popupD
#function:  for double sequence - does the same as create_popup
#input:     (line of the outfile.tab as string, databasepath as string,
#           $something{"unspecific_leading_bases"}, int overflow before and after Sequence)
#output:    HTML formated string
#########################################################################################
sub create_popupD {
      my @line          = @{$_[0]};
      my $report        = "";
      my $target        = "";
      my @match         = split("-", $line[19]);
      my $tmp           = $line[5];
      $tmp              =~s/\s+//ig;
      my @querry        = split("_", $tmp );
      my $leadingbases  = $_[2];
      my $adjustleft    = 0;
      my $start         = 0;
      my $end           = 0;
      my $insert0       = 0;
      my $insert1       = 0;
      my $deletion0     = 0;
      my $deletion1     = 0;
      my $spacer        = "";
      my $match0        = "";
      my $match1        = "";
      
      #adjust the first/second querry-string 
      if ($line[22] eq "rc") {
            $querry[0]   = reverse_comp($querry[0]);
            $match0 = $match[0];
            $match1 = $match[1];
            #count insertions and deletions to adjust the end/start property
            $deletion0   = $match0 =~ tr/D/o/;
            $deletion1   = $match1 =~ tr/D/o/;
            #calc start and end of the target sequence
            $start       = $line[17] - $_[3] - 1 - $deletion1; # start of the match - the wanted overflow on the left side
            $end         = $line[18] + $line[23] + length($match[0]) + $_[3] - ($leadingbases*2) - ($deletion0*2); # end of the match + spacer size + length of the matchstring + the wanted overflow on the right side - the count of unspec leadingbases*2 (1st and 2nd)
      }
      else {
            $querry[1]   = reverse_comp($querry[1]);
            #swap both query and matchstring and change the start and end points ($line[15] is no longer the start of the target cause the first query swaped)
            @querry[0,1] = @querry[1,0];
            @match[0,1]  = @match[1,0];
            $match0 = $match[0];
            $match1 = $match[1];
            #count deletions to adjust the end property
            $deletion0   = $match0 =~ tr/D/o/;
            $deletion1   = $match1 =~ tr/D/o/;
            #calc start and end of the target sequence
            $start       = $line[17] - $line[23] - length($match[0]) - $_[3] + $leadingbases - 1 - $deletion1; # swaped of the if-clause and reverted (-/+), but only the leadingbases of the 2nd string count here
            $end         = $line[18] + $_[3] - $leadingbases - ($deletion0*2); # swaped of the if-clause and reverted (-/+), but only the leadingbases of the 1st string count here
      }
      
      if ($start < 0 ) {
            $adjustleft = abs($start);
            $start = 0;
      }
      
      #count insertions and deletions to adjust the end/start property
      $insert0   = $match0 =~ tr/I/o/;
      $insert1   = $match1 =~ tr/I/o/;
      for (my $i = 0; $i < ($line[23] - $leadingbases*2 - $deletion0 + $deletion1); $i++){
            $spacer .= " ";
      }
      
      #search the target-string in database - substr workaround used, because the Bioperl function does not work proper with large numbers
      my $db             = Bio::DB::Fasta->new( $_[1] . '.all.dna.fa', -makeid => \&make_my_id );
      my @target_name    = split("::", $line[16]);
      my $obj            = $db->get_Seq_by_id($target_name[0]);
      if (!defined $obj) { #have to search in whole chromosom
            $db          = Bio::DB::Fasta->new( $_[1] . '.dna.toplevel.fa', -makeid => \&make_my_id );
            $obj         = $db->get_Seq_by_id($line[16]);
      }
      
      $target  = substr $obj->seq, $start, $end - $start - 1;
      
      #build the table for the popup
      $report .= '<table style="font-size:75%">';
      $report .= '<tr><td>Query:</td><td></td>'.create_popup_string($querry[0], $_[3], 1, 0, $match[0], 0, 0, $adjustleft).'</tr>';
      $report .= '<tr><td>Matchstring:</td><td></td>'.create_popup_string($match[0], $_[3], 1, 1, $match[0], 0, 0, $adjustleft).'</tr>';
      $report .= '<tr><td>Target:</td><td>|'.$start.'*|</td>'.create_popup_string($target, $_[3], 0, 0, $match[0].$spacer.$match[1], 0, 1, $adjustleft).'<td>|'.($end - 1 - $insert1).'*|</td></tr>';
      $report .= '<tr><td>Matchstring:</td><td></td>'.create_popup_string($match[1], $_[3], 1, 1, $match[1], (length($match[0]) + int($line[23]) + $insert0 - ($leadingbases*2) - $deletion0 + $deletion1), 0, $adjustleft).'</tr>';
      $report .= '<tr><td>Query:</td><td></td>'.create_popup_string($querry[1], $_[3], 1, 0, $match[1], (length($match[0]) + int($line[23]) + $insert0 - ($leadingbases*2) - $deletion0 + $deletion1), 0, $adjustleft).'</tr>';
      $report .= '</table>';
      $report .= '<hr /><p style="font-size:50%">*start- and endpoint in the target sequence</p>';
      
      return $report;
}

#########################################################################################
#name:      create_popup_string
#function:  adjusts a sequence for the popup-output
#input:     (string sequence, int overflow from create_popup, boolean if it should be used,
#           boolean if colored or not, string matchstring, 2nd overflow for double seq,
#           boolean if the last letter should  only be popped for 2nd half - 1 = yes,
#           int additional left re-adjustment)
#output:    HTML formated string
#########################################################################################
sub create_popup_string {
      my $string  = "";
      my @letters = split("",$_[0]);
      my @match   = split("",$_[4]);
      my $indel   = 0;
      my $letter  = "";
      my $length  = length($_[4]);
      
      if ($_[2] == 1) {
            for (my $i = 0; $i < $_[1] + $_[5] - $_[7]; $i++){
                  $string .= '<td></td>';
            }
      }
       
      
      for(my $i = 0; $i < scalar(@letters); $i++){
            if ($_[3] == 1) { # matchstring needs color
                  $string .= '<td align="center" valign="middle">'.print_offtarget_string($letters[$i]).'</td>';
            }
            else{
                  #adjust the letters for indel
                  $letter = $letters[$i];
                  if ($_[2] == 1 && $indel < $length) { # for querystring (same size as matchstring)
                        if ($match[$indel] eq "D") {
                              $letter = "_";
                              $i--;
                        }
                  }
                  elsif ($indel - $_[1] >= 0 && ($indel - $_[1]) < $length) { # for targetstring (needs adjustment)
                        if ($match[$indel - $_[1]] eq "I") {
                              $letter = "_";
                              $i--;
                              if ($_[6] == 0) {
                                    pop (@letters); #delete last Element
                              }
                              
                              elsif ($indel - $_[1] > $length/2) {
                                    pop (@letters); #delete last Element only for I's in the 2nd matchstring
                              }
                        }
                  }
                  
                  $string .= '<td align="center" valign="middle">'.$letter.'</td>';
                  $indel++;
            }
      }
      
      return $string;
}

#########################################################################################
#name:      make_temp_fasta_file
#function:  creates a temporary fasta file for the bowtie index and builds a trees
#input:     (given id-Array, tree as referemce, something-Hashreference,
#           enzyme db, temp_dir, 1/0 if file or not)
#output:    N/A
#########################################################################################
sub make_temp_fasta_file {
      if ($_[5] == 0 && scalar(@{$_[0]})>=50) {
            print_error_html( $_[4], "Your input is more than 50 Sequences\n" );
             die;
      }
      if ($_[5] == 1 && scalar(@{$_[0]})>=50) { #@Flo die 500 sind hier Absicht?
            print_error_html( $_[4], "Your input is more than 50 lines with IDs.<br> Please shorten the list, or maybe change to option to FASTA.<br>" );
            die;
      }
      if ( !( $_[2]->{"specific_transcript"} eq "any") && scalar(@{$_[0]}) >1) {
            print_error_html( $_[4], "Transcript specificity is only defined for single gene analyses.\n" );
            die;
      }
      open (my $tempfile, ">", $_[4] . "/tempfile.fasta");
            foreach my $id (@{$_[0]}) { 
                  $id =~ s/\s//ig;
                  my $seq_obj = $_[3]->get_Seq_by_id($id); # get a PrimarySeq obj
                  if ($seq_obj) {
                        my $header = $_[3]->header($id); # get the header, or description line
                        $header =~ m/chrom:([\w\.]+):(\d+)..(\d+)/ig;
                        my $chrom = $1;
                        my $location_offset = $2;
                        my $location_end = $3;
                        if ( !exists $_[1]->{$chrom} ) {
                              $_[1]->{$chrom} = build_tree( "/var/www/E-CRISP/DATABASEFILESNEW/" . $_[2]->{"ref_organism"} . "/" . $chrom . "_indexed" );
                        }
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
#name:      find_and_print_CRISPRS
#function:  find the CRISPRS and build a hash for CRISPRS and the statistics...
#input:     (chrom, location_offset, builded tree as reference, dont_asses_context,
#           global temp_dir, parallel_number, something-Hashreference)
#output:    hash for CRISPRS and hash for statistics
#########################################################################################
sub find_and_print_CRISPRS {
      my $seq_obj_ref               = $_[0];
      my $chrom                     = $_[1];
      my $location_offset           = $_[2];
      my %trees                     = %{ $_[3] };
      my $dont_asses_context        = $_[4];
      my $temp_dir                  = $_[5];
      my $parallel_number           = $_[6];
      my %something                 = %{ $_[7] };
      my $seq_obj                   = $$seq_obj_ref;
      my $whole_seq                 = $seq_obj->seq;
      my $count                     = 0;
      my %finished_CRISPR_hash      = ();
      my $pm                        = Parallel::ForkManager->new($parallel_number);
      my $cutnumber                 = int( length($whole_seq) / int( $parallel_number - 1 ) );
      my $cut                       = 0;
      my @cuts                      = ();
      my %tempstatistics            = ();
      
      while ( $cut <= length($whole_seq) ) {
            push @cuts, $cut;
            $cut = $cut + $cutnumber;
      }
      
      #################################################################################################################################################################################
      # cut the sequence into equal peaces, so that each forked child can work on one part (paralell!)
      #################################################################################################################################################################################
      
      foreach $cut (@cuts) {
            $pm->start and next;
            my $seq = substr( $whole_seq, $cut, $cutnumber );
            my %CRISPR_hash = ();
            my %Gpos = make_pos_index( \$seq, "G" );
            my %Cpos = make_pos_index( \$seq, "C" );
            my %Apos = make_pos_index( \$seq, "A" );
            my %Tpos = make_pos_index( \$seq, "T" );
            my %combined;
                  % {$combined{"G"}}=%Gpos;
                  % {$combined{"A"}}=%Apos;
                  % {$combined{"C"}}=%Cpos;
                  % {$combined{"T"}}=%Tpos;
            my %dont_care_ind=();
            my %dont_care_ind_right=();
            my %PAMindex = ();
            
            ###########################################################################################################################################################################
            # Single Sequence
            ###########################################################################################################################################################################
            
            if ($something{"kind"} eq "single") {
                  
                  #####################################################################################################################################################################
                  # Foward Sequence Calculations
                  #####################################################################################################################################################################
                  
                  if ($something{"preceding"} eq "A") {
                        %dont_care_ind=%Apos;
                  }elsif($something{"preceding"} eq "G"){
                        %dont_care_ind=%Gpos;
                  }elsif($something{"preceding"} eq "C"){
                        %dont_care_ind=%Cpos;
                  }elsif($something{"preceding"} eq "T"){
                        %dont_care_ind=%Tpos;
                  }else{
                        %dont_care_ind=(%Gpos,%Tpos,%Cpos,%Apos);
                  }
                  if ($something{"PAM"} eq "NAG") {
                        %PAMindex=%Apos;
                  } elsif ($something{"PAM"} eq "NGG") {
                        %PAMindex=%Gpos;
                  } else{
                        %PAMindex=(%Apos,%Gpos);
                  }
                  POSLOOP: foreach my $Gposind ( sort( keys(%dont_care_ind) ) ) {
                        LENGTHLOOP: foreach my $length ( ($something{"min_length"}+1) .. ($something{"max_length"}+1) ) {
                              if ( exists $PAMindex{ ( $Gposind + $length ) } && exists $Gpos{ ( $Gposind + $length + 1 ) } ) {
                                    my $taleseq = substr( $seq, $Gposind, $length + 2 );
                                    my @flank_array = find_base_count( $taleseq );
                                    $tempstatistics{"Total number of possible designs"}++;
                                    if (  $something{"min_A"} < $flank_array[0] && $something{"max_A"} > $flank_array[0] &&
                                          $something{"min_C"} < $flank_array[1] && $something{"max_C"} > $flank_array[1] &&
                                          $something{"min_T"} < $flank_array[2] && $something{"max_T"} > $flank_array[2] &&
                                          $something{"min_G"} < $flank_array[3] && $something{"max_G"} > $flank_array[3] &&
                                          !($taleseq=~m/TTTTT/) 
                                    ) {
                                          my $name = $seq_obj->display_id;
                                          #$name =~ s/_/:/ig;
                                          $name .= "_" . $count . "_" . $cut;
                                          my @new_score=(0,0,0);
                                          if(exists($Gpos{$Gposind}) && exists($Gpos{$Gposind+1})){
                                              $new_score[2]++;
                                          }
                                          if(exists($Gpos{$Gposind})){
                                              $new_score[2]++;
                                          }
                                          if(($flank_array[3]+$flank_array[1])>80){
                                              $new_score[2]--;
                                          }
                                          @flank_array = find_base_count( substr( $taleseq, $length-7,6) );
                                          if(($flank_array[3]+$flank_array[1])>70){
                                              $new_score[2]++;
                                          }
                                          
                                         
                                          $new_score[2]=$new_score[2]+score_micro_homolgy(\%combined,30,( $Gposind + $length +2 -5 ),5,\$seq); 
                                          ${ $CRISPR_hash{$name} }{"start"} = ($Gposind) + $cut;
                                          ${ $CRISPR_hash{$name} }{"end"} = ( $Gposind + $length + 2 ) + $cut;
                                          ${ $CRISPR_hash{$name} }{"length"} = $length + 2;
                                          my $start = ${ $CRISPR_hash{$name} }{"start"} + $location_offset;
                                          my $end = ${ $CRISPR_hash{$name} }{"end"} + $location_offset;
                                          my %score = calculate_CRISPR_score(\%trees, \%something, $end-3, $end, $chrom, 1, \@new_score);
                                          
                                          #############################################################################################################################################
                                          #Statistics
                                          #############################################################################################################################################
                                          
                                          if (make_CRISPR_statistics(\%something, \%score, $dont_asses_context, \%tempstatistics) == 1){
                                                delete $CRISPR_hash{$name};
                                                next LENGTHLOOP;
                                          }
                                          
                                          if ( exists $something{"retrieve_recomb_matrix"} ) {
                                                ${ ${ $CRISPR_hash{$name} }{"homology"} }{"left"} = substr( $whole_seq, ( ${ $CRISPR_hash{$name} }{"start"} - $something{"left_homology"} ), ($something{"left_homology"}) );
                                                ${ ${ $CRISPR_hash{$name} }{"homology"} }{"right"} = substr( $whole_seq, ${ $CRISPR_hash{$name} }{"end"}, $something{"right_homology"} );
                                          }
                                          
                                          %{ ${ $CRISPR_hash{$name} }{"context"} } = %score;
                                          ${ $CRISPR_hash{$name} }{"nucseq"} = $taleseq;
                                          ${ $CRISPR_hash{$name} }{"strand"} = "plus";
                                          $count++;
                                          #############################################################################################################################################
                                          
                                    } else {
                                          $tempstatistics{"Number of designs excluded because their nucleotide composition was too invariable or contained TTTTT"}++;
                                          next LENGTHLOOP;
                                    }
                              }
                        }
                  }
                  
                  #####################################################################################################################################################################
                  # Backward Sequence Calculations
                  #####################################################################################################################################################################
                  
                  if ($something{"preceding"} eq "A") {
                        %dont_care_ind=%Tpos;
                  } elsif($something{"preceding"} eq "G"){
                        %dont_care_ind=%Cpos;
                  } elsif($something{"preceding"} eq "C"){
                        %dont_care_ind=%Gpos;
                  } elsif($something{"preceding"} eq "T"){
                        %dont_care_ind=%Apos;
                  } else{
                        %dont_care_ind=(%Gpos,%Tpos,%Cpos,%Apos);
                  }
                  if ($something{"PAM"} eq "NAG") {
                        %PAMindex=%Tpos;
                  } elsif ($something{"PAM"} eq "NGG") {
                        %PAMindex=%Cpos;
                  } else{
                        %PAMindex=(%Tpos,%Cpos);
                  }
                  POSLOOP: foreach my $Cposind ( sort( keys(%Cpos) ) ) {
                        LENGTHLOOP: foreach my $length ( ($something{"min_length"}+1) .. ($something{"max_length"}+1) ) {
                              if ( exists $PAMindex{ ( $Cposind + 1 ) } && exists $dont_care_ind{ ( $Cposind + $length + 1 ) } ) {
                                    my $taleseq = substr( $seq, $Cposind, $length + 2 );
                                    my @flank_array = find_base_count( $taleseq );
                                    $tempstatistics{"Total number of possible designs"}++;
                                    if (  $something{"min_A"} < $flank_array[0] && $something{"max_A"} > $flank_array[0] &&
                                          $something{"min_C"} < $flank_array[1] && $something{"max_C"} > $flank_array[1] &&
                                          $something{"min_T"} < $flank_array[2] && $something{"max_T"} > $flank_array[2] &&
                                          $something{"min_G"} < $flank_array[3] && $something{"max_G"} > $flank_array[3] &&
                                          !($taleseq=~m/AAAAA/) &&
                                          !($taleseq=~m/TTTTT/)
                                    ){
                                          my $name = $seq_obj->display_id;
                                          #$name =~ s/_/:/ig;
                                          $name .= "_" . $count . "_" . $cut;
                                          my @new_score=(0,0,0);
                                          if(exists($Cpos{($Cposind + $length + 1)}) && exists($Cpos{($Cposind + $length)})){
                                                $new_score[2]++;
                                          }
                                          if(exists( $Cpos{($Cposind + $length + 1)}) ){
                                                $new_score[2]++;
                                          }
                                          if(($flank_array[3]+$flank_array[1])>80){
                                                $new_score[2]--;
                                          }
                                          @flank_array = find_base_count( substr( $taleseq, 0,6) );
                                          if(($flank_array[3]+$flank_array[1])>70){
                                                $new_score[2]++;
                                          }
                                          $new_score[2]=$new_score[2]+score_micro_homolgy(\%combined,30,( $Cposind + 5 ),5,\$seq);
                                          ${ $CRISPR_hash{$name} }{"start"} = ($Cposind) + $cut;
                                          ${ $CRISPR_hash{$name} }{"end"} = ( $Cposind + $length + 2 ) + $cut;
                                          ${ $CRISPR_hash{$name} }{"length"} = $length + 2;
                                          my $start = ${ $CRISPR_hash{$name} }{"start"} + $location_offset;
                                          my $end = ${ $CRISPR_hash{$name} }{"end"} + $location_offset;
                                          my %score = calculate_CRISPR_score(\%trees, \%something, $start, $start+3, $chrom, 0,\@new_score);
                                          
                                          #############################################################################################################################################
                                          #Statistics
                                          #############################################################################################################################################
                                          
                                          if (make_CRISPR_statistics(\%something, \%score, $dont_asses_context, \%tempstatistics) == 1){
                                                delete $CRISPR_hash{$name};
                                                next LENGTHLOOP;
                                          }
                                          
                                          if ( exists $something{"retrieve_recomb_matrix"} ) {
                                                ${ ${ $CRISPR_hash{$name} }{"homology"} }{"left"} = substr( $whole_seq, ( ${ $CRISPR_hash{$name} }{"start"} - $something{"left_homology"} ), ($something{"left_homology"}) );
                                                ${ ${ $CRISPR_hash{$name} }{"homology"} }{"right"} = substr( $whole_seq, ${ $CRISPR_hash{$name} }{"end"}, $something{"right_homology"} );
                                          }
                                          
                                          %{ ${ $CRISPR_hash{$name} }{"context"} } = %score;
                                          ${ $CRISPR_hash{$name} }{"nucseq"} = $taleseq;
                                          ${ $CRISPR_hash{$name} }{"strand"} = "minus";
                                          $count++;
                                          
                                          #############################################################################################################################################
                                          
                                    } else {
                                          $tempstatistics{"Number of designs excluded because their nucleotide composition was too invariable or contained TTTTT"}++;
                                          next LENGTHLOOP;
                                    }
                              }
                        }
                  }
            } else{
                  
                  #####################################################################################################################################################################
                  # Double Sequence - only forward calculations needed
                  #####################################################################################################################################################################
                  
                  if ($something{"preceding"} eq "A") {
                        %dont_care_ind=%Tpos;
                        %dont_care_ind_right=%Apos;
                  }elsif($something{"preceding"} eq "G"){
                        %dont_care_ind=%Cpos;
                        %dont_care_ind_right=%Gpos;
                  }elsif($something{"preceding"} eq "C"){
                        %dont_care_ind=%Gpos;
                        %dont_care_ind_right=%Cpos;
                  }elsif($something{"preceding"} eq "T"){
                        %dont_care_ind=%Apos;
                        %dont_care_ind_right=%Tpos;
                  }else{
                        %dont_care_ind=(%Gpos,%Tpos,%Cpos,%Apos);
                        %dont_care_ind_right=(%Gpos,%Tpos,%Cpos,%Apos);
                  }
                  my %PAMindex_right=();
                  if ($something{"PAM"} eq "NAG") {
                        %PAMindex=%Tpos;
                        my %PAMindex_right=%Apos;
                  } elsif ($something{"PAM"} eq "NGG") {
                        %PAMindex=%Cpos;
                        %PAMindex_right=%Gpos;
                  } else{
                        %PAMindex=(%Tpos,%Cpos);                        
                        %PAMindex_right=(%Gpos,%Apos);
                  }
                  POSLOOP: foreach  my $Cposind ( sort( keys(%Cpos) ) ) {
                        SPACERLOOP: foreach my $spacerlength ( ($something{"minspacerlength"}) .. ($something{"maxspacerlength"}) ) {
                              LENGTHLOOP: foreach  my $length ( ($something{"min_length"}+1) .. ($something{"max_length"}+1) ) {
                                    if ( exists $PAMindex{ ( $Cposind + 1 ) } && exists $dont_care_ind{ ( $Cposind + $length + 1 ) } && exists $dont_care_ind_right{ ( $Cposind + $length + 1 + $spacerlength ) } && exists $PAMindex_right{ ( $Cposind + $length +$length + 1 + $spacerlength ) } && exists $Gpos{ ( $Cposind + $length +$length + 2 + $spacerlength ) } ) {
                                          my $left_taleseq = substr( $seq, $Cposind, $length + 2 );
                                          my $right_taleseq = substr( $seq, ( $Cposind + $length + 1 + $spacerlength ), ( $length + 2));
                                          my $completeseq=$left_taleseq.$right_taleseq;
                                          my @flank_array = find_base_count( ($left_taleseq.$right_taleseq) );
                                          $tempstatistics{"Total number of possible designs"}++;
                                          if (  $something{"min_A"} < $flank_array[0] && $something{"max_A"} > $flank_array[0] &&
                                                $something{"min_C"} < $flank_array[1] && $something{"max_C"} > $flank_array[1] &&
                                                $something{"min_T"} < $flank_array[2] && $something{"max_T"} > $flank_array[2] &&
                                                $something{"min_G"} < $flank_array[3] && $something{"max_G"} > $flank_array[3] &&
                                                !($completeseq=~/TTTTT/)
                                          ) {
                                                my $name = $seq_obj->display_id;
                                                #$name =~ s/_/:/ig;
                                                $name .= "_" . $count . "_" . $cut;
                                                my @new_score=(0,0,0);
                                                if(exists($Cpos{( $Cposind + $length + 1 )}) && exists($Cpos{( $Cposind + $length)})){
                                                    $new_score[2]++;
                                                }
                                                if(exists($Cpos{( $Cposind + $length + 1 )})){
                                                    $new_score[2]++;
                                                }
                                                if(exists($Gpos{( $Cposind + $length + 1 + $spacerlength )}) && exists($Gpos{( $Cposind + $length + 2 + $spacerlength )})){
                                                    $new_score[2]++;
                                                }
                                                if(exists($Gpos{( $Cposind + $length + 1 + $spacerlength )})){
                                                    $new_score[2]++;
                                                }
                                                if(($flank_array[3]+$flank_array[1])>80){
                                                    $new_score[2]--;
                                                }
                                                
                                                @flank_array = find_base_count(substr( $left_taleseq, 0,6));
                                                if(($flank_array[3]+$flank_array[1])>70){
                                                    $new_score[2]++;
                                                }
                                                
                                                @flank_array = find_base_count(substr( $right_taleseq, $length-7,6));
                                                if(($flank_array[3]+$flank_array[1])>70){
                                                    $new_score[2]++;
                                                }
                                                $new_score[2]=$new_score[2]+score_micro_homolgy(\%combined,30,( $Cposind + $length +5 ),5,\$seq); 
                                                $new_score[2]=$new_score[2]+score_micro_homolgy(\%combined,30,( $Cposind + $length +2 + $spacerlength -5 ),5,\$seq);
                                                @{${ $CRISPR_hash{$name} }{"lengthcombo"}}=($length,$spacerlength);
                                                ${ $CRISPR_hash{$name} }{"start"} = ($Cposind) + $cut;
                                                ${ $CRISPR_hash{$name} }{"end"} = ( $Cposind + $length+$spacerlength+$length+2 + 2 ) + $cut;
                                                ${ $CRISPR_hash{$name} }{"length"} =  $length+$spacerlength+$length+2 + 2;
                                                my $start = ${ $CRISPR_hash{$name} }{"start"} + $location_offset;
                                                my $end = ${ $CRISPR_hash{$name} }{"end"} + $location_offset;
                                                my %score = calculate_CRISPR_score(\%trees, \%something, $start, $start+$length, $chrom, 0, \@new_score);
                                                
                                                #######################################################################################################################################
                                                #Statistics
                                                #######################################################################################################################################
                                                
                                                if (make_CRISPR_statistics(\%something, \%score, $dont_asses_context, \%tempstatistics) == 1){
                                                      delete $CRISPR_hash{$name};
                                                      next LENGTHLOOP;
                                                }
                                                
                                                if ( exists $something{"retrieve_recomb_matrix"} ) {
                                                      ${ ${ $CRISPR_hash{$name} }{"homology"} }{"left"} = substr( $whole_seq, ( ${ $CRISPR_hash{$name} }{"start"} - $something{"left_homology"} ), ($something{"left_homology"}) );
                                                      ${ ${ $CRISPR_hash{$name} }{"homology"} }{"right"} = substr( $whole_seq, ${ $CRISPR_hash{$name} }{"end"}, $something{"right_homology"} );
                                                }
                                                
                                                %{ ${ $CRISPR_hash{$name} }{"context"} } = %score;
                                                @{${ $CRISPR_hash{$name} }{"nucseq"}} = ($left_taleseq,$right_taleseq);
                                                ${ $CRISPR_hash{$name} }{"strand"} = "plus";
                                                $count++;
                                                
                                                #######################################################################################################################################
                                                
                                          } else {
                                                $tempstatistics{"Number of designs excluded because their nucleotide composition was too invariable or contained TTTTT"}++;
                                                next LENGTHLOOP;
                                          }
                                    }
                              }
                        }
                  }
            }
            
            ##########################################################################################################################################################################
            #store the CRISPR and the Statistics in temporary files - for each child process and end the fork
            ##########################################################################################################################################################################
            
            {
                  my $json = JSON::XS::encode_json(\%CRISPR_hash);
                  write_file( $temp_dir . "/" .$seq_obj->display_id . $cut . '.json', { binmode => ':raw' }, $json );
                  $json = JSON::XS::encode_json(\%tempstatistics);
                  write_file( $temp_dir . "/" . $seq_obj->display_id . $cut . 'stats.json', { binmode => ':raw' }, $json );
            }
            
            $pm->finish();
      }
      
      ##########################################################################################################################################################################
      #parent wait till all children are done and then rebuild the CRISPR and the Statistics out of the temporary files
      ##########################################################################################################################################################################
      
      $pm->wait_all_children();
      foreach  my $cut (@cuts) {
            my $json = read_file( $temp_dir . "/" .$seq_obj->display_id . $cut . '.json', { binmode => ':raw' } );
            %finished_CRISPR_hash = ( %finished_CRISPR_hash, %{ decode_json $json } );
            unlink $temp_dir . "/" . $seq_obj->display_id . $cut . ".json";
            $json = read_file( $temp_dir . "/" . $seq_obj->display_id . $cut . 'stats.json', { binmode => ':raw' } );
            my %sechash=%{ decode_json $json };
            foreach  my $seckey (keys(%sechash)){
                  $tempstatistics{$seckey}=$tempstatistics{$seckey}+$sechash{$seckey};
            }
            unlink $temp_dir . "/" . $seq_obj->display_id . $cut . "stats.json";
      }
      return (\%finished_CRISPR_hash,\%tempstatistics);
}

#########################################################################################
#name:      calculate_CRISPR_score
#function:  helper-function for find_and_print_CRISPRS
#           creates the score for the given CRISPRS
#input:     (builded tree as referemce, something-Hashreference, start, end, chrom
#           1/0 for $score{"CDS"}++)
#output:    the calculated score as hash
#########################################################################################
sub calculate_CRISPR_score {
      my %score = ();
      my @new_score=@{$_[6]};  
      my $expression="[";
      if (($_[1]->{"number_of_CDS"}>0)) {
         foreach my $number(1..$_[1]->{"number_of_CDS"}){
            $expression.=$number;
            }
         $expression.="]";
      }else{
             $expression="[.]";
      }
      $score{"CRISPRi"}=0;
      $score{"CRISPRa"}=0;
      my %transcripts=();
      if ( exists $_[0]->{$_[4]} ) { # check wethere the tree exists
            #search for annotations in the intervall from start (2) to end (3) and count them in score
            my $annotations = $_[0]->{$_[4]}->fetch( int($_[2]), int($_[3]) );
            
            foreach  my $anno ( @{$annotations} ) {
                  if ( $anno =~ m/gene_(\S+)_([0-9]+)_([0-9]+)/ ) {
                       my $gene_annotations = $_[0]->{$_[4]}->fetch( int($2), int($3) );                        
                        foreach  my $gene_anno ( @{$gene_annotations} ) {
                              if ($gene_anno=~m/exon::(\S+)::(\d+)_(\d+)_(\d+)/) {    
                                    ${$transcripts{$1}}{"exon".$2}=$3."_".$4;
                              }
                        }
                        last;
                  }
            }
            foreach  my $anno ( @{$annotations} ) {
                  print $anno."\n";
                  if ( $anno =~ m/gene_(\S+)_[0-9]+_[0-9]+/ ) {
                        $new_score[1]++;
                        ${ $score{"gene"} }{$1}++;                      
                  } elsif ( $anno =~ m/exon::(\S+)::(\d+)\_(\d+)_(\d+)/) {
                        ${ $score{"exon"} }{$2}++;
			$new_score[1]=$new_score[1]+5/$2;
                        if (exists $score{"transcripts"}) {
                              $score{"transcripts"}=$score{"transcripts"}."_".$1;
                        }else{
                              $score{"transcripts"}="_".$1;
                        }
                        ${$score{"transexons"}}{$1."::".$2}++;
                        
                  } elsif ( $anno =~ m/CpG/ ) {
                        $score{"CpG"}++;
			$new_score[1]--;
                  } elsif ( $anno =~ m/CDS::(\S+)::($expression)\_/ ) {
                        $new_score[1]=$new_score[1]+5/$2;
		       if($_[1]->{"specific_transcript"} ne "any"){
			   if ($1 eq $_[1]->{"specific_transcript"}) {
                                    $score{"CDS_1"}++;
                              }
                        } else{
                              $score{"CDS_1"}++;
                        }
                        if ($_[5] == 1) { # only for FORWARD needed @Flo ist das so gewollt?
                              $score{"CDS"}++;
                        }
                  } elsif ( $anno =~ m/CDS/ ) {
                        $score{"CDS"}++;
                        $new_score[1]++;
                  }
            }
            my $strand=1;
            my $first_start=0;
            my $first_end=0;
            my $second_start=0;
            foreach my $trans (keys %transcripts){
                  my %temp=%{$transcripts{$trans}};
                  my @first_temp=split("_",$temp{"exon1"});
                  $first_start=$first_temp[0];
                  $first_end=$first_temp[1];
                  my @second_temp=split("_",$temp{"exon2"});
                  $second_start=$second_temp[0];
                  if ($second_start<$first_start) {
                       $strand=0;
                  }                  
                  if (
                        $strand==1
                        && ($first_start+300)>=int($_[2])
                        && ($first_start-50)<=int($_[2])
                        ) {
                          $score{"CRISPRi"}=1;
                     }elsif(
                           $strand==0
                         && ($first_end-300)<=int($_[2])
                         && ($first_end+50)>=int($_[2])
                     ){
                           $score{"CRISPRi"}=1;
                     }elsif(
                           $strand==1
                         && ($first_start-400)<=int($_[2])
                         && ($first_start-50)>=int($_[2])
                     ){
                           $score{"CRISPRa"}=1;
                     }elsif(
                           $strand==0
                         && ($first_end+400)>=int($_[2])
                         && ($first_end+50)<=int($_[2])
                     ){
                           $score{"CRISPRa"}=1;
                     }
            print $strand."\t".$first_start."\t".int($_[2])."\t".$score{"CRISPRi"}."\t".$score{"CRISPRa"}."\n";    
            }
            
            #search for the start and stop coddon in the intervall from start (2) to end (3) with an up/downstream window
            $annotations = $_[0]->{$_[4]}->fetch( int( $_[2] - $_[1]->{"downstream_window"}), int( $_[3] + $_[1]->{"upstream_window"}) );
            foreach  my $anno ( @{$annotations} ) {
                  if ( $anno =~ m/start_codon/ ) {
                        $score{"start_codon"}++;
                        $new_score[1]++;
                  }elsif ( $anno =~ m/stop_codon/ ) {
                        $score{"stop_codon"}++;
                        $new_score[1]++;
                  }
            }                        
            
      }
      $score{"new_score"}=\@new_score;
      return %score;
}

#########################################################################################
#name:      make_CRISPR_statistics
#function:  helper-function for find_and_print_CRISPRS
#           adjust statistics for the given CRISPRS
#input:     (something-Hashreference, score-Hashreference, dont_asses_context,
#           statistics-hashreference)
#output:    1 if delete and next in loop is needed, 0 if not
#########################################################################################
sub make_CRISPR_statistics {
      if ( exists $_[0]->{"gene_exclusive"} && !exists $_[1]->{"gene"} && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they did not hit any gene"}++;
            return 1;
      }
      if ( exists $_[0]->{"exclude_overlapping_genes"} && scalar(keys(%{$_[1]->{"gene"}}))>1 && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they hit multiple genes"}++;
            return 1;
      }
      if ( exists $_[0]->{"exon_exclusive"} && !exists $_[1]->{"exon"} && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they did not hit any exon"}++;
            return 1;
      }
      if ( exists $_[0]->{"CpG_exclusive"} && exists $_[1]->{"CpG"} && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they were located in an CpG island"}++;
            return 1;
      }
      if ( exists $_[0]->{"purpose_exclusive"} && ($_[0]->{"purpose"} eq "Knockdown/-out") && !(exists $_[1]->{"CDS_1"}) && ($_[2] != 1) ) {
            $_[3]->{"Number of designs excluded because they were not directly behind the ATG of the specified transcript"}++;
            return 1;
      }
      if ( exists $_[0]->{"purpose_exclusive"} && ($_[0]->{"purpose"} eq "N-Terminal tagging") && !exists $_[1]->{"start_codon"} && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they were not located at the start codon"}++;
            return 1;
      }
      if ( exists $_[0]->{"purpose_exclusive"} && ($_[0]->{"purpose"} eq "C-Terminal tagging") && !exists $_[1]->{"stop_codon"} && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they were not located at the stop codon"}++;
            return 1;
      }
      if ( exists $_[0]->{"purpose_exclusive"} && ($_[0]->{"purpose"} eq "CRISPRa") && $_[1]->{"CRISPRa"}!=1 && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they were not amenable for CRISPRa"}++;
            return 1;
      }
      if ( exists $_[0]->{"purpose_exclusive"} && ($_[0]->{"purpose"} eq "CRISPRi") && $_[1]->{"CRISPRi"}!=1 && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they were not amenable for CRISPRi"}++;
            return 1;
      }
      if ( exists $_[0]->{"CDS_only"} && !(exists $_[1]->{"CDS"}) && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they were not located in a coding sequence"}++;
            return 1;
      }
      if ( exists $_[0]->{"exon_exclusive"} && !( $_[0]->{"specific_exon"} eq "any") && $_[2] != 1  ) {
            my $counter = 0;
            foreach  my $exon ( keys (%{ $_[1]->{"exon"} }) ) {
                  if ( $exon == $_[0]->{"specific_exon"} ) {
                        $counter++;
                        last;
                  }
            }
            if ( $counter == 0 ) {
                  $_[3]->{"Number of designs excluded because they did not hit the specific exon"}++;
                  return 1;
            }
      }
      if ( !( $_[0]->{"specific_transcript"} eq "any") && $_[2] != 1 ) {
            my $counter = 0;
            foreach  my $transcript ( split("_",$_[1]->{"transcripts"} ) ) {
                  if ( $transcript eq $_[0]->{"specific_transcript"} ) {
                        $counter++;
                        last;
                  }
            }
            if ( $counter == 0 ) {
                  $_[3]->{"Number of designs excluded because they did not hit the specific transcript"}++;
                  return 1;
            }
      }
}

#########################################################################################
#name:      make_pos_index
#function:  return an index of occurrences of a certain character in any given string
#input:     (string reference,character)
#output:    hash with character positions as keys
#########################################################################################
sub make_pos_index {
      my %pos     = ();
      my $result  = index( ${$_[0]}, $_[1], 0);
      while ( $result != -1 ) {
            $pos{$result}++;
            $result = index( ${$_[0]}, $_[1], ($result + 1) );
      }
      return %pos;
}

#########################################################################################
#name:      find_base_count
#function:  return an array containing the percentage each character in the string
#input:     (string)
#output:    array of percentages
#########################################################################################
sub find_base_count {
      my $seq  = $_[0];
      return      int( ($seq =~ tr/A/x/) * 100 / length($seq) ),
                  int( ($seq =~ tr/C/x/) * 100 / length($seq) ),
                  int( ($seq =~ tr/T/x/) * 100 / length($seq) ),
                  int( ($seq =~ tr/G/x/) * 100 / length($seq) );
}

#########################################################################################
#name:      reverse_comp
#function:  return the reverse complement of a DNA sequence
#input:     (string)
#output:    reverse complement string
#########################################################################################
sub reverse_comp {
      ( my $rev = reverse $_[0] ) =~ tr/ACGTacgt/TGCAtgca/;
      return $rev;
}

#########################################################################################
#name:      mean
#function:  return the mean of an array of number
#input:     (string)
#output:    mean as string
#########################################################################################
sub mean {
      return sum(@_) / @_;
}

#########################################################################################
#name:      variance
#function:  return the variance of an array of number
#input:     (string)
#output:    variance as string
#########################################################################################
sub variance {
      return ( sum( map { ( $_ - mean( @{ $_[0] } ) )**2 } @{ $_[0] } ) / @{ $_[0] } );
}

#########################################################################################
#name:      gene_label
#function:  return the last value in the features  primary tag (Bio::SeqFeature)
#input:     (Bio::SeqFeature::Generic->new)
#output:    label as string
#########################################################################################
sub gene_label {
      my $feature=shift;
      if ($feature->can("primary_tag") ) {
            my $notes=$feature->primary_tag;
            my @notes=split("::",$notes);
            $notes=pop(@notes);
            $notes;
      }else{
            my $notes="";
            $notes;
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
#input:     (float number as string, digits after komma <int>) 
#output:    shortened float
#########################################################################################
sub round_digits{
      if($_[0]=~m/(\d+\.\d{$_[1]}).*/){
            return $1;
      }
      return $_[0];
}

#########################################################################################
#name:      make_mismatch_string
#function:  convert a SAM file format mismatch string into an "Mismatch string" 
#           M for every match X for every mismatch D for every deletion I for every 
#           insertion
#input:     (Sequence as string-reference, $something{"unspecific_leading_bases"},
#           direction) 
#output:    shortened float
#########################################################################################
sub make_mismatch_string{
      my $mismatchstring      = "";
      my $pos                 = 0;
      my @stringarray=split("\t",${$_[0]});
      if(${$_[0]}=~m/MD:Z:(\S+)\s/){
            $mismatchstring=$1;
      }
      my @matchstring=split("",$stringarray[9]);
      my @matches=$stringarray[5]=~m/[0-9]+[MID]/g;
      foreach my $match (@matches){
            $match=~/([0-9]+)([MID])/;
            foreach (1..$1){
                  $matchstring[$pos]=$2;
                  $pos++;
            }
      }
      @matches    = $mismatchstring =~m/[0-9]+|[\^A-Z]+|[0-9]+$/g;
      $pos        = 0;
      foreach my $match (@matches){
            if($match=~/([0-9]+)/){
                  $pos += $1;
            }elsif($match=~/^[A-Z]$/){
                  $matchstring[$pos]="X";
                  $pos++;
            }else{
                  $pos++
            }
      }
      if ($_[2] eq "fw") {
            foreach (1..$_[1]){
                  unshift @matchstring , "n";
            }
      }
      else {
            foreach (1..$_[1]){
                  push @matchstring , "n";
            }
      }
      
      
      return(@matchstring);
}

#########################################################################################
#name:      print_offtarget_string
#function:  convert a "Mismatch string" into html output with different color-highlighting
#           for M,X,D and I
#input:     (Mismatch string) 
#output:    converted string for html output
#########################################################################################
sub print_offtarget_string {
      my $string=$_[0];
      $string=~s/n/<span style="color: black;">n<\/span>/g;
      $string=~s/M/<span style="color: lightgreen;">M<\/span>/g;
      $string=~s/X/<span style="color: red;">X<\/span>/g;
      $string=~s/D/<span style="color: pink;">D<\/span>/g;
      $string=~s/I/<span style="color: pink;">I<\/span>/g;
      return $string;
}
#########################################################################################
#name:      score_micro_homolgy
#function:  calculate an microhomology score between 0 and 12 for a qiven position in a 
# indexed sequence
#input:     (indices hash of hashes, threshold int,length int, sequence ref to string ) 
#output:    micro homolgy score (int)
#########################################################################################
sub score_micro_homolgy {
	my $score=0;
	my $count=0;
	my $lengthhom=1;
	my $seq;
	#my $position=$_[2];
	#my $threshold=$_[1];
	#my $length=$_[3];
	#my $postitions=%{$_[0]};
	my $stuff=0;
	my $outframe=1;
	my $inframe=1;
	my $limit_right=$_[2]+$_[1];
	if(($_[2]+$_[1]+$_[3])>length(${$_[4]})){
	    $limit_right=length(${$_[4]})-$_[3];
	}
	RIGTHSTART:foreach my $rightarmstart ($_[2]..($limit_right)){
		LENGTH:foreach my $length (2..$_[3]){
			if($lengthhom>0){
				$lengthhom=0;
				my @right_seq=split("",substr(${$_[4]},$rightarmstart,$length));
				$stuff=$length+1;
				LEFTARM: while($stuff<$_[1]){
					$stuff++;
					$count=0;
					$seq="";
					foreach my $letter (@right_seq){
						if(exists(${${$_[0]}{$letter}}{($_[2]-$stuff+$count)})){
							$seq.=$letter;
							$count++;
						}else{
							next LEFTARM;
						}
					}
					if($count==scalar(@right_seq)){
						$lengthhom=1;
						my $gaplength=($stuff+($rightarmstart-$_[2]));
						if($gaplength%3 == 0){
							$inframe=$inframe+($count*exp(0.1*(-$gaplength)));
						}else{
							$outframe=$outframe+($count*exp(0.1*(-$gaplength)));
						}
						#print "(@right_seq): $seq ".scalar(@right_seq)." ".$count." $stuff $position $rightarmstart hit"."\n";
						#print substr($string,($position-20),40)."\n";
						#print substr($string,($position-20),(20-$stuff+$count));
						#print '-' x ($stuff+($rightarmstart-$position-$count)); 
						#print substr($string,($rightarmstart),20-($rightarmstart-$position))."\n";
					}else{
						next LEFTARM;
					}
				}
			}else{
				$lengthhom=1;
				next RIGTHSTART;
			}
		}
	}
	return (10*log($outframe/$inframe));
}
#########################################################################################
#name:      hor_bar_chart
#function:  print horizontal bar chart
#input:     (@array(score 1, score 2, score 3))
#output:    N/A
#########################################################################################
sub hor_bar_chart {
	my $outstring='
	 <table>
      <tbody>
        <tr>
          <td style="font-size: 0.6em;">S</td><td class="bar"><div style="width: '.$_[0].'px;" class="item1"></div></td>
           </tr><tr>
          <td style="font-size: 0.6em;">A</td><td class="bar"><div style="width: '.$_[1].'px;" class="item2"></div></td>
           </tr><tr>
          <td style="font-size: 0.6em;">E</td><td class="bar"><div style="width: '.$_[2].'px;" class="item3"></div></td>
        </tr>
      </tbody>
    </table> ';
    return($outstring);
}
