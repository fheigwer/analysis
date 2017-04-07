#!usr/bin/perl
use warnings;
use strict;



use CGI qw(:standard VARS);                             #ist das paket was abgeschickt wird
use CGI::Carp qw ( fatalsToBrowser );                   #kommunikation zwischen perl und html

my $query		= new CGI;                              #CGI eine pfadsprache zwischen perl und html
my @params		= $query->param;
my %something = ();
foreach my $element (@params) {
    $something{$element} = $query->param($element);
    													#in dem hash sind ALL meine inputs des users gespeichert 
    													#(als keys der Name des textfelds ('Targetsequence','PAM' und der zugehörige value )
}
print "Content-type:text/html\r\n\r\n"; 				#muss perl sagen, dass der ooutput wieder eine html seite ist


##########################################################################################
##########################################################################################
my $PAM;
my $PAM1;
my $lengthofPAM= length ($PAM);
my $error1 = "Please correct your PAM, please only use a lenght about 2-30 bp. ";
my $error2 = "Please correct your PAM. ";
my $error3 = "Please correct your Target sequence only use 'ACGTacgt'. ";
my $sequence;

my $sequence1;
my @sequence=();
my $position;
my $design;
my $zahl;
my $N;
my $Count2;
my $Count2rev;
my @positions2=();
my $reversesequence;
my $positionrev;
my @positions2rev=();
my $lengthofdesign;
my @positions3=();
my @positions3rev=();

##########################################################################################
#save entry as strings and uc it
        $PAM= $something{'textfeld'};
        $PAM1= uc$PAM;
        $sequence= $something{'targetsequence'};				#link zur sequence im texterae
        $sequence1= uc$sequence;
        ($reversesequence= $sequence1)=~ tr/ACGTU/TGCAA/; 		#kopiert den string nach$reverse.und ersetzt dort die buchstaben klammersetzung!
     
       
	if ($sequence1 =~ m/[^ACGT]/g){                         	 # test for targetsequence
        print "<p>$error3</p>\n";
        print $query->end_html();
        print "<p><a href='http://e-crisp-test.dkfz.de/E-CRISP/ECRISPpauline.html/'>Go back</a> </p>\n";
           die;
  	}
    elsif($sequence1 eq ""){									#test for empty targetsequence
    print "<p>$error3</p>\n";
        print $query->end_html();
        print "<p><a href='http://e-crisp-test.dkfz.de/E-CRISP/ECRISPpauline.html/'>Go back</a> </p>\n";
           die;				
    }
           
	elsif ((length ($PAM1)) == 0)	{
        goto MATCH2;
	}
    
    elsif ((length ($PAM))< 2) {                         		# test for lenght of entry
            print "<p>$error1 </p>\n";
           print $query->end_html();
            print "<p><a href='http://e-crisp-test.dkfz.de/E-CRISP/ECRISPpauline.html/'>Go back</a> </p>\n";
           die;
    }
    elsif ((length ($PAM))> 30) {                           	# test for lenght of entry
            print "<p>$error1</p>\n";
        print $query->end_html();
        print "<p><a href='http://e-crisp-test.dkfz.de/E-CRISP/ECRISPpauline.html/'>Go back</a> </p>\n";
        die;
    }
	elsif ($PAM1 =~ m/[^ACGTUKMSWRYBDHVN]/g){                   #test for Iupac Code
    print "<p>$error2</p>";
    print $query->end_html();
    print "<p><a href='http://e-crisp-test.dkfz.de/E-CRISP/ECRISPpauline.html/'>Go back</a> </p>\n";
    die;
    }

    elsif ( $PAM1 =~ /[ACGTUKMSWRYBDHVN]/g) {            		 #if entry was correct:convert entry in regular expression
          $PAM1 =~ s/K/[GT]/g;
           $PAM1 =~ s/M/[AC]/g;
           $PAM1 =~ s/S/[CG]/g;
           $PAM1 =~ s/W/[AT]/g;
           $PAM1 =~ s/R/[AG]/g;
           $PAM1 =~ s/Y/[CT]/g;
           $PAM1 =~ s/B/[^A]/g;
           $PAM1 =~ s/D/[^C]/g;
           $PAM1 =~ s/H/[^G]/g;
           $PAM1 =~ s/V/[^T|^U]/g;
           $PAM1 =~ s/N/[^U]/g;
           goto MATCH;
     }
          

  
##########################################################################################
#Search for Matches for target/reversed target sequence for 'textentry'-PAM MATCH


MATCH:

{   my @positions=();
    my $lengthofPAM= length ($PAM);
    my @positionsrev=();
    my $Count=0;
    my $Countrev=0;
    print "</head><body>\n";
    print "<h1>Your search</h1>\n";
    goto TARGET;
    
	TARGET:
    print "<p>Your target sequence (5'-3'): $sequence1.</p>\n";
    print "<p> Your reverse  target sequence (3'-5'): $reversesequence.</p>\n";
    goto PAMentry;
    
	PAMentry:
    $PAM= uc($something{'textfeld'});
    #print "<p>Your PAM sequence: $PAM1.</p>\n"; -->printet mir [][] reguläre ausdrücke will ich als user aber nicht
 	 print "<p>Your PAM sequence (5'-3'): $PAM.</p>\n";
  
	###################################################
	#target sequence
    while ($sequence1 =~ m/($PAM1)/g)  {
    $position = pos ($sequence1);
	$position-= ($lengthofPAM)-1;
            #print "found at position: ", $position-= ($lengthofPAM-1), "\n";
    push (@positions,$position);
    $Count++;
        }
        
	################################################
	#reverse sequence    
    while ($reversesequence =~ m/($PAM1)/g)  {
    $positionrev = pos ($reversesequence);
    $positionrev-= ($lengthofPAM)-1;
            #print "found at position: ", $position-= ($lengthofPAM-1), "\n";
    push (@positionsrev,$positionrev);
    $Countrev++;
        }
        
	################################################ 
	#output of all RESULTS    
        
    print "<p><br> MATCH RESULTS:</br></p>\n";
    print "<p>Target sequence: Your PAM was found $Count times. </p>"; # muss nach der schleife geprintet werden damit das $ count wie es zuletzt gespeichert wurde geprintet wird , würde dieser printbefehl oben vor der while schleife stehen, printet mit program immer $count=o !!!
    print "<p>Position(s):</p>\n";
    print join("<br>",@positions); 								# nicht so wie ich dachte: print join"<p>("\n",@positions)</p>";
    
   	$sequence1=~ s/($PAM1)/<span style=\"background:yellow\">$1<\/span>/g;
    print "<p>$sequence1</p>\n";
    
   	print "<p> Reverse target sequence: Your PAM was found $Countrev times. </p>\n";  
    print "<p> Position(s) :</p>\n";
    print join("<br>",@positionsrev); 
    
   	$reversesequence=~ s/($PAM1)/<span style=\"background:yellow\">$1<\/span>/g;
    print "<p>$reversesequence</p>\n"; 
    print "<p><a href='http://e-crisp-test.dkfz.de/E-CRISP/ECRISPpauline.html/'>Go back</a> </p>\n";
    #achtung wegen der " " in perl!
    
    print "</body></html>\n";
    exit;
}



##########################################################################################
#search for Matches for target/reversed target sequence for
# designed PAM for the case 'length_of_N1'='length_of_N2'

MATCH2:

{	  $N= ""; 
	if ($something{'length_of_N1'} < $something{'length_of_N2'}){ 
        goto MATCH3; 
         }
                
    elsif ($something{'length_of_N1'} > $something{'length_of_N2'}){ 
    	 print "<p> Please correct your designed PAM<p>\n" ;
         print $query->end_html();
    	 print "<p><a href='http://e-crisp-test.dkfz.de/E-CRISP/ECRISPpauline.html/'>Go back</a> </p>\n";
   		 die;
         }
         
    elsif ($something{'length_of_N1'} eq $something{'length_of_N2'}) { 
        $zahl= $something{'length_of_N1'};
        #print $zahl;
       	for(my $i=0; $i < $zahl; $i++){
       		$N.="N";}
        	#print $N;
        	$design= "$something{'ref_PAM2'}"."$N"."$something{'ref_PAM3'}";
        	#print $design;
          } 
              
 		$lengthofdesign= length ($design);
   		#print "<p> $design</p>";
   		#my $sequence= $something{'targetsequence'};  #link zur sequence
    	$Count2=0;
    	$Count2rev=0;
     	print "</head><body>\n";
    	print "<h1>Your search</h1>\n";
   		goto TARGET;
    
	TARGET:
   		print "<p>Your target sequence (5'-3'): $sequence1.</p>\n";
   		print "<p> Your reverse  target sequence (3'-5'): $reversesequence.</p>\n";
    	goto PAMDESIGN;
    
	PAMDESIGN:
		print "<p>Your PAM sequence (5'-3'): $design.</p>\n"; 		#printet mit IUPAc
	
    ################################################
      #translation for matching with regular expression
 		   $design =~ s/K/[GT]/g;
           $design =~ s/M/[AC]/g;
           $design =~ s/S/[CG]/g;
           $design =~ s/W/[AT]/g;
           $design =~ s/R/[AG]/g;
           $design =~ s/Y/[CT]/g;
           $design =~ s/B/[^A]/g;
           $design =~ s/D/[^C]/g;
           $design =~ s/H/[^G]/g;
           $design =~ s/V/[^T|^U]/g;
           $design =~ s/N/[^U]/g;
       				
        
        
   	 ################################################
   	 #target sequence
        while ($sequence1 =~ m/($design)/g)  {
            $position = pos ($sequence1);
            $position-= ($lengthofdesign)-1;
           	# print "found at position: ", $position-= ($lengthofdesign-1), "\n";
            push (@positions2,$position);
           	$Count2++;
           }
      
      
  		#reverse target sequence   
   		while ($reversesequence =~ m/($design)/g)  {
            $positionrev = pos ($reversesequence);
            $positionrev-= ($lengthofdesign)-1;
           	# print "found at position: ", $position-= ($lengthofdesign-1), "\n";
            push (@positions2rev,$positionrev);
           	$Count2rev++;
           
        	}
     	################################################     
    	#Results of MATCH   
 		print "<p>Target sequence: Your PAM was found $Count2 times. </p>\n"; # muss nach der schleife geprintet werden damit das $ count wie es zuletzt gespeichert wurde geprintet wird , würde dieser printbefehl oben vor der while schleife stehen, printet mit program immer $count=o !!!
   		print "<p>Position(s):</p>\n";
   		print join("<br>",@positions2); # nicht so wie ich dachte: print join"<p>("\n",@positions)</p>";
   	 	$sequence1=~ s/($design)/<span style=\"background:yellow\">$1<\/span>/g;
    	print "<p>$sequence1</p>\n";
    	print "<p>Reverse target sequence: Your PAM was found $Count2rev times. </p>\n"; # muss nach der schleife geprintet werden damit das $ count wie es zuletzt gespeichert wurde geprintet wird , würde dieser printbefehl oben vor der while schleife stehen, printet mit program immer $count=o !!!
   		print "<p>Position(s):</p>\n";
    	print join("<br>",@positions2rev); 
    	$reversesequence=~ s/($design)/<span style=\"background:yellow\">$1<\/span>/g;
    	print "<p>$reversesequence</p>\n";
    	print "<p><a href='http://e-crisp-test.dkfz.de/E-CRISP/ECRISPpauline.html/'>Go back</a> </p>\n";
    	print "</body></html>\n";
    	exit;
}






##########################################################################################
# match for intervall-designed pam

MATCH3:
{	my $N;
	my $sequence2 = uc($sequence);
	my @design=();
	my $design2;
	my $zahl1= $something{'length_of_N1'};
	my $zahl2= $something{'length_of_N2'};
	
	my $lengthofdesign2;
	my $Count3=0;
	my $Count3rev=0;
	my @lengths=();
	my @sequences=();
	my $position3;
	my $position3rev;
	
	my $key = $zahl1;	# printen der richtigen pams!! startet ein intervall bei 5 soll das erste gespeicherte pam 5 ns enthalten 
	
	################################################     
	#save all the designed PAM of an defined interval( intervall 2-5 ergibt dann die pams ANCC ANNCC ANNNCC usw bis ANNNNNCC)
	
	
	for (0..$zahl1) {						#start with the first designed PAM
	$N="";
	$design2= "$something{'ref_PAM2'}"."$N"."$something{'ref_PAM3'}";
    }
	
	my %start= ('1'=>"",'2'=>"N",'3'=>"NN",'4'=>"NNN",'5'=>"NNNN",
				'6'=>"NNNNN",'7'=>"NNNNNN", '8'=>"NNNNNNN",'9'=>"NNNNNNNN",
				'10'=>"NNNNNNNNN",'11'=>"NNNNNNNNNN",'12'=>"NNNNNNNNNNN",'13'=>"NNNNNNNNNNNN",
				'14'=>"NNNNNNNNNNNNN",'15'=>"NNNNNNNNNNNNNN",'16'=>"NNNNNNNNNNNNNNN",'17'=>"NNNNNNNNNNNNNNNN",
				'18'=>"NNNNNNNNNNNNNNNNN",'19'=>"NNNNNNNNNNNNNNNNNN");
	$N= $start{$key};
	#print $N;
	
	
	for ($zahl1..($zahl2-1)){				#solange PAMS designen bis zum Ende des Intervals
	#$N= $N{$zahl1};
    #print "$something{'ref_PAM2'}";
    #for(1..$_){
    $N.="N";
    #}
    $design2= "$something{'ref_PAM2'}"."$N"."$something{'ref_PAM3'}";
    push(@design,$design2);					# im array@design speichere ich immer mit zeilenumbruch das jeweilige pam meines intervalls das ich definiert habe
    #$N="";
    }
    
    #################################################     
    #Results of MATCH 
    	
 	print "</head><body>\n";
    print "<h1>Your search</h1>\n";
    print "<p>Your Target sequence (5'-3'): $sequence2.</p>\n";
    print "<p> Your reverse  target sequence (3'-5'): $reversesequence.</p>\n";
	print "<p>Your Pams:</p>\n";

	print join("<br>",@design),"\n"; 			#printet PAMdesign mit IUPACcode


	foreach my $interval(@design){				#hier speichert mir die jeweilige pam länge ab in einem array,dielänge muss hier gespeichert werden, bevor die pams in regul.expression umgewandelt werden
	$lengthofdesign2= length ($interval);
	push (@lengths,$lengthofdesign2);
	}
	#print join ("<br>",@lengths), "\n"; 		 #die längen der iupacpams gespeichert
	
	
	
	#################################################
	# for each interval get match position

	foreach my $interval(@design) {
		$interval =~ s/N/[^U]/g;		#translation in regular expression,ab jetzt ändert sich auch die länge!
	}
     									#print join("<br>",@design); # printet mir die pams mit regulär expressions

		
	##################################################
	#matches target sequence
	
	foreach my $PAM (@design){ 				
 			while ($sequence2 =~ m/($PAM)/g)  {    
    			$position3 = length($`);
    			#print $position3;
    			push (@positions3,$position3);     
    			$Count3++; 	
    		}
    		if( ($sequence2 !~ m/($PAM)/)){
    			my $position3 = undef;	#speichert mir wenn kein treffergefunden wurde
   				# push (@positions3,$position);
   			}
   		    $sequence2=~ s/($PAM)/<span style=\"background:yellow\">$1<\/span>/g;
     }  
    										
   ##################################################
   #matches reverse target sequence 
   
    foreach my $PAM2(@design){
 			while ($reversesequence =~ m/($PAM2)/g)  {
       			$position3rev = length($`);					#$`printet mir sobald match die vordere sequenz, dielänge dieser sequenz die entspricht ja dann der position
        		push (@positions3rev,$position3rev);
        		$Count3rev++;
    		}
     		if( ($sequence2 !~ m/($PAM2)/)){
    			my $position3rev = undef;	
   			}
   		
     		$reversesequence=~ s/($PAM2)/<span style=\"background:yellow\">$1<\/span>/g;								
       }									
       										
   #################################################
   #MATCH RESULTS  
  
    #target sequence matches										
		print "<p> Match position(s) in target sequence:</p>\n";
		
  		print join ("<br>",@positions3),"\n";
  		print "<p>Your PAM was found $Count3 times.</p>\n";
   		#print join ("<br>",@sequences),"n";
   
  	    print "<p>$sequence2</p>\n";
    
    #reverse target sequence matches
    	print "<p> Match position(s) in reverse target sequence:</p>\n";
    	
  		print join ("<br>",@positions3rev),"\n";
  		print "<p>Your PAM was found $Count3rev times.</p>\n";
        print "<p>$reversesequence</p>\n";
    	print "<p><a href='http://e-crisp-test.dkfz.de/E-CRISP/ECRISPpauline.html/'>Go back</a> </p>\n";
  
   		print "</body></html>\n";
   		exit;
 }
 
   
##########################################################################################
##########################################################################################
 
 


















