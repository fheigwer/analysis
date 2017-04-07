/**
 * toggle table onclick
 */
function tableShow(){
	obj = document.getElementById("Advanced_search");
	visible = (obj.style.display != "none");
	if(visible){
		obj.style.display = "none";
	} else {
		obj.style.display = "block";
	}
}

/**
 * change parameter preset to relaxed
 */
function change2relaxed(){
    //change parameters of section design purpose
    document.getElementById("purpose_exclusive").checked=false;
    document.getElementById("min_length").value=20;
    document.getElementById("max_length").value=20;        
    document.getElementById("min_G").value=1;
    document.getElementById("max_G").value=80;
    document.getElementById("min_A").value=1;
    document.getElementById("max_A").value=80;    
    document.getElementById("min_C").value=1;    
    document.getElementById("max_C").value=80;
    document.getElementById("min_T").value=1;
    document.getElementById("max_T").value=80;
    document.getElementById("right_homology").value=500;
    document.getElementById("left_homology").value=500; 
    document.getElementById("downstream_window").value=200; 
    document.getElementById("upstream_window").value=200;
    document.getElementById("number_of_CDS").value=3; 
    document.getElementById("minspacerlength").value=10;
    document.getElementById("maxspacerlength").value=20;
    document.getElementById("preceding").selectedIndex=15;
    document.getElementById("PAM").selectedIndex=7;
    
    //change parameters of section gene annotation filtering
    document.getElementById("gene_exclusive").checked=false;
    document.getElementById("exclude_overlapping_genes").checked=false;
    document.getElementById("exon_exclusive").checked=false;
    document.getElementById("CDS_only").checked=false;
    document.getElementById("CpG_exclusive").checked=false;
    document.getElementById("specific_transcript").value="any";
    document.getElementById("specific_exon").value="any";
    document.getElementById("retrieve_recomb_matrix").checked=false;
    document.getElementById("do_restriction_analysis").checked=false;
    document.getElementById("enzymes").selectedIndex=0;  
  
    //change parameters of section off-target analysis
    document.getElementById("bowtie_version").selectedIndex=1;
    document.getElementById("offtargetdb").selectedIndex=0;
    document.getElementById("off-targets-allowed").value=20;
    document.getElementById("unspecific_leading_bases").value=0;
    document.getElementById("edit_distance_allowed").value=0;
    document.getElementById("bowtie_mode").selectedIndex=0;
    document.getElementById("sec_off_target").checked=false;
    
    //change parameters of section output
    document.getElementById("max_per_exon").value=20;
    document.getElementById("out_gff").checked=false;
    document.getElementById("draw_image").checked=true;
    document.getElementById("print_table").checked=true;
    document.getElementById("match_info").checked=true;
}


/**
 * change parameter preset to medium
 */
function change2medium(){
    //change parameters of section design purpose
    document.getElementById("purpose_exclusive").checked=false;
    document.getElementById("min_length").value=20;
    document.getElementById("max_length").value=20;        
    document.getElementById("min_G").value=1;
    document.getElementById("max_G").value=80;
    document.getElementById("min_A").value=1;
    document.getElementById("max_A").value=80;    
    document.getElementById("min_C").value=1;    
    document.getElementById("max_C").value=80;
    document.getElementById("min_T").value=1;
    document.getElementById("max_T").value=80;
    document.getElementById("right_homology").value=500;
    document.getElementById("left_homology").value=500; 
    document.getElementById("downstream_window").value=100; 
    document.getElementById("upstream_window").value=100;
    document.getElementById("number_of_CDS").value=2; 
    document.getElementById("minspacerlength").value=13;
    document.getElementById("maxspacerlength").value=19;
    document.getElementById("preceding").selectedIndex=0;
    document.getElementById("PAM").selectedIndex=0;
    
    //change parameters of section gene annotation filtering
    document.getElementById("gene_exclusive").checked=false;
    document.getElementById("exclude_overlapping_genes").checked=false;
    document.getElementById("exon_exclusive").checked=true;
    document.getElementById("CDS_only").checked=false;
    document.getElementById("CpG_exclusive").checked=true;
    document.getElementById("specific_transcript").value="any";
    document.getElementById("specific_exon").value="any";
    document.getElementById("retrieve_recomb_matrix").checked=false;
    document.getElementById("do_restriction_analysis").checked=false;
    document.getElementById("enzymes").selectedIndex=0;  
  
    //change parameters of section off-target analysis
    document.getElementById("bowtie_version").selectedIndex=0;
    document.getElementById("offtargetdb").selectedIndex=0;
    document.getElementById("off-targets-allowed").value=10;
    document.getElementById("unspecific_leading_bases").value=1;
    document.getElementById("edit_distance_allowed").value=1;
    document.getElementById("bowtie_mode").selectedIndex=1;
    document.getElementById("sec_off_target").checked=false;
    
    //change parameters of section output
    document.getElementById("max_per_exon").value=10;
    document.getElementById("out_gff").checked=false;
    document.getElementById("draw_image").checked=true;
    document.getElementById("print_table").checked=true;
    document.getElementById("match_info").checked=true;
}


/**
 * change parameter preset to strict
 */
function change2strict(){
    //change parameters of section design purpose
    document.getElementById("purpose_exclusive").checked=true;
    document.getElementById("min_length").value=20;
    document.getElementById("max_length").value=20;        
    document.getElementById("min_G").value=1;
    document.getElementById("max_G").value=80;
    document.getElementById("min_A").value=1;
    document.getElementById("max_A").value=80;    
    document.getElementById("min_C").value=1;    
    document.getElementById("max_C").value=80;
    document.getElementById("min_T").value=1;
    document.getElementById("max_T").value=80;
    document.getElementById("right_homology").value=500;
    document.getElementById("left_homology").value=500; 
    document.getElementById("downstream_window").value=50; 
    document.getElementById("upstream_window").value=50;
    document.getElementById("number_of_CDS").value=1; 
    document.getElementById("minspacerlength").value=15;
    document.getElementById("maxspacerlength").value=17;
    document.getElementById("preceding").selectedIndex=0;
    document.getElementById("PAM").selectedIndex=0;
    
    //change parameters of section gene annotation filtering
    document.getElementById("gene_exclusive").checked=true;
    document.getElementById("exclude_overlapping_genes").checked=true;
    document.getElementById("exon_exclusive").checked=true;
    document.getElementById("CDS_only").checked=true;
    document.getElementById("CpG_exclusive").checked=true;
    document.getElementById("specific_transcript").value="any";
    document.getElementById("specific_exon").value="any";
    document.getElementById("retrieve_recomb_matrix").checked=false;
    document.getElementById("do_restriction_analysis").checked=false;
    document.getElementById("enzymes").selectedIndex=0;  
  
    //change parameters of section off-target analysis
    document.getElementById("bowtie_version").selectedIndex=0;
    document.getElementById("offtargetdb").selectedIndex=0;
    document.getElementById("off-targets-allowed").value=5;
    document.getElementById("unspecific_leading_bases").value=6;
    document.getElementById("edit_distance_allowed").value=2;
    document.getElementById("bowtie_mode").selectedIndex=3;
    document.getElementById("sec_off_target").checked=false;
    
    //change parameters of section output
    document.getElementById("max_per_exon").value=4;
    document.getElementById("out_gff").checked=false;
    document.getElementById("draw_image").checked=true;
    document.getElementById("print_table").checked=true;
    document.getElementById("match_info").checked=true;
}

function change2urldefined(){
	//change organism and gene
	if(getQueryVariable("ref_organism")){document.getElementById("ref_organism").selectedIndex=getQueryVariable("ref_organism");}
	if(getQueryVariable("data_type")){
		if(getQueryVariable("data_type")==0){
			document.getElementById("GENE.SYMB").checked=true;
			document.getElementById("FASTA.SEQ").checked=false;
		}else{
			document.getElementById("GENE.SYMB").checked=false;
			document.getElementById("FASTA.SEQ").checked=true;
		}
		document.getElementById("data_type").selectedIndex=getQueryVariable("data_type");
		}
	if(getQueryVariable("kind")){document.getElementById("kind").selectedIndex=getQueryVariable("kind");}
	if(getQueryVariable("inputtext")){document.getElementById("inputtext").value=getQueryVariable("inputtext");}
    //change parameters of section design purpose
    if(getQueryVariable("purpose_exclusive")){document.getElementById("purpose_exclusive").checked=getQueryVariable("purpose_exclusive");}
    if(getQueryVariable("min_length")){document.getElementById("min_length").value=getQueryVariable("min_length");}
    if(getQueryVariable("max_length")){document.getElementById("max_length").value=getQueryVariable("max_length");}     
    if(getQueryVariable("min_G")){document.getElementById("min_G").value=getQueryVariable("min_G");}
    if(getQueryVariable("max_G")){document.getElementById("max_G").value=getQueryVariable("max_G");}
    if(getQueryVariable("min_A")){document.getElementById("min_A").value=getQueryVariable("min_A");}
    if(getQueryVariable("max_A")){document.getElementById("max_A").value=getQueryVariable("max_A");}    
    if(getQueryVariable("min_C")){document.getElementById("min_C").value=getQueryVariable("min_C");}    
    if(getQueryVariable("max_C")){document.getElementById("max_C").value=getQueryVariable("max_C");}
    if(getQueryVariable("min_T")){document.getElementById("min_T").value=getQueryVariable("min_T");}
    if(getQueryVariable("max_T")){document.getElementById("max_T").value=getQueryVariable("max_T");}
    if(getQueryVariable("right_homology")){document.getElementById("right_homology").value=getQueryVariable("right_homology");}
    if(getQueryVariable("left_homology")){document.getElementById("left_homology").value=getQueryVariable("left_homology"); }
    if(getQueryVariable("downstream_window")){document.getElementById("downstream_window").value=getQueryVariable("downstream_window"); }
    if(getQueryVariable("upstream_window")){document.getElementById("upstream_window").value=getQueryVariable("upstream_window");}
    if(getQueryVariable("number_of_CDS")){document.getElementById("number_of_CDS").value=getQueryVariable("number_of_CDS"); }
    if(getQueryVariable("minspacerlength")){document.getElementById("minspacerlength").value=getQueryVariable("minspacerlength");}
    if(getQueryVariable("maxspacerlength")){document.getElementById("maxspacerlength").value=getQueryVariable("maxspacerlength");}
    if(getQueryVariable("preceding")){document.getElementById("preceding").selectedIndex=getQueryVariable("preceding");}
    if(getQueryVariable("PAM")){document.getElementById("PAM").selectedIndex=getQueryVariable("PAM");}
    
    //change parameters of section gene annotation filtering
    if(getQueryVariable("gene_exclusive")){document.getElementById("gene_exclusive").checked=getQueryVariable("gene_exclusive");}
    if(getQueryVariable("exclude_overlapping_genes")){document.getElementById("exclude_overlapping_genes").checked=getQueryVariable("exclude_overlapping_genes");}
    if(getQueryVariable("exon_exclusive")){document.getElementById("exon_exclusive").checked=getQueryVariable("exon_exclusive");}
    if(getQueryVariable("CDS_only")){document.getElementById("CDS_only").checked=getQueryVariable("CDS_only");}
    if(getQueryVariable("CpG_exclusive")){document.getElementById("CpG_exclusive").checked=getQueryVariable("CpG_exclusive");}
    if(getQueryVariable("specific_transcript")){document.getElementById("specific_transcript").value=getQueryVariable("specific_transcript");}
    if(getQueryVariable("specific_exon")){document.getElementById("specific_exon").value=getQueryVariable("specific_exon");}
    if(getQueryVariable("retrieve_recomb_matrix")){document.getElementById("retrieve_recomb_matrix").checked=getQueryVariable("retrieve_recomb_matrix");}
    if(getQueryVariable("do_restriction_analysis")){document.getElementById("do_restriction_analysis").checked=getQueryVariable("do_restriction_analysis");}
    if(getQueryVariable("enzymes")){document.getElementById("enzymes").selectedIndex=getQueryVariable("enzymes");  }
  
    //change parameters of section off-target analysis
    if(getQueryVariable("bowtie_version")){document.getElementById("bowtie_version").selectedIndex=getQueryVariable("bowtie_version");}
    if(getQueryVariable("offtargetdb")){document.getElementById("offtargetdb").selectedIndex=getQueryVariable("offtargetdb");}
    if(getQueryVariable("off-targets-allowed")){document.getElementById("off-targets-allowed").value=getQueryVariable("off-targets-allowed");}
    if(getQueryVariable("unspecific_leading_bases")){document.getElementById("unspecific_leading_bases").value=getQueryVariable("unspecific_leading_bases");}
    if(getQueryVariable("edit_distance_allowed")){document.getElementById("edit_distance_allowed").value=getQueryVariable("edit_distance_allowed");}
    if(getQueryVariable("bowtie_mode")){document.getElementById("bowtie_mode").selectedIndex=getQueryVariable("bowtie_mode");}
    if(getQueryVariable("sec_off_target")){document.getElementById("sec_off_target").checked=getQueryVariable("sec_off_target");}
    
    //change parameters of section output
    if(getQueryVariable("max_per_exon")){document.getElementById("max_per_exon").value=getQueryVariable("max_per_exon");}
    if(getQueryVariable("out_gff")){document.getElementById("out_gff").checked=getQueryVariable("out_gff");}
    if(getQueryVariable("draw_image")){document.getElementById("draw_image").checked=getQueryVariable("draw_image");}
    if(getQueryVariable("print_table")){document.getElementById("print_table").checked=getQueryVariable("print_table");}
    if(getQueryVariable("match_info")){document.getElementById("match_info").checked=getQueryVariable("match_info");}
}

/**
 * change button color of active button to blue

function btnChange() {
	btn1 = document.getElementById("relaxed");
	btn2 = document.getElementById("medium");
	btn3 = document.getElementById("strict");.checked
	if (document.getElementById("max_per_exon").value==4) {
		btn1.style.color="black";
		btn2.style.color="black";
		btn3.style.color="blue";
	} else if (document.getElementById("max_per_exon").value==10) {
		btn1.style.color="black";
		btn2.style.color="blue";
		btn3.style.color="black";
	} else {
		btn1.style.color="blue";
		btn2.style.color="black";
		btn3.style.color="black";
	}
} */
function btnChange() {
	btn1 = document.getElementById("relaxed");
	btn2 = document.getElementById("medium");
	btn3 = document.getElementById("strict");
	if (document.getElementById("max_per_exon").value==4) {
		btn1.checked=false;
		btn2.checked=false;
		btn3.checked=true;
	} else if (document.getElementById("max_per_exon").value==10) {
		btn1.checked=false;
		btn2.checked=true;
		btn3.checked=false;
	} else {
		btn1.checked=true;
		btn2.checked=false;
		btn3.checked=false;
	}
}

/**
 * change button color of display-advanced-options-button to blue if active
 */	
function btnAdvanced() {
	btn = document.getElementById("advanced");
	if (btn.style.color !="blue") {
		btn.style.color="blue";
	} else {
		btn.style.color="black";
	}
}


/**
 * display which buttons are active
 */
function change() {
	obj = document.getElementById("gene_exclusive");
	btn = document.getElementById("advanced");
	if (obj.checked) {
		if (btn.style.color=="blue") {
			document.getElementById("search").innerHTML="Advanced search strict";
		} else {
			document.getElementById("search").innerHTML="Quick search strict";
		}
	} else {
		if (btn.style.color=="blue") {
			document.getElementById("search").innerHTML="Advanced search relaxed";
		} else {
			document.getElementById("search").innerHTML="Quick search relaxed";
		}		
	}
}


/**
 * read in URL variable from www.stuff.org?x=1&y=2   ....
 */
function getQueryVariable(variable) {
    var query = window.location.search.substring(1);
    var vars = query.split('&');
    for (var i = 0; i < vars.length; i++) {
        var pair = vars[i].split('=');
        if (decodeURIComponent(pair[0]) == variable) {
		if (decodeURIComponent(pair[1])=="false") {
			 return false;
		}else{
			return decodeURIComponent(pair[1]);
		}            
        }
    }
}

/*change checkboxes */
function changestate(variable){
        var obj = document.getElementById(variable);
        if (obj.checked==false) {
                obj.checked=true;
        }else{
        	obj.checked=false
        }
}
/*change checkboxes */
function setchecked(variable){
        var obj = document.getElementById(variable);
        obj.checked=true;
}
/*change checkboxes */
function setunchecked(variable){
        var obj = document.getElementById(variable);
	obj.checked=false
}

/**
 * change text boxes
 */
function changetext(elemid,text){
	document.getElementById(elemid).value=text;
}

/**
 * change selections (dropdown menues)
 */
function changesel(elemid,val){
    var sel = document.getElementById(elemid);
    var opts = sel.options;
    for(var opt, j = 0; opt = opts[j]; j++) {
        if(opt.value == val) {
            sel.selectedIndex = j;
            break;
        }
    }
}

/**
 * change selections of a class of elements
 */
function changeClass() {
	var objArray = document.getElementsByClassName('change');
	for(var i = 0; i < objArray.length; i++) {
		var obj = document.getElementById(objArray[i].id);
		if (obj.checked) {
			obj.checked=false;
		}
	}
}