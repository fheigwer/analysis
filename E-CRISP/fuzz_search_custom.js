function fuzzy_search(stuff){

if (document.getElementById('case').checked == true)
    {
        var cs=1
    }else{
        var cs=0
    }
var key=[];
 if (document.getElementById('ensembl').checked == true)
    {
        key.push('ENS_ID')
    }

if (document.getElementById('gene_symbol').checked == true)
    {
        key.push('GENE_SYMBOL')
    }
    
var options = {
  keys: key ,
  threshold: document.getElementById('search.sensitivity').value ,
  caseSensitive: cs
}

var f = new Fuse(library, options);
var result = f.search(document.getElementById(stuff).value);
var list = document.getElementById('results'); // create UL
list.innerHTML = '';
var i=0;

result.forEach(
	function(entry) {
		if(i<=20){
    		var li = document.createElement("li");
    		var textentry=document.createTextNode(entry.ENS_ID+"    "+entry.GENE_SYMBOL)
		if (/some/.test(entry.GENE_SYMBOL)) {
		    li.setAttribute("onclick","paste_to_chosen(\""+ entry.ENS_ID +"\")");
		}else{
		     li.setAttribute("onclick","paste_to_chosen(\""+ entry.GENE_SYMBOL +"\")");
		}
    		
		li.setAttribute("style","cursor: pointer;");
   		li.appendChild(textentry);
    		list.appendChild(li);
    		i=i+1;
   	 }
	}	
);
}

function paste_to_chosen(chosen){
	var list = document.getElementById('chosen'); // create UL
	var li = document.createElement("li");
    var textentry=document.createTextNode(chosen)
    li.setAttribute("onclick","remove_from_chosen(\""+chosen+"\")");
    li.setAttribute("id",chosen);
    li.setAttribute("style","cursor: pointer;");
    li.appendChild(textentry);
    list.appendChild(li);
    set_as_query('chosen')
}

function remove_from_chosen(chosen){
	document.getElementById(chosen).remove();
	set_as_query('chosen')
}

function set_as_query(chosen){
	var textfeld = document.getElementById('inputtext'); 
	textfeld.value='';
	var nums = document.getElementById("chosen");
	var listItem = nums.getElementsByTagName("li");		
	for (var i=0; i < listItem.length; i++) {
		if(i == 0){
    		textfeld.value=textfeld.value+''+listItem[i].innerHTML; 
    	}else{
    		textfeld.value=textfeld.value+';'+listItem[i].innerHTML; 
    	}
	}
}


Element.prototype.remove = function() {
    this.parentElement.removeChild(this);
}
NodeList.prototype.remove = HTMLCollection.prototype.remove = function() {
    for(var i = 0, len = this.length; i < len; i++) {
        if(this[i] && this[i].parentElement) {
            this[i].parentElement.removeChild(this[i]);
        }
    }
}

function loadjscssfile(filename, filetype){
 if (filetype=="js"){ //if filename is a external JavaScript file
  var fileref=document.createElement('script')
  fileref.setAttribute("type","text/javascript")
  fileref.setAttribute("src", filename)
 }
 else if (filetype=="css"){ //if filename is an external CSS file
  var fileref=document.createElement("link")
  fileref.setAttribute("rel", "stylesheet")
  fileref.setAttribute("type", "text/css")
  fileref.setAttribute("href", filename)
 }
 if (typeof fileref!="undefined")
  document.getElementsByTagName("head")[0].appendChild(fileref)
}

var filesadded="" //list of files already added

function checkloadjscssfile(filename, filetype){
 if (filesadded.indexOf("["+filename+"]")==-1){
  loadjscssfile(filename, filetype)
  filesadded+="["+filename+"]" //List of files added in the form "[filename1],[filename2],etc"
 }
}


function toggle_visibility(id) {
       var e = document.getElementById(id);
       if(e.style.display == 'block')
          e.style.display = 'none';
       else
          e.style.display = 'block';
}
