#!/usr/bin/bash

echo $1;
wget -r -nd --no-parent --reject "index.html*" ftp://ftp.ensembl.org/pub/release-72/genbank/$1/ ;
wget -r -nd --no-parent --reject "index.html*" ftp://ftp.ensembl.org/pub/release-72/gtf/$1/ ;
wget -r -nd --no-parent --reject "index.html*" ftp://ftp.ensembl.org/pub/release-72/fasta/$1/cdna/ ;
echo "done";

