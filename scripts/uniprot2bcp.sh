#!/bin/sh
dbxref2bcp.pl -f dbxref.txt
gzip -cd uniref100.xml.gz | unirefxml2cl.pl > uniref100_clusters.pfa
gzip -cd uniprot_sprot.dat.gz uniprot_trembl.dat.gz | uniprot2bcp.pl -t taxons.tab uniref100_clusters.pfa
