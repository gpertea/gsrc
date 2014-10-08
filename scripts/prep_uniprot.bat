# prepare the taxonomy data/files first
cd /data1/igm3/projects/geo/protdb/ncbi_taxonomy
rm *.bcp *.log *.dmp taxdump.tar.gz
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

tar xvfz taxdump.tar.gz

ncbitax2bcp.pl

mybcpin -T -b common@IGM -I *.bcp

#------------ now we can update uniprot 

cd ~/ann/protdb/new/
rm *.log *.gz *.pfa *.lst dbxref.txt*
wget -i ../wget.urls -nv -o log_wget.log

#>>On a host with database access, generate the taxons.tab file:

echo 'select tax_id, name, c_name from taxon' | pullsql -b common@IGM -o taxons.tab

#--- create the bcp & fasta files
# (warning: this uses >4G RAM)
../uniprot2bcp.sh &

(verify the .bcp files were created)
#---- load the bcp files in the database

mybcpin -b common@IGM -T -I uniprot.bcp.bz2 uniprot_names.bcp.bz2
(add uniprot_xref.bcp.bz2 if you need it - that's big!)

#---- prepare the fasta files
# (warning: these also uses a few GB of memory!)
nrdb -d1 -o uniprot_nr.fa uniprot.fa &

uniprot_nrsort.pl uniprot_nr.fa > uniprot_nr_defsrt.fa

# (this will reorder collapsed entries so the more "informative" entry 
#  comes first; will also remove any entry with less than 10aa)

mv uniprot_nr_defsrt.fa uniprot.fa

cdbfasta -a -n2 uniprot.fa

#-- prepare the mammals protein db:

qtaxon -U -x2 Mammalia > up_mammals_nonpred.lst
# -- evidence level 2 means that there was experimental evidence for protein or mRNA 
# -- use level 3 if you want to include homology-predicted proteins

cut -f2 up_mammals_nonpred.lst | cdbyank uniprot.fa.cidx > uniprot_mammals_np.fa
###
# - remove duplicates (because nrdb collapsed entries would've been pulled multiple times!):
###
cdbfasta uniprot_mammals_np.fa
cdbyank -l uniprot_mammals_np.fa.cidx | sort -u | cdbyank uniprot_mammals_np.fa.cidx > unipr_mammals.fa
cdbfasta -a unipr_mammals.fa
fsize unipr_mammals.fa > unipr_mammals.fsize



# for birds (turkey+chicken family)
qtaxon -U -x3 Phasianidae > up_birds_nonpred.lst
cut -f2 up_birds_nonpred.lst | cdbyank uniprot.fa.cidx > uniprot_birds_np.fa

cdbfasta uniprot_birds_np.fa
cdbyank -l uniprot_birds_np.fa.cidx | sort -u | cdbyank uniprot_birds_np.fa.cidx > unipr_birds.fa
cdbfasta -A unipr_birds.fa
fsize unipr_birds.fa > unipr_birds.fsize
