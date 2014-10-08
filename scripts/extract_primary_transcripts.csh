#!/bin/tcsh -f
#
#use annbygff -w to "compact" the data first:
# annbygff -r in.gff3 -w in_wloci.gff3

grep -v -P '\t(locus|gene|exon|CDS)\t' ncbi_hg19_wloci.gff3 | fgrep 'ID=' | cut -f1,3,4,5,7,9 | \
 perl -ne '@t=split(/\t/);($id)=(m/ID=([^;]+)/);($g)=(m/gene_name=([^;]+)/);($p)=(m/product=([^;]+)/); \
 $l=$t[0].":".$t[2]."-".$t[3].$t[4];print "$l\t$id $t[1] $l";print " gene:$g" if $g;print " product:$p" if $p; \
 print "\n"' | sort -k1,1 -u | seqmanip -f- /fs/szdata/genomes/human_hg19/hg19.fa > hg19_primary_transcripts.fa
 
 #sort -k1,1u | seqmanip -f- /fs/szdata/genomes/human_hg19/hg19.fa > hg19_primary_transcripts.
