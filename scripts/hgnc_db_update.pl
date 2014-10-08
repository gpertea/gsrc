#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q{Usage:
 hgnc_db_update.pl [-f <path/to/hgnc_download.txt>]

 Prepares data from input hgnc data into "bcp" tab delimited
 format for the following tables in database "common":

 hgnc
 hgnc_alias
 hgnc_xref
 
 Expected input columns:
 
 0. HGNC ID
 1. Approved Symbol
 2. Approved Name
 3. Status
 4. Previous Symbols
 5. Aliases
 6. Name Aliases
 7. Chromosome
 8. Accession Numbers
 9. Entrez Gene ID
10. Ensembl Gene ID
11. RefSeq IDs
12. Primary IDs
13. Secondary IDs
14. CCDS IDs
15. Entrez Gene ID (mapped data supplied by NCBI)
16. RefSeq (mapped data supplied by NCBI)
17. UniProt ID (mapped data supplied by UniProt)
18. Ensembl ID (mapped data supplied by Ensembl)
19. UCSC ID (mapped data supplied by UCSC)

Download url:

http://www.genenames.org/cgi-bin/hgnc_downloads.cgi?title=HGNC+output+data&col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_name_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_eg_id&col=gd_pub_ensembl_id&col=gd_pub_refseq_ids&col=gd_primary_ids&col=gd_secondary_ids&col=gd_ccds_ids&col=md_eg_id&col=md_refseq_id&col=md_prot_id&col=md_ensembl_id&col=md_ucsc_id&status=Approved&status=Entry+Withdrawn&status_opt=2&level=pri&=on&where=&order_by=gd_app_sym_sort&limit=&format=text&submit=submit&.cgifields=&.cgifields=level&.cgifields=chr&.cgifields=status&.cgifields=hgnc_dbtag
};
umask 0002;
getopts('f:o:') || die($usage."\n");
my $infile=$Getopt::Std::opt_f || 'hgnc_download.txt';
die "$usage Error: cannot locate input file $infile\n" unless -f $infile;
my $wgeturl='http://www.genenames.org/cgi-bin/hgnc_downloads.cgi?'.
   'title=HGNC+output+data&col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&'.
   'col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_name_aliases&'.
   'col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_eg_id&'.
   'col=gd_pub_ensembl_id&col=gd_pub_refseq_ids&col=gd_primary_ids&'.
   'col=gd_secondary_ids&col=gd_ccds_ids&col=md_eg_id&col=md_refseq_id&'.
   'col=md_prot_id&col=md_ensembl_id&col=md_ucsc_id&status=Approved&'.
   'status=Entry+Withdrawn&status_opt=2&level=pri&=on&where=&'.
   'order_by=gd_app_sym_sort&limit=&format=text&submit=submit&.cgifields=&'.
   '.cgifields=level&.cgifields=chr&.cgifields=status&.cgifields=hgnc_dbtag';
open(INF, $infile) || die("Error opening $infile\n");
open(BCP, '>hgnc.bcp') || die "Error creating file hgnc.bcp\n";
open(BCPA, '>hgnc_alias.bcp') || die "Error creating file hgnc_alias.bcp\n";
open(BCPX, '>hgnc_xref.bcp') || die "Error creating file hgnc_xref.bcp\n";

while (<INF>) {
 next if (m/^HGNC ID\t/); {
 chomp;
 my @t=split(/\t/);
 my $withdrawn=($t[1]=~s/[\~\-_ ]withdrawn$//i);
 $withdrawn=1 if $t[3]=~m/withdrawn$/i);
 my $replacedby=$1 if $withdrawn && $t[2]=~$m/see ([\w\.\-]+)$/;
 my $current=$withdrawn ? 0 : 1;
 # $current=1 if ($t[3]=~m/^approved$/i);
 my $entrez_geneid = $t[9] || $t[15];
 print BCP join("\t",@t[0..2], $t[7], $entrez_geneid, $current, $replacedby)."\n";
 
 }
close(INF);
# write go_term table
# --
#************ Subroutines **************

