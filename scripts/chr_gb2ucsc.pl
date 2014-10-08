#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 chr_gb2ucsc.pl <GRCh_sequence_report.txt> hg19.fa.fai
/;
umask 0002;

getopts('o:') || die($usage."\n");
die("$usage\n") unless @ARGV==2 && -f $ARGV[0] && -f $ARGV[1];

my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my @ucsc_tbl=( #from table at http://genome.ucsc.edu/cgi-bin/hgGateway )
'HSCHR6_MHC_APD_CTG1	chr6_apd_hap1',
'HSCHR6_MHC_COX_CTG1	chr6_cox_hap2',
'HSCHR6_MHC_DBB_CTG1	chr6_dbb_hap3',
'HSCHR6_MHC_MANN_CTG1	chr6_mann_hap4',
'HSCHR6_MHC_MCF_CTG1	chr6_mcf_hap5',
'HSCHR6_MHC_QBL_CTG1	chr6_qbl_hap6',
'HSCHR6_MHC_SSTO_CTG1	chr6_ssto_hap7',
'HSCHR4_1_CTG9	chr4_ctg9_hap1',
'HSCHR17_1_CTG5	chr17_ctg5_hap1'
);
my %ctg2chr; # HSCHR alt name => UCSC chr name

foreach (@ucsc_tbl) {
 my @t=split(/\t/);
 $ctg2chr{$t[0]}=$t[1];
}

## -- we want to prepare this hash (or use the last MT line in the file)
my $Mt_acc = 'NC_012920.1';
my %gb2chr = (  $Mt_acc => 'chrM' ); # GFF accession => UCSC name

my %glacc; # GL##### accession => UCSC chr name
open(CHRNAMES, $ARGV[1]) || die("Error opening file $ARGV[1]\n");
while (<CHRNAMES>) {
  next unless m/^chr/;
  my @t=split(/\t/); # $t[0] = UCSC chr name
  if ($t[0]=~m/_([GJ][LH]\d+)/i) {
    my $acc=uc($1);
    $glacc{$acc}=$t[0];
  }
}
close(CHRNAMES);
open(GBTBL, $ARGV[0]) || die("Error opening assembly report file $ARGV[0]");
while(<GBTBL>) {
  next if m/^#/;
  my @t=split(/\t/);
  if ($t[1] eq 'assembled-molecule') {
    if ($t[0] eq 'MT') { #mouse file has Mt chromosome in there
      $t[0]='M';
      delete($gb2chr{$Mt_acc});
    }
    $gb2chr{$t[6]}='chr'.$t[0];
    next;
  }
  if (my $v=$ctg2chr{$t[0]}) {
    $gb2chr{$t[6]}=$v;
    next;
  }
  if ($t[4]=~m/^([GJ][LH]\d+)\.\d+/) {
    my $acc=$1;
    if (my $v=$glacc{$acc}) {
      $gb2chr{$t[6]}=$v;
    }
  }
}
close(GBTBL);

foreach my $acc (sort keys(%gb2chr)) {
 my $chr=$gb2chr{$acc};
 print "$acc\t$chr\n";
}


# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

