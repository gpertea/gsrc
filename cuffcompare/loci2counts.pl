#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
  loci2counts.pl <cuffcmp.loci> [<cuffcmp.combined.gtf>]

  Outputs a column file:

  <locus_location> <TCONS_count> <sample1_tcount> <sample2_tcount> ...
/;
umask 0002;
getopts('o:') || die($usage."\n");
die ($usage) unless $ARGV[0] && -f $ARGV[0];
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my $lfile=shift(@ARGV);
my %loc; # locusID => [ location, tcons_count, tcount1, tcount2, ...]
my @loclist;
open(LFILE, $lfile) || die ("Error opening $lfile\n");
while (<LFILE>) {
 chomp;
 my @t=split(/\t/);
 my @tcounts;
 for (my $l=3;$l<@t;$l++) {
   $t[$l]=~tr/\- //d;
   if ($t[$l]) {
      my @lst=split(/\,/,$t[$l]);
      push(@tcounts, scalar(@lst));
      }
   else { push(@tcounts, 0); }
   }
 push(@loclist, $t[0]);
 $loc{$t[0]}=[$t[1], 0, @tcounts ];
 }
close(LFILE);
my %tids; #keeping track of transcript_ids
while (<>) {
 chomp;
 my ($g)=(m/gene_id\s+"([^"]+)/);
 my ($t)=(m/transcript_id\s+"([^"]+)/);
 if ($g && $t && not exists($tids{$t})) {
    $tids{$t}=1;
    my $d=$loc{$g};
    die("Error: locus $g not found in loci file!\n") 
      unless $d;
    $d->[1]++;
    }
 }
#print results
foreach my $l (@loclist) {
 my $d=$loc{$l};
 print join("\t",@$d)."\n";
 }
# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************


