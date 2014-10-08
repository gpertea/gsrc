#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 pull_pairs.pl [-o <out_prefix>] <mate_1.fq> <mate_2.fq> \
     <mapped.sam>[ <mapped_2.sam>]
 
 Create files <outprefix_1.fq> and <outprefix_2.fq> with 
 paired reads that have at least one mapped mate found in the 
 input .sam files
/;
umask 0002;
getopts('o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --



# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

