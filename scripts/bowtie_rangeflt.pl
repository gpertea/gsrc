#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
bowtie_rangeflt.pl <chr>[<strand>]:<start_coord>-<end_coord> bowtie_map_data..

Show only those bowtie mappings found in the specified genomic region.
/;
umask 0002;
getopts('o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
die($usage."\n") unless $ARGV[0];
my $range=shift(@ARGV);
my ($chr, $start_coord, $end_coord)=
   ($range=~m/^([^\:]+)\:(\d+)[\.\-]+(\d+)/);
die("$usage\n") unless $start_coord>=1 && $end_coord>$start_coord;
my $strand=$1 if ($chr=~s/([\-\+])$//);
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
while (<>) {
  my $line=$_;
  chomp;
  my @t=split(/\t/);
  my $mlen=length($t[4]);
  next if $strand && $strand ne $t[1];
  next if $chr && $chr ne $t[2];
  next if ($start_coord && $t[3]+$mlen<$start_coord);
  next if ($end_coord && $t[3]>=$end_coord);
  print $line;
  }
# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************


