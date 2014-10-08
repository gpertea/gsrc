#!/usr/bin/perl
use strict;
use Getopt::Std;

my $usage = q/Usage example:
 bowtie2 -a .. | samtophits.pl [max_hits] | ...
 Filters the SAM output of a bowtie2 run with -a option,
 in order to let through only the first max_hits alignments for a read.
 Also discards unmapped reads from the SAM stream.
/;
umask 0002;
getopts('o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $maxHits=shift(@ARGV) || 5;
if ($outfile && $outfile ne '-') {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my ($lread, $lhits);
while (<>) {
  my $line=$_;
  if (m/^@[A-Z][A-Z]\t/) {
    print $_;
    next;
  }
  my @t=split(/\t/);
  next if ($t[2] eq '*' && (int($t[1]) & 4)!=0 );
  if ($lread ne $t[0]) {
    $lread=$t[0];
    $lhits=0;
  }
  ++$lhits;
  if ($lhits<=$maxHits) {
    print $line;
  }
}

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

