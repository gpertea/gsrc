#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 col_calc.pl [-f <col#>] <tab_delimited_data>..

 Compute sum and average values for the specified numeric column
 of the input data.
/;
umask 0002;
getopts('f:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $f=$Getopt::Std::opt_f || 1;
$f-- if $f>0;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my ($n,$t);
while(<>) {
 chomp;
 next unless $_;
 my @d=split(/\t/);
 unless ($d[$f]=~m/^\s*[\d\.\-\+]+\s*$/) {
   print STDERR "Warning: ignored non-numeric field: $d[$f]\n";
   next;
   }
 $t+=$d[$f];
 $n++;
}
print "Sum of $n values is $t, average=".($t/$n)."\n";

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

