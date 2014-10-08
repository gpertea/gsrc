#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 gettop.pl [-f <col#>] [-d <delim>] [-t N] sorted_lines..
 
 Display only the top N (default 1) lines with the same string 
 in column <col#> (default 1), where columns are delimited
 by <delim> (default: tab character).
 ATTENTION: the input lines must be already sorted by column <col#>
 AND the criterion of choice.
/;
umask 0002;
getopts('f:d:t:') || die($usage."\n");

my $fld=$Getopt::Std::opt_f || 1;
$fld--;
my $delim=$Getopt::Std::opt_d || "\t";
my $top=$Getopt::Std::opt_t || 1;

my ($prev, $count);
while (<>) {
 chomp;
 my @f=split($delim);
 if ($f[$fld] eq $prev) { $count++ }
                   else { $count=1 }
  print $_."\n" if $count<=$top;
  $prev=$f[$fld];
 }
