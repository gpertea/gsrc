#!/usr/bin/perl
use strict;
use Getopt::Std;
#use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
  
  cat tab_delim_file | trcol.pl [-c <col#>] [-G] <map_file> 

  Replaces (translates) column <col#> (default 1) from tab delimited
  input with the corresponding string according to the word 
  translation table in <map_file>, which is a simple space delimited,
  two column mapping table of the format:
  
  <string_to_find> <replacement>
  
  Use -G option if the input is GFF3 format and to only output the 
  replaced column and the value of the ID field for the mRNA 
  features only.
  
  
/;
umask 0002;
getopts('Gc:') || die($usage."\n");
my $replcol=$Getopt::Std::opt_c || 1; #column to replace
my $gffmap=$Getopt::Std::opt_G;
$replcol--;
die("$usage\n") if $replcol<0;
my $xtable=$ARGV[0];
die ("$usage\n") unless $xtable && -f $xtable;
my %xw;
open(XTBL, $xtable) || die ("Error opening file $xtable!\n");
while (<XTBL>) {
 chomp; 
 next unless $_;
 my @t=split; 
 $xw{$t[0]}=$t[1];
 }
while (<STDIN>) {
 chomp;
 my @t=split(/\t/); 
 my $tx=$xw{$t[$replcol]};
 $t[$replcol]=$tx if $tx;
 if ($gffmap) {
   next unless m/ID=([^;]+)/;
   my $id=$1;
   $id=~s/\.\w+\d$//;
   print $t[$replcol]."\t$id\n";
   }
  else { 
   print join("\t",@t)."\n";
   }
}
