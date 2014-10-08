#!/usr/bin/perl
use strict;
my $usage=q/
gff_add_descr.pl <cdbyank_file.cidx> <gfflines..>

Adds an extra "Descr" attribute to the main gff lines 
(mRNA or gene) of the input, by querying the ID
from the <cdbyank_file.cidx> file.

/;
my %cache; # ID => descr
my $cdbidx=shift(@ARGV) || die("$usage\nError: no cidx file given\n");
my $skipped=0;
my $skip;
while (<>) {
 if (m/^\s*#/) { print $_;next; }
 my @t=split('\t');
 next unless($t[8]); #can't be a valid gff3 line
 my $descr;
 my $f=lc($t[2]);
 if ($f eq 'mrna' || $f eq 'gene') {
   if ($t[8]=~m/(?:Descr|Info)\s*=\s*"?([^;"]+)/i) {
      $descr=$1;
      if (length($descr)>6) {
        $skipped++;
        print $_;
        next;
        }
    $t[8]=~s/(?:Descr|Info)\s*=\s*"?[^;"]+//i;
    $t[8]=~tr/;/;/s;
    }#has descr/info already
  # -- retrieve the description here..
  my ($id)=($t[8]=~m/\bID\s*=\s*"?([^;" ]+)/i);
  die("Error: no ID found for gff $f line: $_\n") unless $id;
  $id=~s/\.(\d+)$/.m$1/;
  $t[8]=~s/\bID\s*=\s*"?([^;" ]+)/ID=$id/;
  $id=~s/\.[a-z]*\d+$//;
  $descr=fetchDescr($id);
  chomp($t[8]);
  $t[8].=';descr="'.$descr.'"';
  $_=join("\t",@t)."\n";
  } #mRNA/gene line
 else { # 
  s/(Parent=[\w\|]+)\.(\d+)$/$1.m$2/;
  } 
  print $_;
} #while

sub fetchDescr {
my $id=$_[0];
my $def=$cache{$id};
return $def if $def;
#print STDERR "fetching: cdbyank -a '$id' -F $cdbidx\n";
$def=`cdbyank -a '$id' -F $cdbidx`;
chomp($def);
$def=~s/^\S+\s*//; #remove first token (the ID)
return $def;
}
