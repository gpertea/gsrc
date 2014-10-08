#!/usr/bin/perl

use strict;
my ($file,$addcoord)=@ARGV;

open(F,$file) || die ("Error opening file $file\n");
my $p='';
if ($file=~m/\.(\d+)$/) {
  $p=$1;
  # -- only needed for correct naming of gff3 entitied 
  # produced by gmap ran on separate parts

  }
while(<F>) {
    chomp;
    if (substr($_,0,1) eq "#") { print $_."\n";}
    else {
        if ($p) {
          s/\.cds(\d+)/\.cds$p\_$1/g;
          s/\.path(\d+)/\.p$p\_$1/g;
          }
	my @a=split;
        $a[3]+=$addcoord;
        $a[4]+=$addcoord;
	print join("\t", @a)."\n";
	}
}
close(F);
