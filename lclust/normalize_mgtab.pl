#!/usr/bin/perl
#
#--swap qry <-> subj such that they can follow an uniform alphabetical order
use strict;
while (<>) {
 chomp;
 my @t=split(/\t/);
 next unless @t>10;
 if ($t[0] gt $t[4]) {
  #swap:
  if ($t[11] eq '-') {
    #swap coords
    ($t[2],$t[3])=($t[3],$t[2]);
    ($t[6],$t[7])=($t[7],$t[6]);
    }
  @t[0..7]=@t[4..7,0..3];
  ($t[12], $t[13]) = ($t[13], $t[12])
    if (@t>12);
  }
 print join("\t", @t)."\n";
}
