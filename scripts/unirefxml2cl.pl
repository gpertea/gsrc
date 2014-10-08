#!/usr/bin/perl
use strict;
my $uref; #current uniref ID
my $mcount; #member count
my @mlist; #member list

while (<>) {
 if (m/<entry id="([\-\w]+)"/) {
   my $id=$1;
   if ($uref) {
     die("Error: invalid number of members for $uref ($mcount vs ".
          scalar(@mlist).")!\n")
        if $mcount!=scalar(@mlist);
     print STDOUT ">$uref ".join(' ', @mlist)."\n";
     }
   $uref=$id;
   $mcount=0;
   @mlist=();
   undef(@mlist);
   next;
   }
 if (m/<property type="member count"\s+value="(\d+)"/) {
   $mcount=$1;
   next;
   }
 if (m/<dbReference type="[\w ]+" id="([-\w]+)"/) {
   my $mid=$1;
   push(@mlist, $mid);
   next;
   }
 }
 
if ($uref) {
  die("Error: invalid number of members for $uref!\n")
     if ($mcount!=@mlist);
  print STDOUT ">$uref ".join(' ', @mlist)."\n";
  }
