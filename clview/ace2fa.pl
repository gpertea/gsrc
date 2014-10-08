#!/usr/bin/perl
use strict;
$/="\n\n";
while (<>) {
s/^(\s+)//s;
if (length>20 && s/^CO\s+(\S+)[^\n]+//s) { 
 chomp;
 print ">$1\n";
 tr/-*\n\r//d; # remove gaps and line endings
 # print sequence fasta-formatted:
 print join("\n", (unpack('(A72)*',$_)))."\n";
 }
}
