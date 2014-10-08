#!/usr/bin/perl
use strict;
#clustalw file conversion back to multi-fasta
my @sn; #sequence names
my @seqs;
my @block;
while (<>) {
 chomp;
 if (substr($_,0,2) eq '  ' || length($_)<2) {
   #end block encountered
   if (@block<2) { 
      @block=();
      next;
      }
   my $first=(@seqs==0);
   for (my $i=0;$i<@block;$i++) {
     my ($id, $seq)=split(/\s+/, $block[$i],2);
     $seq=~tr/-*//d;
     if ($first) {
       push(@sn,  $id);
       push(@seqs,$seq);
       }
     else {
       die("Error: seq id mismatch ($sn[$i] expected):\n$_\n") 
           unless $id eq $sn[$i];
       $seqs[$i].=$seq;
       }
   }
   @block=();
   next;
 } #block ended
 push(@block, $_);
}

for (my $i=0;$i<@sn;$i++) {
 print ">$sn[$i]\n";
 print join("\n", unpack('(A70)*', $seqs[$i]))."\n";
}
