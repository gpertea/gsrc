#!/usr/bin/perl
use strict;
#use Getopt::Std;
my $usage = q/Usage:
 mgbl_probes.pl < mgblast_D4.output > self-hit-only.lst
/;

# mgblast2 -d mm9_exons.fa -i test_probes.fa -X1 -FF -D4 -W 44 > mgblast.tab

umask 0002;
#getopts('o:') || die($usage."\n");
#my $outlst=$Getopt::Std::opt_o;
my %h; # =>[has_self-hit, has_other-hit];
while (<>) {
 chomp;
 my @t=split("\t");
 next unless $t[2]==1 && $t[3]==$t[1]; # && $t[2]==50 && int($t[7])==100
 my $qid=$t[0];
 my $strand=$t[11];
 my $hid='|'.$t[4].'|';
 my $d=$h{$qid};
 unless ($d) {
   $d=[0,0];
   $h{$qid}=$d;
   }
 if (index($qid,$hid)>0) {
      #print STDERR "looking for $hid in '$qid'\n";
      #$htype='self_hit';
      #$h{$qid}[0]++;
      $d->[0]++;
      }
      else {
      #$htype='other_hit'; 
      #$h{$qid}[1]++;
      $d->[1]++;
      }
   #print join("\t",$qid, $strand, @h, $htype)."\n";
}

  foreach my $q ( keys(%h) ) {
   my $d=$h{$q}; #print if it has a self hit but not other hits
   print "$q\n" if $$d[0]>0 && $$d[1]==0;
   }

# if ($outlst) {
#  open(OUTLST, '>'.$outlst) || die ("Error creating file $outlst\n");
#  foreach my $q ( keys(%h) ) {
#   my $d=$h{$q}; #print if it has a self hit but not other hits
#   print OUTLST "$q\t$$d[0]\n" if $$d[0]>0 && $$d[1]==0;
#   }
#  close(OUTLST);
# }
