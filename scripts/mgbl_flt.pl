#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 mgbl_flt.pl [-p<minpid>] [-l<minovl>[%]] [-v<maxovh>] < mgblast_hits
/;
umask 0002;
getopts('p:l:v:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $minpid=$Getopt::Std::opt_p || 50.00;
my $minovl=$Getopt::Std::opt_l;
my $ovlperc = ($minovl=~s/\%$//);
my $maxovh=$Getopt::Std::opt_v;

if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
while (<>) {
 my $line=$_;
 my @t=split(/\t/);
 next if $t[8]<$minpid;
 my ($ql, $qr)=$t[2]<$t[3] ? ($t[2], $t[3]) : ($t[3], $t[2]);
 my $qlen=$t[1];
 my $hlen=$t[5];
 if ($minovl) {
   if ($ovlperc) {
     my $minlen =  $qlen<$hlen ? $qlen : $hlen;
     my $movl=int(($minlen*$minovl)/100);
     next if $qr-$ql+1<$movl;
   } else {
      next if $qr-$ql+1<$minovl;
   }
 }
 if ($maxovh) {
   my ($hl, $hr) = $t[6]<$t[7] ? ($t[6], $t[7]) : ($t[7], $t[6]);
   my $l_ovh = $ql>$hl ? $hl-1 : $ql-1;
   my $r_ovh = $qlen-$qr > $hlen-$hr ? $hlen-$hr : $qlen-$qr;
   next if $l_ovh>$maxovh || $r_ovh>$maxovh;
 }
 print $line;
}

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

