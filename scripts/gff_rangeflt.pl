#!/usr/bin/perl
use strict;
use Getopt::Std;
getopts('T') || die($usage."\n");

my $usage=q/
 Usage:
 gff_rangeflt.pl [-T] <chr>[<strand>]:<start_coord>-<end_coord> <input_gff..>

 Use the -T option to not check specifically for exon overlap, but any overlap
 of the transcript range (including intron-only overlaps)
/;

my $r=shift(@ARGV) || die("$usage\n");
die("$usage\n") if $r=~m/^\-/;
my ($chr, $cstart, $cend)=($r=~m/^([^:]+):(\d+)[-\.]+(\d+)/);
($cstart, $cend)=($cend,$cstart) if $cstart>$cend;
die("$usage\n") unless $cstart>1 && $cend>$cstart;
my $cstrand=chop($chr);
unless ($cstrand eq '-' || $cstrand eq '+') {
  $chr.=$cstrand; #restore, it wasn't a strand there
  $cstrand='';
  }
my @gffbuf; #buffer with current transcript lines
my $prevtid;
my ($pstart, $pend);
while (<>) {
 my $curline=$_;
 my @t=split(/\t/);
 next unless $chr eq $t[0];
 next if $cstrand && $cstrand ne $t[6];
 my $tid;
 my $gffparent;
 if ($t[8]=~m/transcript_id\s+([^;]+)/) {
   $tid=$1;
   }
  else {
   if ($t[8]=~m/Parent=([^;]+)/) { # GFF3 format, exon/element
    # assumes these are given in order
    $tid=$1;
    }
   elsif ($t[8]=~m/ID=([^;]+)/)  { # GFF3, new parent
    $tid=$1;
    $gffparent=1;
    #($pstart, $pend)=($t[3],$t[4]);
    }
   }
 if ($prevtid && $prevtid ne $tid) {
      #transcript change
      unless ($cstart>$pend || $pstart>$cend) {
         # range overlap
         print join('',@gffbuf);
         }
    @gffbuf=();
    ($pstart, $pend) = $gffparent ? ($t[3],$t[4]) : (0,0);
    }
  $prevtid=$tid;
  push(@gffbuf, $curline);
  ($t[3], $t[4])=($t[4], $t[3]) if $t[3]>$t[4];
  $pstart=$t[3] if ($t[3]<$pstart || $pstart==0);
  $pend=$t[4] if $t[4]>$pend;
  }
 # check for last transcript read  
if ($pstart && @gffbuf>0 && $cstart<=$pend && $pstart<=$cend) {
   print join('',@gffbuf);
   }
