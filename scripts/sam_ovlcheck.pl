#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 sam_ovlcheck.pl <bam_file> <interval_list> 

 Reports SAM formatted mappings from <bam_file> that truly overlap 
 (not through a gap) any of the intervals given in <interval_list>
 
 <interval_list> has the following format:
 
 <chromosome>[<strand>]:<start1>-<end1>[,<start2>-<end2>,...]
 
 Note:
 <bam_file> must have an index file present (<bam_file>.bai)
/;
umask 0002;
getopts('o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my ($fbam, $intv)=@ARGV;
die("Error: no index file for $fbam?\n") unless -f $fbam.'.bai';
my ($chr,$rlst)=split(/\:/,$intv);
die("$usage Incorrect format for the interval list!\n") unless $chr && $rlst;
my $strand=chop($chr);
if ($strand ne '-' && $strand ne '+') {
  $chr.=$strand;
  $strand=undef;
  }
my @rdata=map { [split(/[\-\.]+/)] } (split(/[\,\;\s]+/,$rlst));
foreach my $d (@rdata) {
  ($$d[0], $$d[1])=($$d[1], $$d[0]) if $$d[0]>$$d[1];
  }
my @ex = sort { $a->[0] <=> $b->[0] } @rdata;
my $range=$chr.':'.$ex[0]->[0].'-'.$ex[-1]->[1];
# -- assumes samtools and bam index 
#my $pipecmd="samtools view $fbam $range |";
#print STDERR $pipecmd."\n";
open(SAMPIPE, "samtools view $fbam $range |") || die ("Error opening samtools pipe ($!)\n");
while(<SAMPIPE>) {
 my $samline=$_;
 chomp;
 my ($qname, $flags, $gseq, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual, @extra)=
    split(/\t/);
 if ($strand) {
    my $mstrand= (($flags & 0x10)==0) ? '+' : '-';
    next if $mstrand ne $strand;
    }
 #now extract the CIGAR segments
 my @cigdata=($cigar=~m/(\d+[A-Z,=])/g);
 my ($mstart, $mend);
 my $hasOvl=0;
 my $curpos=$pos;
 $mstart=$pos;
 foreach my $cd (@cigdata) {
   my $code=chop($cd);
   if ($code eq 'N') { #gap
      #process previous interval
      if ($mend && checkOverlap($mstart, $mend, \@ex)) {
         $hasOvl=1;
         last;
         }
      $mstart=$curpos+$cd;
      $mend=undef;
      next;
      }
   $mend=$curpos+$cd-1;
   $curpos+=$cd;
   }
 unless ($hasOvl) { #check the last interval
   if ($mend && checkOverlap($mstart, $mend, \@ex)) {
       $hasOvl=1;
       }
   }
 if ($hasOvl) {
   print $samline;
   }
 } # while <SAMPIPE>


close(SAMPIPE);

sub checkOverlap {
 my ($a, $b, $rx)=@_;
 return 0 if ($a>$$rx[-1]->[1] || $b<$$rx[0]->[0]); # not overlapping the whole exon chain
 foreach my $x (@$rx) {
   #return (start<=d->end && end>=d->start)
   return 1 if ($a<=$$x[1] && $b>=$$x[0]);
   return 0 if $b<$$x[0];
   }
}


# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

