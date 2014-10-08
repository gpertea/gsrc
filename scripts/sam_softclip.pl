#!/usr/bin/perl
use strict;
use Getopt::Std;

my $usage = q/Usage:
  samtools view file.bam | sam_softclip.pl [-l <minlen>] [-S][-A] > file.fq
  
  Unless -S option is given, this script discards the soft-clipped regions
  from the input alignments and outputs a FASTQ file with reads trimmed 
  to only the regions that aligned to the reference  (a FASTA format file 
  will be written if -A option was given).
  Default <minlen> is 24 (if the remaining read sequence is shorter than
  this, the read will be discarded)

  However if the -S option is given the output will consist of just the 
  soft-clipped regions of the reads (which could be useful, say, for adapter 
  analysis), in plain fasta format
/;
umask 0002;
getopts('SAl:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $minlen=$Getopt::Std::opt_l || 24;
my $fasta=$Getopt::Std::opt_A;
my $clippings=$Getopt::Std::opt_S;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --

while (<>) {
  my $line=$_;
  if (m/^@[A-Z][A-Z]\t/) {
    #print $_;
    next;
  }
  my @t=split(/\t/);
  if ($t[2] eq '*' && (int($t[1]) & 4)!=0 ) {
   #print $line if ($clippings);
   next;
  }
 my $isrev = (int($t[1]) & 0x10) != 0;
 my ($rname, $cigarstr, $seq, $quals)=@t[0,5,9,10];
 my @cigar=($cigarstr=~m/(\d+[A-Z])/g);
 my $sign='+';
 if ($isrev) {
   @cigar=reverse(@cigar);
   $seq=reverseComplement($seq);
   $quals=reverse($quals);
   $sign='-';
 }
 my ($clipL, $clipR);
 $clipL=int($cigar[0]) if (substr($cigar[0],-1) eq 'S');
 $clipR=int($cigar[-1]) if (substr($cigar[-1],-1) eq 'S');
 if ($clippings) { #only print clippings
  if ($clipL+$clipR>5) {
    print '>'.$rname.' '.$sign.$cigarstr;
    if ($clipL>$clipR) {
      print " $clipL\n".substr($seq, 0, $clipL)."\n";
    }
    else {
      print " -$clipR\n".substr($seq, -$clipR)."\n";
    }
  }
 }
 else {
  if ($clipL) { substr($seq, 0, $clipL)=''; substr($quals, 0, $clipL)=''; }
  if ($clipR) { substr($seq, -$clipR)=''; substr($quals, -$clipR)='';}
  if (length($seq)>=$minlen) {
    if ($fasta) {
     print '>'.$rname.' '.$sign.$cigarstr."\n";
     print "$seq\n";
    }
    else {
     print '@'.$rname.' '.$sign.$cigarstr."\n";
     print "$seq\n+\n$quals\n";
    }
  }
 }
}

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

sub reverseComplement {
  my $s=reverse($_[0]);
  $s =~ tr/AaCcTtGgUuMmRrWwSsYyKkVvHhDdBb/TtGgAaCcAaKkYyWwSsRrMmBbDdHhVv/;
  return $s;
 }
