#!/usr/bin/perl
use strict;

my ($fastafile,$isofile,$outbase, $count)=@ARGV;

#my ($fastafile,$isofile,$outfile)=@ARGV;

open(F,$fastafile);
my $line=<F>;
chomp($line);
my ($header,$offset)=split(/\s+/,$line);
if ($offset) {
  die ("Invalid fasta defline for this script: a numeric offset is expected after seqID\n")
   unless $offset=~m/^\d+$/;
  }
my $seq='';
while (<F>) {
 chomp;
 last if m/^>/;
 $seq.=$_;
 }
close(F);

open(F,$isofile);
while(<F>) {
    next if m/^\s*#/;
    chomp;
    my @a=split(/\t/);
    my $beg=$a[3];
    my ($len)=($a[8]=~m/len=(\d+)/);
    my $outfile=$outbase.'.'.$count;
    open(FO,">$outfile") || die ("Error creating $outfile!\n");
    $count++;
    my $subseq=substr($seq,$beg-1,$len);
    print FO "$header ".($offset+$beg-1)."\n";
    &printfa($subseq);
    close(FO);
    

}
#select(STDOUT);
close(F);

sub printfa {
 my $seqlen=length($_[0]);
 my $pos=0;
 while ($pos<$seqlen) {
   print FO substr($_[0],$pos,60)."\n";
   $pos+=60;
   }
 }

