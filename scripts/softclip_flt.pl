#!/usr/bin/perl
use strict;
use Getopt::Std;

my $usage = q/Usage:
  samtools view -h file.bam | softclip_flt.pl [-x coord1-coord2] \
  [-g <maxgap>] [-e <max_edit_distance>] [-i <max_indels>] [-s <maxclip>] > file.sam
  
  
  Filter local Bowtie2 alignments by excluding alignments that have
  soft-clipped regions larger than <maxclip> on either side of the read alignment.
  Default <maxclip> value is 10, default <maxgap> value is 2;
  
  If -x option is provided, any clipping which would virtually extend into the region coord1-coord2
  will NOT be taken into account for the filtering.
  
  Note: this script also discards unmapped reads.
/;
umask 0002;
getopts('x:s:o:i:e:g:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $maxclip=$Getopt::Std::opt_s || 10;
my $maxedist=$Getopt::Std::opt_e || 10;
my $maxgap=$Getopt::Std::opt_g || 3;
my $maxindels=$Getopt::Std::opt_i || 2;
my ($xstart, $xend);
my $xreg=$Getopt::Std::opt_x;
if ($xreg) {
  ($xstart, $xend)=($xreg=~m/(\d+)[\.\-]+(\d+)/);
  die("Invalid range specified!\n") if ($xend<=$xstart || $xstart<=0);
}

if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --

while (<>) {
  my $line=$_;
  if (m/^@[A-Z][A-Z]\t/) {
    print $_;
    next;
  }
  my @t=split(/\t/, $_, 12);
  if ($t[2] eq '*' && (int($t[1]) & 4)!=0 ) {
   #unmapped, ignore it
   #print $line if ($clippings);
   next;
  }
 my ($rname, $flags, $pos, $cigarstr, $seq, $quals, $tags)=@t[0,1,4,5,9,10,11];
 my $rlen=length($seq);
 my $isrev = (int($flags) & 0x10) != 0;
 my @cigar=($cigarstr=~m/(\d+[A-Z])/g);
 my $strand='+';
 if ($isrev) {
   @cigar=reverse(@cigar);
   $seq=reverseComplement($seq);
   $quals=reverse($quals);
   $strand='-';
 }
 my $gaptoobig=0;
 my ($numindels, $edist);
 if ($tags=~m/\bNM:i:(\d+)/) {
   next if $1 > $maxedist;
 }
 foreach my $cop (@cigar) {
  my $op=substr($cop, -1);
  if ( $op eq 'D' || $op eq 'I' ) {
    ++$numindels;
    if ($numindels>$maxindels || int($cop)>$maxgap) {
      $gaptoobig=1;
      last;
    }
  }
 }
 next if $gaptoobig;
 #my $fuzz=6; # fuzzy distance requirement from region boundary
 my ($clipL, $clipR);
 $clipL=int($cigar[0]) if (substr($cigar[0],-1) eq 'S');
 $clipR=int($cigar[-1]) if (substr($cigar[-1],-1) eq 'S');
 if ($xstart) {
  if ($pos<$xstart && $xstart-$pos<=$rlen-$clipR+$maxclip) { $clipR=0; } #ignore clipping on the right side
  elsif ($pos+$rlen>$xend && $xend-$pos>=$clipL-$maxclip) { $clipL=0; }
 }
 next if ($clipR>$maxclip || $clipL>$maxclip);
 print $line;
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
