#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 rnaseq_mapflt.pl [-x <chromosome>] [-l <seedlen>] [-n <max_seed_mismatches>] \
                    [-t <max_terminal_mismatches>] <bowtie_map_data>
 
 Filters bowtie map output (sorted with best hits first) such that only
 mappings with <max_mismatches> at 3' end (after <seedlen>) are 
 allowed.
 
 Use option -x to exclude mappings from a specific chromosome or contig
 
 Defaults are : -l16 -n1 -t3 -xchrM
 (exclude mitochondrial chromosome)

 For bowtie output created with -n and --best --strata it will also make sure 
 that only the "best strata" is kept based on the total number of mismatches.
 
/;
umask 0002;
getopts('x:t:l:n:o:') || die($usage."\n");
my $seedlen=$Getopt::Std::opt_l||16;
my $maxmmseed=$Getopt::Std::opt_n || 1;
my $maxmmend=$Getopt::Std::opt_t || 3;
my $nochr=$Getopt::Std::opt_x || 'chrM';
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my ($prevq, $prevmmcount, $prevmmseedcount);
while (<>) {
 my $line=$_;
 chomp;
 next unless $_;
 my ($qry, $strand, $chr, $cstart, $seq, $qv, $d, $mism)=split(/\t/);
 #next if $chr eq $nochr;
 my $len=length($seq);
 next unless $len>8;
 my @mm=map { [split(/\:/)] } (split(/\,/,$mism)) ; 
 #mismatches as a list of [mmcoord, nt_change]
 if ($strand eq '-') { #for reverse complement, reverse and adjust 
   @mm=map { [$len-$$_[0]-1, $$_[1]] } reverse @mm;
   }
 my @mmseed=grep { $$_[0]<$seedlen } @mm;
 next if @mmseed>$maxmmseed;
 #my @mmbefore_end=grep { $$_[0]<$len-$maxmmend-1 } @mm;
 #next if @mmbefore_end>0;
 next if scalar(@mm)-scalar(@mmseed)>$maxmmend;
 next if ($prevq eq $qry && @mm>$prevmmcount);
 $prevq=$qry;
 $prevmmcount=scalar(@mm);
 $prevmmseedcount=scalar(@mmseed);
 #TODO: trim 5' end if any mismatches are grouped there close to that end 
 my $tlen=$len;
 #print STDERR ">> testing $line";
 foreach my $md (reverse @mm) {
     #print STDERR "  .. testing $$md[0] vs $tlen-2\n";
     if ($$md[0]<$tlen-2) { last; }
         else { $tlen=$$md[0]; } #truncate at mismatch
     }
 if ($tlen<$len) {
   #adjust mismatches and read/qual strings
   my $dlen=$len-$tlen;
   @mm=grep { $$_[0]<$tlen } @mm;
   if ($strand eq '-') {
     @mm = map { [$len-$$_[0]-1, $$_[1]] } reverse @mm;
     $seq=substr($seq,$dlen);
     $qv=substr($qv, $dlen);
     $cstart+=$dlen;
     #print STDERR "$qry\n" if $tlen>=$seedlen;
     }
   else { #truncate sequence
     $seq=substr($seq,0,$tlen);
     $qv=substr($qv,0,$tlen);
     #print STDERR "$qry (+)\n" if @mm>0 && $tlen>=$seedlen;
     }
   next if $tlen<$seedlen;
   my @mmstr = map { $$_[0].':'.$$_[1] } @mm;
   $line=join("\t", $qry, $strand, $chr, $cstart, $seq, $qv, $d, join(',',@mmstr))."\n";
   }
 print $line unless $chr eq $nochr;
 }

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

