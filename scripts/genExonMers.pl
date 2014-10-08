#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
  genExonMers.pl [-k <len>] [-X] <fasta_file> <id_tmap>
  
  Generates a multi-fasta file with all exon-contained k-mers
  for each sequence in <fasta_file>, with header augmented with
  data taken from <id_tmap>
  
  Use -X option to allow probes to cross exon bounderies
  (default is: only generate probes fully contained within exons)
  
  Use -M option when only k-mers from multi-match genes 
  should be generated (requires <id_tmap>)
/;
umask 0002;
getopts('XMk:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $klen=$Getopt::Std::opt_k || 50;
my $crossexons=$Getopt::Std::opt_X;
my $multimatch=$Getopt::Std::opt_M;

#================ SUBROUTINES ============
die($usage."\n") unless @ARGV>=1;
my ($fa, $tmap)=@ARGV;
if ($tmap) {
 foreach my $f ($fa, $tmap) {
   die "Error: cannot find file '$f'!\n" unless -f $f;
   die "Error: invalid file size for '$f'!\n" unless -s $f;
   }
 }

my %tm; # cufflinks_id => ref_id, but also ref_id => cufflinks_id
my %hn; # count how many times a cluster pair was found
#        

if ($tmap) {
open(TMAP, $tmap) || die("Error opening file '$tmap'!\n");
while (<TMAP>) {
 chomp;
 my @t=split("\t");
 next unless $t[2];
 $tm{$t[0]}=$t[2];
 $tm{$t[2]}=$t[0];
 if ($t[8]) {
    my $clpair=$t[7].'_'.$t[8];
    $hn{$clpair}++;
    $hn{$t[2]}=$clpair;
    }
 }
close(TMAP);
}
open(FA, $fa) || die("Error opening file '$fa'!\n");
my $c_id;
my @c_exons;
my @c_segs;
my $c_seq;
my $c_loc;
while (<FA>) {
 chomp;
 if (m/^>(\S+)/) {
   my $id=$1;
   genProbes();
   if ($tmap && $multimatch) {
     my $clpair=$hn{$id};
     $c_id='';
     next if $hn{$clpair}<2;
     }
   $c_id=$id;
   ($c_loc)=(m/ loc:(\S+)/);
   my ($exons)=(m/ exons:(\S+)/);
   my ($segs)=(m/ segs:(\S+)/);
   die ("Error parsing line ($_), exons or segs not found!\n") unless $exons && $segs;   
   @c_exons=($exons=~m/(\d+\-\d+)/g);
   @c_exons=map { [split('-')] } @c_exons;
   @c_segs=($segs=~m/(\d+\-\d+)/g);
   @c_segs=map { [split('-')] } @c_segs;
   $c_seq='';
   next;
   }
 $c_seq.=uc($_) if $c_id;
}
close(FA);
genProbes(); # uses c_* global variables

#-----------------

sub genProbes {
 return unless $c_id;
 my ($chr, $chr_pos, $chr_strand)=split('\|',$c_loc);
 my $exnum=0;
 my $excount = @c_segs;
 my $clen=length($c_seq);
 die("Error at $c_id ($c_loc) - exon count mismatch!\n") if ($excount != @c_exons);
 foreach my $seg (@c_segs) {
   my ($start, $end)=@$seg;
   $start--;
   my $seglen=$end-$start;
   $exnum++;
   next if $seglen<$klen && !$crossexons;
   my $r_id=$tm{$c_id} || $c_id;
   my ($xnum, $chrpos);
   if ($chr_strand eq '-') {
      $xnum=$excount-$exnum;
      my $exon=$c_exons[$xnum];
      $chrpos=$$exon[1];
      }
     else {
      $xnum=$exnum-1;
      my $exon=$c_exons[$xnum];
      $chrpos=$$exon[0];
      }
   $xnum++;
   if ($crossexons) {
    for (my $p=$start;$p<$end && $p+$klen<=$clen;$p++) {
        my $xpos=($chr_strand eq '-') ? $chrpos-$klen-($p-$start)+1 : $chrpos+($p-$start);
        #my $exn=$exnum;
        my $exn=$xnum;
        $exn.='_' if ($p+$klen>$end); # crosses into next exon
        my $rstart=$p+1;
        print ">$c_id|$r_id|$chr|$chr_strand|$excount|$exn|$xpos|$rstart\n";
        print substr($c_seq, $p, $klen)."\n";
       }
     }
   else {
   for (my $p=$start;$p+$klen<=$end;$p++) {
       my $xpos=($chr_strand eq '-') ? $chrpos-$klen-($p-$start)+1 : $chrpos+($p-$start);
       my $rstart=$p+1;
       #print ">$c_id|$r_id|$chr|$chr_strand|$excount|$exnum|$xpos|$rstart\n";
       print ">$c_id|$r_id|$chr|$chr_strand|$excount|$xnum|$xpos|$rstart\n";
       print substr($c_seq, $p, $klen)."\n";
      }
    }
   }
 $c_seq='';
 @c_exons=();
 @c_segs=();
 $c_loc='';
 $c_id='';
}

