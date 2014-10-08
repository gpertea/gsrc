#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q{Usage:
 orf_check.pl [-P] [-a <min_aa_len>] <gffread_mrnas.fasta>
 Use -N option to not check for ORFs, but just print the output
 without the last column.
};
umask 0002;
getopts('Pa:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $minaa=$Getopt::Std::opt_a || 30;
my $justPrint=$Getopt::Std::opt_P;
my $minbases=$minaa*3;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
#-- globals
my %codons=(
'AAA'=>'K', 'AAC'=>'N', 'AAG'=>'K', 'AAR'=>'K', 'AAT'=>'N',
'AAY'=>'N', 'ACA'=>'T', 'ACB'=>'T', 'ACC'=>'T', 'ACD'=>'T',
'ACG'=>'T', 'ACH'=>'T', 'ACK'=>'T', 'ACM'=>'T', 'ACN'=>'T',
'ACR'=>'T', 'ACS'=>'T', 'ACT'=>'T', 'ACV'=>'T', 'ACW'=>'T',
'ACY'=>'T', 'AGA'=>'R', 'AGC'=>'S', 'AGG'=>'R', 'AGR'=>'R',
'AGT'=>'S', 'AGY'=>'S', 'ATA'=>'I', 'ATC'=>'I', 'ATG'=>'M',
'ATH'=>'I', 'ATM'=>'I', 'ATT'=>'I', 'ATW'=>'I', 'ATY'=>'I',
'CAA'=>'Q', 'CAC'=>'H', 'CAG'=>'Q', 'CAR'=>'Q', 'CAT'=>'H',
'CAY'=>'H', 'CCA'=>'P', 'CCB'=>'P', 'CCC'=>'P', 'CCD'=>'P',
'CCG'=>'P', 'CCH'=>'P', 'CCK'=>'P', 'CCM'=>'P', 'CCN'=>'P',
'CCR'=>'P', 'CCS'=>'P', 'CCT'=>'P', 'CCV'=>'P', 'CCW'=>'P',
'CCY'=>'P', 'CGA'=>'R', 'CGB'=>'R', 'CGC'=>'R', 'CGD'=>'R',
'CGG'=>'R', 'CGH'=>'R', 'CGK'=>'R', 'CGM'=>'R', 'CGN'=>'R',
'CGR'=>'R', 'CGS'=>'R', 'CGT'=>'R', 'CGV'=>'R', 'CGW'=>'R',
'CGY'=>'R', 'CTA'=>'L', 'CTB'=>'L', 'CTC'=>'L', 'CTD'=>'L',
'CTG'=>'L', 'CTH'=>'L', 'CTK'=>'L', 'CTM'=>'L', 'CTN'=>'L',
'CTR'=>'L', 'CTS'=>'L', 'CTT'=>'L', 'CTV'=>'L', 'CTW'=>'L',
'CTY'=>'L', 'GAA'=>'E', 'GAC'=>'D', 'GAG'=>'E', 'GAR'=>'E',
'GAT'=>'D', 'GAY'=>'D', 'GCA'=>'A', 'GCB'=>'A', 'GCC'=>'A',
'GCD'=>'A', 'GCG'=>'A', 'GCH'=>'A', 'GCK'=>'A', 'GCM'=>'A',
'GCN'=>'A', 'GCR'=>'A', 'GCS'=>'A', 'GCT'=>'A', 'GCV'=>'A',
'GCW'=>'A', 'GCY'=>'A', 'GGA'=>'G', 'GGB'=>'G', 'GGC'=>'G',
'GGD'=>'G', 'GGG'=>'G', 'GGH'=>'G', 'GGK'=>'G', 'GGM'=>'G',
'GGN'=>'G', 'GGR'=>'G', 'GGS'=>'G', 'GGT'=>'G', 'GGV'=>'G',
'GGW'=>'G', 'GGY'=>'G', 'GTA'=>'V', 'GTB'=>'V', 'GTC'=>'V',
'GTD'=>'V', 'GTG'=>'V', 'GTH'=>'V', 'GTK'=>'V', 'GTM'=>'V',
'GTN'=>'V', 'GTR'=>'V', 'GTS'=>'V', 'GTT'=>'V', 'GTV'=>'V',
'GTW'=>'V', 'GTY'=>'V', 'MGA'=>'R', 'MGG'=>'R', 'MGR'=>'R',
'NNN'=>'X', 'RAY'=>'B', 'SAR'=>'Z', 'TAA'=>'.', 'TAC'=>'Y',
'TAG'=>'.', 'TAR'=>'.', 'TAT'=>'Y', 'TAY'=>'Y', 'TCA'=>'S',
'TCB'=>'S', 'TCC'=>'S', 'TCD'=>'S', 'TCG'=>'S', 'TCH'=>'S',
'TCK'=>'S', 'TCM'=>'S', 'TCN'=>'S', 'TCR'=>'S', 'TCS'=>'S',
'TCT'=>'S', 'TCV'=>'S', 'TCW'=>'S', 'TCY'=>'S', 'TGA'=>'.',
'TGC'=>'C', 'TGG'=>'W', 'TGT'=>'C', 'TGY'=>'C', 'TRA'=>'.',
'TTA'=>'L', 'TTC'=>'F', 'TTG'=>'L', 'TTR'=>'L', 'TTT'=>'F',
'TTY'=>'F', 'XXX'=>'X', 'YTA'=>'L', 'YTG'=>'L', 'YTR'=>'L'
);

# --
my ($t, $xloc, $chr, $strand, @exons, @segs, $seq);
while (<>) {
 chomp;
 next unless $_;
 if (m/^>(\S+)/) { #defline
  my ($id, $l)=($1,$_);
  process_seq();
  $t=$id;
  ($xloc)=($l=~m/\bgene=(\w+)/);
  my ($loc)=($l=~m/\bloc:(\S+)/);
  if ($loc) {
    ($chr,$strand)=($loc=~m/(\S+)\|\d+\-\d+\|([\-\+])$/);
    }
  my ($ex)=($l=~m/\bexons:([\d\-\,]+)/);
  if ($ex) {
    @exons= map { [(split(/\-/))] } (split(/\,/,$ex));
    }
  my ($sg)=($l=~m/\bsegs:([\d\-\,]+)/);
  if ($sg) {
    @segs= map { [(split(/\-/))] } (split(/\,/,$sg));
    }
  next;
  }
 #sequence
 $seq.=$_;
}

process_seq();

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************
sub process_seq {
 return unless $seq;
 #print $t.' '.$chr.$strand.' '.join(',', @x)."\n";
 $seq=uc($seq);
 my @orfs; #list of [start_on_mrna, orflen, start_on_genome, stop_on_genome]
 my $pstart=0;
 my $i=-1;
 my ($lasti, @lastend);
 goto SKIPORF if ($justPrint);
 while (($i=index($seq,'ATG',$pstart))>=0) {
   last if $i+3>length($seq)-$minaa*3;
   if ($i<$lastend[ $i % 3 ]) {
      $pstart=$i+3;
      next;
      }
   my $aa=trFrom($seq,$i);
   my $blen=length($aa)*3;
   if ($blen>$minbases) {
     $lasti=$i;
     $lastend[ $i % 3 ]=$i+$blen;
     if ($aa=~m/>trunc$/) { #just discard those without stop codons, or with no 3'UTR 
         $pstart=$i+3;
         next;
         }
     my ($gstart, $gstop)=(pos2genome($i+1), pos2genome($i+$blen));
     push(@orfs, [$i+1,$blen, $gstart, $gstop]);
     #print " [====> orf aa: (ntlen=$blen, at ".($i+1)."..".($i+$blen).")\n  $aa\n";
     }
   $pstart=$i+3;
   }
 SKIPORF:
 if (@orfs>0 || $justPrint) {
   my @o;
   if (@orfs>0) {
      @orfs=sort { $main::b->[1]<=>$main::a->[1] } @orfs;
      @o=map { $_->[0].':'.$_->[1] } @orfs;
      }
   my @x=map { $_->[0].'-'.$_->[1] } @exons;
   print join("\t", $t,$xloc,$chr.$strand, 'exons='.join(',', @x),'orfs='.join(',',@o))."\n";
   #print "     >> max orf has ntlen=$orfs[0]->[1] at $orfs[0]->[0]\n";
   }
 # - end of record processing: clear globals 
 ($t, $xloc, $chr, $strand, $seq)=();
 @exons=();
 @segs=();
}

sub trFrom { #translate seq from offset, returns translation until stop codon
 # ($seq, $ofs) given
 my $ofs=$_[1];
 my $s=substr($_[0], $ofs);
 my @cods = unpack('(A3)*',$s);
 my $r;
 my $i=1;
 foreach my $c (@cods) {
  if (length($c)<3) { # premature ending, this is not a valid ORF
    return $r.'>trunc'; #doesn't have a stop codon
    }
  my $aa=$codons{$c} || 'X';
  last if $aa eq '.';
  $r.=$aa;
  $i++;
  if ($i==@cods) {
     return $r.'>trunc';
     }
  }
 return $r;
}

# map position of a base from transcript to genome
# requires exon coordinates and strand info
sub pos2genome {
 if($strand eq '-') {
  
  }
 # -- positive strand
 
 return 0;
 }
