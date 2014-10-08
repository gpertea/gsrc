#!/usr/bin/perl
use strict;
use Getopt::Std;

my $usage = qq/
mirfind_excise.pl [-o <out_prefix>] <chromosomes.fa.cidx> <bowtie_sorted.map>

Excises potential microRNA precursor sequences from a genome 
using the positions of aligned reads as guidelines.

<chromosomes.fa.cidx> must be the cdbyank index for the multifasta file 
with all chromosomes.
The input file <bowtie_sorted.map> MUST BE SORTED
by chromosome AND then start position. 

Two output files will be created:
  1) <out_prefix>.xc.fa : fasta files with excised regions
  2) <out_prefix>.xc.xmap : a pseudo-fasta file (suitable for indexing) 
                           holding read mappings on to every excised region
                          (with coordinates relative to the excised sequence)

If -o option was not given, <out_prefix> will be set to the input file name.
/;

# -- Note: in miRDeep's script the coordinates in the resulting fasta file
#  are relative to the reverse complemented strand for '-' alignments
# --
# this is NOT the way they are written here: coordinates of the excised 
# ranges are always relative to the forward strand, like in a GFF file


getopts('o:') || die "$usage\n";
my $fgenomecidx=shift(@ARGV) or die $usage;
my $file_map=shift(@ARGV) or die $usage;
die($usage."\n") unless $fgenomecidx=~m/\.cidx$/ && -f $fgenomecidx;
my $fgenome=$fgenomecidx;
$fgenome=~s/\.cidx$//;
die($usage."Error: cannot locate fasta file $fgenome\n") unless -f $fgenome;
die($usage."Error: cannot locate mapping file $file_map\n") unless -f $file_map;

my $outprefix=$Getopt::Std::opt_o;
unless ($outprefix) {
 $outprefix=$file_map;
 $outprefix=~s/\.(\w+)$//;
 }
my ($outfa, $outmap)=($outprefix.'.xc.fa', $outprefix.'.xc.xmap');
open(OUTFA, '>'.$outfa) || die("Error creating file $outfa!\n");
select(OUTFA);
open(OUTMAP, '>'.$outmap) || die("Error creating file $outmap!\n");

my ($gseqname, $gseq, $gseqlen, $gseqcounter);
excise_by_map($file_map);
select(STDOUT);
close(OUTFA);
close(OUTMAP);
exit;
# ------------------- 

#sub parse_file_blast_parsed{

sub excise_by_map { #processes bowtie map output
  my($file)=@_;
  open(FILENAME, $file) or die "Could not open file $file";
  
  my (@r_start, @r_end); 
  #current region start and end coordinates of last region - by strand (0='+', 1='-')
  my @xmap=([],[]); #read mappings data on the selected region, by strand (0='+', 1='-')
  while (<FILENAME>) {
    chomp;
    my @t=split(/\t/);
    next unless (@t>4) && ($t[3]>0);
    my $query=$t[0];
    my $strand=$t[1];
    my $qlen=length($t[4]);
    my $chr=$t[2];
    my $sidx=($strand eq '-')? 1 : 0;
    my $g_start=$t[3]+1;
    my $g_end=$g_start+$qlen-1;
    my $newgseq=($chr ne $gseqname);
    if ($newgseq) {
       foreach my $i (0..1) {
         if ($r_end[$i]) {
            my $rs=($i==0)?'+':'-';
            gExcise($rs, $r_start[$i], $r_end[$i], $xmap[$i]); # assumes $gseqname was set properly
            ($r_start[$i], $r_end[$i])=(0,0);
            }
         }
       @xmap=([],[]);
       loadGSeq($chr); #this will set ($gseqname, $gseq, $gseqlen)
       }
      #elsif ($r_end[$sidx] && $g_start-$r_end[$sidx]>30) {
      elsif ($r_end[$sidx] && $g_start-$r_end[$sidx]>22) {
       # same genomic sequence, but separated by a gap
       gExcise($strand, $r_start[$sidx], $r_end[$sidx], $xmap[$sidx]);
       ($r_start[$sidx], $r_end[$sidx])=(0,0);
       $xmap[$sidx]=[];
       }
    push(@{$xmap[$sidx]}, [$query, $g_start, $g_end, $t[7], $t[8]]);
    $r_start[$sidx]=$g_start if ($r_start[$sidx]==0 || $r_start[$sidx]>$g_start);
    $r_end[$sidx]=$g_end if $r_end[$sidx]<$g_end;
    } # for each bowtie mapping
  
  #process last regions:
  foreach my $i (0..1) {
      if ($r_end[$i]) {
          my $rs=($i==0)?'+':'-';
          gExcise($rs, $r_start[$i], $r_end[$i], $xmap[$i]); # assumes $gseqname was set properly
          ($r_start[$i], $r_end[$i])=(0,0);
          $xmap[$i]=[];
          }
      }
  close(FILENAME);
}

sub loadGSeq {
  my ($subject) = @_;
  return if $gseqname eq $subject && $gseqlen>10;
  print STDERR "processing $subject.. \n";
  # load genomic sequence in $gseq and update ($gseq, $gseqname, $gseqlen)
  my $fofs=`cdbyank -a '$subject' -P $fgenomecidx`;
  chomp($fofs);
  die("Error: cannot cdbyank $subject from $fgenomecidx!\n") unless length($fofs)>0;
  open(FASTA, $fgenome) or die "Error opening $fgenome\n";
  binmode(FASTA);
  seek(FASTA, $fofs, 0);
  my $l=<FASTA>;
  if ($l=~m/^>(\S+)/) {
    my $new=$1;
    die("Error: $new sequence name found, expected $subject\n") unless $new eq $subject;
    $gseqname=$new;
    }
   else {
    die("Error: no FASTA header in $fgenome at offset $fofs (qseq=$subject)\n");
    }
  #now read the sequence
  $gseqlen=0;
  $gseq='';
  while (<FASTA>) {
    last if m/^>/;
    chomp;
    $gseq.=$_;
    $gseqlen+=length($_);
    }
  close(FASTA);
  $gseqcounter=0;
}


sub gExcise{
  my ($gstrand, $gstart, $gend, $xmaps)=@_;
  #print STDERR "gexcise: $gstart..$gend\n";
  return unless $gend>$gstart;
  if ($gend-$gstart>60) {
     #longer alignments should be excised as a single potential precursor, shorter as two
     my $ext=136-($gend-$gstart);
     if ($ext<0) { $ext=0; }
            else {  $ext>>=1; }
     #print STDERR "gexcise ext=$ext\n";
     writeExcision(\$gseqcounter, $gstrand, $gstart-$ext, $gend+$ext, $xmaps);
     }
  else { #short region (<=30)
     writeExcision(\$gseqcounter, $gstrand, $gstart-22, $gend+84, $xmaps);
     writeExcision(\$gseqcounter, $gstrand, $gstart-84, $gend+22, $xmaps);
    }
  }


sub writeExcision {
 my ($k, $strand, $start, $end, $xm)=@_;
 $start=1 if $start<1;
 $end=$gseqlen if $end>$gseqlen;
 my $len=$end-$start+1;
 my $subseq=uc(substr($gseq, $start-1, $len));
 $subseq=revcom($subseq) if $strand eq '-';
 print ">$gseqname\_$$k $strand|$start-$end\n$subseq\n";
 print OUTMAP ">$gseqname\_$$k $strand|$start-$end\n";
 foreach my $m (@$xm) {
   my ($mstart, $mend)=($$m[1]-$start+1, $$m[2]-$start+1);
   # reverse complement coordinates if needed:
   ($mstart, $mend)=($len-$mend+1, $len-$mstart+1) if ($strand eq '-');
   my @wm=@$m;
   @wm[1..2]=($mstart, $mend);
   print OUTMAP join("\t",@wm)."\n";
   }
 $$k++;
 }

sub revcom{
    my $seq=reverse($_[0]);
    $seq=~tr/acgtuACGTU/TGCAATGCAA/;
    return $seq;
}

