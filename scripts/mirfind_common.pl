#!/usr/bin/perl
use strict;
use Getopt::Std;

my $usage = q/Usage:
 mirfind_common.pl [-t total1:total2] <file1.mirfind> <file2.mirfind>
 Takes output files of mirfind_score.pl and reports potential miRNAs
 that have the same location of the mature sequence;
 A coordinate variation of 1 for either ends is allowed.
 
 (input files are assumed sorted by region start coordinate, 
  as they are when they are created by mirfind_score.pl)
 
/;
umask 0002;
getopts('t:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
die ($usage."\n") unless @ARGV==2;
my @alines;
my @blines;
my $totals=$Getopt::Std::opt_t;
my @total=split(/[\:\,\/]/,$totals);
my $ratio=$total[1]/$total[0] if $total[0]>0;
open(FILA, $ARGV[0]) || die ("Error opening file $ARGV[0]!\n");
open(FILB, $ARGV[1]) || die ("Error opening file $ARGV[1]!\n");
while (<FILA>) {
 chomp;
 my @t=split(/\t/);
 push(@alines, [@t,0]);
 }
close(FILA);
@alines = sort { ($$a[1] eq $$b[1]) ? $$a[3]<=>$$b[3] : $$a[1] cmp $$b[1] } @alines;
while (<FILB>) {
 chomp;
 my @t=split(/\t/);
 push(@blines, [@t, 0]);
 }
close(FILB);
@blines = sort { ($$a[1] eq $$b[1]) ? $$a[3]<=>$$b[3] : $$a[1] cmp $$b[1] } @blines;

my ($a_reg, $a_start,$a_chr, $a_strand, $a_end, $a_mstart, $a_mend, $acov,
    $b_reg, $b_start,$b_chr, $b_strand, $b_end, $b_mstart, $b_mend, $bcov);
my $anext;
my @ta; #ta[13] has the intron info, if any
my @tb;
my ($ad, $ai, $bd, $bi);
$bi=-1;
for ($ai=0; $ai<@alines; $ai++) {
  $ad=$alines[$ai];
  @ta=@$ad;
  if ($a_start && $a_chr eq $ta[1] && $a_start>$ta[3]) {
     die("Wrong order for $ARGV[0] lines?! ($a_reg vs $ta[0])!\n");
     }
  #($a_start, $a_end)=($ta[0]=~m/\|(\d+)\-(\d+)$/);
  ($a_start, $a_end)=($ta[3],$ta[4]);
  ($a_reg, $a_chr, $a_strand)=@ta[0..2];
  ($a_mstart, $a_mend)=($ta[8]=~m/m\|(\d+)\-(\d+)/);
  ($acov)=($ta[6]=~m/(\d+)$/);
  $ad->[6]=$acov;
  $ta[6]=$acov;
  die ("Error parsing mstar/mend!\n") unless $a_mend>$a_mstart && $a_mstart>1;
  goto AB_CMP if ($anext);
 NEXT_B:
  $bi++;
  last unless $bi<@blines;
  $bd=$blines[$bi];
  @tb=@$bd;
  if ($b_start && $b_chr eq $tb[1] && $b_start>$tb[3]) {
     die("Wrong order for $ARGV[1] lines?! ($b_reg vs $tb[0])!\n");
     }
  ($b_start, $b_end)=($tb[3],$tb[4]);
  ($b_reg, $b_chr, $b_strand)=@tb[0..2];
  ($b_mstart, $b_mend)=($tb[8]=~m/m\|(\d+)\-(\d+)/);
  ($bcov)=($tb[6]=~m/(\d+)$/);
  $bd->[6]=$bcov;
  $tb[6]=$bcov;
 AB_CMP:
  $anext=0;
  if ($b_chr gt $a_chr) { # a lags behind
    $anext=1;
    next;
    }
  elsif ($a_chr gt $b_chr) { # b lags behind
    goto NEXT_B;
    }
  # same chromosome:  
  if ($b_start>$a_end) { # a lags behind
    $anext=1;
    next;
    }
  if ($a_start>$b_end) { # b lags behind
    goto NEXT_B;
    }
 # - hairpin overlap here, check for mature miRNA overlap
 if (abs($a_mstart-$b_mstart)<5 && abs($a_mend-$b_mend)<5) { # && 
     #abs($a_start-$b_start)<7 && abs($a_end-$b_end)<7) {
    # common miRNA prediction here
    #print join("\t", $ta[0], @ta[5..8], join('|',@ta[10..12]), 
    #                 $tb[0],  @tb[5..8], join('|',@tb[10..12]), $ta[14])."\n";
    my ($hstart, $hend, $mseq, $mcoords) = ($ta[6]>$tb[6]) ? (@ta[3..4], @ta[7..8]) : (@tb[3..4], @tb[7..8]);
    my ($mstart, $mend)=($mcoords=~m/m\|(\d+)\-(\d+)/);
    my $foldchanges='';
    if ($ratio) {
      my ($v1,$v2)=($ta[6]*$ratio, $tb[6]);
      my $fc=($v1>$v2)? -($v1/$v2) : $v2/$v1;
      $foldchanges="\t".sprintf('%.2f',$fc)."\t".sprintf('%.2f',abs($fc));
      
      }
    print join("\t", join('|',@ta[1..2], $hstart.'_'.$hend), $mseq, $mstart, $mend, @ta[5..6], $ta[10].'|'.$ta[11], $ta[12], 
                               @tb[5..6], $tb[10].'|'.$tb[11], $tb[12]).$foldchanges."\n";
    
    $bd->[15]=1;
    $ad->[15]=1;
    }
 if ($b_end>$a_end) {
    $anext=1;
    next;
    }
  else {
    goto NEXT_B;
    }
}

# -- now also add those entries only expressed in one sample but not the other
foreach my $d (@alines) {
 next if $d->[15]==1 || $d->[6]<15;
 my @t=@$d;
 my ($hstart, $hend, $mseq, $mcoords) = (@t[3..4], @t[7..8]);
 my ($mstart, $mend)=($mcoords=~m/m\|(\d+)\-(\d+)/);
 my $foldchanges='';
 if ($ratio) {
      my ($v1,$v2)=($t[6]*$ratio, 1);
      my $fc=($v1>$v2)? -($v1/$v2) : $v2/$v1;
      $foldchanges="\t".sprintf('%.2f',$fc)."\t".sprintf('%.2f',abs($fc));
      }

 print join("\t", join('|',@t[1..2], $hstart.'_'.$hend), $mseq, $mstart, $mend, @t[5..6], $t[10].'|'.$t[11], $t[12],
                               '-','-', '-', '-' ).$foldchanges."\n";

 #print join("\t", $t[0], $t[5], 'x'.$t[6], @t[7..8], join('|',@t[10..12]), 
 #                 '-',  '-','-','-','-', '-', $t[14])."\n";
 }

foreach my $d (@blines) {
 next if $d->[15]==1 || $d->[6]<15;
 my @t=@$d;
 my ($hstart, $hend, $mseq, $mcoords) = (@t[3..4], @t[7..8]);
 my ($mstart, $mend)=($mcoords=~m/m\|(\d+)\-(\d+)/);
 my $foldchanges='';
 if ($ratio) {
      my ($v1,$v2)=(1*$ratio, $t[6]);
      my $fc=($v1>$v2)? -($v1/$v2) : $v2/$v1;
      $foldchanges="\t".sprintf('%.2f',$fc)."\t".sprintf('%.2f',abs($fc));
      }
 print join("\t", join('|',@t[1..2], $hstart.'_'.$hend), $mseq, $mstart, $mend, '-','-', '-', '-', 
                               @t[5..6], $t[10].'|'.$t[11], $t[12]).$foldchanges."\n";
 #print join("\t", $t[0], '-','-','-','-', '-', 
 #                $t[0],  $t[5], 'x'.$t[6], @t[7..8], join('|',@t[10..12]), $t[14])."\n";
 }

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

