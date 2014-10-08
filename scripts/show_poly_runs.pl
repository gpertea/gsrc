#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 show_poly_runs.pl [ -p <minrun>] [ -G | -N ] fasta_file..

Shows location and length of same-base runs - more than <minrun> of 
the same base in a row.
Options:
   -p  : minimum base repeat to consider (default 10)
   -N  : also show poly N runs (gaps)
   -G  : (Gaps-only) show only poly N runs (gaps)

For each poly run, the following tab-delimited columns are shown:
 1. sequence name
 2. base
 3. sequence location of the poly-base run
 4. length of poly-base run
/;
umask 0002;
getopts('NGp:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $gapsOnly=$Getopt::Std::opt_G;
my $gapsAlso=$Getopt::Std::opt_N;
die("{$usage}Error: cannot use both -G and -N\n") if $gapsOnly && $gapsAlso;
my $pLen=$Getopt::Std::opt_p || 10;
my $xL=int($pLen-1);
my @PChars=$gapsOnly? ('N') : ('A','C','G','T');
push(@PChars, 'N') if $gapsAlso;
die("$usage\nValue for -p must be >=1 !\n") unless $xL>0;
if ($xL<1 && !$gapsOnly) {
 die("Error: -p1 is only accepted with -G option\n");
}
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my $seqname;
my $seq;
my $showNs;
while (<>) {
  if (m/^>(\S+)/) {
    my $sn=$1;
    procSeq();
    $seqname=$sn;
  } else {
   tr/\n\r\t //d;
   $seq.=uc($_);
  }
}
procSeq();

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#---------------------
sub procSeq {
  return unless $seqname && $seq;
  #print "[DBG]:>SEQ>$seqname\t".length($seq)."\n";
  my @polys; #array of [$poly, $p_start, $p_end] for each base in @PChars (A,C,G,T,N)
  foreach my $ch (@PChars) {
    my $cseq=$seq;
    my $p=0;
    while (length($cseq)>=$pLen) {
       if ($cseq=~m/(${ch}{$pLen,})/) {
        my ($poly, $p_start, $p_end)=($1, $-[0], $+[0]);
        $p+=$p_start;
        my $pl=$p_end-$p_start;
        #               0     1     2
        push(@polys, [ $ch, $p+1,  $pl]);
        $cseq=substr($cseq, $p_end);
        $p+=$pl;
       }
       else { last; }
    } #while matching a poly
  } #for each poly target base
  # sort poly matches by position
  if (@polys>1) {
    @polys = sort { $$a[1] <=> $$b[1] } @polys;
  }
  foreach my $pd (@polys) {
    #print "[DBG]:$seqname:poly match=$poly($p_start-$p_end), seqlen=".length($seq)."\n";
    my ($ch, $pos, $len)=@$pd;
    print join("\t",$seqname, $ch, $pos, $len)."\n";
  }
  $seq='';
  $seqname='';
}

#
#sub procSeq
#{
    #my $p=0;
    #while (length($cseq)>=$pLen) {
    ##### Perl 5.8.8 bug - crashes on the self-referential regexes below #####
    #my ($poly, $p_start, $p_end);
    #if ($gapsOnly) {
    #  if ($seq=~m/(N{$pLen,})/) { 
    #     ($poly, $p_start, $p_end)=($1, $-[0], $+[0]);
    #  }
    #  else { last; }
    #}
    #else { #A,C,G,T and perhaps N
    #   if ($gapsAlso) {
    #     if ($seq=~m/([ACGTN])\1{$xL,}/) {
    #       ($poly, $p_start, $p_end)=(substr($seq,$-[0], $+[0]-$-[0]), $-[0], $+[0]);
    #     }
    #   }
    #   else {
    #     if ($seq=~m/([ACGT])\1{$xL,}/) {
    #       ($poly, $p_start, $p_end)=(substr($seq,$-[0], $+[0]-$-[0]), $-[0], $+[0]);
    #     }
    #   }
    #  last unless $poly;
    #}
    #    #print "[DBG]:$seqname:poly match=$poly($p_start-$p_end), seqlen=".length($seq)."\n";
    #my $polylen=length($poly);
    #my $b=substr($poly,0,1);
    #$p+=$p_start;
    #print join("\t",$seqname, $b, $p+1, $polylen)."\n";
    #$seq=substr($seq, $p_end);
    #$p+=$polylen;
    ##### end of original code which segfaults on 5.8.8
#}
