#!/usr/bin/perl
use strict;
use Getopt::Std;

my $usage = q/Usage:
 jigsaw2gff.pl [-N] [-p <IDbase>] <jigsaw_output>
 
 Converts jigsaw output to gff3 format. The transcript ID
 is by default built from 1st column (chromosome ID) followed by the suffix
 .jsm.<N> (where <N> is the gene number in the jigsaw output).
 
 Options:
  -p use given base name instead of 1st column for transcript ID base name
     and simply append .jsm<N> to it
  -N use a simple .jsm<N> suffix (without the second '.' before the gene #)
/;
umask 0002;
getopts('Np:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $jbase=$Getopt::Std::opt_p;
my $jsuf='.jsm';
$jsuf.='.' unless ($jbase || $Getopt::Std::opt_N);

my $curmodel; #current model name
my ($curtag, $curstrand, $curchr, $genescore); #chromosome and strand for current model
my @exd; #exons for current model

while (<>) {
 next if m/^\s*#/;
 chomp;
 my ($chr, $jsver, $exontype, $exonstart, $exonend, $jscore, 
     $strand, $frame, $lnum)=split(/\t/);
 #next unless $lnum>0;
 my ($modelnum, $gscore);
 if ($lnum=~m/^(\d+)$/) { #old version
   $modelnum=$1;
   }
  else {
   my @a=split(/\s*;\s*/, $lnum);
   ($modelnum)=($a[0]=~m/\.(\d+)$/);
   ($gscore)=($a[1]=~m/gene_score=([\-\d\.]+)/);
   }
 next unless $modelnum;  
($exonstart, $exonend)=($exonend, $exonstart) if $exonend<$exonstart;
 my $base=$jbase || $chr;
 my $locus=$base.$jsuf.$modelnum;
 if ($locus ne $curmodel) {
    &writeModel() if $curmodel;
    $curmodel=$locus;
    $curtag=$chr.$strand;
    $curchr=$chr;
    $curstrand=$strand;
    $genescore=$gscore;
    @exd=([$exonstart, $exonend, $frame]);
    next;
    }
 push(@exd, [$exonstart, $exonend, $frame]);
}

&writeModel() if $curmodel;

sub writeModel {
 my @ex= sort { $main::a->[0] <=> $main::b->[0] } @exd;
 my ($mstart, $mend)=($ex[0]->[0], $ex[-1]->[1]);
 my $attrs="ID=$curmodel;Name=$curmodel";
 if ($genescore) {
    $attrs.=";gene_score=$genescore";
   }
 print STDOUT join("\t",$curchr, 'jigsaw', 'mRNA', $mstart, $mend, '.', 
             $curstrand, '.', $attrs)."\n";
 my $i=1;            
 foreach my $exon (@ex) {
   my ($estart, $eend, $eframe)=@$exon;
   print STDOUT join("\t",$curchr, 'jigsaw', 'CDS', $estart, $eend, '.', 
             $curstrand, $eframe, "Parent=$curmodel")."\n";
   $i++;
   }
 #my @exw = map { $_->[0].'-'.$_->[1] } @ex; 
 #print STDOUT ">$curmodel $mstart $mend $curtag\n";
 #print STDOUT join(',',@exw)."\n";
 $genescore='';
 }
