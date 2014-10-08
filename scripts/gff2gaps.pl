#!/usr/bin/perl
use strict;
my $usage=q/
gff2gaps.pl <gff_file> <gseq_w_gaps.fa>

It will output a fasta file with the same content in <gseq_w_gaps.fa> 
but with exon nucleotides changed to E and introns to I
Only the first set of exons in the GFF file are considered..
/;

my $gffile=shift(@ARGV);
die("$usage\n") unless $gffile;
my @ex;
my @cd;
open(GFF, $gffile) || die ("Error opening $gffile!\n");
while (<GFF>) {
 next if m/^#/;
 my ($chr, $mt, $ftype, $fstart, $fend, $rest)=split(/\t/, $_, 6);
 if ($ftype eq 'exon' || $ftype eq 'CDS') {
   ($fstart, $fend)=($fend, $fstart) if ($fend<$fstart);
  
   if ($ftype eq 'exon') { push(@ex, [$fstart, $fend]) }
     else { push(@cd, [$fstart, $fend]) }
   }
 else {
  last if @ex>0 || @cd>0; #assumes all exons/CDSs are grouped
  }
}
close(GFF);
my @se = sort { $main::a->[0] <=> $main::b->[0] } @ex;
my @sc = sort { $main::a->[0] <=> $main::b->[0] } @cd;

#print STDERR scalar(@se)." exons found:\n";
#foreach my $e (@se) { print STDERR "  $$e[0]..$$e[1] length: ".($$e[1]-$$e[0]+1)."\n"; }
#exit;
my $ei=0;
my $ci=0;
my $p=0;
my $ap=0;
my $isExon=0;
my $isCoding=0;
my $prevState=0;
while (<>) {
if (m/^>/) {
  s/^>(\S+)/>$1\|refmsk/;
  print $_;
  next;
  }
#sequence
tr/\n\r \t//d;
my $l=$_;
for (my $i=0;$i<length($l);$i++) {
 my $c=substr($l,$i,1);
 my $ich='.';
 $ap++; #allignment position
 goto IREST if ( $c eq '-');
 $isExon=0;
 $isCoding=0;
 $p++;
 goto IREST if ($ei >= @se); #beyond last exon
 if (@sc>0) {
   if ($p<=$sc[$ci]->[1]) { # before cur cdseg end
     $isCoding=1 if ($p>=$sc[$ci]->[0]);
     }
   else { # beyond cur cdseg end
     $ci++;
     }
   } # cds defined
 if ($p<=$se[$ei]->[1]) { # before cur exon end
     $isExon=1 if ($p>=$se[$ei]->[0]);
     }
   else { # beyond cur exon end
     $ei++;
     }
 $ich=$c if ($ei>0 && ($p-$se[$ei-1]->[1]<3 || $se[$ei]->[0]-$p<3));
 IREST:
 if ($prevState != $isExon) {
   my ($sep, $xp) = $prevState==1 ? ("\n",$ap-1):(' - ',$ap);
   print STDERR $xp.$sep;
   $prevState=$isExon;
   }
 substr($l, $i, 1)=$isExon ? ($isCoding?'C':'T') : $ich;
 }
print $l."\n"; 
}
