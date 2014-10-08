#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;
umask 0002;
my $usage=q/
 Usage: 
  twinscan2gff.pl [-p <transcript_suffix>] < twinscan_output.gtf > gff3
/;
getopts('p:') || die($usage."\n");
my $prefix=$Getopt::Std::opt_p || "winscan";
my $curgene; #current gene name
my $curtran; #current transcript
my ($curstrand, $curchr); #chromosome and strand for current model
my @exd; #exons for current model
my %stops; #transcripts with stop_codon info transcript_id=>stop_start

while (<>) {
 next if m/^\s*#/;
 chomp;
 next unless $_;
 my ($chr, $prog, $ftype, $fstart, $fend, $fscore, $strand, $frame, $info)=
  split(/\t/);
 my ($gid)=($info=~m/gene_id\s+\"([\w\.\:\;\'\|\#\-]+)\"/);
 my ($tid)=($info=~m/transcript_id\s+\"([\w\.\:\;\'\|\#\-]+)\"/);
 die ("Error parsing gene/transcript id from $_\n") unless $gid && $tid;
 if ($ftype eq 'stop_codon') {
   $stops{$tid}=($strand eq '-')? $fend : $fstart;
   next;
   }
 next unless $ftype eq 'CDS';
 
 if ($tid ne $curtran) {
   &writeTranscript() if $curtran;
   $curtran=$tid;
   $curchr=$chr;
   $curstrand=$strand;
   @exd=();
   }
  push(@exd, [$fstart, $fend, $frame]);
}

&writeTranscript() if $curtran;


sub writeTranscript {
 my @ex= sort { $main::a->[0] <=> $main::b->[0] } @exd;
 my ($mstart, $mend)=($ex[0]->[0], $ex[-1]->[1]);
 if (my $stpos=$stops{$curtran}) {
  if ($curstrand eq '-') {
    die ("Invalid stop codon position for $curtran! ($stpos vs $mstart)\n")
     unless $stpos==$mstart-1;
    $mstart-=3; #include stop codon (twinscan doesn't..)
    $ex[0]->[0]-=3;
    }
   else {
    die ("Invalid stop codon position for $curtran! ($stpos vs $mend)\n")
     unless $stpos==$mend+1;
    $mend+=3; #include stop codon! 
    $ex[-1]->[1]+=3;
    }
  }
 my ($tnum, $anum)=($curtran=~m/\.(\d+)\.(\d+)$/);
 my $geneid='t'.$prefix.$tnum;
 $curtran=$geneid."tws$anum";
 print join("\t",$curchr, 'twinscan', 'mRNA', $mstart, $mend, '.', 
             $curstrand, '.', "ID=$curtran;Name=$geneid")."\n";
 my $i=1;
 foreach my $exon (@ex) {
   my ($estart, $eend, $eframe)=@$exon;
   print join("\t",$curchr, 'twinscan', 'CDS', $estart, $eend, '.', 
             $curstrand, $eframe, 
             "Parent=$curtran")."\n";
   $i++;
   }
}
