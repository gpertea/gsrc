#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 agp_rebuild_seq.pl [-c <contig_name> [-a <output_agp>]] \
   [-o <output.seq>] [-S | -n <seq_name>] <input_agp>[.cidx] <contigs.fasta.fai>
 
 Options:
 -c only build the sequence (and optionally the AGP record if -a) 
    for the given contig name; this option requires the <input_agp>
    parameter to be a cdbfasta index file (.cidx extension)
 -a outputs the AGP record for a single contig in <input_agp>
    (requires -c option)
 -S write separate scaffold sequences (instead of a single pseudomolecule)
    if input AGP has multiple CONTIG comments
 -n unless -S or -c were given, the name <seq_name> is used in the output 
    FASTA sequence (default: pseudomolecule);
    
/;
umask 0002;
getopts('Sn:a:o:c:') || die($usage."Invalid options!\n");
my $outfile=$Getopt::Std::opt_o;
my $outAGP=$Getopt::Std::opt_a;
my $outSeqName=$Getopt::Std::opt_n || 'pseudomolecule';
my $separateCtgSeqs=$Getopt::Std::opt_S;
my $selContig=$Getopt::Std::opt_c;
my %gapname;
@gapname{('clone', 'contig', 'fragment', 'gap')}=();
#die($usage."\n") if ($outAGP && !$selContig);
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my ($agp, $ctgfa)=@ARGV;
die($usage."Cannot find $ctgfa or $agp\n") unless -f $ctgfa && -f $agp;
die("Error: required samtools fasta index (samtools faidx)!\n") 
 unless ($ctgfa=~s/\.fai$//);
if ($selContig) {
die("Error: for -c option, the AGP file must be indexed with cdbfasta\n") 
 unless ($agp=~m/\.cidx$/);
 open(AGP, "cdbyank -a '$selContig' $agp |") 
    || die("Error at cdbyank from $agp\n");
}
else {
 open(AGP, $agp) || die("Error opening $agp\n");
}
my $agpline=1;
if ($outAGP) {
 open(AGPOUT, '>'.$outAGP) || die("Error creating $outAGP\n");
}

my ($outname, $outpos);
$outpos=1;
my $hWritten=0; #at least 1 header written?
while(<AGP>) {
  if (m/^#[\> \t]*(CONTIG:)\s*(\S+)/) {
    $outname=$1.$2;
    #print STDERR "found contig $outname\n";
    if ((!$selContig || $outSeqName ne 'pseudomolecule') && !$separateCtgSeqs) {
      $outname=$outSeqName;
    }
    if ($separateCtgSeqs || !$hWritten) {
      print "\n" if ($outpos>1 && (($outpos-1) % 72 != 0));
      print ">$outname\n";
      $hWritten=1;
      $outpos=1;
      if ($outAGP) {
        print AGPOUT "#>$outname\n";
      }
    }
  next;
  }
  next if m/^#/;
  chomp;
  my $line=$_;
  my ($pname, $pstart, $pstop, $num, $type, $cname, $cstart, $cstop, $cstrand)
      = split(/\t/);
  next unless $cstop && $pstart<$pstop;
  if ($type eq 'N' && $cname=~m/^\d+$/ ) { # && exists($gaptype{$cstart}) {
    next if $outpos==1; #never start a sequence with a gap
    ##They seem to have skipped these two types of gaps :
    next if $cstart eq 'contig' || $cstart eq 'clone'; ## ??
    my $outpos0=$outpos;
    printSeqGap($outname, \$outpos, $cname);
    if ($outAGP) {
      print AGPOUT join("\t", $outname, $outpos0, $outpos-1, 
           $agpline, $type, $cname, $cstart, $cstop, $cstrand)."\n";
      $agpline++;
    }
    next;
  }
  die("Unrecognized AGP line (type=$type, $cstart-$cstop, $cstrand):\n$line\n") 
      unless $type eq 'D' && $cstop>$cstart && 
         ( $cstrand eq '+' || $cstrand eq '-');
  my $samcmd="samtools faidx $ctgfa '$cname:$cstart-$cstop'";
  my $err="Error retrieving subseq with samtools!($samcmd |)";
  open(CSEQ, "$samcmd |") || die("$err\n");
  my $seq;
  while (<CSEQ>) {
    next if m/^>/;
    chomp;
    tr/\t \n\r//d;
    $seq.=$_;
  }
  close(CSEQ);
  my $slen=length($seq);
  my $elen=($cstop-$cstart+1);
  die("$err\n(length mismatch $slen vs $elen!)\n") unless $slen==$elen;
  if ($cstrand eq '-') {
   $seq=reverseComplement($seq);
  }
  my $outpos0=$outpos;
  printSeq($outname, \$outpos, \$seq);
  if ($outAGP) {
    print AGPOUT join("\t", $outname, $outpos0, $outpos-1, 
         $agpline, $type, $cname, $cstart, $cstop, $cstrand)."\n";
    $agpline++;
  }
} #while <AGP>

# --
close(AGPOUT) if $outAGP;
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

sub printSeqGap {
 my ($name, $rp, $len)=@_;
 my $seq='N'x $len;
 printSeq($name, $rp, \$seq);
}

sub printSeq {
 my ($name, $rp, $rseq)=@_;
 #length of  FASTA line: 72
 my $rest=($$rp-1) % 72;
 my $seqlen=length($$rseq);
 my $sp=0;
 my $sleft=$seqlen;
 my $plen = 72 - $rest; 
 while ($sleft>0) {
  $plen=$sleft if ($plen>$sleft);
  print substr($$rseq, $sp, $plen);
  $sleft -= $plen;
  $sp += $plen;
  $$rp += $plen;
  print "\n" if (($$rp-1) % 72 == 0);
  $plen=72;
 }
}
