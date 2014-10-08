#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
  ctg2Soap.pl [-o ctg2scaf.txt] in_ctgs.fasta out_prefix
  
  out_prefix is the prefix of the output files from SOAPDenovo2
  in_ctgs.fasta is the input contigs given to prepare (finalFusion)
  
/;
umask 0002;
getopts('o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my ($ctgfa, $outprefix)=@ARGV;
die("$usage Error: input contig file $ctgfa not found!\n") unless -f $ctgfa;
system("samtools faidx $ctgfa")==0 
   or die("Error running: samtools faidx $ctgfa \n");
system("samtools faidx $outprefix.contig")==0 
   or die("Error running samtools faidx $outprefix.contig !\n");
my %ictgs; # len => [ ctgid, ...]
my @duplens; # (len1, len2, ...)
my %ictgsb; # ctgid => first-60-bases (only for same length ctgids)
my %octgs; # len => [ outctgid, ...]
my %octgsb;# outctgid => first-60-bases
my @soapctgs; #list of soap contigs
my %inscaff; #list of soap contigs linked into scaffolds
open (INCTG, "$ctgfa.fai") || die ("Error opening $ctgfa.fai\n");
while(<INCTG>) {
  next if m/^#/;
  chomp;
  my ($ctg, $clen, $cofs, $cl1, $cl2)=split;
  if (exists($ictgs{$clen})) {
    my @ids=@{$ictgs{$clen}};
    push(@duplens, $clen) if @ids==1;
    push(@{$ictgs{$clen}}, $ctg);
  } else {
    $ictgs{$clen}=[$ctg];
  }
}
close(INCTG);
#if we have duplicate lengths, get a fingerprint for those contigs
foreach my $clen (@duplens) {
  my @ids=@{$ictgs{$clen}};
  my $tail=$clen-59;
  foreach my $ctg (@ids) {
    my $cmd="samtools faidx $ctgfa $ctg:1-60 $ctg:$tail-$clen";
    my $cseq=`$cmd`;
    $cseq=~s/>[^\n]+\n//gs;
    $cseq=~tr/\n\r \t//d;
    #$ictgsb{$ctg}=$clen.'|'.uc($cseq);
    $ictgsb{$clen.'|'.uc($cseq)}=$ctg;
  }
}

# now process SOAPCTG
open (OCTG, "$outprefix.contig.fai") ||
   die ("Error opening $outprefix.contig.fai\n");
while(<OCTG>) {
  next if m/^#/;
  chomp;
  my ($ctg, $clen, $cofs, $cl1, $cl2)=split;
  next unless $cl2>10 && $cl2>=$cl1;
  push(@soapctgs, $ctg);
  if (exists($octgs{$clen})) {
    my @ids=@{$octgs{$clen}};
    push(@duplens, $clen) if @ids==1;
    push(@{$octgs{$clen}}, $ctg);
  } else {
    $octgs{$clen}=[$ctg];
  }
}
close(OCTG);
#if we have duplicate lengths, get a fingerprint for those contigs
foreach my $clen (@duplens) {
  my @ids=@{$octgs{$clen}};
  my $tail=$clen-59;
  foreach my $ctg (@ids) {
    my $cmd="samtools faidx $outprefix.contig $ctg:1-60 $ctg:$tail-$clen";
    my $cseq=`$cmd`;
    die("Error at: $cmd\n") unless length($cseq)>120;
    $cseq=~s/>[^\n]+\n//gs;
    $cseq=~tr/\n\r \t//d;
    $octgsb{$ctg}=$clen.'|'.uc($cseq);
    #$octgsb{$clen.'|'.uc($cseq)}=$ctg;
  }
}

my %s2c; #soap_ID => in_ctg_name
foreach my $clen (keys(%octgs)) {
 my @ctgs=@{$octgs{$clen}};
 if (@ctgs>1) {
   foreach my $ctg (@ctgs) {
      my $ictg=$ictgsb{$octgsb{$ctg}};
      unless ($ictg) {
       die("Error getting in-ctg for out ctg $ctg\n");
       #my $osb=$octgsb{$ctg};
       #print STDERR "octgsb[$ctg] is '$osb'\n";
       #print STDERR "ictgsb[$osb] is '".$ictgsb{$osb}."'\n";
       #die("\n");
      }
      $s2c{$ctg}= $ictg;
   }
 }
 else  {
   $s2c{$ctgs[0]} = $ictgs{$clen}->[0];
 }
}

#parse the scaffold definitions
open(SCAF, "$outprefix.contigPosInscaff") 
  || die("Error opening $outprefix.contigPosInscaff\n");
my $scaf; #current scaffold
my @scafctgs;
while(<SCAF>) {
 if (m/^>(\S+)/) {
   my $s=$1;
   print $scaf."\t".join("\t", @scafctgs)."\n" if $scaf;
   $scaf=$s;@scafctgs=();
   next;
 }
 my ($sctg, $p1, $strand, $p2)=split;
 #$strand='' if ($strand eq '+');
 #push(@scafctgs, $strand.$s2c{$sctg});
 push(@scafctgs, $s2c{$sctg}.' '.$strand.'['.$p1.':'.$p2.']');
 $inscaff{$sctg}=1;
}
print $scaf."\t".join("\t", @scafctgs)."\n" if $scaf;
close(SCAF);
#print the rest of the ctgs which did not make it into scaffolds
foreach my $sctg (@soapctgs) {
  next if exists($inscaff{$sctg});
  print 'C'.$sctg."\t".$s2c{$sctg}."\n";
}

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

