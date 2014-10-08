#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 fna2fq.pl [-U] [-l <minlen>] [-o <out_fname.fq>] <base_fname>

 Options:
  -l : only write reads with length >=minlen (default 42)
  -U : ignore pairing status (write a single output file);
       default is to split the output into <out_fname>_1.fq,
       <out_fname>_2.fq and <out_fname>_u.fq
/;
umask 0002;
my $unread_line=''; # buffer for "unreading" a line in getFnaQual()
getopts('Ul:o:') || die($usage."\n");
my $iname=shift(@ARGV);
die($usage."\n") unless $iname;
my ($inamefa, $infh, $inqh, $inameq);
if ($iname=~m/\.qual$/) {
  $inameq=$iname;
  $inamefa=$iname;
  $inamefa=~s/\.qual$/.fna/;
}
elsif ($iname=~m/\.fna$/) {
  $inamefa=$iname;
  $inameq=$iname;
  $inameq=~s/\.fna$/.qual/;
}
else {
  $inamefa=$iname.'.fna';
  $inameq.=$iname.'.qual';
}

open($infh, $inamefa) || die ("Error opening file $inamefa\n");
open($inqh, $inameq) || die ("Error opening file $inameq\n");

my $outfile=$Getopt::Std::opt_o || 'fna2fq';
my $minlen=$Getopt::Std::opt_l || 42;
my $rawReads=$Getopt::Std::opt_U;
my $fext='fq';
if ($outfile=~s/\.(f[astq]+)$//i) {
  $fext=$1;
}
my @outf=(undef);
if ($rawReads) {
  open($outf[0], ">$outfile.$fext") || die("Error creating output file $outfile.$fext\n");
}
else {
  @outf=(undef,undef,undef);
  open($outf[0], '>'.$outfile."_u.$fext") ||
      die("Error creating output file ${outfile}_u.$fext\n");;
  open($outf[1], '>'.$outfile."_1.$fext") ||
      die("Error creating output file ${outfile}_1.$fext\n");;
  open($outf[2], '>'.$outfile."_2.$fext") ||
      die("Error creating output file ${outfile}_2.$fext\n");;
}

my @rd;

while ( @rd=getFnaQual($infh, $inqh) ) {
  my $fi=0;
  my $rname=$rd[0];
  unless ($rawReads) {
    if ($rd[0]=~s/_left$/\/1/) {
       $fi=1;
      }
    elsif ($rd[0]=~s/_right$/\/2/) {
       $fi=2;
      }
  }
  if ($fi) {
     die("Error: right read found first ?! ($rd[0])\n") if ($fi==2);
     my $rrname=$rname;
     $rrname=~s/_left$/_right/;
     my @rr=getFnaQual($infh, $inqh);
     die("Error: cannot pair left read $rd[0] (found '$rr[0]' instead)\n")
         if ($rr[0] ne $rrname);
     $rr[0]=~s/_right$/\/2/;
     my ($lenl, $lenr)=(length($rd[2]), length($rr[2]));
     if ($lenl<$minlen) {
       printFastq($outf[0], \@rr) if $lenr>=$minlen;
     }
     elsif ($lenr<$minlen) {
       printFastq($outf[0], \@rd) if $lenl>=$minlen;
     }
     else { #both long enough to be written as a pair
       printFastq($outf[1], \@rd);
       printFastq($outf[2], \@rr);
     }
  }
  elsif (length($rd[2])>=$minlen) {
   printFastq($outf[$fi], \@rd);
  }
}



# -- closing output files
foreach my $fh (@outf) {
 close($fh);
 }

#************ Subroutines **************

sub printFastq {
 my ($fh, $d)=@_;
 my $header=$$d[0];
 $header.=' '.$$d[1] if $$d[1];
 print $fh "\@$header\n$$d[2]\n+\n$$d[3]\n";
}

sub getFastaLine {
 my $fh=$_[0];
 if ($unread_line) {
   $_=$unread_line;
   $unread_line='';
   return $_;
 }
 $_=<$fh>;
 return $_;
}

sub getFnaQual { #read 1 (one) full read record from fasta and qual files
 my ($fh, $qh)=@_; # file handle
 #parses next FASTA record
 #returns ($readname, $readseq, $readquals)
 my ($rname, $qrname, $rdefline, $rseq, $rquals);
 local $/="\n";
 while (getFastaLine($fh)) {
   ($rname)=(m/^>(\S+)/);
   last if $rname;
 } #skip empty lines etc.
 return () unless $rname;
 chomp;
 if (length($_)>length($rname)+1) {
  $rdefline=substr($_,length($rname)+1);
 }
 while (<$qh>) {
   ($qrname)=(m/^>(\S+)/);
   last if $qrname;
 } #skip empty lines etc.
 if ($qrname ne $rname) {
  die("Error: fasta seq names mismatch ($qrname vs $rname)\n");
 }
 while (<$fh>) {
     if (m/^>/) {
        $unread_line=$_;
        last;
     }
     chomp;
     $rseq.=$_;
 }
 # read quality values
 while (<$qh>) {
   chomp;
   $rquals.=join('', map { chr($_ + 33) } split() ) ;
   last if (length($rseq)<=length($rquals));
 }
 if ($rseq) {
    return ($rname, $rdefline, $rseq, $rquals);
 }
 else { return () }
}

