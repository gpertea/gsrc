#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 fq_extract.pl [-o <outprefix> ] insert_names.lst reads_1.fq[.gz] [ reads_2.fq[.gz] ]
 Given a list of read names (actually insert names), extracts the 
 [paired] reads from input fastq files.
/;
umask 0002;
getopts('o:') || die($usage."\n");
my $outprefix=$Getopt::Std::opt_o || 'extracted';
my $flst=shift(@ARGV);
die($usage."\n") unless $flst eq '-' || -f $flst;
die($usage."\n") unless @ARGV>0 && @ARGV<3;

my %h; #hash holding the insert names

open(LST, $flst) || die("Error opening file $flst\n");
while (<LST>) {
 chomp;
 my @l=split;
 next unless $l[0];
 $l[0]=~s/\/[12]$//;
 $h{$l[0]}=1;
}
close(LST) unless $flst eq '-';
my $n=keys(%h);
my @f=@ARGV;
die($usage."\nError: list file ($flst) has no entries !\n") if $n==0;
my $paired=(@f==2);
if ($paired) {
   print STDERR "$n read pairs to be extracted.\n";
   }
  else{
   print STDERR "$n reads to be extracted.\n";
   }


my @fh=(undef, undef);
#determine compression, if any
my @fz=('','');
my @fno=($outprefix.'_1.fq', $outprefix.'_2.fq');
$fno[0]=$outprefix.'.fq' unless $paired;
my @foh=(undef, undef);
for (my $i=0;$i<@f;$i++) {
  last unless $f[$i];
  if ($f[$i]=~m/\.bzi?p?2$/) {
   $fz[$i]='bzip2';
   }
   elsif ($f[$i]=~m/.g?zi?p?$/) {
   $fz[$i]='gzip';
   }
  if ($fz[$i]) {
    open($fh[$i], $fz[$i]." -cd '".$f[$i]."'|") || 
        die("Error creating decompression pipe: $fz[0] -cd '$f[$i]' !\n");
    }
   else {
    open($fh[$i], $f[$i]) || die("Error opening file $f[$i] ! \n");
    }
   open($foh[$i], '>'.$fno[$i]) || die("Error opening file $fno[$i] !\n");
}


while (1) {
 my ($rname, $rseq, $rquals)=getFastq($fh[0]);
 last unless $rname;
 my $iname=$rname;
 $iname=~s/\/[12]$//;
 my ($mname, $mseq, $mquals);
 if ($f[1]) {
     ($mname, $mseq, $mquals)=getFastq($fh[1]);
     my $check=$mname;
     $check=~s/\/[12]$//;
     die("Error: cannot find mate for $rname (found $mname instead)!\n") 
          unless $check eq $iname;
     }
 next unless exists($h{$iname});
 print { $foh[0] }'@'.$rname."\n$rseq\n+\n$rquals\n";
 if ($mname) {
     print { $foh[1] } '@'.$mname."\n$mseq\n+\n$mquals\n";
     }
}


for (my $i=0;$i<@f;$i++) {
  last unless $f[$i];
  close($foh[$i]);
  close($fh[$i]);
}

#************ Subroutines **************

sub getFastq {
 my $fh=$_[0]; # file handle
 #parses next FASTQ record
 #returns ($readname, $readseq, $readquals)
 my ($rname, $rseq, $rquals);
 while (<$fh>) {
   ($rname)=(m/^@(\S+)/);
   last if $rname;
   }
 if ($_) {
   while (<$fh>) {
     last if m/^\+/;
     chomp;
     $rseq.=$_;
     }
   if ($_) {
     while (<$fh>) {
       chomp;
       $rquals.=$_;
       last if (length($rseq)<=length($rquals));
       }
   }
  }
 return ($rname, $rseq, $rquals);
}
