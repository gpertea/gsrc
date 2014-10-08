#!/usr/bin/perl
use strict;

my $usage = q/Usage:
 eqfasta_check.pl <qryfasta.fa> <tgtfasta.fa.cidx>
 Checks if each sequence in qryfasta.fa is the same with
 the sequence with the same ID in <tgtfasta.fa>
/;
umask 0002;


my ($qfile, $fcdb)=@ARGV;
die("$usage\n") unless -f $qfile && -f $fcdb;

open(QFILE, $qfile) || die ("Error opening $qfile\n");
my ($numdiff, $totalchecked)=(0,0);
my ($qid, $qseq);
while (<QFILE>) {
 if (m/^>(\S+)/) {
    my $newqid=$1;
    checkRec() if $qseq;
    $qid=$newqid;
    $qseq='';
    next;
    }
 chomp;
 $qseq.=$_;
}
checkRec() if $qseq;
close(QFILE);

print STDERR "Checked $totalchecked records, "
                  .int($numdiff)." are different.\n";
sub checkRec {
 my $syscmd="cdbyank -a '$qid' $fcdb";
 my $r=`$syscmd`;
 die("Error at $syscmd") if $? || length($r)<10;
 my ($defline, $tseq)=($r=~m/^>([^\n]+)\n(.+)/s);
 $tseq=~tr/\n//d;
 die("Error: $syscmd returned empty sequence!\n") unless $tseq;
 $totalchecked++;
 if (uc($qseq) ne uc($tseq)) {
   print STDERR "Warning: difference found for $qid\n";
   $numdiff++;
   }
}
