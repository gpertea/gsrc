#!/usr/bin/perl
use strict;
while (<>) {
 next unless m/\tmRNA\t/;
 my @t=split("\t");
 my ($id)=(m/ID=([\|\w\.]+)/);
 $id=~s/\.[a-z]+\d+$//;
 my ($score)=(m/score=([\d\.]+)/i);
 my ($cov)=(m/Cov=([\d\.]+)/i);
 my $pid=$t[5];
 my $covpid=int(($pid*$cov)/100);
 #my ($geneid)=(m/(GeneID:\d+)/);
 print join("\t", $id, $score, $pid, $cov, $covpid, # $geneid,
            $t[0], $t[3], $t[4], $t[6])."\n";
}
