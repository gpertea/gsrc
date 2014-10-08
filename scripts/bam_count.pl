#!/usr/bin/perl
use strict;
# usage: samtools view your.bam | bam_count.pl
my $non_uniq; #found with a mapping count of at least 2
my %rh; 
my $list=0; 
if (lc($ARGV[0]) eq '-l') {
  $list=1;
  shift(@ARGV);
  }
while (<>) {
 my @t=split(/\t/);
 my $name=$t[0];
 my $num=0;
 $num+=1 if ($t[1] & 0x40) == 0x40;
 $num+=2 if ($t[1] & 0x80) == 0x80;
 if ($num>2) {
   print STDERR "Warning: $name without mate info:\n$_";
   }
  else { $num = ($num<1) ? '' : '/'.$num ; }
 $non_uniq++ if (++$rh{$name.$num}) > 1;
}

my $mapped=scalar(keys %rh);
if ($list) {
  while (my @hd=each %rh) {
    print join("\t", @hd)."\n";
    }
  }
printf STDERR "Reads mapped: %9d\n", $mapped;
my $uniq=$mapped - $non_uniq;
printf STDERR "Uniq. mapped: %9d\n", $uniq;
