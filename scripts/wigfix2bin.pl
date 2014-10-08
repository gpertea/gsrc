#!/usr/bin/perl
use strict;
my $usage=q/Usage:
 wigfix2bin.pl <wigFixstream>
 Converts wigFix to a fixed binary format:
 
 bytes 00-08 :  'WIGB'<32bit-start-offset>
 bytes 09... :  2-byte int = value*1000 for each base
 
/;
$_=<STDIN>;
die("Error: wigFix header not recognized!\n") 
  unless (m/^fixedStep\s/);
my ($chr, $cpos)=(m/chrom=(\S+)\s+start=(\d+)/); 
die("Error: wigFix header not recognized!\n") unless $chr && $cpos;

open(BIN, ">$chr.wigfixbin") || die("Error creating file $chr.wigfixbin!\n");
binmode(BIN);
print STDERR "Writing $chr.wigfixbin ..\n";
print BIN "WIGB";
print BIN pack('I', $cpos);
while (<STDIN>) {
 if (m/start=(\d+)/) {
    my $skipto=$1;
    my $skiplen=$skipto-$cpos;
    for (my $i=0;$i<$skiplen;$i++) {
       print BIN pack('s',0);
       }
    $cpos=$skipto;
    next;
    }
 my ($v)=(m/^([\d\.\-]+)$/);
 die("Error parsing value for line $. : $_\n") unless length($v)>0;
 my $n=int($v*1000);
 if ($n>32767 || $n<-32767) {
    die("Error: value out of range ($n)!\n");
    }
 print BIN pack('s',$n);
 $cpos++;
 }
close(BIN);
print STDERR "    ..done.\n";
#print STDERR "Waiting for 4s..\n";
#sleep(4);
#print STDERR "Now reading data back..\n";

# read back (testing):
# open(BIN, "$chr.wigfixbin");
# my ($tag, $r);
# read(BIN, $tag, 4);
# read(BIN, $r, 4);
# my @v=unpack('I',$r);
# print STDERR "Header tag=$tag, offset=$v[0]\n";
# while (read(BIN,$r,2)==2) {
#   @v=unpack('s',$r);
#   my $n=sprintf('%.3f', $v[0]/1000);
#   print STDERR "read value: $n\n";
#   }
# close(BIN);
