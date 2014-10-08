#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 samtools view <in.bam> | sam2reads.pl 
 sam2reads.pl <in.sam>

 Rebuild the raw read names (for paired reads) from the input 
 SAM\/BAM file and list them.
/;
umask 0002;
getopts('o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --

while (<>) {
 my @t=split(/\t/);
 my $name=$t[0];
 my $num=0;
 $num+=1 if ($t[1] & 0x40) == 0x40;
 $num+=2 if ($t[1] & 0x80) == 0x80;
 if ($num<1 || $num>2) {
   print STDERR "Warning: $name without mate info:\n$_";
   }
 print $name.'/'.$num."\n";
}

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************
