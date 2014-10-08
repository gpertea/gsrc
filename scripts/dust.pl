#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 dust.pl [-c <repeat#>] <fasta_file>
 
 Masks all repeats of unit length 1 or greater that are repeated at 
 least 4 times.
 Options:
 -c    sets the repeat count to <repeat#> (default 4)
/;
umask 0002;
getopts('c:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $c=$Getopt::Std::opt_c || 4;

{
 local $/="\n>";
 while (<>) {
   s/^>//;
   chomp;
   my ($header, $seq)=(m/^([^\n]+)\n(.+)/s);
   $seq =~ tr/\t \n\r//d; 
   
   $seq =~ s/((\w{2,}?)\2{$c,})/'N' x length $1/oeg;
   # for poly-nucleotide
   #$seq =~ s/((\w)\2{6,})/'N' x length $1/oeg;
   $seq =~ s/((\w)\2{5,})/'N' x length $1/eg;
   #print $seq."\n";
   print ">$header\n";
   print join("\n", unpack('(A70)*', $seq))."\n";
   }
}
