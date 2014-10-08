#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 pullctghits <mgtab_hits> <acefile> <contigname> [-o <mgtab_file>]
 by default outputfile <mgtab_file> is <contigname>.mgtab
 
/;
die "$usage\n" unless @ARGV>=3;
my $inhits=shift @ARGV;
my $acefile=shift @ARGV;
my $ctgname=shift @ARGV;
die "$usage\nCannot locate input hits file $inhits\n" unless -f $inhits;
die "$usage\nCannot open ACE file $acefile\n" unless -f $acefile;
umask 0002;
getopts('o:') || die($usage."\n");
my $fctgname=$ctgname;
$fctgname=~tr/|/_/;
my $outfile=$Getopt::Std::opt_o || $fctgname.'.mgtab';
my $flst=$fctgname.'.rlst';
my $cmd="acetool $acefile -c '$ctgname' | grep '^AF ' | cut -f2 -d' '";
system("$cmd > $flst") && die "Error $! at command $cmd\n";
if (-s $flst) {
 my $pre='cat';
 if ($inhits=~m/\.gz$/i || $inhits=~m/\.Z$/i) {
    $pre='gzip -cd';
    }
   elsif ($inhits=~m/\.bZ$/i || $inhits=~m/\.bz2/) {
     $pre='bzip2 -cd'
     }
 
 $cmd="$pre $inhits | tclust PID=80 -o t.cls -r $flst -f $outfile";
 print STDERR "Running $cmd:\n";
 system($cmd);
 print STDERR "done.\n";
 }
else { #zero length
 die "Error: no such contig in $acefile?\n";
 } 
