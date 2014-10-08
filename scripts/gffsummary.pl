#!/usr/bin/perl
use strict;
use Getopt::Std;

my $usage = q/Usage:
 gffstat [-a] [-t <ftype>] [-f <subfeature_name>] input.gffs..
 
 Returns statistics about the number, length and composition of 
 'mRNA' features and their 'CDS' segments. 
 Options:
    -q   omit header line
    -a   treat all input files as a single input stream, 
         so only the overall counts (totals) are shown
    -t   only consider entries with <ftype> in the 
         second column (e.g. "jigsaw")
    -f   use <subfeature_name> instead of default 'CDS'
         (needed for exon statistics)
/;
getopts('qat:f:') || die($usage."\n");

my ($readall, $tflt)=($Getopt::Std::opt_a, $Getopt::Std::opt_t);
my $subfname= $Getopt::Std::opt_f || 'CDS';
$subfname='CDS' unless $subfname;
print join("\t",'     ','G_count', 'G_cov', 'Ex_count', 'Ex_cov')."\n"
   unless ($Getopt::Std::opt_q);
if ($readall || @ARGV==0) {
  processInput();
  }
 else {
  foreach my $f (@ARGV) {
    processInput($f);
    }
 }

sub processInput {
 my $inf=$_[0];
 my ($gcount, $gextent, $exoncount, $exoncov);
 my $rlinefunc; 
 if ($inf) {
   open(INFILE, $inf)
      || die("Error opening file $inf!\n");
   $rlinefunc = \&readfile;
   } else {
   $rlinefunc=\&readstdin;
   }
  
while (&$rlinefunc()) {
 next if m/^\s*#/;
 my ($gobj, $ftype, $fname, $fstart, $fend, $fscore, $fstrand, $frame, $fdesrc)=split(/\t/); 
 next if ($tflt && $tflt ne $ftype);
 ($fstart, $fend)=($fend, $fstart) if $fstart>$fend;
 if ($fname eq 'mRNA') { 
    $gcount++;
    $gextent+= ($fend-$fstart+1);
    next;
    }
 if ($fname eq $subfname) {
    $exoncount++;
    $exoncov+=($fend-$fstart+1);
    }
 }#while parsing
 my $label='Total';
 if ($inf) {
   close(INFILE);
   ($label)=($inf=~m/([^\/]+)$/);
   }
 print join("\t", $label.':', $gcount, $gextent, $exoncount, $exoncov)."\n"; 
}

sub readfile {
 $_=<INFILE>;
 return $_;
}

sub readstdin {
 $_=<>;
 return $_;
}
