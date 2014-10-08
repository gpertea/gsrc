#!/usr/bin/perl
use strict;
#use Getopt::Std;
#use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 jigsaw2itt.pl <input_jigsaw_gff> <output_iit_file>
 
/;
umask 0002;
my ($infile, $outfile) = @ARGV;
die("$usage\nOutput file name required!\n") unless $outfile;

my $inh;

if ($infile eq '-') {
  $inh=\*STDIN;
  }
 else {
  open(INFILE, $infile) || die("Error opening input file $infile!\n");
  $inh=\*INFILE;
  }

my $curmodel; #current model name
my $curtag; #chromosome and strand for current model
my @exd; #exons for current model

open(TOIIT, '| iit_store -o '.$outfile) 
 || die("Error opening pipe to iit_store ($!)!\n");
while (<$inh>) {
  next if m/^\s*#/;
  chomp;
  my ($chr, $jsver, $exontype, $exonstart, $exonend, $jscore, 
      $strand, $frame, $lnum)=split(/\t/);
  next unless $lnum>0;
  ($exonstart, $exonend)=($exonend, $exonstart) if $exonend<$exonstart;
  my $locus=$chr.'.jsm.'.$lnum;
  if ($locus ne $curmodel) {
    &writeModel() if $curmodel;
    $curmodel=$locus;
    $curtag=$chr.$strand;
    @exd=([$exonstart, $exonend]);
    next;
    }
  push(@exd, [$exonstart, $exonend]);
  }

writeModel() if $curmodel;

close(INFILE) unless $infile eq '-';
close(TOIIT);

sub writeModel {
 my @ex= sort { $main::a->[0] <=> $main::b->[0] } @exd;
 my ($mstart, $mend)=($ex[0]->[0], $ex[-1]->[1]);
 my @exw = map { $_->[0].'-'.$_->[1] } @ex;
 print TOIIT ">$curmodel $mstart $mend $curtag\n";
 print TOIIT join(',',@exw)."\n";
 }
