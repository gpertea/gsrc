#!/usr/bin/perl
use strict;
use Bio::DB::Fasta;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;
my $gdefault='/fs/szdata/genomes/human_hg19/all_chrs.fa';
my $usage = qq/Usage:
 chrange.pl [-g <genome_fasta>] [-C][-F][-U] <chromosome>:<coord_start>-<coord_end>
 Options: 
  -C reverse complement the extracted sequence
  -U convert sequence to uppercase
  -F output fasta formatted sequence
  -g provide the custom genome file instead of the default:
     $gdefault
/;
umask 0002;
getopts('BCUFg:o:') || die($usage."\n");
die($usage."\n") unless @ARGV>0;
my $outfile=$Getopt::Std::opt_o;
my $gfile=$Getopt::Std::opt_g || $gdefault;
die("Error: no such file: $gfile\n") unless -f $gfile;
my ($fasta, $upcase, $revcompl) = ($Getopt::Std::opt_F, $Getopt::Std::opt_U, $Getopt::Std::opt_C);

my %opts = (-maxopen=>1);
$opts{-reindex}=1 if $Getopt::Std::opt_B;
my $db = Bio::DB::Fasta->new($gfile, %opts)
  || die("Error creating Bio::DB::Fasta($gfile)\n");
#print STDERR "    ..fasta file opened.\n";
foreach my $gloc (@ARGV) {
 my ($chr, $start, $rangetype, $stop)=($gloc=~m/^([\w\-\|]+)\:(\d+)([\.\-_\,\:]+)(\d+)/);
 die($usage."\n") unless $chr && $start>0 && $stop>0;
 $stop=$start+$stop-1 if ($rangetype eq ':');
 my $seq=$db->seq($chr.':'.$start.','.$stop);
 die("No sequence extracted at $chr:$start-$stop\n") unless $seq;
 $seq=uc($seq) if $upcase;
 $seq=revCompl($seq) if $revcompl;
 if ($fasta) {
   print ">$chr|$start"."_$stop\n".
       join("\n", unpack('(A70)*', $seq))."\n";
   }
  else {
   print "$seq\n";
   }
 } # for each range given

sub revCompl {
  my $s=reverse($_[0]);
  $s =~ tr/AaCcTtGgUuMmRrWwSsYyKkVvHhDdBb/TtGgAaCcAaKkYyWwSsRrMmBbDdHhVv/;
  return $s;
 }
