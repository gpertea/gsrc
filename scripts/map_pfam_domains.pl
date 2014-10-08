#!/usr/bin/perl
use strict;
use Getopt::Std;
my $usage = q{
 Usage:
  map_pfam_domains.pl [-s <suffix>] <pfam_search_results>
  
  Assumes tran6frames.pl with -t option was used to
  generate the sequence names in the <pfam_search_results>, so it 
  separates the input data into <tag>-named files appropriately.

  Use -s option to specify a custom suffix instead of the 
  default "pfam" one.
  
 };

getopts('s:') || die($usage."\n");
 
my $suffix=$Getopt::Std::opt_s || 'pdom';
while(<>) {
    chomp;
    my @t=split;
    my @ndata=split(/\|/,$t[0]);
    my $coords=pop(@ndata);
    my ($beg, $end)=split(/_/,$coords);
    my $ftag=pop(@ndata);
    $beg=1 if $beg eq '0';
    $end=1 if $end eq '0';
    die(qq/Parsing error: invalid seq ID $t[0]. It should have the format:
        <seqname>|<filetag>|<coord1>_<coord2>/."\n") unless $beg && $end && $ftag;
    my $seqid=join('|',@ndata); #remaining '|' fields
    my ($newbeg, $newend) = ($beg<$end) ? 
               ($beg+($t[1]-1)*3, $beg+($t[2]-1)*3+2) :
               ($beg-($t[1]-1)*3, $beg-($t[2]-1)*3-2);
   my $newname=$ftag;
   $newname=~s/(\.\d+)$/\.$suffix$1/;
   open(OF, '>>'.$newname) || die("Error appending to file $newname!\n");
   print OF join(' ',$seqid, $newbeg, $newend, @t[3..9])."\n";
   close(OF);
   #printf("%s %20d %10d %20s %5d %5d %5s %10s %s\n",
   #    $seqid,$newbeg,$newend,@t[3..9]);
}

