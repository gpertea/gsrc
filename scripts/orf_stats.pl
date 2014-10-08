#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 orf_stats.pl <orf_check.tab> <ref_mrnas.iit>

 It will try to report various stats about the ORFs 
 and junctions of transfrags given in <orf_check.tab>

 <orf_check.tab> is the output of orf_check.pl, with the columns:
 <transcript_id>,<xlocus>,<chromosome><strand>,<exons>,<orfs>
 
 ORF coordinates are local - i.e. relative to the transcript sequence
/;
umask 0002;
getopts('o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my ($forfs, $fiit)=@ARGV;
die($usage."\n") unless -f $forfs; # && -f $fiit;
my %ci; #cufflinks isoforms data
# transcript_id => [$chr, $strand, $cstart, $cend, [ exondata..], [ orfdata.. ]]
my $fiitqry=$forfs.'.iitqry';

open(QDATA, '>'.$fiitqry) || die("Error creating iit query file $fiitqry\n");

open(FORFS, $forfs) || die("$usage\n");
while(<FORFS>) {
 chomp;
 next unless $_;
 my ($t, $xloc, $chr, $exonstr, $orfstr)=split(/\t/);
 my $strand=chop($chr);
 $exonstr=~s/^exons=//;
 $orfstr=~s/^orfs=//;
 my @ex= map { [split(/\-/)] } split(/\,/,$exonstr);
 my @orfs;
 if ($orfstr) {
   @orfs=map { [split(/\:/)] } split(/\,/,$orfstr);
   }
 my @t2g;
 my @gorfs;
 my $x=0; #position on transcript sequence
 if ($strand eq '-') {
   foreach my $d (reverse @ex) {
     foreach my $gx (reverse($$d[0] .. $$d[1])) { $t2g[++$x]=$gx }
     }
   #@gorfs= map { [$t2g[$$_[0]+$$_[1]]-2, $t2g[$$_[0]]] } @orfs;
   # keep the order so we can more easily match GFF CDS end coordinate
   if ($orfstr) {
     @gorfs= map { [$t2g[$$_[0]], $t2g[$$_[0]+$$_[1]]-2] } @orfs;
     }
   }
 else {  
   foreach my $d (@ex) {
     foreach my $gx ($$d[0] .. $$d[1]) { $t2g[++$x]=$gx }
     }
   if ($orfstr) {
     @gorfs= map { [$t2g[$$_[0]], $t2g[$$_[0]+$$_[1]]+2] } @orfs;
     }
   }
 # - debug print:
 #print join("\t", $t, $chr.$strand, join(',', (map { $$_[0].'-'.$$_[1] } @ex)),
 #  join(',', (map { $$_[0].'|'.$$_[1] } @gorfs)))."\n";
 print QDATA join(' ',$ex[0]->[0], $ex[-1]->[1], $chr.$strand, '<'.$t)."\n";
 $ci{$t}=[$chr, $strand, $ex[0]->[0], $ex[-1]->[1], [@ex], [@gorfs]];
}
close(QDATA);
close(FORFS);
# -- now query the iit_file with $fiitqry
goto ENDP unless ($fiit && -f $fiit);

my $cmd="iit_get $fiit < $fiitqry";
open(IITGET, "$cmd |") || die ("Error opening pipe $cmd |\n");
my @refs; #refs overlapping query transfrag -- forming this cluster
my $tf; # current transfrag;
my $strand; #strand for current cluster
my $lref; #last ref id
my $lchrstrand; #last ref chr and strand
my %accs; #all acceptor positions for this cluster
my %dons; #all donor positions for this cluster
my %cdstops; #all CDS stop locations
my $jrecomb; # counter for transfrags that have all junctions present in the cluster
my $samestop; #counter for transfrags that have an orf stopping 
              # at a known stop in the cluster (e.g. in %cdstops)
my $haveORFs;
while (<IITGET>) {
 chomp;
 next unless $_;
 #print ">> $_\n";
 if (m/^# Query: /) {
    ($strand, $tf)=(m/([\+\-]) <([^<]+)$/);
    die("Error parsing current transfrag in: $_\n") 
        unless $strand && $tf && exists($ci{$tf});
    $haveORFs=(@{$ci{$tf}->[5]}>0);
    next;
    }
 if (m/>(\S+) \d+ \d+ (\S+)/) {
    #header for a ref interval
    ($lref, $lchrstrand)=($1, $2);
    }
    elsif (m/^\d+\-\d+/) {
     #exons line:
     if ($strand eq '-') {
       my @acc=(m/\-(\d+)\,/g);
       @accs{@acc}=();
       my @don=(m/\,(\d+)\-/g);
       @dons{@don}=();
       }
      else {
       my @acc=(m/\,(\d+)\-/g);
       @accs{@acc}=();
       my @don=(m/\-(\d+)\,/g);
       @dons{@don}=();
       }
     }
    elsif($haveORFs && m/^C:\d+\-/) {
     # CDS line
     my $cdstop;
     if ($strand eq '-') {
       ($cdstop)=(m/^C:(\d+)/);
       #print STDERR ">>>> storing stop $cdstop\n";
       }
      else {
       ($cdstop)=(m/(\d+)$/);
       }
     $cdstops{$cdstop}=1;
     }
 if (m/^# End$/) {
    # collect stats here
    my ($tchr, $tstrand, $tmin, $tmax, $texs, $torfs)=@{$ci{$tf}};
    die ("ChrStrand mismatch!\n") unless $tstrand eq $strand;
    my @acc;
    my @don;
    if ($tstrand eq '-') {
      @acc= map { $$_[1] } @$texs;
      pop(@acc);
      @don= map { $$_[0] } @$texs;
      shift(@don);
      }
     else {
      @don= map { $$_[1] } @$texs;
      pop(@don);
      @acc= map { $$_[0] } @$texs;
      shift(@acc);
      }
    # check if all acceptors and all donors are in the cluster
    my $alljmatch=0;
    my $print=0;
    my @pcol=('','');
    foreach my $a (@acc) { 
       goto JNOTFOUND if !exists($accs{$a});
       }
    foreach my $d (@don) { 
       goto JNOTFOUND if !exists($dons{$d});
       }
    $alljmatch=1;
    JNOTFOUND:
    if ($alljmatch) {
       #print "$tr\tjrecomb\n";
       $jrecomb++;
       $print=1;
       $pcol[0]='jrecomb';
       }
    #now check if same stop codon
    foreach my $od (@$torfs) {
      #print STDERR ">>>> >> checking stop $$od[1]\n" if $strand eq '-';
      if (exists($cdstops{$$od[1]})) {
        $print=1;
        $pcol[1]='samestop';
        $samestop++;
        last;
        }
      }
    print join("\t",$tf,@pcol)."\n" if $print;
    # -- clear for next run  
    ($strand, $tf)=();
    %accs=();
    %dons=();
    %cdstops=();
    @refs=();
    }
 }
close(IITGET);
print STDERR "$jrecomb transfrags are recombinations of known junctions.\n";
print STDERR "$samestop transfrags have an ORF that match a known stop codon.\n";

# --
ENDP:
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************
