#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 mimap_asm.pl [-o <outfile>] [-s <minscore>] <sorted_bowtie.map> 
  Takes sorted bowtie .map data from miRNA mappings and outputs "assembled"
  transcripts\/loci with depth of coverage information

  The bowtie input MUST have been sorted by chromosome and coordinate,
  with a command like this:
      sort -k3,3 -k4,4n bowtie_output.map > bowtie_sorted.map
/;
umask 0002;
getopts('Ds:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my ($maxasmlen, $maxasm);
my $debug=$Getopt::Std::opt_D;
if ($outfile) {
  open(FOUT, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(FOUT);
  }
my $minscore=$Getopt::Std::opt_s || 1.8;

my ($locchr, $locstrand, $locstart, $locend, $loccov, $locseq);
my @locreads; # list of [readname, chr_alnstart, chr_alnend, num_mismatches, strand]
my $lcount=0;

while (<>) {
 next if m/^#/;
 chomp;
 my ($read, $strand, $chr, $chrstart, $rseq, $rquals, $num, $mism)=split();
 next unless $rseq && length($rseq)==length($rquals);
 $chrstart++; # because bowtie coordinate is 0-based
 my ($rcov)=($read=~m/_x(\d+)$/);
 my ($rtrimseq, $tstart, $tend, $num_mism)=mapTrim($strand, $chrstart, $rseq, $mism);
 # DEBUG display:
 my $debug_line;
 if ($debug) {
  my $lpad=' 'x($tstart-$chrstart);
  
  $debug_line="DBG: ".join("\t", $chr.$strand, $chrstart, $rseq,        length($rseq), $mism." \t\t   locpos=[$locstart,$locend]\n", 
                               '     ',$tstart, $lpad.$rtrimseq, $tend, "($num_mism mismatches left)  tpos=[$tstart,$tend]")."\n";
  }
 if (($locchr ne $chr) || ($locstart==0) || ($tstart+12>$locend)) {
     #new chromosome, or different locus
     writeLoc() if $locchr; #flush previous locus, if any
     # now just create the new locus:
     ($locchr, $locstrand, $locstart, $locend, $loccov, $locseq)=
     ($chr,    $strand,    $tstart,   $tend,   $rcov,   $rtrimseq);
     @locreads=([$read, $tstart, $tend, $num_mism, $strand]);
     if ($debug) {
        print STDOUT $debug_line;
        }
     next;
     }
 if ($debug) {
   print STDOUT $debug_line;
  #if ($locstrand && $locstrand ne $strand) {
 #      print STDERR "Warning: merging opposite strand match ($read) into locus at $locchr:$locstart)\n";
 #      }
   }
 # -- merge into current locus:
 if ($tstart<$locstart) { 
    # - shouldn't happen because the mappings are sorted !
    #copy the part of $rtrimseq that's sticking out on the left
    $locseq=substr($rtrimseq, 0, $locstart-$tstart).$locseq;
    $locstart=$tstart;
    print STDOUT "DBG:  pre: $locseq\n" if $debug;
    }
 if ($tend>$locend) {
    #append the part of rtrimseq that's sticking out on the right
    $locseq.=substr($rtrimseq, -($tend-$locend));
    print STDOUT "DBG:  post: $locseq\n" if $debug;
    $locend=$tend;
    }
 $loccov+=$rcov;
 push(@locreads, [$read, $tstart, $tend, $num_mism, $strand]);
 # -- merging into current locus
 }
 
writeLoc(); 

print STDERR "Done. Maximum assembly length=$maxasmlen ($maxasm).\n";
 
# END of script here:
if ($outfile) {
  select(STDOUT);
  close(FOUT);
  }

#------- SUBROUTINES -----------

sub mapTrim {
 my ($strand, $cstart, $seq, $mism)=@_;
 my $len=length($seq);
 my $cend=$cstart+$len-1;
 return ($seq, $cstart, $cend, 0) unless $mism;
 my @b=unpack('(A)*',$seq);
 my @mm=map { $_=[ (m/^(\d+)\:(\w)\>(\w)/) ] } (split(/\,/, $mism));
 my ($ltrimstart, $ltrimlen, $rtrimstart, $rtrimlen)=(-1,0,0,0);
 my $prevtrim=1000000;
 my $halflen = $len>>1;
 #mismatches are sorted
 @mm = reverse(@mm) if ($strand eq '-');
 my @rmm;
 my ($rmmdel, $lmmdel); #trimmed mismatches
 foreach my $m (@mm) {
   #my ($x, $gn, $rn)=($m=~m/^(\d+)\:(\w)\>(\w)/);
   my ($x,$gn, $rn) = @$m; # $x values are always increasing
   $x=$len-$x-1 if ($strand eq '-');  # $ x values are always increasing
   die("Error: nucleotide correction failure at pos $x (len $len, strand $strand)!\n(seq=$seq with $mism, $rn != $b[$x])\n") 
       unless $rn eq $b[$x];
   $b[$x]=$gn; #restore the genomic sequence
   if ($x>$halflen) {
      unshift(@rmm, [$x, $gn, $rn]);
      next;
      }
   # -- trimming left end
   if ($x<2 && $ltrimlen==0) { $ltrimstart=$x; $lmmdel++; $ltrimlen++; $prevtrim=$x; next; }
   if (($x-$prevtrim)==1 || ($ltrimlen>1 && $x-$prevtrim==2)) {
          $ltrimlen+=($x-$prevtrim);
          $lmmdel++;
          $prevtrim=$x;
          }
   }#for each mismatch
 #try trimming the other end
 $prevtrim=1000000;
 foreach my $m (@rmm) {
   my ($x,$gn, $rn) = @$m; # $x values are always decreasing
   if ($x>$len-3 && $rtrimlen==0) { $rtrimstart=$x; $rmmdel++; $rtrimlen++; $prevtrim=$x; next; }
   if ($prevtrim-$x==1 || ($rtrimlen>1 && $prevtrim-$x==2)) {
     $rtrimlen+=($prevtrim-$x);
     $rmmdel++;
     $prevtrim=$x;
     }
   }

 # 0-based coordinates (within the read) of the "non-trimmed" part
 my $num_mism=@mm;
 my $tstart= ($ltrimstart==0 || $ltrimlen>1) ? $ltrimstart+$ltrimlen : 0 ;
 $num_mism-=$lmmdel if $tstart>0;
 my $tend = ($rtrimstart==$len-1 || $rtrimlen>1) ? $rtrimstart-$rtrimlen : $len-1 ;
 $num_mism-=$rmmdel if $tend<$len-1;
 return (join('',@b[$tstart .. $tend]), $cstart+$tstart, $cstart+$tend, $num_mism);
}


sub writeLoc {
return unless $locchr && $locstart;
$lcount++;
my @rdata;
my ($mapscore, $ftotal, $rtotal);
#also recompute strand, in case we have both strand alignments
foreach my $d (@locreads) {
  my ($read, $rmstart, $rmend, $rmism, $rstrand)=@$d;
  push(@rdata, $read.$rstrand.':'.($rmstart-$locstart).'-'.($rmend-$locstart));
  my ($rcov)=($read=~m/_x(\d+)$/);
  my $rmscore=(($rmend-$rmstart-$rmism)*$rcov)/($rmend-$rmstart);
  if ($rstrand eq '-') { $rtotal+=$rmscore; }
                  else { $ftotal+=$rmscore; }
  $mapscore+=$rmscore;
  }

$locstrand = ($rtotal>$ftotal)?'-':'+';
$locseq = reverseComplement($locseq) if $locstrand eq '-';
#format:
#
# asm_id, chr, strand, chr_start, chr_end, coverage, mapping_score, annotation, sequence, reads
my $nscore=$mapscore/$loccov; 
if ($mapscore>$minscore && (@rdata>1 || $nscore>0.96)) {
  print join("\t", sprintf('miasm_%06d', $lcount), $locchr,
       $locstrand, $locstart, $locend, $loccov.'x', sprintf('%.2f',$mapscore), '.', $locseq)."\t".
       join(',', @rdata)."\n";
  if (length($locseq)>$maxasmlen) {
    $maxasmlen=length($locseq);
    $maxasm=sprintf('miasm_%06d', $lcount);
    }
  }
#DEBUG:
print "-----------------------------------------\n" if $debug; 
@locreads=();
($locchr, $locstrand, $locstart, $locend, $locseq, $loccov)=();
}

sub reverseComplement {
  my $s=reverse($_[0]);
  $s =~ tr/AaCcTtGgUuMmRrWwSsYyKkVvHhDdBb/TtGgAaCcAaKkYyWwSsRrMmBbDdHhVv/;
  return $s;
 }
