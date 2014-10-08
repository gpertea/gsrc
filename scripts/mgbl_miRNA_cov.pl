#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 mgbl_miRNA_cov.pl [-P] [-v <max_mms>] [-g <maxgaps>] [-f <valid_overlaps>] \
       [-o <output_cov_data>] <mgblast_output.tab>
 
Filters the overlaps resulted from a mgblast search against a database of 
known miRNAs and estimates expression level (for mature miRNAs).

 Options:
  -P precursor miRNAs were used for the mgblast search
  -g discard overlaps with more than <maxgaps> gaps (default 2)
  -v discard overlaps with more than <max_mms> mismatched bases
     at the 3' end of the read (default 4 at the 3' end of the read);
  -U unique mappings of reads to miRNAs: assuming the hits are sorted
     by score, only use the first hit for a read and discard the rest
  -R append column with all mapped reads for each miRNA
  -B consider aligments on both strands (normally only Plus strand 
     alignments are kept)
  Other implied restrictions:
  * only 1 base mismatch is allowed at the 5' end of the read;
  * in case of mature miRNAs, an overlap should involve at least 90% of the 
    read length and at least 70% of the mature miRNA length
  * in case of precursor miRNAs, an overlap should involve at least 90% of 
    the read length
/;
umask 0002;
getopts('BRUPf:g:v:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $uniqmap=$Getopt::Std::opt_U;
my $rshow=$Getopt::Std::opt_R;
my $bothstrands=$Getopt::Std::opt_B;

my $fflt=$Getopt::Std::opt_f;
my $maxgaps=$Getopt::Std::opt_g || 2;
my $maxovh=$Getopt::Std::opt_v || 4;
my $premiRNA=$Getopt::Std::opt_P;
my $minmcov=70;
my $minrcov=90;
#my $totalreads=$Getopt::Std::opt_n;
#die($usage."Error: total number of reads in the sample is required! (-n parameter)\n")
#   unless $totalreads>1;

if ($outfile) {
 open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
 select(OUTF);
 }
if ($fflt) {
 open(FLT, '>'.$fflt) || die("Error creating output file $fflt\n");
 }
 
my %mdata; # mirna => [$cov, $read1, $read2,...]
my %rd; # reads assigned uniquely to miRNAs - first come first served, 
        # assuming first hit is the best one!
while (<>) {
 my $mgline=$_;
 chomp;
 my ($read, $rlen, $rstart, $rend, $mirna, $mlen, $mstart, $mend,
     $pid, $bitscore, $pvalue, $strand, $rgaps, $mgaps)=split(/\t/);
     
 next unless $read && $mirna;
 next if $uniqmap && exists($rd{$read}); #read already assigned to a mature mRNA
 next if ($strand eq '-' && !$bothstrands);
 my ($read_mult)=($read=~m/_x(\d+)$/); # read multiplicity
 print STDERR ">checking $mgline";
 print STDERR "read_multiplier=$read_mult, pid=$pid\n";
 next if ($read_mult>=3 && $pid<100);
 ($rstart, $rend)=($rend, $rstart) if $rend<$rstart;
 next if ($rend-$rstart+1)*100<$minrcov*$rlen;
 unless ($premiRNA) {
  next if ($mend-$mstart+1)*100<$minmcov*$mlen;
  }
 my ($rovh3, $rovh5)=($strand eq '-') ? ($rstart-1, $rlen-$rend) : ($rlen-$rend, $rstart-1);
 my ($movh3, $movh5)=($mlen-$mend, $mstart-1);
 my $mm3=($rovh3<$movh3) ? $rovh3 : $movh3;
 my $mm5=($rovh5<$movh5) ? $rovh5 : $movh5;
 my $maxmm=$maxovh;
 $maxmm-- if $pid<100;
 next if $mm3>$maxmm || $mm5>1;
 my $numgaps=0;
 if ($rgaps) {
    my @rg=split(/\,/, $rgaps);
    foreach my $g (@rg) {
      $numgaps+= ($g=~m/\+(\d+)$/) ? $1 : 1;
      }
    next if $numgaps>2;
    }
 if ($mgaps) {
    my @rg=split(/\,/, $mgaps);
    foreach my $g (@rg) {
      $numgaps+= ($g=~m/\+(\d+)$/) ? $1 : 1;
      }
    next if $numgaps>2;
    }
 next if ($read_mult>=3 && $numgaps>0);
 # -- here we are, valid overlap, keep it:
 print FLT $mgline if $fflt;
 $rd{$read}=$mirna;
 my ($rmul)=($read=~m/_x(\d+)$/); #read multiplicity
 my $md=$mdata{$mirna};
 unless ($md) {
   my @mc=(0) x $mlen;
   @mc[($mstart-1)..($mend-1)]=($rmul) x ($mend-$mstart+1);
   # -------------   0      1        2       3     4      5
   $mdata{$mirna}=[$rmul, $mlen, $mstart, $mend, [@mc], $read];
   next;
   }
 # add to existing miRNA
 $md->[0]+=$rmul;
 $$md[2]=$mstart if ($$md[2]>$mstart);
 $$md[3]=$mend   if ($$md[3]<$mend);
 my $mc=$$md[4];
 for (my $i=$mstart-1;$i<$mend;$i++) {
    $$mc[$i]+=$rmul;
    }
 push(@$md, $read);
 } #while input lines
# report coverage per known miRNA 
my $k=keys(%mdata);
my ($maxmirna, $maxcov);
while (my ($mirna, $md)=each(%mdata)) {
 my @r=@$md;
 my $mcov=shift(@r);
 my $mlen=shift(@r);
 my $mcovstart=shift(@r);
 my $mcovend=shift(@r);
 my $covdata=shift(@r); #coverage per base
 # find maximum coverage for this miRNA
 my $maxc=$$covdata[0];
 foreach my $c (@$covdata) { $maxc=$c if $c>$maxc; }
 my $mlencov = ($mcovstart==1 && $mcovend==$mlen)? 'full' : sprintf('%.2f',($mcovend-$mcovstart)/$mlen);
 ($maxcov, $maxmirna)=($mcov, $mirna) if $maxcov<$mcov;
 #print join("\t",$mirna, sprintf('%.4f',$mcov/$totalreads), $mcov.'x', 
 my $rlist=$rshow ? "\t".join(',',@r) : '';
 print join("\t",$mirna, $mcov, $maxc, $mlencov)."$rlist\n";
 }
print STDERR "Done. $k miRNAs covered.\n".
  "Maximum coverage was found for $maxmirna at ${maxcov}x.\n";
#- script ending here:
close(FLT) if $fflt;
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#======== Subroutines ============#
