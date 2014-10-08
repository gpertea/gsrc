#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 miRNA_cmp.pl [-o <output.tab>] [-d diffexp.tab] \
   <sample1_cov.micov>[:<totalreads1>] <sample2_cov.micov>[:<totalreads2>]
 
 Compares the expression levels of known miRNAs in two experiments
 based on the the coverage files produced by mgbl_miRNA_cov.pl
 If -d option is given, that file will have entries with a significant 
 p-value according to Fisher's exact test with Bonferroni correction.
 
 If total reads counts are not given, a "context total" is computed by 
 summing up existing counts (coverage values) in each sample
/;
umask 0002;
print STDERR "#Command line was:\n".$FindBin::Script.
              ' '.join(' ',@ARGV)."\n";
getopts('d:o:') || die($usage."\n");
die($usage."\n") unless @ARGV==2;
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my $dfile=$Getopt::Std::opt_d;
if ($dfile) {
  open(FDIFF, '>'.$dfile) || die("Error creating output file $dfile\n");
  }
my ($fd1, $fd2)=@ARGV;
my ($f1, $total1)=split(/[\:\,~]+/, $fd1);
die("$usage File $f1 not found!\n") unless -f $f1;
#die("$usage Total reads incorrect!\n") unless $total1>10;
my ($f2, $total2)=split(/[\:\,~]+/, $fd2);
die("$usage File $f2 not found!\n") unless -f $f2;
#die("$usage Total reads incorrect!\n") unless $total2>10;
if ($total1>0 && $total2>0) {
print STDERR "Total reads for $f1: $total1\n";
print STDERR "Total reads for $f2: $total2\n";
}
# -- Fisher's exact test - utility constants
my $PI=atan2(0,-1);
my $twoPI=2*$PI;
my $pi_3=$PI/3;
my @lnfcache=(0,0);
#---

my %mi; # miRNA => [ cov_sample_1, cov_sample_2 ]
open(FIN, $f1) || die("Error opening file $f1\n");
while (<FIN>) {
  chomp;
  my @t=split(/\t/);
  #            n1    n2  foldchange  normalized  Fisher p
  $t[0]=lc($t[0]);
  $mi{$t[0]}=[$t[2], 0,      0,       0,      0];
  }
close(FIN);

open(FIN, $f2) || die("Error opening file $f2\n");
while(<FIN>) {
  chomp;
  my @t=split(/\t/);
  $t[0]=lc($t[0]);
  my $d=$mi{$t[0]};
  if ($d) {
     $$d[1]=$t[2];
     }
    else {
     $mi{$t[0]}=[0, $t[2], 0, 0];
     }
  }
close(FIN);
#print '#'.join("\t", 'microRNA', 'count1', 'count2', 'fold_change', 'nfold_change', 'FET_p')."\n";
if ($total1<10 && $total2<10) { #compute contextual totals here
 $total1=0;
 $total2=0;
 while (my ($mirna, $md)=each(%mi)) {
   #      0   1      2          3         4       
   my ($n1, $n2, $fchange, $nfchange, $fisher)=@$md;
   next if $n1<2 && $n2<2;
   $total1+=$n1;
   $total2+=$n2;
   }
 print STDERR "Total reads given for sample 1: $total1\n";
 print STDERR "Total reads given for sample 2: $total2\n";
 } # compute totals if not given
my $num_miRNAs=keys(%mi);
print STDERR "$num_miRNAs total miRNAs seen.\n";
my $p_corrected=0.05/$num_miRNAs;
print STDERR "Bonferroni-corrected p-value is 0.05/$num_miRNAs = $p_corrected\n";
print join("\t", '#totals:', 'placebo|'.$total1, 'statin|'.$total2, $p_corrected)."\n";
my $numexpr;
my $numdiffp; # count of miRNAs with p<$p_corrected
my $numdiffp0; # count of miRNAs with p<=0.001
my $numdiff; # count of miRNAs with p<=0.001 AND abs(foldchange)>1.5
my $t_ratio = $total1/$total2;

while (my ($mirna, $md)=each(%mi)) {
#      0   1      2          3         4       
 my ($n1, $n2, $fchange, $nfchange, $fisher)=@$md;
 next if $n1<2 && $n2<2;
 my ($n1real, $n2real)=(int($n1),int($n2));
 $n1=1 if $n1==0;
 $n2=1 if $n2==0;
 my ($nn1, $nn2)=($n1, $n2*$t_ratio);
 $fchange =($n1<$n2) ? '+'.sprintf('%.2fx',$n2/$n1) : '-'.sprintf('%.2fx', $n1/$n2);
 $nfchange = ($nn1<$nn2) ? '+'.sprintf('%.2fx',$nn2/$nn1) : '-'.sprintf('%.2fx', $nn1/$nn2);
 $$md[2]=$fchange;
 $$md[3]=$nfchange;
 my $nfc=$nfchange;
 $nfc=~tr/x+-//d;
 my ($a, $b, $c, $d)=($n1, $total1-$n1, $n2, $total2-$n2);
 $fisher=calcFET($a, $b, $c, $d);
 $$md[4]=$fisher;
 $numexpr++;
 my $diffprinted=0;
 #print join("\t", $mirna, $n1, $n2, $fchange, $nfchange, $fisher, $nfc)."\n";
 
 print join("\t", $mirna, $n1real, $n2real, $nfchange, $fisher)."\n";
 if ($fisher<=0.001 && $nfc>=1.5) {
  print FDIFF join("\t", $mirna, $n1real, $n2real, $fchange, $nfchange, $fisher, $nfc)."\n";
  $diffprinted=1;
  $numdiff++;
  }
 $numdiffp0++ if ($fisher<=0.001);
 if ($fisher<=$p_corrected) {
     $numdiffp++;
     print FDIFF join("\t", $mirna, $n1real, $n2real, $fchange, $nfchange, $fisher, $nfc)."\n"
       if ($dfile && $diffprinted==0);
     }
}

print STDERR "$numexpr miRNAs with read count>=2 in either sample.\n";
print STDERR "$numdiffp0 with FET p-value<=0.001\n";
print STDERR "$numdiffp with FET p-value<=$p_corrected\n";
print STDERR "$numdiff with nfold change >1.5 and p-value<=0.001\n";

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }
if ($dfile) {
 close(FDIFF);
 }

sub LnFactorial {
 my $n=shift;
 die "Bad args to log_factorial: $n" if $n < 0;
 my $ln_fact;
 if ($n<900) {
    return $lnfcache[$n] if defined $lnfcache[$n];
    for (my $i = scalar(@lnfcache); $i <= $n; $i++) {
       $lnfcache[$i] = $lnfcache[$i-1] + log($i);
       }
    return $lnfcache[$n];
    }
   else { 
    # Gosper's approximation; from
    # http://mathworld.wolfram.com/StirlingsApproximation.html
    $ln_fact = 0.5 * log($twoPI*$n + $pi_3) + $n*log($n) - $n;
    }
 return $ln_fact;
 }


# Compute the probability of getting this exact table
# using the hypergeometric distribution
sub ProbOneTable {
  my ($a , $b , $c, $d) = @_;
  my $n = $a + $b + $c + $d;
  my $LnNumerator     = LnFactorial($a+$b)+
                        LnFactorial($c+$d)+
                        LnFactorial($a+$c)+
                        LnFactorial($b+$d);

  my $LnDenominator   = LnFactorial($a) +
                        LnFactorial($b) +
                        LnFactorial($c) +
                        LnFactorial($d) +
                        LnFactorial($n);

  my $LnP = $LnNumerator - $LnDenominator;
  return exp($LnP);
}

# Compute the cumulative probability by adding up individual
# probabilities
sub calcFET {
  my ($a, $b, $c, $d) = @_;

  my $min;

  my $n = $a + $b + $c + $d;

  my $p = 0;
  $p += ProbOneTable($a, $b, $c, $d);
  if( ($a * $d) >= ($b * $c) ) {
    $min = ($c < $b) ? $c : $b;
    for(my $i = 0; $i < $min; $i++) {
      $p += ProbOneTable(++$a, --$b, --$c, ++$d);
    }
  }

  if ( ($a * $d) < ($b * $c) ) {
    $min = ($a < $d) ? $a : $d;
    for(my $i = 0; $i < $min; $i++) {
      $p += ProbOneTable(--$a, ++$b, ++$c, --$d);
    }
  }
  return $p;
}
