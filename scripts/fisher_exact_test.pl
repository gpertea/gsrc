#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;
use Fisher qw(fishers_exact log_factorial);


my $usage = q/Usage:
 fisher_exact_test.pl [-T] a1 b1 a2 b2
 Computes Fisher's exact test p-values for 2x2 matrix:
   
   a1   b1
          
   a2   b2
 
 -T option will make the input values to be taken in terms of totals:
 
   t1 t2 n1 n2
   
  ..where n1 is a subset of t1 and n2 a subset of t2, i.e. 
    t1=a1+b1 and t2=a2+b2
 
/;
umask 0002;
getopts('To:') || die($usage."\n");
die($usage." Error: exactly 4 values are needed (totals first)!\n") unless @ARGV==4;
my $totals=$Getopt::Std::opt_T;
my @fncache=(0,0);

# -- needed by Audic-Calverie: calc_audic()
my $MAXIT = 500;
my $EPS   = 3.0E-30;
my $FPMIN = 1.0E-30;
my @COF = ( 76.18009172947146, -86.50532032941677,
            24.01409824083091, -1.231739572450155,
            0.1208650973866179E-2,-0.5395239384953E-5 );



{
 my ($a, $b , $c, $d)=@ARGV;
 if ($totals) {
   my ($t1, $t2, $n1,$n2)=@ARGV;
   die("Error: totals should be larger than parts!\n") 
    unless $n1<$t1 && $n2<$t2;
   ($a, $b , $c, $d)=($n1, $t1-$n1, $n2, $t2-$n2);
   print STDERR " a b c d = $a $b $c $d\n";
   #my $ap=audic($t1, $t2, $n1, $n2);
   my ($ap, $apdir)=calc_audic($n1, $n2, $t1, $t2, 1);
 print ">>>>   Audic & Claverie  : $ap ($apdir)\n";
   }
 my $FisherProb = ProbCTable($a, $b , $c, $d);
 print "-----> Fisher Probability: $FisherProb\n";
 my $fet1=fishers_exact($a,$b,$c,$d);
 my $fet2=fishers_exact($a,$b,$c,$d, 1); 
 print "Fisher.pm :   one-tailed:  $fet1\n";
 print "Fisher.pm :   two-tailed:  $fet2\n";
 exit;
 #--- rest is optional:
 my $a_b = $a + $b;
 my $c_d = $c + $d;
 my $a_c = $a + $c;
 my $b_d = $b + $d;
 my $total = $a + $b + $c + $d;
 my $p1 = $a/$a_b;
 my $p2 = $c/$c_d;
 my $deltap = abs($p1 - $p2);
 my $p  = $a_c / $total;
 my $sd = sqrt($p * (1-$p) * (1/$a_b + 1/$c_d));
 my $z  = abs($p1 - $p2) / $sd;
 my $zp = abs(ltqnorm($FisherProb));
 my $z95 = abs(ltqnorm(0.95));

 my $z_ratio = $z / $zp;
 my $deltap95 = $deltap * ($z95 / $zp);
 my $ratio95 = $z95 / $zp;
 print "Deltap:$deltap\tProjectedDeltap95:$deltap95\tRatio:$ratio95\n\n";

 print "$a\t$b\t$a_b\n$c\t$d\t$c_d\n".
   "$a_c\t$b_d\t$total\nFisherProb: $FisherProb <----\n" .
   "z:$z\nzp:$zp\nratio:$z_ratio\n$p1\t$p2\t$p\n\n";

  my $p_this_table = ProbOneTable($a , $b , $c , $d);
  print "This one table p:$p_this_table\n";

  my $mid_p = $FisherProb - $p_this_table/2;
  print "Mid-p:$mid_p\n\n";
}

#================ SUBROUTINES ============

sub LnFactorial {
    my $n = shift;
    return $fncache[$n] if defined $fncache[$n];
    return undef if $n < 0;
    for (my $i = scalar(@fncache); $i <= $n; $i++) {
      $fncache[$i] = $fncache[$i-1] + log($i);
    }

    return $fncache[$n];
}
sub factorial {
 my $r=1;
 $r *= $_ foreach 2..$_[0];
 return $r;
}

#compute Audic & Claverie probability
sub audic { 
 my ($n1, $n2, $x, $y)=@_;
 my $v=$y*log($n2/$n1)+log_factorial($x+$y)-log_factorial($x)-log_factorial($y)-
          ($x+$y+1)*log(1+$n2/$n1);
 return exp($v);
 #my $nratio=$n2/$n1;
 #my $v=($nratio ** $y)*factorial($x+$y)/(factorial($x)*factorial($y)*((1+$nratio)**($x+$y+1)));
 
}

sub calc_audic {
=pod

=head2 calc_audic $x, $y, $Nx, $Ny, <$signedValue>

Determines the statistical significance of the difference
in tag count (expression) between two libraries.  This
function uses the method described by Audic and 
Claverie (1997).  This method can be called on an
instantiated object, as well as statically.

B<Arguments>

I<$x,$y>

  The number of tags in the x- and y-axis 
  libraries, respectively.

I<$Nx,$Ny>

  The total number of tags in the x- and y-axis
  libraries, respectively.

I<$signedValue> (optional)

  A boolean value (>=1 is FALSE).  If this flag is
  set to TRUE, downregulated comparisons will return
  both a p-value and either +1, -1, or 0 to indicate
  up/down/same-regulation (i.e. -1 if the expression 
  ratio of tags in the x-axis library(s) is greater 
  than that of the y-axis library(s)).  This flag
  is FALSE by default.

B<Returns>

  The p-value for the observation.  A lower number is
  more significant.  Typically, observations with
  p <= 0.05 are considered statistically significant.

  If $signedValue is set to TRUE, the function also
  returns a 0, -1 or +1 to indicate same/down/up-regulation.

B<Usage>

  # the function is static, so it can be accessed directly
  my $p = calc_audic( 3, 10, 50000, 60000 );

  # or:
  my ( $p, $sign ) = calc_audic( 3, 10, 50000, 60000, 1 );
  if( $p <= 0.05 ) {
    if( $sign == +1 ) { print "Significantly upregulated.\n"; }
    if( $sign == -1 ) { print "Significantly downregulated.\n"; }
    if( $sign == 0 ) { die( "Same expression should never be significant!" ); }
  }

=cut
   my $x = shift;
    if( !defined( $x ) ) { die( " calc_audic no arguments provided." ); }
    my $y = shift;
    my $M = shift; # cf n1
    my $N = shift; # cf n2
    my $bSign = shift || 0;

    my $p = $M / ( $M+$N );

    my $thisproba = __betai( $x+1, $y+1, $p );
    my $thisproba2 = __betai( $y+1, $x+1, 1.0-$p );

    if( $bSign >= 1 ) {
        my $ratio1 = $x / $M;
        my $ratio2 = $y / $N;
        my $sign = 0;
        if( $ratio1 > $ratio2 ) { $sign = -1; }
        if( $ratio1 < $ratio2 ) { $sign = +1; }
        #return ( ( $thisproba < $thisproba2 ? ( 2*$thisproba ) : ( 2*$thisproba2 ) ), $sign );
        return ( ($thisproba < $thisproba2) ? $thisproba : $thisproba2, $sign );
    }

    # return ( $thisproba < $thisproba2 ? ( 2*$thisproba ) : ( 2*$thisproba2 ) );
    return ( $thisproba < $thisproba2) ?  $thisproba  :  $thisproba2 ;

}

###################################
# Audic and Claverie C->Perl Port #
###################################

sub __gammln {

    my $xx = shift;

    my $x = $xx;
    my $y = $xx;

    my $tmp = $x + 5.5;
    $tmp -= ( $x + 0.5 ) * log( $tmp );
    my $ser = 1.000000000190015;
    for( my $j = 0; $j <= 5; $j++ ) { $ser += $COF[$j] / ++$y; }

    return -$tmp + log( 2.5066282746310005 * $ser / $x );

}

sub __betai {
    my $a = shift;
    my $b = shift;
    my $x = shift;

    if( $x < 0.0 || $x > 1.0 ) {
        die( "Bad x in routine betai." );
    }

    my $bt;

    if( $x == 0.0 || $x == 1.0 ) { 
        $bt = 0.0; 
    } else {
        $bt = exp( __gammln( $a+$b ) - __gammln( $a ) - __gammln( $b ) + $a*log( $x ) + $b*log( 1.0-$x ) );
    }

    if( $x < ( $a+1.0 )/( $a+$b+2.0 ) ) {
        return $bt * __betacf( $a, $b, $x ) / $a;
    }

    return 1.0 - $bt * __betacf( $b, $a, 1.0-$x ) / $b;

}

sub __fabs {
    my $x = shift;
    return ( $x < 0 ? -$x : $x );
}

sub __betacf {

    my $a = shift;
    my $b = shift;
    my $x = shift;

    my $qab = $a + $b;
    my $qap = $a + 1.0;
    my $qam = $a - 1.0;
    my $c = 1.0;
    my $d = 1.0 - $qab * $x / $qap;

    if( __fabs( $d ) < $FPMIN ) { $d = $FPMIN; }

    $d = 1.0 / $d; # inverse d
    my $h = $d;
    my $m;
    for( $m = 1; $m <= $MAXIT; $m++ ) {
        my $m2 = 2 * $m;
        my $aa = $m * ( $b-$m ) * $x / ( ( $qam + $m2 ) * ( $a + $m2 ) );

        $d = 1.0 + $aa*$d;
        if( __fabs( $d ) < $FPMIN ) { $d = $FPMIN; }

        $c = 1.0 + $aa/$c;
        if( __fabs( $c ) < $FPMIN ) { $c = $FPMIN; }

        $d = 1.0 / $d;  # inverse d

        $h *= $d*$c;

        $aa = -($a+$m)*($qab+$m)*$x/(($a+$m2)*($qap+$m2));

        $d = 1.0 + $aa * $d;
        if( __fabs( $d ) < $FPMIN ) { $d = $FPMIN; }

        $c = 1.0 + $aa / $c;
        if( __fabs( $c ) < $FPMIN ) { $c = $FPMIN; }

        $d = 1.0 / $d; # inverse d;

        my $del = $d*$c;
        $h *= $del;
        if( __fabs( $del-1.0 ) < $EPS ) { last; }
    }

    if( $m > $MAXIT ) {
       # die( "a or b too big, or MAXIT too small in __betacf" );
       print STDERR "a($a) or b($b) too big, or MAXIT($MAXIT) too small in __betacf!\n";
    }
    return $h;
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
sub ProbCTable {
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

# Lower tail quantile for standard normal distribution function.
#
# This function returns an approximation of the inverse cumulative
# standard normal distribution function.  I.e., given P, it returns
# an approximation to the X satisfying P = Pr{Z <= X} where Z is a
# random variable from the standard normal distribution.
#
# The algorithm uses a minimax approximation by rational functions
# and the result has a relative error whose absolute value is less
# than 1.15e-9.
#
# Author:      Peter J. Acklam
# Time-stamp:  2000-07-19 18:26:14
# E-mail:      pjacklam@online.no
# WWW URL:     http://home.online.no/~pjacklam
sub ltqnorm  {
  my $p = shift;
  die "input argument must be in (0,1)\n" unless 0 < $p && $p < 1;

  # Coefficients in rational approximations.
  my @a = (-3.969683028665376e+01,  2.209460984245205e+02,
           -2.759285104469687e+02,  1.383577518672690e+02,
           -3.066479806614716e+01,  2.506628277459239e+00);
  my @b = (-5.447609879822406e+01,  1.615858368580409e+02,
           -1.556989798598866e+02,  6.680131188771972e+01,
           -1.328068155288572e+01 );
  my @c = (-7.784894002430293e-03, -3.223964580411365e-01,
           -2.400758277161838e+00, -2.549732539343734e+00,
            4.374664141464968e+00,  2.938163982698783e+00);
  my @d = ( 7.784695709041462e-03,  3.224671290700398e-01,
            2.445134137142996e+00,  3.754408661907416e+00);

  # Define break-points.
  my $plow  = 0.02425;
  my $phigh = 1 - $plow;

  # Rational approximation for lower region:
  if ( $p < $plow ) {
    my $q  = sqrt(-2*log($p));
    return ((((($c[0]*$q+$c[1])*$q+$c[2])*$q+$c[3])*$q+$c[4])*$q+$c[5]) /
           (((($d[0]*$q+$d[1])*$q+$d[2])*$q+$d[3])*$q+1);
  }

  # Rational approximation for upper region:
  if ( $phigh < $p ) {
    my $q  = sqrt(-2*log(1-$p));
    return -((((($c[0]*$q+$c[1])*$q+$c[2])*$q+$c[3])*$q+$c[4])*$q+$c[5]) /
            (((($d[0]*$q+$d[1])*$q+$d[2])*$q+$d[3])*$q+1);
  }

  # Rational approximation for central region:
  my $q = $p - 0.5;
  my $r = $q*$q;
  return ((((($a[0]*$r+$a[1])*$r+$a[2])*$r+$a[3])*$r+$a[4])*$r+$a[5])*$q /
         ((((($b[0]*$r+$b[1])*$r+$b[2])*$r+$b[3])*$r+$b[4])*$r+1);
}

