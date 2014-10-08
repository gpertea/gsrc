#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 cuffdiff_loc_check.pl [-c <min_cov>] [-f <min_fpkm>] [-l <labels>] [-b <bamfiles>] \
   <gene_exp.diff> <cuffcompare.loci> <cuffcompare.tracking> <transfrags1.ifa> <transfrags2.ifa>
 
 (assume that transfrags?.ifa files were created using gff2iit -XTA )
 If BAM files are provided for specificity validation, they must be indexed.
/;
umask 0002;
getopts('b:f:c:l:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $lnames=$Getopt::Std::opt_l;
my $bamfiles=$Getopt::Std::opt_b;
my $MIN_COV=$Getopt::Std::opt_c || 2;
my $MIN_FPKM=$Getopt::Std::opt_f || 2;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
  
my @labels=split(/[\,\s\;]/,$lnames) if $lnames;
my @bam=split(/[\,\s\;]/,$bamfiles) if $bamfiles;
foreach my $bfile (@bam) {
 die("Error: file $bfile not found!\n") unless -f $bfile;
 die("Error: bam index $bfile.bai not found!\n") unless -f $bfile.'.bai';
 }
unless ($labels[1]) {
  @labels=('Exp1','Exp2');
  }
# --
my ($gfexp, $lfile, $kfile, $tfile1, $tfile2)=@ARGV;
die($usage) unless $tfile2;
#                           0       1      2      3       4       5          6
my %locdiff; # XLOC_id => [chr, c_start, c_end, fpkm_1, fpkm_2, p-value, is_significant]
             # cuffdiff loci data

          #               0      1       2     3    4        5       6
my %loci; # locus_id => [chr, strand, start, end, [refs], [@ts1], [@ts2]]
my %t1; # t_id => [$loc, ref, code, FPKM, cov, covlen, $exd]
my %t2; # t_id => [$loc, ref, code, FPKM, cov, covlen, $exd]
        #           0    1      2    3      4      5     6
open(DFILE, $gfexp) || die("Error opening file $gfexp\n");

while (<DFILE>) {
 next unless m/^XLOC/;
 chomp;
 my @t=split(/\t/);
 #next unless $t[-1] eq 'yes';
 my $is_significant = ($t[-1] eq 'yes');
 my ($chr, $cstart, $cend)=($t[2]=~m/^([^\:]+):(\d+)\-(\d+)/);
 ($cstart, $cend)=($cend,$cstart) if $cend<$cstart;
 die("Error: no chromosome coordinates parsed ($_)!\n") unless $cstart && $cend;
 $locdiff{$t[0]}=[$chr, $cstart, $cend, $t[6], $t[7], $t[10], $is_significant];
 }
close(DFILE);
die("Error: no gene diff data loaded from $lfile\n") unless keys(%locdiff)>0;
print STDERR "Loading $gfexp file..\n";

open(LFILE, $lfile) || die("Error opening file $lfile!\n");
print STDERR "Loading loci file..\n";
while (<LFILE>) {
  chomp;
  my @t=split(/\t/);
  my ($chr, $strand, $cstart, $cend)=
    ($t[1]=~m/^([^\[]+)\[([\+\-\.])\](\d+)\-(\d+)/);
  my @refs=split(/\,/,$t[2]) if $t[2] ne '-';
  my @ts1=split(/\,/,$t[3]) if $t[3] ne '-';
  my @ts2=split(/\,/,$t[4]) if $t[4] ne '-';
  $loci{$t[0]}=[$chr, $strand, $cstart, $cend, [@refs], [@ts1], [@ts2]];
  }
close(LFILE);
print STDERR "Loading tracking file..\n";
open(KFILE, $kfile) || die("Error opening file $kfile!\n");
while (<KFILE>) {
  chomp;
  my @t=split(/\t/);
  #next unless ($t[2] eq '-' && $t[3] ne 'r') || $t[3]=~tr/.p\-ucj=eoi//; 
  next if $t[3]=~tr/xs//; # tr/xsr//
  $t[2]='' if $t[2]=~m/^\s*\-\s*$/; # reference transcript
  unless ($t[4] eq '-') {
     my @q=split(/\|/,$t[4]);
     $t1{$q[1]}=[$t[1], $t[2], $t[3], $q[3], $q[6], 0, []];
     }
  unless ($t[5] eq '-') {
     my @q=split(/\|/,$t[5]);
     $t2{$q[1]}=[$t[1], $t[2], $t[3], $q[3], $q[6], 0, []];
     }
  if ($t[2]) {
    # update ref list in loci as sometimes it misses some o: and i: genes (!)
    my $ld=$loci{$t[1]};
    die("Error: locus $t[1] not found in the loci file but found in tracking file!\n") unless $ld;
    my %g;
    if (@{$ld->[4]}>0) {
       @g{@{$ld->[4]}}=();
       }
    $g{$t[2]}=1;
    $ld->[4]=[(keys(%g))];
    }
  }
close(KFILE);

print STDERR "Loading transcript file 1 ..\n";
open(TFILE, $tfile1) || die("Error opening file $tfile1!\n");
my $t; #current transcript
while (<TFILE>) {
  chomp;
  if (m/^>(\S+)\s+([^:]+):(\d+)\.\.(\d+)/) {
    my ($tid, $chr, $cstart, $cend)=($1,$2,$3,$4);
    my $td=$t1{$tid};
    next unless $td;
    my ($covlen, $numexons, @ex);
    while (<TFILE>) {
      chomp;
      if (m/^\d+\-/) {
        @ex=map { [split(/\-/)] } (split(/\,/));
        $numexons=scalar(@ex);
        map { $covlen+=$_->[1]-$_->[0]+1 } @ex;
        next;
        }
      if (m/^i:/) {
        my ($fpkm)=(m/FPKM=([^;]+)/);
        my ($cov)=(m/cov=([^;]+)/);
        ($$td[3], $$td[4], $$td[5], $$td[6])=
         ($fpkm,   $cov,   $covlen, [@ex]);
        # -- debug: 
        #print STDERR "..saving ".scalar(@ex)." exons \n";
        #my ($loc, $ref, $code, $fpk, $cc, $covl, $x)=@$td;
        #my $numex=scalar(@$x);
        #die "$tid: $loc $ref\[$code\] cov=$cc covlen=$covl numexons=$numex\n";
        #--
        last;
        }
      } #while record body
    } #header found
  }
close(TFILE);
print STDERR "Loading transcript file 2 ..\n";
open(TFILE, $tfile2) || die("Error opening file $tfile2!\n");
while (<TFILE>) {
  if (m/^>(\S+)\s+([^:]+):(\d+)\.\.(\d+)/) {
    my ($tid, $chr, $cstart, $cend)=($1,$2,$3,$4);
    my $td=$t2{$tid};
    next unless $td;
    my ($covlen, $numexons, @ex);
    while (<TFILE>) {
      chomp;
      if (m/^\d+\-/) {
        @ex=map { [split(/\-/)] } (split(/\,/));
        $numexons=scalar(@ex);
        map { $covlen+=$_->[1]-$_[0]+1 } @ex;
        next;
        }
      if (m/^i:/) {
        my ($fpkm)=(m/FPKM=([^;]+)/);
        my ($cov)=(m/cov=([^;]+)/);
        ($$td[3], $$td[4], $$td[5], $$td[6])=
         ($fpkm,   $cov,   $covlen, [@ex]);
        last;
        }
      } #while record body
    } #header found
  }
close(TFILE);
print STDERR "Now compute coverage and print results.\n";
my $maxcoverage=0;
my $max_fpkm=0;
my ($maxcov_loc, $maxfpkm_loc);
# now finally check each locus
while (my ($l,$ld)=each(%loci)) {
   my ($chr, $strand, $cstart, $cend, $rrefs, $rt1, $rt2)=@$ld;
   my @refs=@$rrefs;
   # check if this "prediction" is confirmed by cuffdiff
   my $tlist1 = (@$rt1>0) ? join(',',@$rt1) : '-';
   my $tlist2 = (@$rt2>0) ? join(',',@$rt2) : '-';
   #merge all exons of @ts1
   #my $debug=("$chr:$cstart-$cend" eq 'chr3:41288638-41289701');
   my ($mcov1, $mcovlen1, $maxcovlen1, $maxexons1, $maxfpkm1, $dcode1)=getAvgCov($rt1, \%t1);
   my ($mcov2, $mcovlen2, $maxcovlen2, $maxexons2, $maxfpkm2, $dcode2)=getAvgCov($rt2, \%t2);
   ($maxcoverage, $maxcov_loc)=($mcov1,$l) if $mcov1>$maxcoverage;
   ($maxcoverage, $maxcov_loc)=($mcov2,$l) if $mcov2>$maxcoverage;
   ($max_fpkm, $maxfpkm_loc)=($maxfpkm1,$l) if $maxfpkm1>$max_fpkm;
   ($max_fpkm, $maxfpkm_loc)=($maxfpkm2,$l) if $maxfpkm2>$max_fpkm;
   next if ($maxcovlen1==0 && $maxcovlen2==0);
   my $dcode=$dcode1.$dcode2;
   my $refcode;
   if (length($dcode)==2 && $dcode1 eq $dcode2 && $dcode1 ne '.') {
        $refcode=$dcode1.':'; # prefix for reference column
        }
       elsif (length($dcode)==1 && $dcode ne '.') {
        $refcode=$dcode.':';
        }
   #print STDERR "$l has uniform refcode: $refcode\n" if $refcode;
   my $reflist = (@refs>0) ? join(',',@refs) : '-';
   #-- discard single exon "novel" genes for now 
   #  (though they could be some relevant non-coding RNAs..)
   next if ($reflist eq '-' || $strand eq '.') && $maxexons1<2 && $maxexons2<2;
   my $status;
   my $pnovel=($reflist eq '-' || $refcode eq 'i:' || $refcode eq 'p:'); # has potential to be a novel gene
   my $isnovel=0;
   #my ($fpkm1, $fpkm2, $pval)=('-','-','-');
   #if ($pnovel && ($mcov1+$mcov2>=$MIN_COV)) {
   #   $status='novel';
   #   $isnovel=1;
   #   }

   # if not "novel", discard it unless
   # cuffdiff found a significant change in expression level
   # in a region that is consistent with this locus
   my $diffd=$locdiff{$l};
   my ($dchr, $dstart, $dend, $fpkm1, $fpkm2, $pval, $is_significant);
   unless ($diffd) {
      #print STDERR "Warning: $l not found in $gfexp file, skipping..\n";
      #next;
      ($fpkm1, $fpkm2)=($maxfpkm1, $maxfpkm2); # may not be reliable..
      next if ($fpkm1+$fpkm2<$MIN_FPKM && $mcov1+$mcov2<$MIN_COV);
      if ($pnovel) {
          $status='novel';
          $isnovel=1;
          }
      }
    else {
    #next unless ($isnovel || $diffd); 
    #if ($diffd) {
     ($dchr, $dstart, $dend, $fpkm1, $fpkm2, $pval, $is_significant)=@$diffd;
       #($fpkm1, $fpkm2, $pval)=($v1,$v2, $dp);
     if ($chr ne $dchr) {
         print STDERR "Warning: $l not found on the same chromosome in $gfexp file?\n";
         next;
         }
     if ($pnovel && ($is_significant || ($fpkm1+$fpkm2>=$MIN_FPKM && $mcov1+$mcov2>=$MIN_COV))) {
      $status='novel';
      $isnovel=1;
      }
     next unless ($isnovel || $is_significant);
     my $loclen=($cend-$cstart+1);
     my $dloclen=($dend-$dstart+1);
     # check if the quantification interval can be trusted
     if (!$isnovel && ($dloclen-$loclen>400 || $dend-$cend>200 || $cstart-$dstart>200)) {
         # bad quantification interval
         $is_significant=0;
         #next if ($dloclen-$loclen>400 || $dend-$cend>200 || $cstart-$dstart>200);
         }
     }
   my $libtag=$isnovel ? '-novel' : '-specific';
   if ($tlist1 eq '-' && $fpkm1==0 && ($is_significant || ($fpkm2>=$MIN_FPKM && $mcov2>=$MIN_COV))) {
        if ($bamfiles && !bamHasOvl($chr, $cstart, $cend, $rt2, \%t2, $bam[0])) {
            $status=$labels[1].$libtag;
            }
        }
       elsif ($tlist2 eq '-' && $fpkm2==0 && ($is_significant || ($fpkm1>=$MIN_FPKM && $mcov1>=$MIN_COV))) {
        if ($bamfiles && !bamHasOvl($chr, $cstart, $cend, $rt1, \%t1, $bam[1])) {
           $status=$labels[0].$libtag;
           }
        }
   unless ($status) {
      next if $refcode eq 'o:'; #discard dubious overlaps
      if ($mcov1+$mcov2>0 && $fpkm1+$fpkm2>0) {
        if ($is_significant) {
          $status = ($fpkm1>$fpkm2) ? $labels[0].'-inc_cdiff' : $labels[1].'-inc_cdiff';
          }
        #else {
          my $covratio=signedRatio($mcov1, $mcov2); #decrease ratio
          my $covlenratio=signedRatio($mcovlen1,$mcovlen2); #decrease ratio, negative if there is an increase
          if (($covlenratio>=2 && $covratio>=3) || (abs($covlenratio)<2 && $covratio>4)) {
            updStat(\$status, $labels[0],'inc_cov') if $fpkm1>$fpkm2;
            }
          elsif (($covlenratio<=-2 && $covratio<=-3) || (abs($covlenratio)<2 && $covratio<-4)) {
            updStat(\$status, $labels[1],'inc_cov') if $fpkm2>$fpkm1;
            }
          elsif (abs($covlenratio)<2) {
            if ($fpkm1>2*$fpkm2 && $mcov1>$mcov2*2) { updStat(\$status, $labels[0],'inc_cov') }
             elsif ($fpkm2>2*$fpkm1 && $mcov2>$mcov1*2) { updStat(\$status, $labels[1],'inc_cov')}
            }
        #  }
        }
      }
   next unless $status;
   
   print join("\t",$chr, $strand, $cstart, $cend, "$chr\:$cstart-$cend", $l, $status,
           $refcode.$reflist, join('|', $mcov1, $maxexons1, $maxfpkm1, $tlist1), 
                    join('|', $mcov2, $maxexons2, $maxfpkm2, $tlist2), $fpkm1, $fpkm2, $pval)."\n";
   } # for each locus
print STDERR "Max. coverage found: $maxcoverage for locus $maxcov_loc\n";
print STDERR "Max. FPKM found: $max_fpkm for locus $maxfpkm_loc\n";
# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************
sub updStat {
 my ($rstatus, $lib, $tag)=@_;
 if (length($$rstatus)==0) {
    $$rstatus=$lib.'-'.$tag;
    }
   else {
    $$rstatus.='-'.$tag;
    }
}

sub getAvgCov {
 my ($tlst, $th)=@_;
 my $maxcovlen=0;
 my $tcov=0;
 my $tcovlen=0;
 my $tcount=0;
 my $maxexons=0;
 my %codes;
 my $domcode='';
 my $maxcov=0;
 my $maxfpkm=0;
 foreach my $t (@$tlst) {
   my $thd=$th->{$t};
   next unless ($thd);
   my ($loc, $ref, $code, $fpkm, $cov, $covlen, $ex)=@$thd;
   $domcode='.' if $ref;
   $codes{$code}++ unless $code=~tr/u.\-r//;
   my $numexons=scalar(@$ex);
   #die "$t: $loc $ref\[$code\] covlen=$covlen, exons=$numexons\n";
   
   $maxcovlen=$covlen if ($maxcovlen<$covlen);
   $maxexons=$numexons if $maxexons<$numexons;
   $tcovlen+=$covlen;
   $tcount++;
   $maxcov=$cov if ($cov>$maxcov);
   $maxfpkm=$fpkm if $fpkm>$maxfpkm;
   $tcov+=$cov;
   }
 my @classes=keys(%codes);
 $domcode=$classes[0] if (@classes==1);
 #my ($mcov, $mcovlen)=($tcount==0) ? ( 0, 0 ) :
 #       (sprintf('%.3f', $tcov/$tcount), sprintf('%d', $tcovlen/$tcount));
 my ($mcov, $mcovlen)=($tcount==0) ? ( 0, 0 ) : (sprintf('%.2f',$maxcov), $maxcovlen);
 return ($mcov, $mcovlen, $maxcovlen, $maxexons, sprintf('%.2f',$maxfpkm), $domcode);
}

sub signedRatio { 
 # decrease ratio from $a to $b (negative if there is an increase)
 my ($a,$b)=@_;
 if ($a>$b) { return ($b==0)?1000 : $a/$b }
   elsif ($a==$b) {
     return ($a==0) ? 0 : 1
     }
     else { return ($a==0) ? -1000 : -$b/$a }
}

sub bamHasOvl {
 my ($chr, $cstart, $cend, $tlst, $th, $fbam)=@_;
 
 my $range=$chr.':'.$cstart.'-'.$cend;
 #my $debug=($range eq 'chr3:41288638-41289701');
 my @exons;
 foreach my $t (@$tlst) {
   my $thd=$th->{$t};
   next unless $thd;
   my ($loc, $ref, $code, $fpkm, $cov, $covlen, $ex)=@$thd;
   push(@exons, @$ex);
   } #gather exons from all transcripts in this locus
 return 0 unless @exons>0;
 my @sorted_exons= sort { $a->[0]<=>$b->[0] } @exons;
 open(SAMPIPE, "samtools view $fbam $range |") || die ("Error opening samtools pipe ($!)\n");
 while(<SAMPIPE>) {
   my $samline=$_;
   chomp;
   my ($qname, $flags, $gseq, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual, @extra)=
      split(/\t/);
  #    if ($strand) {
  #       my $mstrand= (($flags & 0x10)==0) ? '+' : '-';
  #       next if $mstrand ne $strand;
  #       }
   #now extract the CIGAR segments
   my @cigdata=($cigar=~m/(\d+[A-Z,=])/g);
   my ($mstart, $mend);
   my $hasOvl=0;
   my $curpos=$pos;
   $mstart=$pos;
   foreach my $cd (@cigdata) {
     my $code=chop($cd);
     if ($code eq 'N') { #gap
        #process previous interval
        if ($mend && checkOverlap($mstart, $mend, \@sorted_exons)) {
           $hasOvl=1;
           last;
           }
        $mstart=$curpos+$cd;
        $mend=undef;
        next;
        }
     $mend=$curpos+$cd-1;
     $curpos+=$cd;
     }
   unless ($hasOvl) { #check the last interval
     if ($mend && checkOverlap($mstart, $mend, \@sorted_exons)) {
         $hasOvl=1;
         }
     }
   if ($hasOvl) {
     close(SAMPIPE);
     return 1;
     }
   } # while <SAMPIPE>
  close(SAMPIPE);
  return 0; 
}



sub checkOverlap {
 my ($a, $b, $rx)=@_;
 return 0 if ($a>$$rx[-1]->[1] || $b<$$rx[0]->[0]); # not overlapping the whole exon chain
 foreach my $x (@$rx) {
   #return (start<=d->end && end>=d->start)
   return 1 if ($a<=$$x[1] && $b>=$$x[0]);
   return 0 if $b<$$x[0]; #@$rx is sorted
   }
}
