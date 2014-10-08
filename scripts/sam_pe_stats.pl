#!/usr/bin/perl
use strict;
use Getopt::Std;

my $usage = q/Report SAM stats for PE reads. Usage:
 sam_pe_stats.pl [-r <ctg_regions_file>]  [-f <min_frag_size>] [-p <minpid>] [-l <minovl>[%] ] \
  [-m <maxclip>] [-o <stats.txt>] < file.sam
 
  Options:
   -N                 : no coverage calculation, assumes input SAM is sorted by read name
   -m <maxclip>       : ignore any alignments with soft clipping larger than <maxclip>
   -l <minovl>[%]     : ignore any alignments where matching region gets below <minovl>
                       (if % is given, this is min. percentage from read length)
   -p <minpid>        : discard alignments with percent identity below <minpid>
   -f <min_frag_size> : discard PE alignments where reads are closer than <min_frag_size>
   -v <max_ovhang>    : discard alignments with overhangs larger than <max_ovhang>

/;
umask 0002;
getopts('Nf:v:r:o:l:p:m:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $minovl=$Getopt::Std::opt_l; #could be percentage if it ends with '%'
my $minfraglen=$Getopt::Std::opt_f;
my $regfile=$Getopt::Std::opt_r;
my $movlperc=0; #percentage was used
$movlperc=1 if ($minovl=~s/\%$//);
my $maxclip=$Getopt::Std::opt_m;
my $maxovh=$Getopt::Std::opt_v;
my $minpid=$Getopt::Std::opt_p; #computed from edit distance
my $NoCov = $Getopt::Std::opt_C;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}

# --
my @reflst;
my %refs; # refname => [ reflen, [@regions], [ base-cov ], [ ins-cov] ]
#my ($qname, $hitcount, $qSeq, $qQual); #data for current read
if ($regfile) {
 open(RF, $regfile) || die("Error opening ctg region file $regfile\n");
 while (<RF>) {
   chomp;
   my @t=split(/\t/);
   #my @r=split(/\,/, $t[3]);
   $t[3]=~s/ \[[^\]]+\]//g;
   $t[3]=~s/ \(\d+\)//g;
   #print STDERR ">$t[0]\t$t[3]\n";
   my @regs =  map { [ (split(/\-/)) ] } split(/\,/, $t[3]);
   checkRegs($t[0], $t[2], \@regs);
   #                 len   merges     base cov     insert cov  
   $refs{$t[0]} = [$t[2], [@regs], [ (0)x$t[2] ], [ (0)x$t[2] ]  ];
 }
 close(RF);
}
#####
#exit(1); #debug
#####
my ($tname, $qname, $prevtname, $hitcount, $qSeq);
           #              0        1       2          3
my @pdata; #pair data: ($tname, $tlen, [@r1data], [@r2data])
           # where @rdata=($refname, $rpos, $strand, $ovlen, $pid,  $rlen, $hitindex, $nhits)
           #                   0       1       2        3      4      5       6         7
my ($fsum, $fsumcount, $favg, $favgsum, $favgcount); 
while (<>) {
  my $line=$_;
  #check for header:
  if (m/^@[A-Z][A-Z]\t/) {
    #print $_;
    #keep refseq length
    if (m/^\@SQ\tSN:(\S+)\tLN:(\d+)/) {
      my ($refname, $reflen)=($1, $2);
      my $rd=$refs{$refname};
      if ($rd) { #previously loaded info (-r )
        die("Error: mismatch between \@SQ ref data and previous info ($refname)\n")
         if ($rd->[0] != $reflen);
      }
      else { 
         if ($NoCov) { $refs{$refname}=[$reflen, [], [ ], [ ] ]; }
                else { $refs{$refname}=[$reflen, [], [ (0)x$reflen ], [ (0)x$reflen ]]; }
      }
      push(@reflst, $refname);
    }
    next;
  }
 chomp;
 my ($rname, $flags, $refname, $pos, $mapq, $cigarstr, 
     $rnext, $pnext, $tlen, $seq, $quals, $tags)=split(/\t/, $_, 12);
 next unless length($tags)>5;
 my $newpair=($rname ne $prevtname);
 if ($newpair) {
   #process pair data
   
   @pdata=();
 }
 
 $prevtname=$rname;
 $tname=$rname;
  #if ($t[2] eq '*' && 
 $flags=int($flags);
 my $sflag = $flags & 0xc0;
 my $first;
 my $firstseen;
 if ($sflag == 0x40) {
    $rname.='/1';
    $first = 1;
 } elsif ($sflag == 0x80) {
    $rname.='/2';
 }
 my $isrev = ($flags & 0x10) != 0;
 my $materev= ($flags & 0x20) != 0;

 if ($rname ne $qname) {
  $hitcount=0;
  $qname = $rname;
  #$qSeq = $isrev ? reverseComplement($seq) : $seq;
  #($qSeq, $qQual) = $isrev ? (reverseComplement($seq), reverse($quals)) : ($seq, $quals);
 }
 
 #if ($maxhits && $hitcount>=$maxhits) {
 #  next;
 #}

 my $unmapped = (($flags & 4)!=0);

 next if $unmapped;

 my ($edist) = ( $tags=~m/\bNM:i:(\d+)/ );
 if ($edist==0 && $tags=~m/\bnM:i:(\d+)/) {
   $edist=$1;
   $edist+=$_ foreach ($cigarstr=~m/(\d+)[ID]/g);
 }
 my ($numhits) = ( $tags=~m/\bNH:i:(\d+)/ );
 my ($hitindex) = ( $tags=~m/\bHI:i:(\d+)/ );
 my ($PEstatus) = ( $tags=~m/\bYT:Z:(\w\w)/ );
 my $rlen=length($qSeq); #too short?
 die("Error (sequence too short) at line:\n$line\n") if $rlen<8;
 #if ($rlen<$minlen) {
 #  next;
 #}
 
 my @cigar=($cigarstr=~m/(\d+[A-Z])/g);
 my $refovlen=0;
 my $strand=$isrev ? '-' : '+';

 my ($clipL, $clipR)=(0,0);
 $clipL=int($cigar[0]) if (substr($cigar[0],-1) eq 'S');
 $clipR=int($cigar[-1]) if (substr($cigar[-1],-1) eq 'S');
 # ovhL and ovhR - L and R refer to the location of aln on reference 
 my ($ovhL, $ovhR) = ($clipL, $clipR);
 #but make clipL and clipR to refer to the fwd strand of the read
 ($clipL, $clipR) = ($clipR, $clipL) if $isrev;

 chomp($line);
 my $refdata=$refs{$refname}
    || die("Error: couldn't get reference info for $refname\n$line\n");
 my $reflen=$refdata->[0];

 $refovlen+=$_ foreach (grep(/\d+[NDXM=]/, @cigar) ); #aln length on ref
 my $ovlen = $rlen-$clipL-$clipR; #aln length on read
 # 0 based coords !!
 # rL and rR are aln coordinates on the forward strand of the read
 my $rL = $clipL;
 my $rR = $rL+$ovlen-1;
 my $refL = $pos-1;
 my $refR = $refL + $refovlen-1;
 my $ovperc=($ovlen*100.00)/$rlen;
 
 #my ($pseq, $pquals)= ($clipL>$clipR) ? (substr($qSeq, 0, $clipL),substr($qQual, 0, $clipL)) : 
 #                                   ( $clipR ? (substr($qSeq, -$clipR), substr($qQual, -$clipR)) :
 #                                     ('', '') ) ;

 my $pid=(($ovlen-$edist)*100.00)/$ovlen;
 
 if ($minovl) {
  if ($movlperc) {     
     next if $ovperc<$minovl;
   }
   else {
    next if $ovlen<$minovl;
   }
 }
 next if $maxclip && ($clipL>$maxclip || $clipR>$maxclip);
 next if ($minpid && $pid<$minpid);
 next if $minfraglen && $tlen && abs($tlen)<$minfraglen;
 if ($maxovh) {
   $ovhL=$refL if $ovhL>$refL; # the read extends outside the left end of reference
   $ovhR = $reflen-$refR if $ovhR>$reflen-$refR; # the read extending outside the right end of ref
   next if ($ovhR>$maxovh || $ovhL>$maxovh);
 }
 #########----------
 #print $line;
 #print "\tYI:Z:".sprintf("%.2f",$pid) unless $line=~m/\bYI:Z:/;
 # print "\n";
 #########----------
 #print STDERR join("\t", $rname, $rlen, $isrev?'-':'+', $rL, $rR, $refname.'-@-'.$pos, $reflen, $refL, $refR, 
 #      $cigarstr, "tlen:$tlen", "ploc:$rnext-\@-$pnext")."\n";
 
 ##- update base coverage
 unless ($NoCov) {
   my $covinc=1/$numhits;
   map { $_+=$covinc } @{$refdata->[2]}[$pos .. $refR+1];
 }
 $hitcount++;
 #if ($outClean && !$qTrash) {
   #see how can we trim this read based on this alignment
 #$qTrash='V' if ($ovperc>=80) || (length($pseq)<$minlen);
   #next;
 #}
} #for each SAM line
if ($favgcount==0) { $favgsum=$fsum; $favgcount=$fsumcount;}
my $favglen=int($favgsum/$favgcount);
print STDERR "Average fragment len: $favglen\n";
# --
unless ($NoCov) {
  foreach my $refname (@reflst) {
    my $refdata=$refs{$refname} || die("Error retrieving data for $refname\n");

    my @regs=@{$refdata->[1]};
    my ($reflen, $bcov, $icov) = ($refdata->[0], $refdata->[2], $refdata->[3]);
    my $mbcov=avgCov($bcov, $favglen, $reflen-$favglen+1);
    my $micov=avgCov($icov, $favglen, $reflen-$favglen+1);
    print ">$refname $reflen avg_bcov=".sprintf('%.2f',$mbcov).
       ' avg_icov='.sprintf('%.2f',$micov)."\n";
    #printRegs($refdata->[1]);  
    if (@regs>0) {
      foreach my $reg (@regs) {
        my $reglen=$reg->[1]-$reg->[0]+1;
        print "$reg->[0]-$reg->[1] ($reglen)\t";
        my ($rbcov, $xbcov) = avgRCovs($reflen, $favglen, $bcov, $reg->[0], $reg->[1]);
        my ($ricov, $xicov) = avgRCovs($reflen, $favglen, $icov, $reg->[0], $reg->[1]);
        print 'bcov='.sprintf('%.2f',$rbcov).' (vs'.sprintf('%.2f',$xbcov).")\t";
        print 'icov='.sprintf('%.2f',$ricov).' (vs'.sprintf('%.2f',$xicov).') '.
             sprintf('[%.2f  %.2f]',$icov->[$reg->[0]-1], $icov->[$reg->[1]-1])."\n";
      }
    }
    else {
    # reportCovSpikes($bcov, $icov); #possible collapsed repeats
    # reportCovDrops($bcov, $icov);  #possible misassembly
    }
  }
}

if ($outfile) {
 select(STDOUT);
 close(OUTF);
}

#************ Subroutines **************

sub reverseComplement {
  my $s=reverse($_[0]);
  $s =~ tr/AaCcTtGgUuMmRrWwSsYyKkVvHhDdBb/TtGgAaCcAaKkYyWwSsRrMmBbDdHhVv/;
  return $s;
 }

sub processPair {
 #uses global @pdata;
 my ($fname, $tlen, $rd1, $rd2) = @pdata;
 return unless $rd1 && $rd2;
 my ($r1ref, $r1pos, $r1strand, $r1ovlen, $r1pid, $r1len,$r1hi, $r1nh)=@$rd1;
 my ($r2ref, $r2pos, $r2strand, $r2ovlen, $r2pid, $r2len,$r2hi, $r2nh)=@$rd2;
 #update insert coverage data
 if (!$NoCov && $tlen!=0 && $r1ref eq $r2ref) {
  my $refdata=$refs{$r1ref}
    || die("Error: couldn't get reference info for $r1ref\n");
   my $covinc=1/$r1nh;
   my $hlen=int($r1len/2);
   my ($u0,$u1)=($r2pos<$r1pos) ? ($r2pos, $r1pos) : ($r1pos, $r2pos);
   $u0+=$hlen+6; $u1+=$hlen-6;
   map { $_+=$covinc } @{$refdata->[3]}[$u0 .. $u1];
   #for (my $i=$u0;$i<$u1;$i++)
 if ($r1hi<2) {
  #account for fragment average (only for first hit)
  $fsum+=abs($tlen);
  $fsumcount++;
  if ($fsumcount>=100000) {
    $favg=$fsum/$fsumcount;
    $favgsum+=$favg;
    $favgcount++;
    $fsumcount=0;
    $fsum=0;
  }
 }
 }
}


sub checkRegs {
 my ($a, $b, $regs)=@_;
 if (@$regs<1) {
  die("Error: no merge regions found for $a ($b)\n");
 }
 my @rs;
 foreach my $r (@$regs) {
   push(@rs, $r->[0].'-'.$r->[1]);
   if ($$r[0]==0 || $$r[1]==0 || $$r[0]>=$$r[1]) {
      die("Error: invalid region ($$r[0]-$$r[1]) found for $a ($b)\n");
   }
 }
 #print STDERR "$a\t$b\t".join(',',@rs)."\n";
}

sub avgCov {
 my ($cov, $a, $b)=@_;
 my $fsum=0;
 map { $fsum+=$_ } @{$cov}[$a-1 .. $b-1];
 return $fsum/($b-$a+1);
}

sub avgRCovs {
 my ($rlen, $flen, $cov, $a, $b)=@_;
 my $x0=$a-$flen;
 $x0=$flen if $x0<$flen;
 my $x1=$b+$flen;
 $x1=$rlen-$flen if $x1>$rlen-$flen;
 my $cov0=avgCov($cov, $x0, $a-1) if ($a-$x0>100);
 my $cov1=avgCov($cov, $b+1, $x1) if ($x1-$b>100);
 my $xcov=$cov0+$cov1;
 $xcov/=2 if $cov0 && $cov1;
 my $rcov=avgCov($cov, $a, $b);
 return ($rcov, $xcov);
}
