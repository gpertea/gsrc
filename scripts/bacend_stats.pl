#!/usr/bin/perl
use strict;
use Getopt::Std;

my $usage = q/Report SAM stats for BAC end reads, per contig.
Usage:
 bacend_stats.pl [-p <minpid>] [-l <minovl>[%] ] [-m <maxclip>] [-o <stats.txt>] < file.sam
  
  Options:
   -m <maxclip>       : ignore any alignments with soft clipping larger than <maxclip>
   -l <minovl>[%]     : ignore any alignments where matching region gets below <minovl>
                       (if % is given, this is min. percentage from read length)
   -p <minpid>        : ignore any alignments with percent identity below <minpid>
   -A                 : show all pairs (default is to filter out same-contig ~100k BEs)
   -D                 : no distance check (don't show RF pairs on the same scaffold and opposite strands)
   -S                 : ignore any RF pairs on the same scaffold (no distance or strand check)

/;
umask 0002;
getopts('SDAf:r:o:l:p:m:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $discardPaired = !$Getopt::Std::opt_A;
my $noDistanceCheck = $Getopt::Std::opt_D;
my $noSameTarget = $Getopt::Std::opt_S;
($noDistanceCheck, $discardPaired)=(1,1) if $noSameTarget;
my $minovl=$Getopt::Std::opt_l; #could be percentage if it ends with '%'
#my $minfraglen=$Getopt::Std::opt_f;
#my $regfile=$Getopt::Std::opt_r;
my $movlperc=0; #percentage was used
$movlperc=1 if ($minovl=~s/\%$//);
my $maxclip=$Getopt::Std::opt_m || 200;
#my $maxovh=$Getopt::Std::opt_v;
my $minpid=$Getopt::Std::opt_p; #computed from edit distance

if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}
my %be; #BE clone => [ $fdata, $rdata ]
# where $[fr]data has this format: [ $aln1, $aln2, $aln3 .. ]
# where $alnX is alignment data in this format:
#          [$ctg, $ctgpos, $strand, $pid, $score, $clipL, $clipR]

# --
my @reflst;
my %refs; # refname => [ reflen ] ## [ base-cov ], [ ins-cov] ]
#my ($qname, $hitcount, $qSeq, $qQual); #data for current read

#####
#exit(1); #debug
#####

my ($qname, $pname, $hitcount, $qSeq);
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
      else { $refs{$refname}=[ $reflen ]; }
      push(@reflst, $refname);
    }
    next;
  }
 chomp;
 my ($rname, $flags, $refname, $pos, $mapq, $cigarstr, 
     $rnext, $pnext, $tlen, $seq, $quals, $tags)=split(/\t/, $_, 12);
 next unless length($tags)>5;
 my $clone=$rname;
 my ($s1, $s2)=('a', 'b');
 if ($clone=~s/^([^_RF]+)([RF])_/$1_/) {
  $s1=$2;
  }
 if ($clone=~s/([RF])(_[^_]+)$/$2/) {
  $s2=$1;
  }
 if ($s1 ne $s2 || $s1!~m/^[RF]/ || $s2!~m/[RF]/) {
    $s1=$s2 if $s1!~m/^[RF]/;
    print STDERR "Warning: problem determining clone ID for $rname, assuming $s1\n";
    $s2 = $s1;
    #next;
 }
    
 my $cdata_idx = ($s1 eq 'F') ? 0 : 1 ; #'F' , 'R'
 #my $newpair=($rname ne $pname);
 #$pname=$rname;
  #if ($t[2] eq '*' && 
 $flags=int($flags);
 my $sflag = $flags & 0xc0;
 #my $first;
 #if ($sflag == 0x40) {
 #   $rname.='/1';
 #   $first = 1;
 #} elsif ($sflag == 0x80) {
 #   $rname.='/2';
 #}
 my $isrev = ($flags & 0x10) != 0;
 my $isprimary = ($flags & 0x100) == 0;
 next unless $isprimary;
 if ($rname ne $qname) {
  $hitcount=0;
  $qname = $rname;
  $qSeq = $isrev ? reverseComplement($seq) : $seq;
  #($qSeq, $qQual) = $isrev ? (reverseComplement($seq), reverse($quals)) : ($seq, $quals);
 }

 #if ($maxhits && $hitcount>=$maxhits) {
 #  next;
 #}

 my $unmapped = (($flags & 4)!=0);

 next if $unmapped;

 my $edist = 0; ;
 if ($tags=~m/\bNM:i:(\d+)/) {
  $edist=$1;
  }
 elsif ($tags=~m/\bnM:i:(\d+)/) {
  $edist=$1;
  $edist+=$_ foreach ($cigarstr=~m/(\d+)[ID]/g);
 }
 my ($numhits) = ( $tags=~m/\bNH:i:(\d+)/ );
 my ($hitindex) = ( $tags=~m/\bHI:i:(\d+)/ );
 my ($score) = ( $tags=~m/AS:i:(\d+)/ );
 my $rlen=length($qSeq); #too short?
 die("Error (sequence too short) at line:\n$line\n") if $rlen<10;
 #if ($rlen<$minlen) {
 #  next;
 #}
 
 my @cigar=($cigarstr=~m/(\d+[A-Z])/g);
 my $refovl=0;
 
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

 $refovl+=$_ foreach (grep(/\d+[NDXM=]/, @cigar) ); #aln length on ref
 my $ovlen = $rlen-$clipL-$clipR; #aln length on read
 # 0 based coords !!
 # rL and rR are aln coordinates on the forward strand of the read
 my $rL = $clipL;
 my $rR = $rL+$ovlen-1;
 my $refL = $pos-1;
 my $refR = $refL + $refovl-1;
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
# next if $minfraglen && abs($tlen)<$minfraglen;
 $pid=sprintf('%.1f',$pid);
 $pid=~s/\.0$//;
 my $cdata=$be{$clone};
 if (!$cdata) {
  $be{$clone}=[ [], [] ];
  $cdata=$be{$clone};
 }
 
 push(@{$cdata->[$cdata_idx]}, 
  [$refname, $reflen, $pos, $strand, $pid, $score, $clipL, $clipR, $ovlen] );
 
 #if ($maxovh) {
 #  $ovhL=$refL if $ovhL>$refL; # the read extends outside the left end of reference
 #  $ovhR = $reflen-$refR if $ovhR>$reflen-$refR; # the read extending outside the right end of ref
 #  next if ($ovhR>$maxovh || $ovhL>$maxovh);
 #}
 #########----------
 #print $line;
 #print "\tYI:Z:".sprintf("%.2f",$pid) unless $line=~m/\bYI:Z:/;
 # print "\n";
 #########----------
 #print STDERR join("\t", $rname, $rlen, $isrev?'-':'+', $rL, $rR, $refname.'-@-'.$pos, $reflen, $refL, $refR, 
 #      $cigarstr, "tlen:$tlen", "ploc:$rnext-\@-$pnext")."\n";

 $hitcount++;
 #if ($outClean && !$qTrash) {
   #see how can we trim this read based on this alignment
 #$qTrash='V' if ($ovperc>=80) || (length($pseq)<$minlen);
   #next;
 #}
} #for each SAM line

foreach my $clone (sort keys(%be)) {
 my $cdata=$be{$clone};
 my ($falns, $ralns)=@$cdata;
 next unless (@$falns>0 && @$ralns>0);
 #collapse equivalent alignments
 $falns=collapse($falns) if @$falns>1;
 $ralns=collapse($ralns) if @$ralns>1;
 my $haveAdjacent=0;
 if ($discardPaired) {
   my $properlyPaired=0;
   foreach my $faln (@$falns) {
     foreach my $raln (@$ralns) {
       if ($$faln[0] eq $$raln[0] && 
           $$faln[3] ne $$raln[3]) {
           if ($noDistanceCheck) {
               $properlyPaired=1;
               last;
           }
           else {
             my $dist=abs($$faln[2]-$$raln[2]);
             if ($dist>30000 && $dist<300000) {
               $properlyPaired=1;
               last;
             }
           }
       }
       elsif ($$faln[0] ne $$raln[0]) {
        my ($forder)=($$faln[0]=~m/_(\d+)$/);
        my ($rorder)=($$raln[0]=~m/_(\d+)$/);
        $haveAdjacent=1 if ($forder && $rorder && 
                   abs(int($forder)-int($rorder))==1);
       }
     }
     last if $properlyPaired;
   }
   next if $properlyPaired;
 }
 if ($noSameTarget && @$falns==1 && @$ralns==1) {
   next if $$falns[0]->[0] eq $$ralns[0]->[0];
 }
 my $adj=$haveAdjacent ? ' adj' : '';
 print ">$clone\tF:".scalar(@$falns)."\tR:".scalar(@$ralns)."$adj\n";
 foreach my $a (@$falns) {
    #     0      1       2      3       4      5      6        7       8
    my ($ctg, $ctglen, $pos, $strand, $pid, $score, $clipL, $clipR, $ovlen)=@$a;
    my $refR=$pos+$ovlen-1;
    print "F\t$ctg($ctglen)\[$strand\]:$pos-$refR~$pid\% $score $clipL>$ovlen<$clipR\n";
   }
 foreach my $a (@$ralns) {
    my ($ctg, $ctglen, $pos, $strand, $pid, $score, $clipL, $clipR, $ovlen)=@$a;
    my $refR=$pos+$ovlen-1;
    print "R\t$ctg($ctglen)\[$strand\]:$pos-$refR~$pid\% $score $clipL>$ovlen<$clipR\n";
   }
}

if ($outfile) {
 select(STDOUT);
 close(OUTF);
}

#************ Subroutines **************
sub collapse {
 my $r=$_[0];
 my %h;
 for (my $i=0; $i<@$r ; ++$i) {
   #     0       1      2       3     4       5       6       7       8
   my ($ctg, $ctglen, $pos, $strand, $pid, $score, $clipL, $clipR, $ovlen)=@{$$r[$i]};
   $h{$ctg.'_'.$pos.'_'.$strand.'_'.$score}=sprintf('%06d',$i);
 }
 #my $res=[];
 my @a;
 foreach my $i (sort values(%h)) {
   push(@a, $$r[int($i)]);
 }
 # also collapse mappings on the same contig & strand in about 30k distance from each-other
 @a=sort { $a->[0].sprintf('%09d', $a->[2]) cmp $b->[0].sprintf('%09d', $b->[2]) } @a;
 my $changed=0;
 do {
  $changed=0;
  for (my $i=1;$i<@a;$i++) {
   if ($a[$i-1]->[0] eq $a[$i]->[0] && 
        $a[$i-1]->[3] eq $a[$i]->[3] &&
         abs($a[$i-1]->[2] - $a[$i]->[2])<40000) {
      #keep the one with a higher pid?
      my $d = ($a[$i-1]->[5] > $a[$i]->[5]) ? $i-1 :  
               ($a[$i-1]->[5] < $a[$i]->[5]) ? $i : 
                ($a[$i-1]->[8] > $a[$i]->[8]) ? $i-1 : $i ;
      $changed=1;
      splice(@a, $d, 1);
      last;
   }
  }
 } while ($changed);
 return [@a];
}


sub reverseComplement {
  my $s=reverse($_[0]);
  $s =~ tr/AaCcTtGgUuMmRrWwSsYyKkVvHhDdBb/TtGgAaCcAaKkYyWwSsRrMmBbDdHhVv/;
  return $s;
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
