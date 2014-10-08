#!/usr/bin/perl
use strict;
use Getopt::Std;

my $usage = q/Perform vector\/contaminant cleaning by trimming or discarding 
reads based on local alignments (assumes Bowtie was run using --local against a 
vector reference database)

  Usage:
  bowtie2 --local .. | sam_vtrim.pl [-r <minreadlen>] [-p <minpid>] [-l <minovl>[%] ] \
     [-o <clean&trim.fq>] [-i <trim_info.txt>] [-c <max_clip>] [-v <max_ovhang>] \

Options:
 -i also write a "cleaning" info\/report as <trim_info.txt>
 -l ignore alignments of length below <minovl> bases (or shorter than this
    percent of read length if the '%' character suffix is used)
 -o the clean & trimmed reads will be written in this file (default: stdout)
 -p ignore alignments with percent identity lower than <minpid>
 -r discard reads shorter than <minreadlen> (default 30)
 -c ignore alignments with soft clipping greater than <max_clip>
 -v ignore alignemts with unmatched overhangs longer than <max_ovhang> (like <max_clip> 
    but without counting bases falling outside the reference boundaries)
/;
umask 0002;
getopts('o:i:l:o:p:r:c:') || die($usage."\n");
#my $outfile=$Getopt::Std::opt_o;
my $minlen=$Getopt::Std::opt_r || 30;

my $minovl=$Getopt::Std::opt_l; #could be percentage if it ends with '%'
my $movlperc=0; #percentage was used ?
$movlperc=1 if ($minovl=~s/\%$//);
my $maxclip=$Getopt::Std::opt_c;
my $maxovh=$Getopt::Std::opt_v;
my $minpid=$Getopt::Std::opt_p; #computed from edit distance
my $outClean=$Getopt::Std::opt_o; # trimming/cleaning functionality
my $clnInfo=$Getopt::Std::opt_i;
#my $keepUnmapped=$Getopt::Std::opt_U;

if ($outClean) {
  open(OUTF, '>'.$outClean) || die("Error creating output file $outClean\n");
  select(OUTF);
}
if ($clnInfo) {
  open(CLNFO, '>'.$clnInfo) || die("Error creating output file $clnInfo\n");
}

# --
my %refs; # refname => reflen
my ($qname, $hitcount, $qSeq, $qQual); #data for current read
my @qTrims; # for -c, list of [trimSeq, trimQuals, hit_info] mappings 
              # that passed the filters for the current read and qualifies for trimming/cleaning
my $qTrash;
my $prev_alnscore;

while (<>) {
  my $line=$_;
  #check for header:
  if (m/^@[A-Z][A-Z]\t/) {
    print $_;
    #keep refseq length
    if (m/^\@SQ\tSN:(\S+)\tLN:(\d+)/) {
     $refs{$1}=$2;
     }
    next;
  }
 chomp;
 my ($rname, $flags, $refname, $pos, $mapq, $cigarstr, 
     $rnext, $pnext, $tlen, $seq, $quals, $tags)=split(/\t/, $_, 12);
 next unless length($tags)>5;
 $flags=int($flags);
 my $sflag = $flags & 0xc0;
 my $isrev = ($flags & 0x10) != 0;
 my $materev = ($flags & 0x20) != 0;
 if ($rname ne $qname) {
  $hitcount=0;
  $prev_alnscore=0;
  if ($qname) {
    writeClean($qname, $qTrash, \@qTrims);
    @qTrims=();
    }
  $qname = $rname;
  ($qSeq, $qQual) = $isrev ? (reverseComplement($seq), reverse($quals)) : ($seq, $quals);
  $qTrash='';
 }

 my $frag_idx=0;
 if ($sflag == 0x40) {
    $rname.='/1';
    $frag_idx=1;
 } elsif ($sflag == 0x80) {
    $rname.='/2';
    $frag_idx=2;
 }

 my $unmapped = (($flags & 4)!=0);

 if ($unmapped) {
   #writeClean() will write it
   #print $line;
   push(@qTrims, [$qSeq, $qQual, '']); #passes as is
   next;
 }
 
 my ($edist) = ( $tags=~m/\bNM:i:(\d+)/ );
 if ($edist==0 && $tags=~m/\bnM:i:(\d+)/) {
   $edist=$1;
   $edist+=$_ foreach ($cigarstr=~m/(\d+)[ID]/g);
 }
 
 my $rlen=length($qSeq); #too short?
 if ($rlen<$minlen) {
   push(@qTrims, [$qSeq, $qQual, '']); #trashed
   $qTrash='L';
   next;
 }
 
 my @cigar=($cigarstr=~m/(\d+[A-Z])/g);
 my $refovl=0;
 
 my $strand=$isrev ? '-' : '+';

 my ($clipL, $clipR)=(0,0);
 $clipL=int($cigar[0]) if (substr($cigar[0],-1) eq 'S');
 $clipR=int($cigar[-1]) if (substr($cigar[-1],-1) eq 'S');
 # ovhL and ovhR - L and R refer to the location of aln on *reference*
 # these will be adjusted later for alignments extending beyond ref boundary
 my ($ovhL, $ovhR) = ($clipL, $clipR); 
 #but clipL and clipR to refer to the fwd strand of the *read*
 ($clipL, $clipR) = ($clipR, $clipL) if $isrev;

 chomp($line);

 my $reflen=$refs{$refname} 
   || die("Error: couldn't get reference length for $refname\n$line\n");

 $refovl+=$_ foreach (grep(/\d+[NDXM=]/, @cigar) ); #aln length on ref
 my $ovlen = $rlen-$clipL-$clipR; #aln length on read
 # Note: 0 based coords for rL, rR, refL and refR !
 # rL and rR are aln coordinates on the forward strand of the read
 my $rL = $clipL;
 my $rR = $rL+$ovlen-1;
 #refL - refR :alignment coords on the reference (fwd strand)
 my $refL = $pos-1;
 my $refR = $refL + $refovl-1;
 my $ovperc=($ovlen*100.00)/$rlen;
 my $sovperc=sprintf("%.2f",$ovperc);
 $sovperc=~s/\.?0+$//;
 my ($pseq, $pquals)= ($clipL>$clipR) ? (substr($qSeq, 0, $clipL),substr($qQual, 0, $clipL)) : 
                                    ( $clipR ? (substr($qSeq, -$clipR), substr($qQual, -$clipR)) :
                                      ('', '') ) ;


 my $pid=(($ovlen-$edist)*100.00)/$ovlen;
 

 push(@qTrims, [$pseq, $pquals, 
      sprintf('%d-%d(%s):%s(%d-%d)|%.2f',$rL+1, $rR+1, $strand, $refname, $refL+1, $refR+1, $pid)]);

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
 if ($maxovh) {
   #ovhL/R are really non-matching regions at the end of the read
   #exclude parts that fall off the edge of the reference
   $ovhL=$refL if ($ovhL>$refL);
   $ovhR = $reflen-$refR if ($ovhR>$reflen-$refR);
   
   next if ($ovhR>$maxovh || $ovhL>$maxovh);
 }
 
 #if ($outClean && !$qTrash) {
   #see how can we trim this read based on this alignment
 $qTrash='V' if ($ovperc>=80) || (length($pseq)<$minlen);
   #next;
 #}
} #for each SAM line
# --
if ($qname) {
    writeClean($qname, $qTrash, \@qTrims);
    @qTrims=();
    }
if ($outClean) {
  select(STDOUT);
  close(OUTF);
}
if ($clnInfo) {
  close(CLNFO);
}


#************ Subroutines **************

sub reverseComplement {
  my $s=reverse($_[0]);
  $s =~ tr/AaCcTtGgUuMmRrWwSsYyKkVvHhDdBb/TtGgAaCcAaKkYyWwSsRrMmBbDdHhVv/;
  return $s;
 }

sub writeClean {
 my ($rname, $trash, $rd)=@_;
 my $tcode = $trash || '.';
 die("Error at processing $rname, no trim data found!\n") unless @$rd>0;
 my @rdata;
 if (@$rd > 1) {
    @rdata = sort { length($main::a->[0])<=>length($main::b->[0]) } @$rd;
 }
 else {
   @rdata = @$rd;
 }
 my ($seq, $quals, $info)=@{$rdata[0]};
 if ($clnInfo && ($trash || $info)) {
    $info='.' unless $info;
    print CLNFO join("\t",$rname, $tcode, length($seq), $info)."\n";
 }
 if ($trash) {
   return;
 }
 print '@'.$rname.' '.$info."\n";
 print "$seq\n+\n$quals\n";
}
