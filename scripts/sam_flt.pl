#!/usr/bin/perl
use strict;
use POSIX qw(strftime);
use Getopt::Std;

my $usage = q/Filter SAM records by various conditions.
Usage:
  bowtie2 --reorder --local --no-unal .. | sam_flt.pl [-U] [-r <minreadlen>] \
     [-p <minpid>] [-a <minovl>[%] ] [-f <min_frag_len>[,<max_frag_len>][RF|FR]] \
     [-c <clean&trim.fq>] [-I] [-c <max_clip>] [-e <max_dist>[FR|RF]] \
     [-l <libname>:<min_len>,<max_len>[RF|FR]] [-b <bambus_links.xml>]
     [-v <max_overhang>] [-S] [-w <contig_stats.txt>] filtered.sam

Options:
 -a alignment length must be at least <minovl> bases (or this percent of
    read length if '%' character is used)
 -e for paired reads mapped on different references, discard mappings farther 
    than <max_dist> from either ends of the reference; if <max_dist> is followed 
    by FR or RF, the orientation of mates is also considered and <max_dist> is 
    measured from the end or start of the read, respectively
 -f for paired reads mapped on the same reference, discard mates mapped closer 
    than <min_frag_len> and, if given, farther than <max_frag_len>;
    if followed by FR or RF, also discards pairs mapped in the "wrong" orientation
 -S for paired reads, discard alignments if both mates aligned to the same reference
 -B for paired reads, discard single mate alignments (both mates are required
    to have valid alignments)
 -M for paired reads, only keep alignments if both mates have valid alignments on
    the same reference
 -b (requires -l) for paired reads, write Bambus evidence file with only the 
     mappings of paired reads linking different contigs
 -l Bambus library info for -b option; also accepts RF or FR suffix like the
    -f option but does not filter the alignment output for same ref mates
 -c discard alignments with soft clipping greater than <max_clip>
 -p alignments with percent identity lower than <minpid> are discarded
 -r alignments of reads shorter than <minreadlen> are discarded
    (default 30) 
 -t multiple alignments ranking up to <#hits> for each read are kept
    (if -t1 is provided, only alignments with the same best score are kept);
    this assumes alignments are already sorted by decreasing score per read;
 -U keep SAM records of unmapped reads (by default these are discarded)
 -P keep unmapped mates for paired reads with one mate mapped properly
 -v maximum length of unmatched overhangs for each read (like <max_clip> 
    but without counting bases falling outside the reference boundaries)
 -w collect and report stats about PE reads mappings on each contig
    (uses pair orientation and distance; to improve this info,
    the use of -e, -f with FR or RF suffixes is recommended)
 -X fix read names for 454 PE reads which may still have
    \/1 and \/2 or .left\/.right suffixes in their name
 IMPORTANT: the input SAM records are expected to be grouped together by read\/pair, 
   (use --reorder option of Bowtie2, or sort by read name).
/;


umask 0002;
my $cmdargs=join(' ', @ARGV);
getopts('UPSBXMa:f:b:t:r:l:p:e:c:v:w:') || die($usage."\n");

#my $outfile=$Getopt::Std::opt_o;
my $minlen=$Getopt::Std::opt_r || 30;
my $minfraglen=$Getopt::Std::opt_f;
my $fixReadNames=$Getopt::Std::opt_X;

my $frag_ori;
if ($minfraglen && $minfraglen=~s/(\d+)([FfRr][RrFf])$/$1/) {
   $frag_ori=lc($2);
   die("Error at -f option value, only RF or FR suffixes are recognized!\n")
       if ($frag_ori ne 'fr' && $frag_ori ne 'rf');
}

my $maxfraglen;
if ($minfraglen=~m/^(\d+)[\,\;\~\-\.](\d+)$/) {
 ($minfraglen, $maxfraglen)=($1,$2);
 die("Error at -f option (invalid range)!\n") if ($maxfraglen<$minfraglen);
 print STDERR "req: fraglen= $minfraglen .. $maxfraglen, ori=$frag_ori\n";
}
my $blibinfo=$Getopt::Std::opt_l;
my ($blibname, $bminlen, $bmaxlen, $b_ori);
if ($blibinfo) {
  ($blibname, $bminlen)=split(/:/,$blibinfo);
  die("Error getting library name from library info string ($blibinfo)\n") unless $blibname;
  if ($bminlen && $bminlen=~s/(\d+)([FfRr][RrFf])$/$1/) {
    $b_ori=lc($2);
    die("Error at -l option value, only RF or FR suffixes are recognized!\n")
       if ($b_ori ne 'fr' && $b_ori ne 'rf');
  }
  else { $b_ori='fr'; }
  if ($bminlen=~m/^(\d+)[\,\;\~\-\.](\d+)$/) {
     ($bminlen, $bmaxlen)=($1,$2);
     die("Error at -l option (invalid range)!\n") if ($bmaxlen<$bminlen);
  }
  #print STDERR "blibname=$blibname '$bminlen' - '$bmaxlen' ($b_ori)\n";
  #exit;
}

my $minovl=$Getopt::Std::opt_a; #could be percentage if it ends with '%'
my $movlperc=0; #percentage was used
$movlperc=1 if ($minovl=~s/\%$//);
my $maxclip=$Getopt::Std::opt_c;
my $maxovh=$Getopt::Std::opt_v;
my $minpid=$Getopt::Std::opt_p; #computed from edit distance
my $maxhits=$Getopt::Std::opt_t;
my $bambusOut=$Getopt::Std::opt_b;
my $max_edist=$Getopt::Std::opt_e;
my $noSameRefMates = $Getopt::Std::opt_S;
my $noSingleMates = $Getopt::Std::opt_B;
my $onlySameRefMates = $Getopt::Std::opt_M;
$noSingleMates = 1 if $onlySameRefMates;
my $ctgPEstats = $Getopt::Std::opt_w;
my $max_edist_ori;
if ($max_edist && $max_edist=~s/(\d+)([FfRr][RrFf])$/$1/) {
  $max_edist_ori=lc($2);
  die("Error at -e option value, only RF or FR suffixes are recognized!\n")
       if ($max_edist_ori ne 'fr' && $max_edist_ori ne 'rf');
}
my $keepUnmapped=$Getopt::Std::opt_U;
my $keepUnmappedMate=$Getopt::Std::opt_P;

my $pairProcessing = $noSameRefMates || $noSingleMates || $bambusOut || $minfraglen || $maxfraglen || $keepUnmappedMate;
if ($bambusOut) {
 die("Error: -b option required -l !\n") unless $blibname;
 #open(BAMBUS, '>'.$bambusOut.'.tmp') || die("Error creating file $bambusOut.tmp !");
}

# --
my @bfrags; # list of linking fragments [tname, r_id, r_name, mate_id, mate_name]

my %refs; # refname => [ reflen, num_links, [read1_dta], [read2_dta], ] 
          #   where read_dta is [ r_id, strand, ref_l, ref_r ]
my @reflst; #simple enumeration of reference sequences as found in the SAM header
my ($tname, $tlen, $qname, $qfrag, $hitcount, $qSeq, $qQual); #data for current read

my @pdata; #pair data: ($tname, $tlen, $sflag, [@rdata1], [@rdata2])
    # where @rdata=($refname, $pos, $strand, $ovlen, $pid,  $rlen, $hitindex, $nhits, $alnscore, $samline)
    #                   0       1       2       3      4      5       6         7        8          9
my $checkPair;
while (<>) {
  if ($fixReadNames) {
    # fix some badly named 454 PE reads
    s/^(\S+?)\/[12]\t/$1\t/; 
    s/^(\S+?)\.(left|right)\t/$1\t/;
  }
  my $line=$_;
  #check for header:
  if (m/^@[A-Z][A-Z]\t/) {
    print $_;
    #keep refseq length
    if (m/^\@SQ\tSN:(\S+)\tLN:(\d+)/) {
     my ($refname, $reflen)=($1, $2);
     $refs{$refname}=[ $reflen, 0 ];
     push(@reflst, $refname);
     }
    next;
  }
 chomp;
 my ($rname, $flags, $refname, $pos, $mapq, $cigarstr, 
     $rnext, $pnext, $tlen, $seq, $quals, $tags)=split(/\t/, $_, 12);
 my @mytags; #list of novel SAM tags to be added to the SAM record
 next unless length($tags)>5;
 $flags=int($flags);
 my $sflag = $flags & 0xc0; #either mate in a pair
 my $isrev = ($flags & 0x10) != 0;
 my $materev = ($flags & 0x20) != 0;
 my $new_pair=0;
 if ($rname ne $qfrag) { #new pair!
  $new_pair=1;
  if ($checkPair) {
     processFragment();
     @pdata=();
  }
  $tname=$rname;
  $checkPair = $pairProcessing && $sflag;
  if ($checkPair) {
    #          0       1       2     3     4
    @pdata=($rname, $tlen,  $sflag,  [],   []);
  }
  $hitcount=0;
#  $prev_alnscore=0;
  $qfrag=$rname;
  ($qSeq, $qQual) = $isrev ? (reverseComplement($seq), reverse($quals)) : ($seq, $quals);
 }
 
# my $frag_idx=0; # unpaired
 my $dix=3; #index of aln data in @pdata
 if ($sflag == 0x40) {
    $rname.='/1';
 #   $frag_idx=1;
    $dix=3;
 } elsif ($sflag == 0x80) {
    $rname.='/2';
 #   $frag_idx=2;
    $dix=4;
 }

 if ($new_pair==0 && $rname ne $qname) {
   #next read in the pair
   ($qSeq, $qQual) = $isrev ? (reverseComplement($seq), reverse($quals)) : ($seq, $quals);
 }
 $qname=$rname;
 if ($checkPair) {
     #initialize as if it's unmapped
     #this data will be updated if the read eventually passes the filters
     $pdata[$dix]->[1]=0;
     my $uflags=$flags;
     $uflags &= ~0x04; # remove the mapped flag
     $uflags &= ~0x10; # remove the reverse mapped flag
     $pdata[$dix]->[9]=join("\t", $tname, $uflags, '*', 0, 0, '*', 
                               $rnext, $pnext, 0, $qSeq, $qQual); #any tags?
 }

 if ($maxhits && $hitcount>=$maxhits) {
   next;
 }

 my $unmapped = (($flags & 4)!=0);
 #my $unmapped_mate = (($flags & 8)!=0);
 #my $mapped_both=$frag_idx && !$unmapped && !$unmapped_mate;
 #$rnext=$rname if $rnext eq '=';
 #my $mapped_pair = ($mapped_both && ($rnext eq $rname));

 #next if $mapped_pair && ( ($minfraglen && abs($tlen)<$minfraglen) 
 #    || ($maxfraglen && abs($tlen)>$maxfraglen) ) ;

 chomp($line);

 if ($unmapped) {
   next unless $keepUnmapped;  #|| $keepUnmappedMate);
   #keep unmapped:
   print "$line\n" unless ($checkPair);
   next;
 }
 
 #if ($frag_ori && $mapped_pair) {
 #  next unless ($isrev xor $materev); #must be mapped on opposite strands
 #  if ($frag_ori eq 'fr') { #FR
 #    next if ($isrev ? ($pnext>$pos) : ($pnext<$pos) );
 #  }
 #  else { # RF
 #    next if ($isrev ? ($pnext<$pos) : ($pnext>$pos) ); 
 #  }
 #}
 my ($alnscore) = ( $tags=~m/\bAS:i:(\-?\d+)/ );
 my ($hitindex) = ( $tags=~m/\bHI:i:(\-?\d+)/ );
 my ($nhits)    = ( $tags=~m/\bNH:i:(\-?\d+)/ );
 my ($edist) = ( $tags=~m/\bNM:i:(\d+)/ );
 if ($edist==0 && $tags=~m/\bnM:i:(\d+)/) {
   $edist=$1;
   $edist+=$_ foreach ($cigarstr=~m/(\d+)[ID]/g);
 }
 
 my $rlen=length($qSeq); #too short?
 if ($rlen<$minlen) {
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

 my $reflen=$refs{$refname}->[0]
   || die("Error: couldn't get reference length for $refname\n$line\n");

 $refovl+=$_ foreach (grep(/\d+[NDXM=]/, @cigar) ); #aln length on ref
 my $ovlen = $rlen-$clipL-$clipR; #aln length on read
 if ($ovlen==0) {
  print STDERR "Error ovlen==0 at:\n$line\n";
  print STDERR join("\n","rlen=$rlen", "clipL,R=$clipL,$clipR")."\n";
  die(1);
 }
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
 push(@mytags, 'YV:Z:'.$sovperc.'%') unless $line=~m/\bYV:Z:/;
 #my ($pseq, $pquals)= ($clipL>$clipR) ? (substr($qSeq, 0, $clipL),substr($qQual, 0, $clipL)) : 
 #                                   ( $clipR ? (substr($qSeq, -$clipR), substr($qQual, -$clipR)) :
 #                                     ('', '') ) ;

 my $pid=(($ovlen-$edist)*100.00)/$ovlen;
 
 if ($max_edist) {
   if ($max_edist_ori) {
     my $ori_dist;
     if ($max_edist_ori eq 'fr') { #FR
        #check only distance from the end of the read to the ref end on that side
        $ori_dist = $isrev ? $refL : $reflen-$refR;
      }
     else { # RF
        #check only distance from the start of the read to the ref end on that side
        $ori_dist = $isrev ? $reflen-$refR : $refL;
     }
     next if $ori_dist>$max_edist;
   }
   else { #no orientation requirement provided
     next if ($refL>$max_edist && $reflen-$refR>$max_edist);
   }
 }

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
 my $spid=sprintf("%.2f",$pid);
 $spid=~s/\.?0+$//;
 push(@mytags, 'YI:Z:'.$spid) unless $line=~m/\bYI:Z:/;
 if ($maxovh) {
   #ovhL/R are really non-matching regions at the end of the read
   #exclude parts that fall off the edge of the reference
   if ($ovhL>$refL) {
     # the read hangs over the left end of reference
     push(@mytags, 'YH:Z:XL'.($ovhL-$refL));
     $ovhL=$refL;
   }
   if ($ovhR>$reflen-$refR) {
     # the read hangs over the right end of ref
     push(@mytags, 'YH:Z:XR'.($ovhR-($reflen-$refR)));
     $ovhR = $reflen-$refR;
   }
   
   next if ($ovhR>$maxovh || $ovhL>$maxovh);
 }
 my $samline=$line;
 $samline.="\t".join("\t",@mytags) if @mytags>0;
 if ($checkPair) {
   next if (defined($pdata[$dix]->[8]) && $pdata[$dix]->[8]>$alnscore); #already got best hit
   my $newline=$line;
   $newline=~s/\tZT:Z:[CDPU][CDPU]\b//;
   $line.="\t".join("\t",@mytags) if @mytags>0;
                #    0       1      2        3      4      5        6        7       8          9
   $pdata[$dix]=[$refname, $pos, $strand, $ovlen, $pid, $rlen, $hitindex, $nhits, $alnscore, $newline];
 }
 else {
   print $line;
   print "\t".join("\t",@mytags) if @mytags>0;
   print "\n";
 }
 #$hitcount++ if ($prev_alnscore != $alnscore);
 
 $hitcount++;
 
 #$prev_alnscore=$alnscore;
 #$qTrash='V' if ($ovperc>=80) || (length($pseq)<$minlen);
 #next;
} #for each SAM line

processFragment() if $checkPair; #for the last fragment/read


if ($bambusOut) {
  my $date=US_Date();
  open(BAMBUS, '>'.$bambusOut) || die("Error creating file $bambusOut!\n");
  my $evid=$blibname.'_1';
  print BAMBUS qq{<EVIDENCE ID="$evid" DATE="$date" PROJECT="sam_$blibname"
   PARAMETERS="$cmdargs" >}."\n";
  print BAMBUS qq{<LIBRARY ID="lib_1" NAME="$blibname" MIN="$bminlen" MAX="$bmaxlen" >}."\n";
  my $i=0;
  foreach my $ins (@bfrags) {
    my ($tname, $r_id, $r_name, $m_id, $m_name)=@$ins;
    print BAMBUS qq{<INSERT ID="ins_$i" NAME="$tname">}."\n";
    $i++;
    print BAMBUS qq{<SEQUENCE ID="$r_id" NAME="$r_name" />}."\n";
    print BAMBUS qq{<SEQUENCE ID="$m_id" NAME="$m_name" />}."\n";
    print BAMBUS "</INSERT>\n";
  }
  print BAMBUS "</LIBRARY>\n";
  
  foreach my $ctg (@reflst) {
    my $refdata=$refs{$ctg} || die("Error getting ref data for contig $ctg\n");
    next if $$refdata[1]==0;
    my @refd=@$refdata;
    my $ctglen=shift(@refd);
    my $numlinks=shift(@refd);
    print BAMBUS qq{<CONTIG ID="$ctg" NAME="$ctg" LEN="$ctglen">}."\n";
    foreach my $rd (@refd) {
      my ($rid, $strand, $al, $ar)=@$rd;
      my $ori;
      if ($b_ori eq 'rf') {
        #swap strand
        $ori = ($strand eq '-') ? 'BE' : 'EB';
      }
      else {
        $ori = ($strand eq '+') ? 'BE' : 'EB';
      }
      print BAMBUS qq{<SEQUENCE ID="$rid" ORI="$ori" ASM_LEND="$al" ASM_REND="$ar" />}."\n";
    }
    print BAMBUS "</CONTIG>\n";
  }
  print BAMBUS "</EVIDENCE>\n";

  close(BAMBUS);
}
# --

#************ Subroutines **************
sub US_Date {
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
  return sprintf('%02d/%02d/%02d', $mon+1, $mday, $year-100);
}

sub timeStamp{
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    my $nice_timestamp = sprintf ( "[%04d.%02d.%02d %02d:%02d:%02d]\t",
                                   $year+1900,$mon+1,$mday,$hour,$min,$sec);
    return $nice_timestamp;
}

sub reverseComplement {
  my $s=reverse($_[0]);
  $s =~ tr/AaCcTtGgUuMmRrWwSsYyKkVvHhDdBb/TtGgAaCcAaKkYyWwSsRrMmBbDdHhVv/;
  return $s;
 }

sub processFragment {
 my $bi=3; #base index
 return unless defined($pdata[$bi]) && defined($pdata[$bi+1]);
 my ($tname, $tlen, $sfrag)=@pdata[0..2];
 my @rd1=@{$pdata[$bi]};
 my @rd2=@{$pdata[$bi+1]};
 if ($rd1[1]>0 && $rd2[1]>0) { #both mates are mapped
   if ($rd1[0] eq $rd2[0]) {
     #mapped on same contig
     return if $noSameRefMates;
     if ($frag_ori) {
        my $isrev=$rd1[2] eq '-';
        my $materev=$rd2[2] eq '-';
        return unless ($isrev xor $materev); #must be mapped on opposite strands
        if ($frag_ori eq 'fr') { #FR
           next if ($isrev ? ($rd2[1]>$rd1[1]) : ($rd2[1]<$rd1[1]) );
        }
        else { # RF
           return if ($isrev ? ($rd2[1]<$rd1[1]) : ($rd2[1]>$rd1[1]) ); 
        }
     }
     return if ( ($minfraglen && abs($tlen)<$minfraglen) 
         || ($maxfraglen && abs($tlen)>$maxfraglen) );
     
     print $rd1[9]."\tYT:Z:CP\n";
     print $rd2[9]."\tYT:Z:CP\n";
   }
   else { 
     return if $onlySameRefMates;
     #mapped on different contigs
     my $lnkid=scalar(@bfrags);
     if ($bambusOut) {
       #write this linkage info for Bambus
       push(@bfrags, [$tname, "r$lnkid".'_1', $tname.'_1', "r$lnkid".'_2', $tname.'_2']);
       my $refdata1=$refs{$rd1[0]} || die("Error getting ref data for $tname/1\n");
       my $refdata2=$refs{$rd2[0]} || die("Error getting ref data for $tname/2\n");
       $refdata1->[1]++; #inc number of links for these contigs
       $refdata2->[1]++;
       push( @{$refdata1}, [ "r$lnkid".'_1', $rd1[2], $rd1[1], $rd1[1]+$rd1[3] ] );
       push( @{$refdata2}, [ "r$lnkid".'_2', $rd2[2], $rd2[1], $rd2[1]+$rd2[3] ] );
     }
     print $rd1[9]."\tYT:Z:DP\n";
     print $rd2[9]."\tYT:Z:DP\n";
   }
 }
 else {
  #one or both are NOT mapped
  return if $noSingleMates;
  if ($rd1[1]>0) {
    print $rd1[9]."\tYT:Z:UP\n";
    print $rd2[9]."\tYT:Z:UP\n" if $keepUnmappedMate;
  }
  else {
    print $rd1[9]."\tYT:Z:UP\n" if $keepUnmappedMate;
    print $rd2[9]."\tYT:Z:UP\n";
  }
 }
}
