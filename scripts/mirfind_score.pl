#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;
use IO::Handle qw( );
use RNA;

my $usage = q/Usage:
  mirfind_score.pl [-c <cutoff>] [-D] <regions.fa> <regions.xmap>
 Options
  -c set the score cut-off (default 1)
  -A print all entries that qualify structurally, even with bad score
  -R use Randfold
  -D consider Drosha processing (number of base pairings in the lower stems)
/;
umask 0002;
getopts('VRDAc:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $showall=$Getopt::Std::opt_A;
my $verbose=$Getopt::Std::opt_V;
my ($use_randfold, $use_drosha)=($Getopt::Std::opt_R, $Getopt::Std::opt_D);
my $score_min=$Getopt::Std::opt_c || 1;
die("$usage\n") unless @ARGV==2;
my ($frna,$fxmaps)=@ARGV;
die("$usage\n") unless -f $frna && -f $fxmaps;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# -- setup RNA folding parameters:
$RNA::dangles=0;
#parameters
my $nucleus_lng=7;
my $score_star=3.9;
my $score_star_not=-1.3;
my $score_nucleus=3;
my $score_nucleus_not=-0.6;
my $score_randfold=1.6;
my $score_randfold_not=-2.2;
my $score_intercept=0.3;
my @scores_stem=(-3.1,-2.3,-2.2,-1.6,-1.5,0.1,0.6,0.8,0.9,0.9,0);
my $e=2.718281828;


 
# --                         0            1          2             3            4          5
my @xmaps; #   list of   [readID_xN, read_freq, read_alnstart, read_alnend, mismatches, intronInfo]
my %qmaps; # $read_name=>[readID_xN, read_freq, read_alnstart, read_alnend, mismatches, intronInfo]
#                            0            1          2             3           4            5
#TODO: while parsing this read data, collapse reads that are included in other reads (summing up the frequencies)
my $reg_id; #current region id
#                    0             1            2              3
my @rnafold; # (region_descr, region_seq, region_struct, region_mfe)
my %hash_bp; # base pairing data for current region
my $qmaxfreq; # read with maximum coverage in current region
my @qmaxfreqdata; # list of (readID_xN, read_freq, read_alnstart, read_alnend, mismatches, intronInfo) for the read with maximum coverage in current region
my @prevmir; # -- to avoid reporting duplicates
            # [ 0=mature_loc, 1=hairpin_loc, 2=final_score, 3=mfe_score, 4= report line ]
#------
open(INRNA, $frna) || die ("Error opening file $frna\n");
open(INMAPS, $fxmaps) || die ("Error opening file $fxmaps\n");
my $prevchr;
while (($reg_id=RNAfold(*INRNA, \@rnafold))) {
 #print STDERR "DBG: processing region: $reg_id\n" if $verbose;
 my $x_id=readXmaps(*INMAPS, \@xmaps, \%qmaps);
 die("Error: couldn't match RNAfold id ($reg_id) with mappings id ($x_id)\n") unless $x_id eq $reg_id;
 my $chr=$x_id;
 $chr=~s/_\d+$//;
 if ($chr ne $prevchr) {
   print STDERR "processing regions on $chr..\n";
   $prevchr=$chr;
   }
 processRegion(\@rnafold, \@xmaps);
 }
print $prevmir[4] if $prevmir[4];
close(INMAPS);
close(INRNA);
print STDERR "Done.\n";
# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }
#************ Subroutines **************

sub processRegion {
 my ($rfold, $rmaps)=@_;
 my ($g_id, $g_descr, $g_seq, $g_struct, $g_mfe)=@$rfold;
 my $g_seqlen=length($g_seq);
 #print STDERR ">greg: $g_id mfe=$g_mfe\nseq: $g_seq\nstr: $g_struct\n";
 # --> fill_structure()
 #reads the dot bracket structure into the 'bp' hash where each key and value are basepaired
 %hash_bp=();
 my $lng=length($g_struct); 
 #local stack for keeping track of basepairings
 my @bps;
 for(my $i=1;$i<=$lng;$i++){
      #my $pstruct=excise_str($g_struct,$i,$i);
      my $ps=substr($g_struct,$i-1,1);
      if ($ps eq "(") {
          push(@bps,$i);
          next;
          }
      if ($ps eq ")") {
          my $i_prev=pop(@bps);
          $hash_bp{$i_prev}=$i;
          $hash_bp{$i}=$i_prev;
          }
      }
 # --> fill_pri() 
 #fills basic specifics on the precursor into the 'comp' hash
 # we already have the values in $g_id, $g_descr, $g_seq, $g_struct, $g_mfe and $g_seqlen
 
 # --> fill_mature()
    my $mature_query=$qmaxfreq;
    #my($mature_beg,$mature_end)=get_qloc($mature_query);
    my ($mature_beg, $mature_end)=($qmaxfreqdata[2], $qmaxfreqdata[3]);
    my $mature_seq=excise_str($g_seq, $mature_beg, $mature_end);
    my $mature_struct=excise_str($g_struct,$mature_beg,$mature_end);
    my $mature_arm=arm_mature($mature_struct); # 1, 2 or 0 
#     $hash_comp{"mature_query"}=$mature_query;
#     $hash_comp{"mature_beg"}=$mature_beg;
#     $hash_comp{"mature_end"}=$mature_end;
#     $hash_comp{"mature_strand"}=$mature_strand;
#     $hash_comp{"mature_struct"}=$mature_struct;
#     $hash_comp{"mature_seq"}=$mature_seq;
#     $hash_comp{"mature_arm"}=$mature_arm;

  #  --> fill_star{
  #fills specifics on the expected star strand into 'comp' hash ('component' hash)
  #if the mature sequence is not plausible, don't look for the star arm
    my ($star_arm, $star_beg, $star_end, $star_seq, $star_struct);

    return unless ($mature_arm>0); # EXIT point if no mature arm determined
    ($star_beg,$star_end)=find_star($mature_beg, $mature_end);
    ($star_arm, $star_struct)=arm_star($star_beg, $star_end, $g_seqlen, $mature_arm, $mature_beg, $mature_end, $g_struct);
    return unless $star_struct; #EXIT point if no star sequence

    #excise expected star sequence and structure
    $star_seq=excise_str($g_seq,$star_beg,$star_end);
  #     $hash_comp{"star_beg"}=$star_beg;
  #     $hash_comp{"star_end"}=$star_end;
  #     $hash_comp{"star_seq"}=$star_seq;
  #     $hash_comp{"star_struct"}=$star_struct;
  #     $hash_comp{"star_arm"}=$star_arm;


    return if ($mature_arm+$star_arm!=3); #bail out if no valid mature and star combo

 # --> fill_loop{
    #fills specifics on the loop sequence into the 'comp' hash
    #unless both mature and star sequences are plausible, do not look for the loop
    my ($loop_beg, $loop_end, $loop_seq, $loop_struct);
    #defining the begin and end positions of the loop from the mature and star positions
    #excision depends on whether the mature or star sequence is 5' of the loop ('first')
    ($loop_beg, $loop_end)=($mature_arm==1) ? ($mature_end+1, $star_beg-1) : ($star_end+1, $mature_beg-1);
    #unless the positions are plausible, do not fill hash
    return unless ($loop_beg<$loop_end && $loop_end<=$g_seqlen && $loop_beg>0);

    $loop_seq=excise_str($g_seq,$loop_beg,$loop_end);
    $loop_struct=excise_str($g_struct,$loop_beg,$loop_end);
 # --> fill_lower_flanks()
  #fills specifics on the lower flanks and unpaired strands into the 'comp' hash
    my ($flank_first_end,$flank_second_beg)=($mature_arm==1) ? ($mature_beg-1, $star_end+1) : 
                      ($star_beg-1, $mature_end+1);
   return unless ($flank_first_end<=$flank_second_beg && $flank_first_end>0 && $flank_second_beg<=$g_seqlen);
   my $flank_first_seq=excise_str($g_seq, 1, $flank_first_end);
   my $flank_second_seq=excise_str($g_seq, $flank_second_beg, $g_seqlen);
   my $flank_first_struct=excise_str($g_struct, 1, $flank_first_end);
   my $flank_second_struct=excise_str($g_struct, $flank_second_beg, $g_seqlen);
# -- fill stems Drosha   
   my ($stem_first, $stem_second, $stem_bp_first, $stem_bp_second, $stem_bp);
   if ($use_drosha) {
     $stem_first=substr($flank_first_struct, -10);
     $stem_second=substr($flank_second_struct,0,10);
     $stem_bp_first=($stem_first=~tr/(/(/);
     $stem_bp_second=($stem_second=~tr/)/)/);
     $stem_bp=($stem_bp_first<$stem_bp_second)?$stem_bp_first : $stem_bp_second;
     }
# -- final checks:
   my ($pre_struct, $pre_seq, $mstar_struct) = ($mature_arm==1) ? 
                   ($mature_struct.$loop_struct.$star_struct, $mature_seq.$loop_seq.$star_seq, $mature_struct.$star_struct) :
                   ($star_struct.$loop_struct.$mature_struct, $star_seq.$loop_seq.$mature_seq, $star_struct.$mature_struct);
  #simple pattern matching checks for bifurcations
  my $unloop=($loop_struct=~tr/././); #how many unpaired bases in loop
  #print STDERR ">>>>> $g_id : loop_struct_len=".length($loop_struct)." unpaired bases: $unloop | $loop_struct\n";
  unless (($unloop<<1)>length($loop_struct)) { #too many 
    print STDERR "[rejected]\t$g_id\t$g_descr($g_seqlen)\ttoo many pairings found in loop structure\n" if $verbose;
    }
  unless ($mstar_struct=~m/^[\.\(]+[\.\)]+$/) {
    print STDERR "[rejected]\t$g_id\t$g_descr($g_seqlen)\tbifurcation detected in precursor\n" if $verbose;
    return;
    }
  #minimum 14 base pairings required in duplex
  my $duplex_bps=($mature_struct=~tr/()/()/);
  unless ($duplex_bps>=14) {
    print STDERR "[rejected]\t$g_id\t$g_descr($g_seqlen)\ttoo few base pairings in duplex\n" if $verbose;
    return;
    }
  #no more than 6 nt difference between mature and star length
  if (abs(length($mature_seq)-length($star_seq))>6) {
    print STDERR "[rejected]\t$g_id\t$g_descr($g_seqlen)\ttoo big difference between mature and star lengths\n" if $verbose;
    return;
    }
 # ---> pass_filtering_signature() 
  #number of reads that map in consistence with Dicer processing
  my $consistent=0;
  #number of reads that map inconsistent with Dicer processing
  my $inconsistent=0;
  # number of potential star reads map in good consistence with Drosha/Dicer processing
  #  (3' overhangs relative to mature product)
  my $star_perfect=0;
  # number of potential star reads that do not map in good consistence with 3' overhang
  my $star_fuzzy=0;
  #sort queries (deep sequences) by their position on the hairpin
  my @xmapos = sort { $a->[2]<=>$b->[2] } @xmaps;
  foreach my $qd (@xmapos) {
    my ($qid, $qfreq, $qbeg, $qend,$qmism, $qinfo)=@$qd;
    #test which Dicer product (if any) the read corresponds to
    #the deep sequence is allowed to stretch 2 nt beyond the expected 5' end
    my $fuzz_beg=2;
    #the deep sequence is allowed to stretch 5 nt beyond the expected 3' end
    my $fuzz_end=5;
    my $product=0; # 1= mature, 2=star, 3=loop
    if (contained($qbeg, $qend, $mature_beg-$fuzz_beg, $mature_end+$fuzz_end)) {
        $product=1;
        }
      elsif (contained($qbeg, $qend, $star_beg-$fuzz_beg, $star_end+$fuzz_end)) {
        $product=2;
        if (abs($qbeg-$star_beg)<2) { $star_perfect+=$qfreq; }
                               else { $star_fuzzy+=$qfreq; }
        }
      elsif (contained($qbeg, $qend, $loop_beg-$fuzz_beg, $loop_end+$fuzz_end)) {
        $product=3;
      }
    #if the deep sequence corresponds to a Dicer product, add to the 'consistent' variable  
    if ($product) {  $consistent+=$qfreq;  }
    #if the deep sequence do not correspond to a Dicer product, add to the 'inconsistent' variable
             else {  $inconsistent+=$qfreq; }
    }
# if the majority of potential star sequences map in good accordance with 3' overhang
#  score for the presence of star evidence
  my $star_read=1 if ($star_perfect>$star_fuzzy);
   #total number of reads mapping to the hairpin
  my $freq=$consistent+$inconsistent;
   #$hash_comp{"freq"}=$freq;
  unless($freq>1) {
       print STDERR "[rejected]\t$g_id\t$g_descr($g_seqlen)\tread frequency too low\n" if $verbose;
       return;
       }
  #unless >90% of the reads map in consistence with Dicer processing, the hairpin is discarded
  my $cfraction=$inconsistent/($inconsistent+$consistent);
  unless($cfraction<=0.1) { 
     print STDERR "[rejected]\t$g_id\t$g_descr($g_seqlen)\tpoor mapping frequency: inconsistent:$inconsistent, consistent:$consistent\n"
          if $verbose;
     return;
     }
# --> pass_threshold_score() 
  #minimum free energy of the potential precursor
  my $score_mfe=score_mfe($g_mfe);
  #count of reads that map in accordance with Dicer processing
  my $score_freq=score_freq($freq);
  #basic score
  my $score=$score_mfe+$score_freq;
  $score+= $star_read ? $score_star : $score_star_not; #penalty if there are no "perfect" star reads
  #score lower stems for potential for Drosha recognition
  if($use_drosha){
    $score+=$scores_stem[$stem_bp];
    }
  $score+=$score_intercept; # ?? why
  if ($use_randfold) {
   # only calculate randfold value if it can make the difference between the potential precursor
   # being accepted or discarded
     if($score+$score_randfold>=$score_min and $score+$score_randfold_not<=$score_min) {
        #randfold value<0.05
       if(test_randfold($g_seq)) { $score+=$score_randfold; }
         #randfold value>0.05
                          else{ $score+=$score_randfold_not; }
       }
     } #use randfold test
  #round off values to one decimal
  my $round_mfe=sprintf('%.1f', $score_mfe);
  my $round_freq=sprintf('%.1f',$score_freq);
  my $round=sprintf('%.1f', $score);
  #print scores
  my $passed=($score>=$score_min) ? '[passed]':'';
  # print STDERR "score_mfe\t$round_mfe\nscore_freq\t$round_freq\nscore\t$round\t$passed\n";
  #return 1 if the potential precursor is accepted, return 0 if discarded
  #--- print all that made it up to here
  my $g_chr=$g_id;
  $g_chr=~s/_\d+$//;
  my ($g_strand, $g_range)=split(/\|/,$g_descr);
  my ($g_beg, $g_end)=split(/[\-\._]+/,$g_range);
  my ($mature_gbeg, $mature_gend)=coords2genome($g_strand, $g_beg, $g_end, $mature_beg, $mature_end);
  my ($star_gbeg, $star_gend)=coords2genome($g_strand, $g_beg, $g_end, $star_beg, $star_end);
  #TODO: perhaps we can extend $flank_first_end and $flank_second_beg based on actual read coverage and RNA structure?
  my ($hairpin_gbeg, $hairpin_gend)=coords2genome($g_strand, $g_beg, $g_end, $flank_first_end, $flank_second_beg);
  my $m_loc_id=$g_chr.$g_strand.$mature_gbeg.'-'.$mature_gend;
  my $h_loc_id=$g_chr.$g_strand.$hairpin_gbeg.'-'.$hairpin_gend;
  
  my $pline=join("\t", $g_id.'|'.$g_descr, $g_chr, $g_strand, $hairpin_gbeg, $hairpin_gend, $mature_query, 'x'.$qmaxfreqdata[1],
                $mature_seq, 'm|'.$mature_gbeg.'-'.$mature_gend, '*|'.$star_gbeg.'-'.$star_gend, $round_mfe, 
                $round_freq, $round, $passed, $qmaxfreqdata[5])."\n";
  #print $pline;
  if ($passed) {
     if ($prevmir[0] eq $m_loc_id) { #  && $prevmir[1] eq $h_loc_id) {
      # same miRNA as before
      if ($prevmir[2]>$round || ($prevmir[2]==$round && $prevmir[3]>$round_mfe)) {
         #prev miRNA was the same but had a better score, print that
         print $prevmir[4];
         }
        else { #prev didn't have a better score, print this
         print $pline;
         }
       @prevmir=();
      }
     else { #different miRNA than before
      print $prevmir[4] if $prevmir[4];
      @prevmir=($m_loc_id, $h_loc_id, $round, $round_mfe, $pline);
      #print $pline;
      }
     } # good prediction
    else { # optionally, report also the potentially valid structures that failed the scoring test
     print $prevmir[4] if $prevmir[4];
     print $pline if $showall;
     @prevmir=();
     }
 STDOUT->flush();
}

sub coords2genome { #convert local coordinates to genome coordinates
 my ($strand, $g1, $g2, $l1, $l2)=@_;
 return ($strand eq '-') ? ($g2-$l2+1, $g2-$l1+1) : ($g1+$l1-1, $g1+$l2-1);
}

sub contained{
    #is c1-c2 interval contained within a1-a2 ?
    my($c1,$c2,$a1,$a2)=@_;
    #testbeginend($beg1,$end1,$beg2,$end2);
    return ($c1>=$a1 and $c2<=$a2);
}

sub arm_mature {
   #tests whether the mature sequence is in the 5' ('first') or 3' ('second') arm of the potential precursor
  my ($struct)=@_; 
  #-- there should be no bifurcations and minimum one base pairing at that respective end
  #my $struct=excise_seq($hash_comp{"pri_struct"},$beg,$end,$strand);
  my $r1=index($struct,'(');
  my $r2=index($struct,')');
  if ($r2<0 && $r1>=0) {   #struct=~m/^[\(\.]+$/ && $struct=~m/\(/) {
       return 1;
       }
  if ($r1<0 && $r2>=0) {   # ($struct=~m/^[\)\.]+$/ && $struct=~m/\)/){
       return 2;
       }
  return 0;
}

sub find_star {
    #uses the 'bp' hash to find the expected star begin and end positions from the mature positions
    my ($mature_beg, $mature_end)=@_;
    $mature_end-=2; #the -2 is for the overhang
    my $mature_lng=$mature_end-$mature_beg+1;
    #in some cases, the last nucleotide of the mature sequence does not form a base pair,
    #and therefore does not basepair with the first nucleotide of the star sequence.
    #In this case, the algorithm searches for the last nucleotide of the mature sequence
    #to form a base pair. The offset is the number of nucleotides searched through.
    my $offset_star_beg=0;
    my $offset_beg=0;
    #the offset should not be longer than the length of the mature sequence, then it
    #means that the mature sequence does not form any base pairs
    while(!$offset_star_beg and $offset_beg<$mature_lng){
       if($hash_bp{$mature_end-$offset_beg}){
           $offset_star_beg=$hash_bp{$mature_end-$offset_beg};
        }else{
           $offset_beg++;
          }
      }
    #when defining the beginning of the star sequence, compensate for the offset
    my $star_beg=$offset_star_beg-$offset_beg;
    #same as above
    my $offset_star_end=0;
    my $offset_end=0;
    while (!$offset_star_end and $offset_end<$mature_lng){
       if($hash_bp{$mature_beg+$offset_end}){
         $offset_star_end=$hash_bp{$mature_beg+$offset_end};
        }else{
          $offset_end++;
        }
      }
    #the +2 is for the overhang
    my $star_end=$offset_star_end+$offset_end+2;
    return ($star_beg, $star_end);
}

sub arm_star{
  #tests whether the star sequence is in the 5' ('first') or 3' ('second') arm of the potential precursor
  my ($beg, $end, $pri_end, $mature_arm, $mature_beg, $mature_end, $g_struct)=@_;
  #unless the begin and end positions are plausible, test negative
  return (0,'') unless ($beg>0 and $end>=$beg and $end<=$pri_end);
  #must be no overlap between the mature and the star sequence
  return (0,'') if ($mature_arm==1 && $mature_end>=$beg); #check ?
  return (0,'') if ($mature_arm==2 && $mature_beg<=$end);
  #there should be no bifurcations and minimum one base pairing

  my $struct=excise_str($g_struct,$beg,$end);
  my $r1=index($struct,'(');
  my $r2=index($struct,')');
  if ($r2<0 && $r1>=0) {   #($struct=~/^(\(|\.)+$/ and $struct=~/\(/)
       return (1,$struct);
       }
  if ($r1<0 && $r2>=0) {   # ($struct=~/^(\)|\.)+$/ and $struct=~/\)/)
       return (2, $struct);
       }
  return (0, $struct);
}



sub get_qloc {
 my $qd=$qmaps{$_[0]} || die ("Error: couldn't retrieve read $_[0] info (on $reg_id)\n");
 return ($$qd[2], $$qd[3]);
}


sub excise_str {
    #excise sub structure
    my($struct,$beg,$end)=@_;
    my $lng=length($struct);
    #begin can be equal to end if only one nucleotide is excised
    die("Error: invalid begin..end coordinates given to excise_str ($reg_id\:$beg..$end) lng=$lng\n(struct: $struct)\n") 
        if ($beg>$end || $beg>$lng || $end>$lng);
    #rarely, permuted combinations of signature and structure cause out of bound excision errors.
    #this happens once appr. every two thousand combinations
    #return 0 unless ($beg<=length($struct));
    
    #if excising relative to minus strand, positions are reversed
    #if($strand eq "-"){($beg,$end)=rev_pos($beg,$end,$lng);}
    return substr($struct,$beg-1,$end-$beg+1);
}


sub RNAfold {
 my ($fh, $rfd)=@_;
 #             0             1            2              3           4
 # @$rfd=(region_id,  region_descr, region_seq, region_struct, region_mfe)
 
 @$rfd=();
 my $rgid;
 $_=<$fh>; #header line
 return '' unless $_;
 my $rgdescr;
 ($rgid, $rgdescr)=(m/^>(\S+)\s+(\S+)/);
 push(@$rfd, $rgid, $rgdescr);
 my $seq=<$fh>; # sequence line
 return '' unless $seq; #terminated prematurely?
 $seq=~tr/\n\r \t//d;
 my ($structure, $mfe)=RNA::fold($seq);
 if (length($structure)!=length($seq)) {
   die("Error: RNA::fold returned a structure with different length!\n".
     "seq: $seq\nstr: $structure\n");
   }
 push(@$rfd, $seq, $structure, $mfe);
 return $rgid;
 }
 
 
sub readXmaps {
 my ($fh, $rxmaps, $hxmaps)=@_;
 @$rxmaps=(); #clear the array
 %$hxmaps=();
 $qmaxfreq='';
 @qmaxfreqdata=();
 #my $freqmax=0;
 my $rgid;
 my @rmaps;
 while (<$fh>) {
   if (m/^>(\S+)/) {
     if ($rgid) {
       seek($fh, -length($_), 1);
       last; # --> jump to processing
       }
     $rgid=$1;
     next;
     }
   chomp;
   my @t=split(/\t/);
   die("Error: invalid xmaps line format ($_)!\n") 
       unless ($t[1]>0 && $t[2]>$t[1]);
   my $qname=shift(@t);
   my ($qfreq)=($qname=~m/_x(\d+)$/);
   #insert read freq and store the whole thing
   my $rdata=[$qname, $qfreq, @t];
   #      0          1           2              3            4         5
   # [readID_xN, read_freq, read_alnstart, read_alnend, mismatches, intronInfo]
   #push(@$rxmaps, $rdata);
   push(@rmaps, $rdata);
   #$hxmaps->{$qname}=$rdata;
   }
 die("Error: no mappings found for region $rgid !\n") unless @rmaps>0;
 #collapse near-duplicate mappings (and sum up the frequencies)
 # -- sort by start position
 @rmaps = sort { ($a->[2]==$b->[2])? $a->[3]<=>$b->[3] : $a->[2] <=> $b->[2] } @rmaps;
 
my $i=0;
while ($i+1 < @rmaps) {
 my $pdel=0;
 for (my $j=$i+1;$j<@rmaps;$j++) {
      my $pd=$rmaps[$i];
      my $cd=$rmaps[$j];
      my ($fr, $idel, $ikeep) = ($$cd[1]>$$pd[1])? ($$cd[1]/$$pd[1], $i, $j) : ($$pd[1]/$$cd[1], $j, $i); #freq ratio
      my $plen=$$pd[3]-$$pd[2]+1;
      my $clen=$$cd[3]-$$cd[2]+1;
      my $pmm=$rmaps[$i]->[4];
      my $cmm=$rmaps[$j]->[4];
      my $dbg=0;
      #if ($rgid eq 'chrX_5' && !$cmm && !$pmm) {
      #   print STDERR "checking merge between $$pd[0] ($$pd[2]-$$pd[3]) and $$cd[0] ($$cd[2]-$$cd[3])\n";
      #   $dbg=1;
      #   }
      $pmm=~s/\d+\://g;
      $cmm=~s/\d+\://g;
      last if ($$cd[2]-$$pd[2]>16);
      $fr=200 if ($$cd[1]<20 && $$pd[1]<20); 
      if ( ($$cd[2]-$$pd[2]<3 && abs($$cd[3]-$$pd[3])<3) && abs($plen-$clen)<5 && $pmm eq $cmm && $fr<3.6) { #collapse very close, similar mappings
         #if ($dbg) {
         #  print STDERR "..merging $rmaps[$idel]->[0] ($rmaps[$idel]->[1]\:$rmaps[$idel]->[2]-$rmaps[$idel]->[3]) into ".
         #  "$rmaps[$ikeep]->[0] ($rmaps[$ikeep]->[1]\:$rmaps[$ikeep]->[2]-$rmaps[$ikeep]->[3])\n";
           #     if $rmaps[$ikeep]->[0] eq 'E2s5T00116060_x5438';
         #  }
         $rmaps[$ikeep]->[1] += $rmaps[$idel]->[1];
         $rmaps[$ikeep]->[2] = $pd->[2];
         $rmaps[$ikeep]->[3] = ($$cd[3]>$$pd[3])? $$cd[3] : $$pd[3];
         $pdel=1 if $idel==$i;
         splice(@rmaps, $idel, 1);
         last if $pdel;
         }
      } #for each $j (read mapping after $i)
    $i++ unless $pdel;
   } #for each read mapping $i
 #also populate the hash
 my $freqmax=0;
 foreach my $rd (@rmaps) {
   push (@$rxmaps, $rd);
   $hxmaps->{$$rd[0]}=$rd;
   if ($$rd[1]>$freqmax) {
      $qmaxfreq=$$rd[0];
      $freqmax=$$rd[1];
      }
   }
 #print STDERR " DBG: $rgid : qmaxfreq=$qmaxfreq\n" if $verbose;
 @qmaxfreqdata=@{$hxmaps->{$qmaxfreq}};
 return $rgid;
 }


sub test_randfold {
    my ($rseq)=@_;
    #print sequence to temporary file, test randfold value, return 1 or 0
    #print_file("pri_seq.fa",">pri_seq\n".$hash_comp{"pri_seq"});
    my $p_value=`echo ">\\n$rseq" | randfold -s - 200 | cut -f3`;
    chomp($p_value);
    if($p_value<=0.05){return 1;}
    return 0;
}

sub score_mfe{
#   scores the minimum free energy in kCal/mol of the potential precursor
#   Assumes Gumbel distribution as described in methods section of miRDeep manuscript
    my $mfe=shift;
    #numerical value, minimum 1
    my $mfe_adj = (-$mfe>1) ? -$mfe : 1 ;
    #parameters of known precursors and background hairpins, scale and location
    my $prob_test=prob_gumbel_discretized($mfe_adj,5.5,32);
    my $prob_background=prob_gumbel_discretized($mfe_adj,4.8,23);
    my $odds=$prob_test/$prob_background;
    my $log_odds=log($odds);
    return $log_odds;
}

sub prob_gumbel_discretized{
#   discretized Gumbel distribution, probabilities within windows of 1 kCal/mol
#   uses the subroutine that calculates the cdf to find the probabilities
    my ($var,$scale,$location)=@_;
    my $bound_lower=$var-0.5;
    my $bound_upper=$var+0.5;
    my $cdf_lower=cdf_gumbel($bound_lower,$scale,$location);
    my $cdf_upper=cdf_gumbel($bound_upper,$scale,$location);
    my $prob=$cdf_upper-$cdf_lower;
    return $prob;
}

sub cdf_gumbel{
#   calculates the cumulative distribution function of the Gumbel distribution
    my ($var,$scale,$location)=@_;
    my $cdf=$e**(-($e**(-($var-$location)/$scale)));
    return $cdf;
}

sub score_freq{
#   scores the count of reads that map to the potential precursor
#   Assumes geometric distribution as described in methods section of manuscript
    my $freq=shift;
    #parameters of known precursors and background hairpins
    my $parameter_test=0.999;
    my $parameter_control=0.6;
    #log_odds calculated directly to avoid underflow
    my $intercept=log((1-$parameter_test)/(1-$parameter_control));
    my $slope=log($parameter_test/$parameter_control);
    my $log_odds=$slope*$freq+$intercept;
    #if no strong evidence for 3' overhangs, limit the score contribution to 0
    #unless($options{x} or $hash_comp{"star_read"}){$log_odds=min2($log_odds,0);}
    return $log_odds;
}
