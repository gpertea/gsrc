#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;
my $wigbindir='/fs/szannotation/human_rnaseq/hg19_phastCons';
#my $fmid='phyloP46way'; #middle part of the file name

my $usage = qq/Usage:
 hg19_getcons.pl  {-r <chr>:<interval_list> | -g <gff>} \
     [-d <wigfixbin_dir>]
 Returns conservation average and sample standard deviation across 
 bases in given intervals (exons)
 
 Requires a directory with .wigfixbin files for each chromosome.
 (default: $wigbindir)
/;

umask 0002;
my %ch; # file access cache: chr->[filepath, start_cpos]
    # the value for start_cpos is always at offset 9 in the file!
    
getopts('g:r:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }

my $intv=$Getopt::Std::opt_r;
if ($intv) {
 my ($chr,$rlst)=split(/\:/,$intv);
 die("$usage Incorrect format for the interval list!\n") unless $chr && $rlst;
 my $strand=chop($chr);
 if ($strand ne '-' && $strand ne '+') {
   $chr.=$strand;
   $strand=undef;
   }
 my @rdata=map { [split(/[\-\.]+/)] } (split(/[\,\;\s]+/,$rlst));
 foreach my $d (@rdata) {
  ($$d[0], $$d[1])=($$d[1], $$d[0]) if $$d[0]>$$d[1];
  }
 my @ex = sort { $a->[0] <=> $b->[0] } @rdata;
 my ($avg, $std)=fetchCons($chr, \@ex);
 print STDOUT "$intv\t$avg\t$std\n";
 exit;
 }
#else { 
# GFF/GTF input
my $fgff=$Getopt::Std::opt_g;
my %gffrecs;
# ---------------------  0     1         2         3       4      5       6       7        8   
#           recID =>  [ chr, strand, feat_type,  gname, tdescr, fstart, fend, [@exons], [@cds] ]

# --
die("$usage Error: -r or -g required!\n") unless $fgff;
my $gffh;
if ($fgff eq '-') {
  open($gffh, "<&=STDIN") || die("Error: couldn't alias STDIN. $!\n");
  }
 else {
  open($gffh, $fgff) || die("Error: failed to open $fgff. $!\n");
  } 

loadGff($gffh, \%gffrecs);
# -- sort features by chromosome:
#print STDERR "GFF data loaded. Sorting by chromosome location..\n";
my @sorted_features=sort sortByLoc keys(%gffrecs);
# print STDERR "Writing BED file..\n";
processGffRecs(\%gffrecs, \@sorted_features);
# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************
sub loadGff {
 my ($fh, $recs)=@_; # $hr=\%gffrecs 
 # ---------------------  0     1         2         3       4      5       6       7        8   
 #           recID =>  [ chr, strand, feat_type,  gname, tdescr, fstart, fend, [@exons], [@cds] ]
 while (<$fh>) {
   next if m/^\s*#/;
   chomp;
   my ($chr, $track, $f, $fstart, $fend, $fscore, $strand, $frame, $lnum)=split(/\t/);
   next unless $fstart>1 && $lnum;
   next if $f eq 'gene' || $f eq 'locus'; # Warning: skipping any 'gene' or 'locus' features, unconditionally
   my $gff3_ID;
   my $gff3_Parent;
   my ($gname, $tdescr);
   ($fstart, $fend)=($fend, $fstart) if $fend<$fstart;
   ($gff3_ID)=($lnum=~m/\bID=([^;]+)/);
   ($gff3_Parent)=($lnum=~m/\bParent=([^;]+)/);
   if ($gff3_ID || $gff3_Parent) { # GFF format
      $gff3_ID=~tr/"//d; #"
      $gff3_Parent=~tr/"//d; #"
      $gff3_Parent='' if ($f eq 'mRNA');
      if ($gff3_ID && !$gff3_Parent) { #top level feature
         if ($f=~m/RNA/i || $f=~/gene/) {
           # try to parse the description, if any
           $tdescr='';
           $gname='';
           if ($lnum=~m/\b(?:descr|tophit|info|product)\s*=\s*"?([^;"]+)/i) {
             $tdescr=$1;
             }
            elsif ($lnum=~m/Name\s*=\s*"?([^;"]+)/) {
             $tdescr=$1;
             }
           if ($lnum=~m/\bgene_name[\s=]+"?([^;"]+)/i) {
             $gname=$1;
             }
            elsif ($lnum=~m/Name\s*=\s*"?([^;"]+)/) {
             $gname=$1;
             }
           $tdescr='' if ($tdescr eq $gname);
           $gname='' if $gname eq $gff3_ID;
           }
         die("Error: duplicate feature $gff3_ID on $chr\n") if (exists($$recs{"$chr|$gff3_ID"}));
         my $recID="$chr|$gff3_ID";
         $$recs{$recID} = [$chr, $strand, $f, $gname, $tdescr, $fstart, $fend, [], [] ];
         next;
         } # parent/top-level feature
      } #GFF
     else { #GTF format
      next if ($f eq 'transcript'); #exception: GTF with parent 'transcript' feature
      }
   # -------------- exon/CDS line here:
   my $recID;
   if ($track=~/^jigsaw/ && $lnum=~m/^\d+$/) {
     $recID=$chr.'.jsm.'.$lnum;
     }
    elsif ($lnum=~m/Parent=(['"\:\w\|\-\.]+)/) {
     $recID=$1;
     $recID=~tr/"//d; #"
     }
    elsif ($lnum=~m/transcript_id[= ]+(['"\:\w\.\|\-]+)/) {
     $recID=$1;
     $recID=~tr/"//d; #"
     }
    else {
     die("Error: cannot parse locus/transcript name from input line:\n$_\n");
     }
   if (!$gname && $lnum=~m/gene_id[= ]+(['"\:\w\.\|\-]+)/) {
      $gname=$1;
      $gname=~tr/"//d; #"
      }
   $tdescr='' if index($recID, $tdescr)>=0;
   $gname='' if index($recID, $gname)>=0;
   $recID=$chr.'|'.$recID;
   my $ld = $$recs{$recID};
   if ($ld) { #existing entry
     my $i=($f eq 'CDS') ? 8 : 7;
     my ($lstart, $lend)=($$ld[5], $$ld[6]);
     $$ld[5]=$fstart if $fstart<$lstart;
     $$ld[6]=$fend if $fend>$lend;
     push(@{$$ld[$i]}, [$fstart, $fend, $fscore]);
     }
    else { # first time seeing this locus/gene
     $$recs{$recID} = ($f eq 'CDS') ? 
           [$chr, $strand, $f, $gname, $tdescr, $fstart, $fend, [], [[$fstart, $fend, $fscore]] ] :
           [$chr, $strand, $f, $gname, $tdescr, $fstart, $fend, [[$fstart, $fend, $fscore]], [] ] ;
     }     # 0       1      2     3       4        5       6     7 (exons)                   8 (CDS)
 } #while <$fh>
} # loadGff

sub sortByLoc {
 my $da=$gffrecs{$a};
 my $db=$gffrecs{$b};
 if ($$da[0] eq $$db[0]) {
    return ($$da[5]==$$db[5]) ? $$da[6] <=> $$db[6] : $$da[5] <=> $$db[5] ;
    }
  else { return $$da[0] cmp $$db[1] ; }
}

sub processGffRecs {
 #return if keys(%recs)==0;
 my ($recs, $rlist)=@_;
 my @recs_keys;
 unless ($rlist) {
    @recs_keys=keys(%$recs);
    $rlist=\@recs_keys;
    }
 foreach my $gffid (@$rlist) {
    my $ld=$$recs{$gffid};
   # my $ld=$recs{$l} || die ("Error: locus $l found in list but not in hash!\n");
   #     0       1       2       3       4      5         6    7    8
   my ($chr, $strand, $ftype, $gname, $descr, $lstart, $lend, $er, $cr) = @$ld;
   # my ($mstart,$mend)=($lstart, $lend);
   my $CDexons=0;
   my @ex;
   my @cds;
   #some records might lack exons, but have only CDS segments (e.g. mitochondrial genes)
   if (@$er<1 && @$cr>0) {
     @ex = sort { $a->[0] <=> $b->[0] } @$cr;
     $CDexons=1;
     }
    else {
     @ex = sort { $a->[0] <=> $b->[0] } @$er;
     #if (@$cr>0) { # sort cds segments too
     #  @cds= sort { $a->[0] <=> $b->[0] } @$cr;
     #  }
     }
   # get the more accurate version of the start-end coords for the feature
   #($mstart, $mend) = ($ex[0]->[0], $ex[-1]->[1]);
   my @exonlst = map { $_->[0].'-'.$_->[1] } @ex;
   my ($avg, $std)=fetchCons($chr, \@ex);
   substr($gffid, 0, length($chr)+1)='';
   print join("\t", $gffid, $chr.$strand.':'.join(',',@exonlst),$avg, $std)."\n";
  } #for each stored transcript
}


# -- fetch conservation data for a sorted set of exons
sub fetchCons {
 my ($chr, $exr)=@_; # $exr  is \@exons (sorted)
 my ($tstart, $tend)=($$exr[0]->[0], $$exr[-1]->[1]);
 my $chd=$ch{$chr};
 my ($cfname, $cfstart);
 my $fh;
 if ($chd) {
    ($cfname, $cfstart)=@$chd;
    open($fh, $cfname) || die("Error opening file $cfname!\n");
    binmode($fh);
    }
   else {
    $cfname="$wigbindir/$chr.wigfixbin";
    open($fh, $cfname) || die("Error opening file $cfname!\n");
    binmode($fh);
    my ($tag, $r);
    read($fh, $tag, 4);
    die("Error: $cfname is not a wigbin file (tag=$tag)!\n") unless $tag eq 'WIGB';
    read($fh, $r, 4);
    my @v=unpack('I',$r);
    $cfstart=$v[0];
    #print STDERR "start=$cfstart\n";
    $ch{$chr}=[$cfname, $cfstart];
    }
 my $tspan=($tend-$tstart+1);
 my @rdata = (0) x $tspan;
 if ($tend>=$cfstart) {
     my ($rstart, $rspan)=($tstart, $tspan); # real start and read values
     if ($rstart<$cfstart) {
        $rstart=$cfstart;
        $rspan=($tend-$rstart+1);
        }
     #print STDERR "..seeking to ".(8+2*($rstart-$cfstart))."\n";
     seek($fh, 8+2*($rstart-$cfstart), 0);
     my $rdata;
#     my $rb;
#     my $nr=0;
#     my @v;
#      for (my $i=0;$i<$rspan;$i++) {
#        $rb=read($fh,$rdata,2);
#        last if $rb<2;
#        my @cv=unpack('s',$rdata);
#        push(@v,$cv[0]);
#        $nr++;
#        }
     my $rb=read($fh, $rdata, 2*$rspan);
     $rb>>=1; # number of values read (each value is 2 bytes)
     my @v=unpack('s'.$rb, $rdata);
#
#     print STDERR "rspan=$rspan; unpacking s$rb values: (".join(',',@v).")\n";
     @rdata[$rstart-$tstart .. $rstart-$tstart+$#v]=@v;
     }
 #now compute average&std for exons only
 my ($sum, $vcount);
 foreach my $e (@$exr) {
    $vcount+=$$e[1]-$$e[0]+1;
    map { $sum+=($_/1000) } @rdata[$$e[0]-$tstart .. $$e[1]-$tstart];
    }
 my $avg=$sum/$vcount;
 $sum=0;
 foreach my $e (@$exr) {
    map { $sum+=((($_/1000)-$avg)**2) } @rdata[$$e[0]-$tstart .. $$e[1]-$tstart];
    }
 my $std=sqrt($sum/($vcount-1));
 return(sprintf('%.3f',$avg), sprintf('%.3f',$std));
}
