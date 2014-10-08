#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Reconcile jigsaw results on separate strands.
Usage:
 recon_jigsaw.pl [-P <gseq.fa>] [-o <output.gff3>] <both_strands_jigsaw.jgff> \
                 <forward_strand.jgff> <reverse_strand.jgff>
 
 Requires the raw jigsaw output from the 3 jigsaw runs on the same
 genomic sequence.
 Outputs a gff3 file with the resulting "reconciliation".
 
 Use the -P option to discard any partial gene predictions made by jigsaw.
 <gseq.fa> is a fasta providing the genomic sequence required for checking
 start&stop codons. 
/;
umask 0002;
getopts('P:Do:') || die($usage."\n");
my $seqfile=$Getopt::Std::opt_P; #remove partial predictions
my $seq='';
my $outfile=$Getopt::Std::opt_o;
my $debug=$Getopt::Std::opt_D;
die($usage."\n") unless @ARGV==3;

my ($bfile, $ffile, $rfile)=@ARGV;
my %stopCodons=('TAA'=>1,'TAG'=>1,'TAR'=>1, 'TGA'=>1, 'TRA'=>1);
my $startCodon='ATG';

#foreach ($bfile, $ffile, $rfile) {
die ($usage."File $_ not found or empty! \n") unless -s $bfile; 
# }

if ($seqfile) { #load sequence in memory
 open(SEQFILE, $seqfile) || die ("Error opening fasta file '$seqfile'!\n");
 while (<SEQFILE>) {
  if (m/^>/) {
    last if length($seq)>10;
    next;
    }
  tr/\n\r \t//d;
  $seq.=$_;
  }
 close(SEQFILE);
 }

my (@bgenes, @fgenes, @rgenes);
## 
##  @genes= (generef1, generef2, ... )
##    where a generef := [ predset,  strand, start, end, \@exons, cluster, total_exonlen, score]
##                           0        1       2      3      4       5         6            7
##        where @exons := list of [exonstart, exend, frameshift]
##        and  predset := jigsaw prediction set code ( represented as one of 
##                     the letters 'b', 'f', 'r', meaning one of: 
##                     both strands at once, forward-only, reverse-only
#
## a "cluster" is a reference to a list of coordinates followed by overlapping loci:
##   $cluster := [ $cstart, $cend, $generef1, $generef2 ... ]

my %clusters; #just a set of references to formed clusters

my $genomicseq;

$genomicseq=&readJigsaw($bfile, \@bgenes, 'b');
&readJigsaw($ffile, \@fgenes, 'f');
&readJigsaw($rfile, \@rgenes, 'r');



## -- now cluster the loci and check them for conflicts:
## for now only consider predictions which overlap
## the "both-strands" prediction set
for (my $bl=0;$bl<@bgenes;$bl++) {
 for (my $fl=0;$fl<@fgenes;$fl++) {
   last if $fgenes[$fl]->[2]>$bgenes[$bl]->[3];
   # check Gene overlap and update cluster if so
   &processGovl($bgenes[$bl], $fgenes[$fl]);
   }
 for (my $rl=0;$rl<@rgenes;$rl++) {
   last if $rgenes[$rl]->[2]>$bgenes[$bl]->[3];
   # check Gene overlap and update cluster if so
   &processGovl($bgenes[$bl], $rgenes[$rl]);
   }
}


# now %clusters contains the loci clusters (except for 'f' or 'r' singletons we don't care about)
my @cls=values(%clusters);
@cls=sort { $main::a->[0] <=> $main::b->[0]  } @cls;
my $clno=1;
my $gno=1;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die ("Error creating file $outfile !\n");
  select(OUTF);
  }
foreach my $c (@cls) {
  $c=&processCluster($c);

  my @cl=@$c;
  my $clstart=shift(@cl);
  my $clend=shift(@cl);
  #only generefs in @cl now:
  
  print STDERR "Locus Cluster #$clno ($clstart-$clend)\n" if $debug;
  foreach my $gr (@cl) {
   print STDERR "gene #$gno:\t".join("\t",@{$gr}[0..3], $$gr[7])."\n" if $debug;
   #if ($$gr[7]>0) {
   printGff($genomicseq, $gno, $gr);
   $gno++;
   #  }
   }
  $clno++;
 }

if ($outfile) {
 close(OUTF);
 select(STDOUT);
 }

##### ============== SUBROUTINES ============ #####


sub printGff {
 my ($srcseq, $gno, $g)=@_;
 my $curmodel="tjsm.$gno";
 my @xflags;
 push(@xflags, 'NoStartCodon') if $$g[8];
 push(@xflags, 'NoStopCodon') if $$g[9];
 my @xsplit;
 if ($$g[10]) {
   push(@xflags, 'ExonGap='.$$g[10]);
   @xsplit=split(/:/,$$g[10]);
   }
 my $xfl='';
 if (@xflags) {
   $xfl=';'.join(';',@xflags);   
   }
 print join("\t", $srcseq, 'jigsaw', 'mRNA', $$g[2], $$g[3], '.', $$g[1], '.',
            "ID=$curmodel;Name=$curmodel".$xfl)."\n";

 my $exons=$g->[4];
 my $i=1;
 foreach my $exon (@$exons) {
   my ($estart, $eend, $eframe, $splitframe)=@$exon;
   my $xattr="Parent=$curmodel";
   $xattr.=';ExonSplit' if ($estart==$xsplit[1]);
   $xattr.=';SplitFrame' if $splitframe;
   print join("\t",$srcseq, 'jigsaw', 'CDS', $estart, $eend, '.', 
             $$g[1], $eframe, $xattr)."\n";
   $i++;
   }

}


sub processCluster {
 my $c=shift(@_);
  my @cl=@$c;
  my $clstart=shift(@cl);
  my $clend=shift(@cl);
  # check for exon collisions
  #my %removal; # list of genes to remove due to exon conflicts
  #my %score; # $generef => score based on overlap with other model on the same strand
  foreach my $g (@cl) {
   $$g[7]+=$$g[6] if length($$g[0])>1;
   }
  my %removal; # list of genes to remove due to exon conflicts
  for (my $i=0;$i<@cl;$i++) {
    my $gx=$cl[$i];
    for (my $j=$i+1;$j<@cl;$j++) {
      my $gy=$cl[$j];
      # CDS lengths and exon overlap bps
      my $exovl=&exonOverlap($gx, $gy);
      if ($exovl>=0 && ($gx->[1] eq $gy->[1])) { #overlapping, on the same strand: must choose one!
        if ($$gx[6]+$$gx[7] > $$gy[6] +$$gy[7]) {
           # gx "better" than gy
           $$gx[7]+=$exovl;
           $removal{$gy}=$gy;
           }
          else {
           $$gy[7]+=$exovl;
           $removal{$gx}=$gx;
           }
        } #overlapping on the same strand
      } #inner loop
    } #outer loop
   # second pass, to remove conflicting models:  
    for (my $i=0;$i<@cl;$i++) {
    my $gx=$cl[$i];
    next if $removal{$gx};
    for (my $j=$i+1;$j<@cl;$j++) {
      my $gy=$cl[$j];
      next if $removal{$gy};
      # CDS lengths and exon overlap bps
      if ($gx->[1] ne $gy->[1] && &exonOverlap($gx, $gy)) { #opposite strands
        #opposite strands -- one of them should be thrown away
        my $grm = ($$gx[6]+$$gx[7]>$$gy[6]+$$gy[7]) ? $gy : $gx;
        $removal{$grm}=$grm;
        } #opposite strands
      } #inner for loop
    } #outer for loop
  my @newcl= grep { ($_->[7]>0 || ($_->[7]==0 && $_->[0] eq 'b')) && not exists($removal{$_}) } @cl;
  #my @newcl= grep { ($_->[7]>0 && not exists($removal{$_}) } @cl;
  #DBG#:
  foreach my $g (values(%removal)) {
   if ($$g[0] =~ m/b/) {
     print STDERR "W: eliminated base model at $$g[2]-$$g[3]($$g[1]) from $bfile\n";
     }
   }
  return [$clstart, $clend, @newcl]
}


sub exonOverlap { #returns the total length of overlapping exonic regions
 my ($ga, $gb)=@_;
 return -1 unless ( $$ga[2]<=$$gb[2] ? $$gb[2]<=$$ga[3] : $$ga[2] <= $$gb[3] );
 my ($ea, $eb)=($ga->[4], $gb->[4]);
 my $ovl=0;
 foreach my $a (@$ea) {
   foreach my $b (@$eb) {
    if ($$a[0]<=$$b[0]) {
     if ($$b[0]<=$$a[1]) {
       $ovl += ( $$b[1]<$$a[1] ? $$b[1]-$$b[0]+1 : $$a[1]-$$b[0]+1 );
       }
     }
    else {
     if ($$a[0]<=$$b[1]) {
       $ovl += ( $$b[1]<$$a[1] ? $$b[1]-$$a[0]+1 :  $$a[1]-$$a[0]+1 );
       }
     }
   } #for b
  } # for a
 return $ovl;
}

sub processGovl {
 my ($ga, $gb)=@_;
 my ($ovl, $min, $max);
 #in fact, check if their cluster ends overlap
 my ($ca, $cb)=($ga->[5], $gb->[5]);
 return if ($ca eq $cb); # already in the same cluster
 if ($$ca[0]<=$$cb[0]) {
    $min =  $$ca[0] if ($ovl=($$cb[0]<=$$ca[1]));
    }
   else {
    $min =  $$cb[0] if ($ovl = ($$ca[0]<=$$cb[1]));
    }
 return unless $ovl;
 #overlap: update cluster data
 $max = ($$ca[1]>$$cb[1]) ? $$ca[1] : $$cb[1];
 $ca->[0]=$min;
 $ca->[1]=$max;
 
 #we know that one of them is SURELY 'b'
 #collapse them if that is the case

 # add the $cb's loci here to this cluster and update their cluster pointer
 my @bg=@$cb;shift(@bg);shift(@bg);
 my $sameLocus=&sameLocus($ga, $gb);
 $ga->[0].=$gb->[0] if $sameLocus; 
 delete($clusters{$cb}); #-- only 'b'-based clusters were stored, but they could get merged
 foreach my $g (@bg) {
   unless ($sameLocus && $g eq $gb) {
     push(@$ca, $g);
     $g->[5] = $ca;
     }
   }
}


sub sameLocus {
 my ($ga, $gb)=@_; 
 return 0 if $$ga[1]!=$$gb[1] || $$ga[2]!=$$gb[2] || $$ga[3]!=$$gb[3];
 my $warn="WARNING ($bfile): $$ga[0] and $$gb[0] predictions at $$ga[2]-$$ga[3] do not have the same exon structure!\n";
 my ($ea, $eb)=($$ga[4], $$gb[4]);
 if (@$ea!=@$eb) {
   print STDERR $warn;
   return 0;
   }
 for (my $i=1;$i<@$ea;$i++) {
   if ($$ea[$i]->[0]!=$$eb[$i]->[0]  || 
       $$ea[$i]->[1]!=$$eb[$i]->[1]) {
         print STDERR $warn;
         return 0
         }
   }
 return 1;
}

sub readJigsaw {
 my ($file, $aref, $predset)=@_;
 return unless -s $file;
 local *FIN;
 open(FIN, $file) || die ("Error opening file $file\n");
 
 my @exd;
 my ($curgnum, $cstrand);
 my $curchr;
 while (<FIN>) {
  next if m/^\s*#/;
  chomp;
  next unless length($_)>4;
  my ($chr, $jsver, $ftype, $exonstart, $exonend, $jscore, 
      $strand, $frame, $lnum)=split(/\t/);

  my $gnum;
  if ($lnum=~m/^(\d+)$/) { #old version
    $gnum=$1;
    }
   else {
    my @a=split(/\s*;\s*/, $lnum);
    ($gnum)=($a[0]=~m/\.(\d+)$/);
    }

  die("Invalid strand ($strand) found for expected '$predset' type locus")
   if ($predset eq 'f' && $strand ne '+') || ($predset eq 'r' && $strand ne '-');
  ($exonstart, $exonend) = ($exonend, $exonstart) if $exonstart>$exonend;
  die("Error at parsing jigsaw line: $_\n") unless ($jsver=~/^jigsaw/ && $gnum);
  if ($curgnum ne $gnum) {
    &storeGene($aref, \@exd, $predset, $cstrand, $curgnum, $file) if $curgnum;
    $curgnum=$gnum;
    $cstrand=$strand;
    $curchr=$chr;
    @exd=();
    }
  push(@exd, [$exonstart, $exonend, $frame]);
  }#while reading jgff lines
 close(FIN);
 &storeGene($aref, \@exd, $predset, $cstrand, $curgnum, $file) if $curgnum;
 # sort genes by coordinates 
 @$aref = sort { $main::a->[2] <=> $main::b->[2] } @$aref;
 return $curchr;
}


sub storeGene {
 my ($gr, $er, $pset, $strand, $mgnum, $mfile)=@_;
 my @ex= sort { $main::a->[0] <=> $main::b->[0] } @$er;
 my ($nostart, $nostop);
 if ($seq) { #validate that it has start and stop codons
   my ($cstart, $cend) = ($strand eq '-') ? 
        ( uc(revCompl(substr($seq, $ex[-1]->[1]-3,3))), uc(revCompl(substr($seq, $ex[0]->[0]-1,3))) ) :
        ( uc(substr($seq, $ex[0]->[0]-1,3)), uc(substr($seq, $ex[-1]->[1]-3,3)) );
   print STDERR "$mgnum ($mfile) : start codon: $cstart  | stop codon: $cend\n" if $debug;
   if ($cstart ne 'ATG') {
     my $discard =  (@ex<3) ? 'discarded ':''; 
     print STDERR "..${discard}partial gene $mgnum ($mfile): lacking START codon ($cstart)\n";
     return if $discard;
     $nostart=1;
     }
   unless ($stopCodons{$cend}) {
     my $discard =  (@ex<3) ? 'discarded ':''; 
     print STDERR "..${discard}partial gene $mgnum ($mfile): lacking STOP codon ($cend)\n";
     return if $discard;
     $nostop=1;
     }
   }
 my $exlen=0;
 
 my ($xgap, $nex) = exonGapCheck(\@ex, $strand, $mgnum, $mfile) if $seq;
 if ($xgap) {
   @ex=@$nex;
   }
 map { $exlen += $_->[1]-$_->[0]+1 } @ex;
 
 #           0       1          2           3            4      5      6    7      8         9      10
 my $gref=[$pset, $strand, $ex[0]->[0], $ex[-1]->[1], [@ex], undef, $exlen, 0, $nostart, $nostop, $xgap];
 my $selfcl=[$ex[0]->[0], $ex[-1]->[1], $gref];
 $gref->[5]=$selfcl;
 $clusters{$selfcl}=$selfcl if $pset eq 'b'; #only store "both-strands" clusters
 push(@$gr, $gref);
 
}


sub exonGapCheck {
 my ($er, $strand, $mgnum, $mfile)=@_;
 #check each exon for N-gaps;
 my @xce; #exons, possibly split, will be stored here
 my @gapdata; # gapstart_gapend
 for (my $i=0;$i<@$er;$i++) {
   my ($x0, $x1, $xfr) = @{$$er[$i]};
   my $xseq=uc(substr($seq, $x0-1, $x1-$x0+1));
   my $gapstart=index($xseq,'NN');
   my $gapend=rindex($xseq, 'NN');
   if ($gapstart>=0) { # Houston, we have a problem: exon gap
     my ($x0e, $x1b)=($x0+$gapstart-1,$x0+$gapend+2);
     push(@gapdata, ($x0e+1).':'.($x1b-1));
     #fix the frame:
     my $xf;
     if ($strand eq '-') {
        #$xf=(($x1-$x0e)+$xfr)%3;
        $xf=(($x1-$x1b+1)+$xfr)%3;
        # print STDERR "exon gap for ($x0,$x1,$xfr) at $x0e .. $x1b\n";
        if ($x0e>$x0) {
           push(@xce, [$x0, $x0e, $xf, 1]);
           #print STDERR "    added split exon: ($x0, $x0e, $xf)\n";
           }
        if ($x1>$x1b) {
           push(@xce, [$x1b, $x1, $xfr]);
           #print STDERR "    added split exon: ($x1b, $x1, $xfr)\n";
           }
        }
       else {        
        $xf=(($x1b-$x0+1)+$xfr)%3;
        #print STDERR "exon gap for ($x0,$x1,$xfr) at $x0e .. $x1b (new frame: $xf)\n";
        push(@xce, [$x0, $x0e, $xfr]) if $x0e>$x0;
        push(@xce, [$x1b, $x1, $xf, 1]) if $x1>$x1b;
        }
     next; 
     }
   push(@xce, [$x0, $x1, $xfr]);
   }
 my $gapds=$gapdata[0]; #should probably be: join(',', @gapdata)
 if (@gapdata>1) {
   print STDERR "WARNING: multiple exons with sequencing gaps found for $mgnum ($mfile)!\n";
   }
 #print STDERR " >> gapdata returned: $gapds\n" if $gapds;
 return ($gapds, \@xce);
}

sub revCompl {
  my $s=reverse($_[0]);
  $s =~ tr/AaCcTtGgUuMmRrWwSsYyKkVvHhDdBb/TtGgAaCcAaKkYyWwSsRrMmBbDdHhVv/;
  return $s;
 }


