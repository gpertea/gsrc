#!/usr/bin/perl
use strict;
use Getopt::Std;
my $usage=q/Usage:
gffann.pl [-o <outpuf.gff3>] [-t <fprefix>] [-R] [-m <mrna_gmap.iit>]
  -p <uniref_pmap.iit> [-P <uniref.fa.cidx> <models.gff>]
 
Create a gff file with the same entries like <models.gff> but with
an annotation attribute ("TopHit=...") added to the last field, based
on the best pmap alignment with UniRef proteins.

Prior to running this, the pmap output should be converted with 
gff2iit to the required input file <uniref_pmap.iit>.

If -t option was given, the top 5 protein mappings will be written in gff3 
format into the file <fprefix>.top5pmap.gff3; if -m option was also given
an equivalent gff file with the best mRNA\/EST hits will be created as
<fprefix>.top5gmap.gff3


/;
umask 0002;
getopts('Dp:m:P:o:t:') || die($usage."\n");

my $ftophits=$Getopt::Std::opt_t;
my $outfile=$Getopt::Std::opt_o;
my $debug=$Getopt::Std::opt_D;
my $iitfile=$Getopt::Std::opt_p;
my $miitfile=$Getopt::Std::opt_m;
my $mrefseq=$Getopt::Std::opt_M; #use the -m file as refseq (strong evidence) for annotation
my $protcdb=$Getopt::Std::opt_P;
my $gff=shift(@ARGV) || die("$usage\n");

#my %dup; #to avoid duplicate pmap/gmap entries

#-------- local array holdin uniformative annotation patterns:
my @uninformative=(
'\bunknown\b',
#'unknown protein\b',
'\bhypothetical\b',
'[A-Z]+\d+ cDNA sequence',
'cDNA sequence [A-Z]+\d+',
'uncharacterized protein',
'unnamed protein product',
'open reading frame',
'\borf\b',
'\bputative\b',
'\bhomologue\b',
'\bsimilar to',
'^expressed sequence \S+$',
'\bHA\d{4}\b',
'\bDKFZP\S+\b',
'PROTEIN FOR MGC:\d+',
'PROTEIN FOR IMAGE:\d+',
'\bR\d{5}\_\d\b',
'\bPRO\d{4}\b',
# 'KIAA\d+ GENE PRODUCT',
#'KIAA\d+ PROTEIN',
#'\bKIAA\d+\b',
'\bHSPC\d+\b',
# HSPC\d+ PROTEIN
#'\bC\d+ORF\d+\b',
'FLJ\d+ PROTEIN',
'\bDJ\d+[A-Z]\d+(\.\d+)*',
'NOVEL PROTEIN',
'CG\d+ PROTEIN',
'CG\d+ GENE PRODUCT',
'^\s*CG\d+\s*$',
'CGI\-\d+ PROTEIN',
'CGI\-\d+',
'CDNA:? FLJ\d+ FIS, CLONE \w+',
'BA\d+[A-Z]\d+[A-Z]?\.\d(\.\d)?',
#'\bRIKEN CDNA .{10} GENE\b',
'\bRIKEN.+?CDNA\b',
'MRNA, COMPLETE CDS, CLONE:\d+(\+\d[A-Z])?\-\d+',
'MRNA, COMPLETE CDS, CLONE:SMAP\d+\-\w+',
'BRAIN CDNA, CLONE MNCB-\d+',
'.{10}RIK PROTEIN',
'^MY\d{3}\s*$',
'MY\d{3} PROTEIN^',
'^probable\b',
'BRAIN MY\d{3}$',
'NPD\d{3} PROTEIN',
'[A-Z]\d{2}[A-Z0-9]+\.\d+ PROTEIN',
'WUGSC:H_\w+\.\w+ PROTEIN',
#'DNA SEGMENT, CHR [0-9XY]+, WAYNE STATE UNIVERSITY \d+, EXPRESSED',
#'DNA SEGMENT, CHR [0-9XY]+, KL MOHLKE \d+',
#'DNA SEGMENT, CHR [0-9XY]+, BAYLOR \d+',
'\bDNA SEGMENT\b',
'PROTEIN HSPC\d+',
#'HYPOTHETICAL [\.\d]+\s*KDA PROTEIN \S+ IN CHROMOSOME \S+',
'EG:[0-9A-Z\.]+ PROTEIN',
'GENOMIC DNA, CHROMOSOME \d+, P1 CLONE:\S+',
'[^,]+, RIKEN FULL-LENGTH ENRICHED LIBRARY, CLONE:.{10}, FULL INSERT SEQUENCE',
'ZK\d+\.\d+ PROTEIN',
'\bEST \w+',
'B2 ELEMENT'
);

open(GFF, $gff) || die ("Error: cannot open $gff\n");

#expects the models to be given in the proper order
# 'CDS' exons are the only ones used;
# the first 'mRNA' entry is going to be annotated
my @linebuf; #current gff gene/mRNA line buffer
my ($gchr, $gstrand, @m_cds, $m_cdslen) ;  #current predicted mRNA info
my $curmodel;

if ($outfile) {
  open(OUTF, '>'.$outfile) || die ("Error creating $outfile\n");
  select OUTF;
  }

if ($ftophits) {
  open(TOPHIT, ">$ftophits.top5pmap.gff3") || die ("Error creating $ftophits.top5pmap.gff3\n");
  if ($miitfile) {
    open(MTOPHIT, ">$ftophits.top5gmap.gff3") || die ("Error creating $ftophits.top5gmap.gff3\n");
    }
  }
#print STDERR "gffann.pl processing $gff ..\n";
my $gffline=1;
while ($gffline) {
 $gffline=<GFF>;
 if (defined($gffline) && ($gffline=~m/^\s*#/ || $gffline=~m/^\s+$/ || length($gffline)<4)) {
   #print $_;
   next;
   };
 chomp($gffline);
 my ($chr, $v, $f, $fstart, $fend, $fscore, $strand, $frame, $info)=split(/\t/, $gffline);
 next if ($f eq 'gene');# skip useless 'gene' entries and other nonsensical lines
 if ($f eq 'mRNA' || !defined($gffline)) { #new record starting, or EOF
   my $line=$gffline;
   if (@m_cds>0) { #previous model's exons are loaded
     @m_cds= sort { $main::a->[0] <=> $main::b->[0] } @m_cds;
     my @phits;
     my ($model_id)=($linebuf[0]=~m/ID=([^;]+)/);
     my ($pann, $pgeneid, $p_mcov, $p_xgaps, $alt_pann, $alt_pgeneid, $alt_pcov, $alt_pxgaps)
                    =annModel($iitfile, \@phits); #protein annotations
     my @mhits;
     my ($mann, $mgeneid, $m_mcov, $m_xgaps, $alt_mann, $alt_mgeneid, $alt_mcov, $alt_mxgaps)
                    =annModel($miitfile, \@mhits, 1) if $miitfile; #mrna/refseq annotations
     if ($p_mcov<2) { 
         if ($m_mcov<2) {
            print STDERR "Warning (gffann.pl $gff): model $model_id has no protein or mRNA coverage!\n";
            }
        else {
            print STDERR "Warning (gffann.pl $gff): model $model_id has no protein coverage!\n";
            }
        }
     
     #--
     my ($ann, $geneid,   $mcov,   $qxgaps,  $altann,   $altgid,      $altcov,  $altxgaps) = 
                ($mgeneid && ($m_mcov>=90 || $p_mcov-$m_mcov<10)) ? 
         # if we have a decent RefSeq annotation, use it
        ($mann, $mgeneid, $m_mcov, $m_xgaps, $alt_mann, $alt_mgeneid, $alt_mcov, $alt_mxgaps) : 
        ($pann, $pgeneid, $p_mcov, $p_xgaps, $alt_pann, $alt_pgeneid, $alt_pcov, $alt_pxgaps);
     unless ($altann) {
       ($altann, $altgid, $altcov, $altxgaps) = $alt_mann ? ($alt_mann, $alt_mgeneid, $alt_mcov, $alt_mxgaps) : 
                                                            ($alt_pann, $alt_pgeneid, $alt_pcov, $alt_pxgaps);
       }
#      if ($mcov<95) {
#        if ($mcov<$m_mcov) {
#         $ann=$mann; $mcov=$m_mcov; $qxgaps=$m_xgaps;
#         }
#        elsif ($mcov<$p_mcov) {
#         $ann=$pann;$mcov=$p_mcov; $qxgaps=$p_xgaps;
#         }
        print STDERR "Warning (gffann.pl $gff): possible merge? model $model_id coverage by best hit is only $mcov\%!\n" 
           if $mcov<80 && $mcov>1;
#       }
     if ($mcov>10) {
       $linebuf[0].=";mcov=$mcov";
       $linebuf[0].=";geneId=$geneid" if $geneid;
       $linebuf[0].=";qxgap=$qxgaps" if $qxgaps;
       if ($ann) {
         my @s=split(/\x01/, $ann);
         $ann=$s[0] if @s>1;
         $linebuf[0].=";tophit=\"$ann\"";
         }
       $linebuf[0].=";altGeneId=$altgid" if $altgid;
       if ($altann) {
         my @s=split(/\x01/, $altann);
         $altann=$s[0] if @s>1;
         $linebuf[0].=";altTopHit=\"$altann\"";
         $linebuf[0].=";altCov=\"$altcov\"";
         $linebuf[0].=";altQXGap=\"$altxgaps\"" if $altxgaps;
         }
       }
     # -- check the model exons values for exons with 0 evidence count
     foreach my $d (@m_cds) {
       #if ($$d[3]==0) {
      $linebuf[$$d[4]].=';mappingEvCount='.int($$d[3]);
      #if ($$d[3]==0) {
      # print STDERR "Warning: exon $$d[0]-$$d[1] of model $model_id ($geneid) has no mapping evidence registered!\n";
      # }
     }
     print join("\n",@linebuf)."\n";
     
     @m_cds=();$m_cdslen=0;
     @linebuf=();
     @phits=();
     @mhits=();
     $gstrand=$strand;
     $gchr=$chr;
     }
   last unless $gffline;
   chomp($line);
   push(@linebuf, $line);
   next;
   }
 chomp($gffline);  
 if ($f eq 'CDS') { # model exon
   ($curmodel)=($info=~m/Parent=(['"\:\w\|\-\.]+)/);
   $gstrand=$strand;
   $gchr=$chr;
   ($fstart, $fend)=($fend, $fstart) if ($fstart>$fend);
   push(@m_cds, [$fstart, $fend, $frame, 0, scalar(@linebuf)]);
   $m_cdslen+=$fend-$fstart+1;
   }
 push(@linebuf, $gffline);
} #while <GFF>
 
close(GFF);

# if (@m_cds>0) { #final record
#  my ($ann, $geneid, $mcov)=annModel($iitfile);
#  annModel($miitfile, 1) if $miitfile; #discard mrna annotation, just build the top5 file
#  $linebuf[0].=";GeneId=\"$geneid\"" if $geneid;
#  $linebuf[0].=";TopHit=\"$ann\"" if $ann; 
#  print join("\n",@linebuf)."\n";
#  }


sub annModel {
 my ($iitdb, $hits, $mrna)=@_;
 return '' unless (@m_cds>0);
 my ($mstart, $mend)=($m_cds[0]->[0], $m_cds[-1]->[1]);
 #get all overlapping features from iit database
 print STDERR "iit_get $iitdb $mstart $mend\n" if $debug;
 open(IITGET, "iit_get $iitdb $mstart $mend |") ||
  die("Error opening pipe from: iit_get $iitdb $mstart $mend \n");
 #my $maxscore=0;
 #map { $maxscore+= $_->[1]-$_->[0]+1 } @m_cds;
 my ($l_id, $l_exon, $l_cds, $l_info, $l_cov, $qgaps, $track);
 my $line=1;
 while ($line) {
   $line=<IITGET>;
   if (!$line || $line=~m/^>\S+\s/) {
     # header line -- new record
     
     if ($l_id) { ### process the previous record's lines
        $l_exon=$l_cds unless $l_exon;
        my @intv=split(/\,/, $l_exon);
        my @xdata;
        my ($ovlscore, $ovlen)=checkOverlap(\@intv, $l_id, \@xdata, $l_exon);
        if ($ovlscore>0) { #has exon overlap
           my $pid=sprintf('%4.1f',($ovlscore*100.0)/$ovlen);
           my $mcov=sprintf('%d', ($ovlen*100)/$m_cdslen); #model coverage by this hit
           my @cdata; #cds overlap data, if any
           #if ($l_cds) {
           # my @cintv=split(/\,/,$l_cds);
           # my ($cdsovlscore, $cdsovlen)=checkOverlap(\@cintv, $l_id, \@cdata, $l_cds);
           # @cdata=() unless $cdsovlen>10;
           # }
           my $geneid;
           if ($l_info=~m/(?:gid|Gene)\:(\S+)/) {
               $geneid=$1;
               $geneid=~s/\]$//;
               }
              else {
               my $n=$l_id;
               $n=~s/\.[pmrnaexobltsi]+\d+$//;
               #$n=~s/\.mrna\d+$//;
               $geneid=$1 if $n=~m/\|gid\|(.+)$/;
               }
           push(@$hits, [$l_id, $ovlscore, [@xdata], [@cdata], $l_cov, $l_info, $pid, $geneid, $mcov, $track, $qgaps]);
           print STDERR "..adding hit: $l_id, $ovlscore, $l_cov, $l_info, $pid, $geneid, $mcov, $track, $qgaps\n" if $debug;
           }
        #-----
        }
     last unless $line; #end of <IITGET> ?
     ($l_exon, $l_cds, $l_info, $l_cov, $track, $qgaps)=();
     my @d=split(/\s+/, $line);
     if (substr($d[3], -1) ne $gstrand) {
      $l_id=undef;
      next; #not a valid match (opposite strand)
      }
     # get the name of this
     $track=$d[-1];
     $track=~s/gblat/blat/;
     ($l_id)=($line=~m/^>(\S+)/);
     next;
     } #>header line
 next unless $l_id; #skip record
 print STDERR ">parsing data for hit $l_id :\n" if $debug;
 print STDERR "       $line\n" if $debug;
 if ($line=~m/^([iCvx])\:(.+)/) {
    my ($t, $d)=($1,$2);
    if ($t eq 'C') { $l_cds=$d; }
      elsif ($t eq 'v') { $l_cov=$d; }
      elsif ($t eq 'i') { $l_info=$d; 
         $l_info=~s/(Rec|Sub)Name: Full=//;
         $l_info=~s/\. Flags: [^\{]+/. /;
         $l_info=~s/ EC=[\d\-\._]+//g;
         $l_info=~s/\. AltName: [^\{]+/. /;
         $l_info=~s/\. Short=[^\{]+/. /;
         }
      elsif ($t eq 'x') { $qgaps=$d; # print STDERR "added qxgaps: $qgaps\n" if $debug;}
                            }
    }
   elsif ($line=~m/^\d+\-\d+/){
    $l_exon=$line;
    }
 } #while <IITGET>

 close(IITGET);
 return getBestHits($hits, $mrna);
}

if ($outfile) {
 select STDOUT;
 close(OUTF);
 }

if ($ftophits) {
 close(TOPHIT);
 close(MTOPHIT) if $miitfile;
 }

sub checkOverlap {
 my ($intv, $id, $xdata, $l)=@_;
 my ($ovlscore, $ovlen)=(0,0);
 foreach my $iv (@$intv) {
    my ($xstart, $xend, $xscore)=($iv=~m/^(\d+)\-(\d+)\:(\d+)/);
    #$xscore='89' if $xscore eq '.';
    die ("Error parsing $id exon line: $l\n") unless $xstart>0 && $xstart<=$xend;
     my ($ovlsc, $ovll)=&checkExonOvl($xstart, $xend, $xscore);
     push(@$xdata, [$xstart, $xend, $xscore, $ovll]);
     $ovlscore+=$ovlsc;
     $ovlen+=$ovll;
     }
 return ($ovlscore, $ovlen);
}

#---
sub checkExonOvl {
 my ($cstart, $cend, $cscore)=@_;
 # uses global variable @cds for the current model's exons to scan
 # returns the score sum of (overlap_length * $cscore/100 for each exon) 
 my $ovlscore=0;
 my $covlen=0; #cumulative bases coverage
 foreach my $d (@m_cds) { #for each model's exon
   last if $$d[0]>$cend;
   my $covl;
   #get the overlap size:
   if ($$d[0]<$cstart) {
      if ($cstart<$$d[1]) {
        $covl=($cend>$$d[1])? ($$d[1]-$cstart+1) : ($cend-$cstart+1);
        }
      }
     else {
      if ($$d[0]<$cend) {
        $covl=($cend<$$d[1])? $cend-$$d[0]+1 : $$d[1]-$$d[0]+1;
        }
      }
   if ($covl) { #there was overlap with this exon
     $$d[3]++; # increase evidence overlap counter for this model's exon
     $ovlscore += (($covl*$cscore)/100.00);
     $covlen+=$covl;
     }
   }# for each exon of the model
 return ($ovlscore, $covlen);
}


sub ovlTest {
 my ($a1, $a2, $b1, $b2)=@_;
 return ($a1<$b1) ? $b1<$a2 : $a1<$b2;
}

sub getBestHits {
 my ($href, $mrna)=@_;
 #sort by decreasing score
 my @hits = sort { $main::b->[1] <=> $main::a->[1] } @$href;
 $href=[@hits];
 #print STDERR ">$curmodel $gstrand ($mstart .. $mend) :\n" if $debug && !$mrna;
 my ($firsthid, $firstdescr, $firstgid, $first_mcov, $first_qxgaps);
 # uses cdbyank -F -a $pID $protcdb to get the defline of the protein $$href[0]->[0]
 # skips to next best hit until isInformative confirms it
 my ($besthid, $besthit, $best_gid, $best_mcov, $best_qxgaps, 
     $alt_hid, $alt_hit, $alt_gid, $alt_mcov, $alt_qxgaps);
 my $hcount=0;
 my @fregion; #first region covered by the best hit
 foreach my $h (@hits) {
    #my ($prot, $pscore, $pmapr,  $cov, $descr, $mcov) = @$h;
    my ($hid, $ovlscore, $xmapr, $cmapr, $cov, $descr, $pid, $geneid, $mcov, $track, $qgaps)=@$h;
    my $qid=$hid;
    $qid=~s/\.[mrnapbltexosi]+\d+$//; 
    #$qid=~s/\.p[mp]\d+$//;
    my @region=($xmapr->[0]->[0], $xmapr->[-1]->[1]); #start-end coordinates for this hit
    if (@fregion && $firstdescr) { # is this another region?
      unless (&ovlTest(@region, @fregion)) {
        #alternate partial annotation
        $alt_hid=$qid;
        $alt_hit=$descr;
        $alt_gid=$geneid;
        $alt_mcov=$mcov;
        $alt_qxgaps=$qgaps;
        }
      }
     else {
      @fregion=@region;
      }
   #----
   if ($mrna) { #mrna hit
      unless ($firstdescr) {
        $firsthid=$qid;
        $firstdescr=$descr; 
        $firstgid=$geneid;
        $first_mcov=$mcov;
        $first_qxgaps=$qgaps;
        }
      unless ($besthit) {
        if (&isInformative($descr)) {
             $besthid=$qid;
             $besthit=$descr;
             $best_gid=$geneid;
             $best_mcov=$mcov;
             $best_qxgaps=$qgaps;
             }
        }
      if ($ftophits && $hcount<5) {
        my @pex=@$xmapr;
        print MTOPHIT join("\t", $gchr, $track, 'mRNA', $pex[0]->[0], $pex[-1]->[1], '.',
                $gstrand, '.', "ID=$hid;Name=$hid");
        print MTOPHIT ";cov=$cov" if ($cov);
        print MTOPHIT ";gene=$geneid" if $geneid; # shouldn't this be Name instead?
        print MTOPHIT ";qxgap=$qgaps" if $qgaps;
        print MTOPHIT ";info=\"$descr\"" if $descr;
        
        print MTOPHIT "\n";
        my $c=0;
        foreach my $px (@pex) {
         print MTOPHIT join("\t", $gchr, $track, 'exon', $px->[0], $px->[1], $px->[2],
                $gstrand, '.', "Parent=$hid")."\n";
         $c++;
         }
        $hcount++; 
        }
      } #mRNA
     else { #protein hit
      my $defline = $descr;
      if ($protcdb && !$defline) {
        $defline=`cdbyank -F -a '$qid' $protcdb`;
        printf STDERR "  (%04d) %s", $ovlscore, $defline if $debug;
        chomp($defline);
        $defline=~s/^\S+\s+//;
        $defline=~s/^\[\d+\]\s+//;#my old Uniref entries
        $defline=~s/( \- )/ \[$hid\] /;
        }
      unless ($firstdescr) {
        $firsthid=$qid;
        $firstdescr=$defline; 
        $firstgid=$geneid;
        $first_mcov=$mcov;
        $first_qxgaps=$qgaps;
        }
      unless ($besthit) {
        if (&isInformative($defline)) {
             $besthid=$qid;
             $besthit=$defline;
             $best_gid=$geneid;
             $best_mcov=$mcov;
             $best_qxgaps=$qgaps;
             }
        }
      if ($ftophits && $hcount<5) {
        my @pex=@$xmapr;
        print TOPHIT join("\t", $gchr, $track, 'mRNA', $pex[0]->[0], $pex[-1]->[1], '.',
                $gstrand, '.', "ID=$hid;Name=$hid");

        print TOPHIT ";cov=$cov" if $cov;        
        print TOPHIT ";gene=$geneid" if $geneid;
        print TOPHIT ";qxgap=$qgaps" if $qgaps;
        print TOPHIT ";info=\"$defline\"" if $defline;
        print TOPHIT "\n";
        my $c=0;
        foreach my $px (@pex) {
         print TOPHIT join("\t", $gchr, $track, 'exon', $px->[0], $px->[1], $px->[2],
               # $gstrand, '.', "ID=$prot.e$c;Parent=$prot")."\n";
                $gstrand, '.', "Parent=$hid")."\n"; #optimized for argo applet
         $c++;
         }
        $hcount++; 
        }
       }
 } #foreach hit
 unless ($besthit) {
  $besthit=$firstdescr;
  $besthid=$firsthid;
  $best_gid=$firstgid;
  $best_mcov=$first_mcov;
  $best_qxgaps=$first_qxgaps;
  }
 $alt_hit=$alt_hid.' '.$alt_hit if $alt_hit;
 $besthit=$besthid.' '.$besthit if $besthit;
 return ($besthit, $best_gid, $best_mcov, $best_qxgaps,
        $alt_hit, $alt_gid, $alt_mcov, $alt_qxgaps);
}



#===============================================
# bool isInformative($description) 
# expects only the descripts - not the accession
#===============================================
sub isInformative {
 local $_=$_[0];
 s/^\s+//g;s/\s+$//g;
 return 0 if length($_)<2;
 foreach my $pat (@uninformative) {
   if (m/$pat/i) {
     #&flog("uninformative by /$pat/i : '$_'") if ($debug);
     return 0;
     }
   }
return 1;
}
