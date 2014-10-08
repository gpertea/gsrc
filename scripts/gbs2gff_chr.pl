#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 gbs2gff_chr.pl -p <prot_id2transcript_id-file> [-m <seq_contig.md>] [-M][-T] \
      [-f <mapflt>] [-c genome.fa.fai] [-t <track>][-A][-a <attr1,attr2,..>] \
      <gbs_stream>
 Creates chromosome gff annotation files from NCBI GB data
 
 Takes the contig mappings (onto chromosomes) from <seq_contig.md> (optional),
 filtering it by group_label if requested (e.g. -f Primary_Assembly)
 and converts the input GenBank formatted data into a gff  file.

 Note: only these features are parsed: 
    gene, mRNA, CDS, misc_RNA, ncRNA, rRNA, tRNA
 Currently these are discarded:
   C_region, V_segment
 Options: 
   -A write all attributes
   -a write these attributes (comma delimited list)
   -M parse and write mRNA annotations only
   -T do not write parent "gene" records 
   -K for -c option, skip genome sequences for which the sequence
      length cannot be verified
/;
#   -S write fasta files with the chromosome\/contig sequences

umask 0002;
getopts('p:c:m:f:o:t:a:KTMA') || die($usage."\n");
my $ctgmapfile=$Getopt::Std::opt_m;
my $f_out=$Getopt::Std::opt_o;
my $skip_no_len_check=$Getopt::Std::opt_K;
my $writeseq=$Getopt::Std::opt_S;
my $addattrs=$Getopt::Std::opt_a;
my $mapflt=$Getopt::Std::opt_f;
my $p2nfile=$Getopt::Std::opt_p;
my $allattrs=$Getopt::Std::opt_A;
my $t_only=$Getopt::Std::opt_T; # no gene parent, print transcripts only
my $mrnaonly=$Getopt::Std::opt_M;
my $chrlen_test=$Getopt::Std::opt_c; #must be a SAMTools faidx file
my %chrsizes;
if ($chrlen_test && -f $chrlen_test) {
 open(FAI, $chrlen_test) || die("Error opening $chrlen_test!\n");
 while (<FAI>) {
  chomp;
  next unless $_;
  next if m/^\s*#/;
  my @t=split;
  $chrsizes{$t[0]}=$t[1];
  }
 close(FAI);
 }
die("Error: -p is required -- to map protein_id to transcript_id\n")
   unless $p2nfile && -f $p2nfile;
my %pmap; # protein_id => transcript_id
open(PROTMAP, $p2nfile) || die("Error opening $p2nfile!\n");
while (<PROTMAP>) {
 chomp;
 my ($p, $n)=split;
 next unless $n;
 $p=~s/\.\d+$//;
 $n=~s/\.\d+$//; # discard versioning
 $pmap{$p}=$n;
}
close(PROTMAP);

my %haplotypes = (
'GL000250' => 'chr6_apd_hap1',
'GL000251' => 'chr6_cox_hap2',
'GL000252' => 'chr6_dbb_hap3',
'GL000253' => 'chr6_mann_hap4',
'GL000254' => 'chr6_mcf_hap5',
'GL000255' => 'chr6_qbl_hap6',
'GL000256' => 'chr6_ssto_hap7',
'GL000257' => 'chr4_ctg9_hap1',
'GL000258' => 'chr17_ctg5_hap1'
 );

#                       0      1         2        3         4          5            6
my %ctgmap; # $ctg=> [chr, strand, chr_start, chr_end, ctg_gl_acc,  chr_name,  seq_definition]
if ($ctgmapfile) {
	open(CTGMAP, $ctgmapfile) || die("Error opening $ctgmapfile\n");
	while (<CTGMAP>) {
	 next if m/^#/;
	 my @t=split(/\t/);
	 next if $t[6] eq 'na';
	 next if $mapflt && $t[8]!~/$mapflt/i;
	 my $ctgacc=$t[5];
	 $ctgacc=~s/\.\d+$//; # remove version
	 my $chrname='chr'.$t[1];
	 $chrname=~s/\.\d+$//;
	 my ($chrstart, $chrend)=($t[2], $t[3]);
	 ($chrstart, $chrend)=($chrend, $chrstart) if $chrend<$chrstart;
	 my $chrstrand=($t[4] eq '-') ? '-' : '+';
	 my $ctgdata=[$chrname, $chrstrand, $chrstart, $chrend, '', ''];
	 $ctgmap{$t[5]}=$ctgdata;
	 $ctgmap{$ctgacc}=$ctgdata;
	}
	close(CTGMAP);
 }
#key attributes to always keep :
my %kattrs; @kattrs{'gene', 'product', 'transcript_id', 'codon_start'}=();
my %xattrs; #extra attributes to keep
if ($addattrs) {
 $addattrs=~tr/ //d;
 my @add=split(/\,/,$addattrs);
 @xattrs{@add}=();
 }
my $moreattrs = ( $allattrs || $addattrs );
# ---------------
my %gids; #ensure unique gene IDs
          #allow for duplication $gids{$gene}=[ [gstart1, gend1], [gstart2, gend2] ]
my %tids; # and transcript IDs
#----
# -- per gbase: 

my %genes; # genes{$gene}=[$gbase, $strand, $gstart, $gend, [@gxrefs], \%gattrs,  \%rnas,  isPseudo, [@genesyms] ]
           #                  0        1        2      3        4         5          6        7          8
# %rnas hash holds all mRNA data for a single gene (i.e. all isoforms)
#                             0         1      2         3        4          5         6           7       
# $rnas{transcript_id}=[ $featname, [@exon], [@cds], \%tattrs, [@txrefs], $partial3, $partial5, [@genesyms] ]
#current values:
my ($curctg, $curctglen, $curcmt); # per contig
my $cursection;
my ($curfeat, $curfstrand, @fsegs, %fattrs,  $fpartial5, $fpartial3); # per feature
# $gpseudo,  $gstart, $gend, $partial5, $partial3);
my $last_line;
my $lastattr;
my $track=$Getopt::Std::opt_t || 'ncbi';
#only parse these sections:
my %sparse;@sparse{'FEATURES', 'COMMENT'}=(1,1); #sections for which to parse multiple rows

#only parse these features (WARNING: check this everyonce in a while, NCBI may change these tags
# gzip -cd *.gbs.gz | perl -ne 'print $1."\n" if m/^\s{4,8}(\S+)\s{2,}\S+/' | sort -u
#
my %fparse;
@fparse{'gene','misc_RNA','mRNA','ncRNA', 'CDS', 'rRNA',
        'tRNA', 'snoRNA','miRNA','microRNA','micro_RNA', 'exon',
        'C_region', 'D_segment', 'J_segment', 'V_segment'}=();
my %fskip; @fskip{'source', 'gap', 'repeat'}=();
my $skipctg;
if ($f_out) {
  open(F_OUT, '>'.$f_out) || die("Error creating output file $f_out\n");
  select(F_OUT);
  }
while (<>) {
 LSKIP:
 $last_line=$_;
 if (m/^ORIGIN\s*$/) {
   # read sequence here
   my $ctgdata=$ctgmap{$curctg};
   my ($c_chr,$c_acc, $c_def)=($ctgdata->[0], $ctgdata->[4], $ctgdata->[6]);
   my $defline=">$c_chr $c_acc";
   $defline.=" $c_def" if $c_def;
   my $fname=$c_chr.'.fa';
   flushGSeq();
   if ($writeseq) {
     open(FASEQ, '>'.$fname) || die ("Error creating fasta file $fname!\n");
     print FASEQ $defline."\n";
     }
   my $l;
   while (($l=<>)) {
    $l=~s/[\s\n\r]+$//s;
    last if $l eq '//';
    my ($coord, $seq)=($l=~m/^\s*(\d+)\s+(\w[ \w]*)/);
    $seq=~tr/ \t//d;
    if ($writeseq) {
      print FASEQ uc($seq)."\n";
      }
    }
   close(FASEQ) if $writeseq; 
   next;
   }
 if (m/^([A-Z,0-9]+)\s+(\S.+)/) {
   ($cursection, my $rest)=($1,$2);
   if ($cursection eq 'LOCUS') {
     flushGSeq(); # flush contig features accumulated so far
     my ($ctgacc, $ctglen)=($rest=~m/^(\w+)\s+(\d+)\s+bp/);
     die("Error: cannot parse the LOCUS line:\n$last_line\n") unless $ctglen>1;
     my $cd;
     if ($ctgmapfile) {
        $cd=$ctgmap{$ctgacc};
        unless ($cd) {
          print STDERR "Skipping contig $ctgacc (no chromosome entry)\n";
          $skipctg=1;
          next;
          }
        die("Error: mismatch between contig $ctgacc length and its chromosome positioning!\n")
          unless $ctglen == $cd->[3]-$cd->[2]+1;
        }
     $skipctg=0;
     $curctg=$ctgacc;
     $curctglen=$ctglen;
     next;
     }
    next if $skipctg;
    if ($cursection eq 'DEFINITION') {
     unless ($ctgmapfile) {
       my ($chrnum)=($rest=~m/chromosome\s+([\dXYM]+)/i);
       unless ($chrnum) {
         if (m/ mitochondrion/) {
             $chrnum='M';
             }
           elsif (m/unplaced genomic contig/) {
             $chrnum='Un';
             }
         }
       die("Error parsing chromosome name! ($rest)\n") unless $chrnum;
       $rest=~s/\.$//;
       my $chrid='chr'.$chrnum;
       $ctgmap{$curctg}=[$chrid, '+', 1, $curctglen, '', '', $rest];
       }
     next;
     }
    if ($cursection eq 'VERSION') {
     unless ($ctgmapfile) {
       my ($accver)=($rest=~m/^\s*\w+(\.\d+)/);
       die("Error parsing accession.version! ($rest)\n") unless $accver;
       $ctgmap{$curctg}->[4]=$curctg.$accver;
       }
     next;
     }
    if ($cursection eq 'COMMENT') {
     $curcmt=$rest;
     }
    next;
    } #big section start
 next if $skipctg;
 next unless exists $sparse{$cursection};
 if ($cursection eq 'COMMENT') {
    chomp;
    s/^\s+/ /;
    $curcmt.=$_;
    }
  elsif ($cursection eq 'FEATURES') {
    # current contig information has been parsed by now
    if ($curctg && $curcmt) {
       my $ctgdata=$ctgmap{$curctg};
       die("Error: no chromosome mapping for contig $curctg!\n") unless $ctgdata;
       my ($glacc)=($curcmt=~m/sequence is identical to\s+(GL\d+)/s);
       my $chrname=$$ctgdata[0];
       if ($glacc) {
          $$ctgdata[4]=$glacc;
          my $chr=$$ctgdata[0];
          if ($chr=~m/^(chr\w+)/) {
             my $chrbase=$1;
            # try to match UCSC naming scheme
             if ($chrbase eq 'chrUn') {
                 $chrname=$chrbase.'_'.lc($glacc);
                 }
               else {
                 if (exists($haplotypes{$glacc})) {
                    $chrname=$haplotypes{$glacc};
                    }
                  else {
                    $chrname=$chrbase.'_'.lc($glacc).'_random';
                    }
                 }
             }
          }
        $$ctgdata[5]=$chrname;
        if ($chrlen_test) {
           my $oldchrlen=$chrsizes{$chrname};
           if ($oldchrlen) {
              die("Error: $chrname length mismatch (expected '$chrsizes{$chrname}' but found $curctglen) !\n")
                 unless ($oldchrlen==$curctglen);
              }
           else {
             my $msg="Warning: couldn't check length of $chrname.";
             if ($skip_no_len_check) {
                $skipctg=1;
                delete $ctgmap{$curctg};
                GMessage($msg." Skipped.");
                next;
                }
             GMessage($msg);
             }
           }
       }
    #big FEATURES block
    if (m/^\s{3,8}(['\-\w]+)\s{2,}(\S+)/) { # a feature starts here (gene/*RNA/CDS/exon
       my ($newfeat, $rest)=($1,$2);
       flushFeature() if $curfeat;
       $lastattr='';
       $curfeat=$newfeat;
       
       # next unless exists($fparse{$curfeat});
       next if exists($fskip{$curfeat});
       unless (exists($fparse{$curfeat})) {
         GMessage("Warning: skipping unknown feature: $curfeat!");
         next;
         }
       # check for coordinate definition
       $_=$rest;
       @fsegs=();
       $curfstrand=(s/complement\(//)?'-':'+';
       my $multi=(s/join\(//);
       $fpartial5=1 if m/^\<\d+/;
       s/\s+$//;
       my $readnext;
       do {
         my @s=(m/([\d+\.\>]+)/g);
         $readnext=m/\,$/;
         foreach my $seg (@s) {
           my ($x1,$x2)=($seg=~m/(\d+)/g);
           unless ($x2) { # microexons may be given as only one coordinate
               $x2=$x1;
               }
             else {
              ($x1,$x2)=($x2,$x1) if $x1>$x2;
              }
           push(@fsegs, [$x1,$x2]);
           }
         $_=<> if $readnext;
         } while ($readnext); #parse feature segments (exons)
       $fpartial3=1 if m/^\.\.\>\d+\)*$/;
       die("Error: couldn't parse segment coordinates for feature $curfeat on $curctg\n$_\n") unless @fsegs>0;
       next;
       } # coordinate line for a gene/mRNA/CDS/exon/*RNA feature
     #next unless exists($fparse{$curfeat}); #only parse recognizable features
     next if exists($fskip{$curfeat});
     unless (exists($fparse{$curfeat})) {
         GMessage("Warning: skipping unknown feature: $curfeat!");
         next;
         }
     
     # -- here we are within a feature block (after coordinates were parsed)
     #-- attribute parsing
     chomp;
     if (m/^\s{9,}\/([^=]+)(.*)/) { # attribute start
       my ($attr, $value)=($1,$2);
       if ($value) {
         $value=~s/^=//;
         $value=~tr/;/,/d; #"
         $value=~s/^"//;
         $value=~s/"$//;
         }
        else { $value=1 };
       $attr=~s/^gene_syn\w+/gene_syn/;
       if ($attr eq 'db_xref' || $attr eq 'gene_syn') {
         $value=~tr/ \n\t\r//d;
         my @vmulti=split(/[\;\,]/, $value);
         if ($attr eq 'gene_syn' && @vmulti>1) {
           push(@{$fattrs{$attr}}, @vmulti);
           }
         else {
           push(@{$fattrs{$attr}}, $value);
           }
         }
        else {
         $fattrs{$attr}=$value;
         }
       $lastattr=$attr;
       } #attribute start line
     else { #value continuation line for last attribute:
      die("Error: invalid line encountered when attribute continuation was expected!\n".
           "$last_line\ncurctg=$curctg, curfeat=$curfeat, lastattr=$lastattr\n")
         unless $lastattr && exists $fattrs{$lastattr};
      s/^\s+/ /;
      s/"$//;
      my $v=$_;
      if ($lastattr eq 'gene_syn') {
        $v=~tr/ \n\t\r//d;
        my @vmulti=split(/[\;\,]/,$v);
        push(@{$fattrs{$lastattr}}, @vmulti);
        }
       else {
        $fattrs{$lastattr}.=$v;
        }
      } #attribute continuation
      #if ($lastattr eq 'note') {
      #  my $v=$fattrs{$lastattr};
      #  if ($v=~s/Derived by automated(.+?)prediction method: ([^\.]+)\..*$/predicted by: $2/) 
      #    { $fattrs{$lastattr}=$v }
      #  }
 
    }# we are in FEATURES block
} #for each input line

flushGSeq();

if ($f_out) {
 close(F_OUT);
 }


#---------------------------------------------------------
sub flushGSeq {
 return unless $curctg;
 flushFeature(); # close last feature read
 return unless keys(%genes);
 
 my $gcount=keys(%genes);
 while (my ($g,$gd)=each(%genes)) {
    writeGene($g,$gd);
 }
 %genes=();
 $curctg='';
 $curctglen=0;
 $curcmt='';
 $cursection='';
 #-- debug only:
 #exit;
}


sub makeUniqID {
 my ($idref, $href)=@_;
 my $id=$$idref;
 $id=~s/~\d+$//;
 my $c = ++$$href{$id};
 if ($c>1) {
   $c--;
   $$idref=$id.'.m'.$c;
   }
 return $$idref;  
}

sub createGeneEntry {
 my ($gene, $gstart, $gend)=@_;
 push(@{$gids{$gene}},[$curctg, $gstart, $gend]);
 my $gno=scalar(@{$gids{$gene}})-1;
 $gene.='~'.$gno if $gno>0;
 if (exists($genes{$gene})) { #should never happen, we just made sure it's unique!
       my $prev_gstart=$genes{$gene}->[2];
       my $prev_gend=$genes{$gene}->[3];
      # if (abs($prev_gstart-$gstart)>30000) {
         die("Error -- gene $gene is duplicated ($prev_gstart vs $gstart) !?\n");
      #   }
      # else { #update min..max coordinates
      #   $genes{$gene}->[2]=$gstart if $gstart<$prev_gstart;
      #   $genes{$gene}->[3]=$gend if $gend>$prev_gend;
      #   ${$genes{$gene}->[5]}{'fragmented'}=1;
      #   }
      }
    else 
    { # create a new gene entry here
     my $pseudo= delete($fattrs{'pseudo'}) ? 1 : 0;
     #                 0        1           2       3     4    5   6      7 
     $genes{$gene}=[$curctg, $curfstrand, $gstart, $gend, [],  {}, {}, $pseudo, []];
     }
 delete $fattrs{'gene'};
 if ($moreattrs) {
   if (keepAttr('db_xref', \%fattrs)) {
      $genes{$gene}->[4]=[@{$fattrs{'db_xref'}}];
      delete $fattrs{'db_xref'};
      }
   #another special case: gene_syn
   if (keepAttr('gene_syn', \%fattrs)) {
      $genes{$gene}->[8]=[@{$fattrs{'gene_syn'}}];
      delete $fattrs{'gene_syn'};
      }
   my $garef=$genes{$gene}->[5];
   if ($allattrs) {
      my @k=keys(%fattrs);
      @{$garef}{@k}=(values(%fattrs));
      } else{
      foreach my $k (keys(%fattrs)) {
        $garef->{$k}=$fattrs{$k} if exists($xattrs{$k});
        }
      }
   } #more attributes are kept
 return $gene;
}


sub flushFeature { #called after a feature attributes block has been parsed into %fattrs
 return unless $curfeat && keys(%fattrs)>0;
 my $v=$fattrs{'note'};
 if ($v && 
      $v=~s/Derived by automated(.+?)prediction method: ([^\.]+)\..*$/predicted by: $2/) {
   $fattrs{'note'}=$v;
   }
 my $gene=$fattrs{'gene'};
 my @exon= sort { $main::a->[0] <=> $main::b->[0] } @fsegs;
 my ($xstart, $xend)=($exon[0]->[0], $exon[-1]->[1]);
 if (!$gene && $curfeat=~m/[^m]RNA/) {
   my $t=$fattrs{'product'};
   if (length($t)<14) {
     $t=~tr/- ,/___/s;
     $gene=$t;
     }
    else {
     $gene=$curfeat;
     }
   $gene='NOT_GENE|'.$gene;  
   # clearly this is an orphan RNA entry (with no gene)
   # so we create a dummy gene entry AND a transcript entry here
   createGeneEntry($gene, $xstart, $xend);   
   }
   
 die("Error: no gene attribute found for current feature $curfeat\n")
    unless $gene;
 if ($curfeat eq 'gene') { # gene attributes processing
   #make sure it's unique
   die("Error: gene feature ($gene) has multiple segments?\n") if @fsegs>1;
   die("Error: gene feature ($gene) has no coordinates?\n") if @fsegs==0;
   createGeneEntry($gene, $xstart, $xend);
   @fsegs=();
   %fattrs=(); 
   ($fpartial3, $fpartial5)=();
   return; #gene data added
   } # stored gene attributes
# else 
# - storing RNA attributes
# now make sure we place this RNA within the right gene - in case of gene duplications
# locate the matching $gene
my $gxdata=$gids{$gene};
unless ($gxdata) {
  # strange, this looks like another orphan xRNA entry
  createGeneEntry($gene, $xstart, $xend);
  $gxdata=$gids{$gene};
  }

die("Error: gene $gene not found in \%gids !\n") unless $gxdata; #should never happen
if (@$gxdata>1) {
  my $i=0;
  my $found;
  foreach my $gd (@$gxdata) {
    if ($gd->[0] eq $curctg && $gd->[1]<=$xstart && $gd->[2]>=$xend) {
      $found=1;
      last;
      }
    $i++;
    }
  die("Error: cannot find gene enclosing RNA data ".
      $fattrs{'transcript_id'}." ".$fattrs{'protein_id'}." (on $curctg, $xstart..$xend)\n")
       unless $found;
  $gene.='~'.$i if $i>0; 
  }

my $geneid=$gene;
my $gdata=$genes{$gene};
die("Error: \$genes data for gene '$gene' cannot be found!\n") unless $gdata;
my $rnas=$gdata->[6];
my $gpseudo=$gdata->[7];
$geneid=~tr/~/_/;
my ($tid, $pid);
my $isMito=($ctgmap{$curctg}->[0] eq 'chrM'); #different format for mitochondrial entries
$pid=$fattrs{'protein_id'};
my $isCDS = ($curfeat eq 'CDS');
if ($isCDS && $pid) { # just finished a CDS parsing
  #unless ($pid) {
  #   print STDERR "Warning: no protein_id found for CDS feature of gene $geneid (on $curctg).\n";
  #   @fsegs=();
  #   %fattrs=();
  #   ($fpartial3, $fpartial5)=();
  #   return;
  #   }
  my @t=split(/\|/,$pid);
  $pid=$t[1] if (@t>1);
  $pid=~s/\.\d+$//; #remove version
  $tid=$pmap{$pid};
  die("Error: cannot get transcript for protein id $pid (gene $geneid)!\n") unless $tid;
  unless (exists($rnas->{$tid})) { #this should never happen - CDS without mRNA definition?
                                  # --> yeah, it does happen for Mitochondrial chromosomes!
  if ($isMito) {
       $tid=$pid;
       #                0        1       2      3   4   5  6  7
       #$rnas->{$tid}=['mRNA', [@exon], [],     {}, [], 0, 0, []];
       $rnas->{$tid}=['mRNA', [@exon], [@exon], {}, [], 0, 0, []];
       }
     else {  
      $tid=~s/\.\d+$//;
      die("Error: no mRNA entry for CDS/protein $pid (reported as coming from transcript $tid)\n")
       unless exists($rnas->{$tid});
       }
    }
   # #just add the CDS
  push(@{$rnas->{$tid}->[2]}, @exon);
  } #CDS
 else { # exons, other RNA features or CDS without protein_id
  if ($isCDS) {
       my $msg="Warning: no protein_id for CDS feature of gene $geneid (on $curctg)";
       if ($gpseudo) {
        # discard CDS segments for pseudogenes
        @fsegs=();
        %fattrs=();
        ($fpartial3, $fpartial5)=();
        $msg=~s/ gene/ pseudogene/;
        GMessage("$msg. Discarded.");
        return;
        }
       GMessage("$msg.");
       }
 $tid=$fattrs{'transcript_id'};
 if ($tid) {
    my @t=split(/\|/,$tid);
    $tid=$t[1] if (@t>1);
    #$tid=~s/^[a-z]+\|//i;
    #$tid=~s/\|\w*$//; # some entries have a suffix - e.g. ref|NM_006824.1|P40
    $tid=~s/\.\d+$//; #remove version!
    }
   else {
    $tid=$geneid.'.t'; #some special cases with features parented directly by [pseudo] genes
    }
  unless (exists($rnas->{$tid})) {
    my $frna=$curfeat;
    if ($frna!~m/RNA/) {
      $frna= $gpseudo? 'pseudoRNA' : 'misc_RNA';
      }
    if ($isCDS) {
       #               0         1       2    3   4   5  6    7
       $rnas->{$tid}=[$frna, [], [@exon], {},  [], 0, 0,  [] ];
       }
     else {
       #               0           1     2   3   4   5  6    7
       $rnas->{$tid}=[$frna, [@exon], [], {}, [], 0, 0,  [] ];
       }
    }
   else { #just add the exons
    my $segidx = $isCDS ? 2 : 1;
    push(@{$rnas->{$tid}->[$segidx]}, @exon);
    }
  } #exon segments
my $tdata=$rnas->{$tid};
$tdata->[5]=1 if $fpartial5;
$tdata->[6]=1 if $fpartial3;
# add attributes
delete @fattrs{'gene', 'transcript_id'};
delete $fattrs{'protein_id'} unless exists($xattrs{'protein_id'});
my $product = delete $fattrs{'product'};
my $ahref=$tdata->[3];
$ahref->{'product'}=$product if $product;
if (my $ncRNAclass=delete $fattrs{'ncRNA_class'}) {
  $ahref->{'ncRNA_class'}=$ncRNAclass;
  }
if ($moreattrs) {
 my %xrefs;
 my @cx;
 @cx=@{$fattrs{'db_xref'}}  if keepAttr('db_xref', \%fattrs); #if exists $fattrs{'db_xref'};
 push(@cx, @{$tdata->[4]});
 @xrefs{@cx}=();
 @cx=keys(%xrefs);
 $tdata->[4]=[@cx];
 delete $fattrs{'db_xref'};
 #if (exists $fattrs{'gene_syn'}) {
 if ($isCDS && (my $v=$fattrs{'note'})) {
   if ($v=~s/isoform .+? encoded by transcript [\w \-\,\.]+//) {
     $v=~s/^\;?\s+//;
     $v=~s/\;?\s+$//;
     $v=~s/\;\s*/, /;
     if (length($v)>2) { $fattrs{'note'}=$v }
                  else { delete $fattrs{'note'} }
     }
   }
 if (keepAttr('gene_syn', \%fattrs)) {
   %xrefs=();
   @cx=@{$fattrs{'gene_syn'}};
   push(@cx, @{$tdata->[7]});
   @xrefs{@cx}=();
   @cx=keys(%xrefs);
   $tdata->[7]=[@cx];
   delete $fattrs{'gene_syn'};
   }
 #add new attributes
 if ($allattrs) {
   foreach my $k (keys %fattrs) {
     $ahref->{$k}=$fattrs{$k};
     }
   }
   else { # $addattrs
   foreach my $k (keys %fattrs) {
     $ahref->{$k}=$fattrs{$k} if exists($xattrs{$k});
     }
   }
 }  
#--> clean up
@fsegs=();
%fattrs=();
($fpartial3, $fpartial5)=();
$lastattr='';
$curfeat='';
}

sub writeGene {
 my ($gene, $gdata)=@_;
 my $notgene=($gene=~s/^NOT_GENE\|//);
 my $geneid=$gene;
 my ($gene_name)=($geneid=~m/^([^~]+)/);
 $geneid=~tr/~/_/;
 my ($cbase,       $cstrand,      $gstart,     $gend,      $gpseudo)=
  ($gdata->[0], $gdata->[1], $gdata->[2], $gdata->[3], $gdata->[7]);
 my $gfeature=($gpseudo)?'pseudogene' : 'gene';
 #makeUniqID(\$gene, \%gids); # $gene will be updated if necessary
 #my $geneid=$gene;
 #$geneid=~tr/~/_/;
 my $attrs="ID=$geneid;Name=$gene_name";
 my @gxrefs=@{$gdata->[4]};
 $attrs.=';xrefs='.join(',', @gxrefs) if (@gxrefs>0); 
 my @gsyn=@{$gdata->[8]};
 $attrs.=';gene_syn='.join(',', @gsyn) if (@gsyn>0); 
 my $gattrs=$gdata->[5];
 my $rnas=$gdata->[6];
 my $gcoords;
 my ($gbase, $strand, $gchrlen);
 ($gbase, $strand, $gcoords, $gchrlen)=contig2chr($cbase, $cstrand, [[$gstart, $gend]]);
 ($gstart, $gend)=($$gcoords[0]->[0], $$gcoords[0]->[1]);
 my $showpseudo = $gpseudo && (keys(%$rnas)==0);
 unless (!$showpseudo && ($mrnaonly || $t_only || $notgene)) {
   foreach my $k (keys(%$gattrs)) {
       $attrs.=";$k=".$gattrs->{$k};
       }
   print join("\t", $gbase, $track, $gfeature, $gstart, $gend, '.', $strand, '.', $attrs)."\n";
   }
 
 foreach my $tid (keys(%$rnas)) {
    #die("Error: no transcript_id found at write_RNA for gene $gene ($geneid)\n$last_line\n") unless $tid;
    my $tdata=$rnas->{$tid};
    my $rnafeat=$tdata->[0];
    if ($rnafeat eq 'RNA_exon' || $rnafeat eq 'misc_RNA') {
       $rnafeat = $gpseudo ? 'pseudoRNA' : 'misc_RNA';
       }
    if ($mrnaonly) {
      $rnafeat=~s/misc_RNA/mRNA/;
      }
      
    next if $mrnaonly && $rnafeat ne 'mRNA';
    # my $texons=$tdata->[1];
    my $texons;
    ($gbase, $strand, $texons)=contig2chr($cbase, $cstrand, $tdata->[1]);
    my $tcds=$tdata->[2];
    if ($tcds && @$tcds>0) {
      ($gbase, $strand, $tcds)=contig2chr($cbase, $cstrand, $tdata->[2]);
      }
    else {
      $tcds='';
      }
    my $tahref=$tdata->[3]; #attributes
    my $txar=$tdata->[4]; # db_xrefs
    my $txgsyn=$tdata->[7]; # gene_syns
    my $ncRNA;
    if (my $rnaclass=$tahref->{'ncRNA_class'}) {
      $rnafeat=$rnaclass;
      delete($$tahref{'ncRNA_class'});
      $ncRNA=1;
      }
    my ($part3, $part5)=($tdata->[5], $tdata->[6]);
    makeUniqID(\$tid, \%tids);
    my @exon=sort { $main::a->[0] <=> $main::b->[0] } @$texons;
    my ($tstart, $tend)=($exon[0]->[0], $exon[-1]->[1]);
    $tid=~s/^NOT_GENE\|//;
    my $attrs="ID=$tid";
    my $gffname=$tid;
    if ($notgene) {
      $attrs.=";Name=$tid";
      }
     else {
      $attrs.=";Parent=$geneid" unless ($mrnaonly || $t_only);
      $gffname=$gene_name;
      $attrs.=";Name=$gene_name;gene_name=$gene_name";
      }
  
    my $product=$tahref->{'product'};
    $ncRNA=1 if $rnafeat eq 'tRNA' || $gffname=~m/NCRNA\d+/;
    $attrs.=';ncRNA=1' if $ncRNA;
    $attrs.=';partial=1' if ($part5 || $part3);
    if ($product) {
      $product=~tr/"/'/; #"
      $product=~tr/;/./;
      # $attrs.=';product="'.$product.'"';
      $attrs.=";product=$product";
      }
    if ($moreattrs) {
      $attrs.=';xrefs='.join(',', @$txar) if ($txar && @$txar>0);
      $attrs.=';gene_syn='.join(',', @$txgsyn) if ($txgsyn && @$txgsyn>0);
      foreach my $k (keys(%$tahref)) {
        next if exists($kattrs{$k});
        $attrs.=";$k=".$tahref->{$k};
        }
      }
    print join("\t", $gbase, $track, $rnafeat, $tstart, $tend, '.', $strand, '.', $attrs)."\n";
    unless($ncRNA && @exon==1) {
      foreach my $ex (@exon) {
        print join("\t", $gbase, $track, 'exon', $$ex[0], $$ex[1], '.', $strand, '.', "Parent=$tid")."\n";
        }
     }
    @exon=();
    #write_CDS();
    if ($tcds) { #have CDS, write it
      my @cds= sort { $main::a->[0] <=> $main::b->[0] } @$tcds;
      my $codonstart=$tahref->{'codon_start'};
      if ($codonstart>1) {
        $codonstart--;
        if ($strand eq '-') {
          $cds[-1]->[1]-=$codonstart;
          }
         else {
          $cds[0]->[0]+=$codonstart;
          }
        }
       else { $codonstart='0'; }

      foreach my $cd (@cds) { #TODO: calculate phase ?
        print join("\t", $gbase, $track, 'CDS', $$cd[0], $$cd[1], '.', $strand, '.', "Parent=$tid")."\n";        
        }
      }
  } #for each transcript
}

sub contig2chr {
my ($ctg, $ctgstrand, $ctgsegs)=@_;
my $cd=$ctgmap{$ctg};
die("Error: cannot find chromosome mapping for contig $ctg!\n") unless $ctg;
my ($chr, $chrstrand, $chrstart, $chrend, $glacc, $chrname)=@$cd;
unless ($ctgmapfile) {
  return ($chrname, $ctgstrand, [@$ctgsegs], $ctgmap{$ctg}->[3]);
  }
my $clen=$chrend-$chrstart+1;
my @chrsegs;
foreach my $c (@$ctgsegs) {
 my ($start, $end)=($chrstrand eq '-') ? ($clen+$chrstart-$$c[1],$clen+$chrstart-$$c[0]) :
               ($$c[0]+$chrstart-1, $$c[1]+$chrstart-1); 
 push(@chrsegs, [$start ,$end]);
 }
my $nstrand=($ctgstrand eq $chrstrand) ? '+' : '-' ;
return ($chr, $nstrand, [@chrsegs], $clen, $chrname);
}


sub keepAttr {
 my ($attr, $href)=@_;
 
 if ($allattrs) {
   if ($href) { return exists($href->{$attr}); }
         else { return 1; }
   }
  elsif ($addattrs)  {
   if (exists($xattrs{$attr})) {
       if ($href) { return exists($href->{$attr}); }
             else { return 1; }
       }
     else { return 0; }
   }
  else { return 0; }
}

sub GMessage {
 print STDERR join("\n",@_)."\n";
}

