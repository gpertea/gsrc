#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;
use dbSession;

my $usage = q/Usage:
 prep4sequin.pl -f <fasta_file> [-a <acc2version_table>] [-A] [-t <track>] 
    [-o <out_prefix>] [-c 'organism_code'] [-n '<organism_name>'] <gff_file>


 Default <track> used for features is 'jigsaw'.

 Creates a .fsa and .tbl files to be used by tbl2asn
 
/;
umask 0002;
getopts('n:c:Aa:f:t:o:') || die($usage."\n");
my $fastafile=$Getopt::Std::opt_f || die("$usage A fasta file is required!\n");
my $defpre=$fastafile;
$defpre=~s/\.\w+$//;
my $orgcode=$Getopt::Std::opt_c;
unless ($orgcode) { 
 ($orgcode)=($defpre=~m/([a-z,_]+)\d+$/i);
 $orgcode=~tr/-_//d;
 $orgcode=uc($orgcode);
 $orgcode='BB' if $orgcode eq 'BS';
 }
my $accfile=$Getopt::Std::opt_a;
my $accparse=$Getopt::Std::opt_A;

my $outprefix=$Getopt::Std::opt_o || $defpre;
my $orgname=$Getopt::Std::opt_n;
my $track=$Getopt::Std::opt_t || 'jigsaw';
# die("$usage A feature track must be specified!\n");
my $gff=shift(@ARGV) || die("$usage An input gff3 file must be given!\n");
die("$gff: nothing to do.") if (-e $gff && (-s $gff<10));

my $ds = dbSession->new('geanno@NEOSYBASE');
my $sth=$ds->prep('select xref_data from uniprot_xref where xrefdb in ("GenBank", "EMBL") and up_id=?');
my %acc;
if ($accfile) {
 open(ACCFILE, $accfile) || die ("Error opening $accfile!\n");
 chomp;
 while (<ACCFILE>) {
  my @a=split;
  $acc{$a[0]}=$a[1];
  }
 close(ACCFILE);
}
open(FSA, '>'.$outprefix.'.fsa') 
  || die ("Error creating file $outprefix.fsa !\n");
open(FA, $fastafile)||die("Error opening $fastafile !\n");
my ($seqid, $seqlen);
while (<FA>) {
 if (m/^>(\S+)\s*(.*)/) {
   last if $seqid;#only one sequence !
   ($seqid, my $rest)=($1, $2);
   if ($accfile) {
     my $accver=$acc{$seqid};
     $seqid=$accver if $accver;
     }
   my ($enc)=($rest=~m/EN([mr]\d+)/);
   $orgcode.='_'.$enc if ($enc);
   if ($accparse) {
     my ($short, $acc)=split(/\|/,$seqid);
     if ($acc=~m/(\w+)v(\d+)$/) {
       $seqid=$1.'.'.$2;
       }
     }
   print FSA '>'.$seqid;
   print FSA " [organism=$orgname]" if $orgname;
   print FSA " [gcode=1] [primary=$seqid] $rest\n";
   next;
   } 
 print FSA $_;
 tr/\n\r\t //d;
 $seqlen+=length($_);
}
close(FSA);

open(FTBL, '>'.$outprefix.'.tbl') 
  || die("Error creating $outprefix.tbl file!\n");
open(GFF, $gff) || die("Error opening the input gff file: $gff !\n");
print FTBL ">Feature\t$seqid\n";
#print FTBL join("\t", 1, $seqlen, 'misc_feature')."\n";
my $credits='Predicted annotation generated by JIGSAW and other methods, and provided'.
' by Geo Pertea, Mihaela Pertea, Jonathan Allen and Steven Salzberg, Center for'.
' Bioinformatics and Computational Biology, University of Maryland';
#print FTBL join("\t", '', '', '','note', $credits)."\n";

my $cparent; #current parent - ASSUMES lines are ordered properly (children following parent)
my $cf; #current parent feature name (e.g. 'mRNA' or 'gene')
my $cxf; # child feature type: 'exon' or 'CDS'

my $cstrand;
my ($cstart, $cend);
#TODO: also use @ex for mRNA fkey (when UTRs are available)
my @cx;      #children intervals
my ($cdescr, $gene_name, $mcov);
my %gn; #unique gene names
while (<GFF>) {
 next if m/^\s*#/;
 chomp;
 my ($chr, $ftrack, $f, $fstart, $fend, $fscore, $strand, $frame, $attrs)=split(/\t/);
 next unless lc($track) eq lc($ftrack);
 ($fstart, $fend)=($fend, $fstart) if $fend<$fstart;
 $f='exon' if ($f =~m/\-exon$/);
 my ($fid)=($attrs=~m/ID\s*=\s*"?([^;"]+)/);
 my ($fp)=($attrs=~m/Parent\s*=\s*"?([^;"]+)/);
 unless ($fid || $fp) {
   print STDERR "Warning: feature line $_ has neither ID nor Parent! skipping..\n";
   next;
   }
 if ($fp) { #child feature
   if ($cparent && $cparent ne $fp) {
     die("Error: invalid order of input lines (unknown parent $fp while still processing $cparent\n$_\n");
     }
   if ($cxf && $cxf ne $f) {
     # TODO: will permit both 'exon' and 'CDS' subfeatures later
     die("Error: multiple subfeatures are not accepted yet (parent $fp, $cf has subfeature $cxf already)\n");
     }
   $cxf=$f;
   push(@cx, [$fstart, $fend]);
   }
 else { #parent
   if ($cparent ne $fid) {
      writeFeature() if @cx>0;
      $cparent=$fid;
      $cstrand=$strand;
      $cf=$f;
      $cxf='';
      ($cstart, $cend)=($fstart, $fend);
      @cx=();
      $cdescr='';
      }
     else {
      die("Error: duplicate parent entry for $fid\n");
      }
     # try to parse the description/annotation, if any
     ($mcov)=($attrs=~m/MCov=([\d\.]+)/);
     ($gene_name)=($attrs=~m/GeneId=([^\;]+)/);
     if ($attrs=~m/TopHit\s*=\s*"?([^;"]+)/i) {
          $cdescr=$1;
          }
         elsif ($attrs=~m/Descr\s*=\s*"?([^;"]+)/) {
          $cdescr=$1;
          }
         elsif ($attrs=~m/Info\s*=\s*"?([^;"]+)/) {
          $cdescr=$1;
          }
         elsif ($attrs=~m/BestHit\s*=\s*"?([^;"]+)/i) {
          $cdescr=$1;
          }
         elsif ($attrs=~m/Name\s*=\s*"?([^;"]+)/) {
          $cdescr=$1;
          }
     if ($cdescr) {
       $cdescr=~s/^\s+//;$cdescr=~s/\s+$//;
       $cdescr='' if index(lc($fid),lc($cdescr))>=0;
       }
   #     if ($cdescr=~m/hypothetical/ && $cdescr=~m/\bLOC(\d+)/) {
   #        
   #        }
   } #parent
  
}

writeFeature() if @cx>0;

close(GFF);
close(FTBL);

$sth->finish();
$ds->logout();

sub writeFeature {
  #TODO: use @ex to decide mRNA and CDS 
  $cparent=~tr/.,:/___/s;
  #print STDERR join("\t",$cparent, $cstrand, $cstart, $cend)."\n";
  #my ($pl, $pr)=('<','>');
  #$cstart=$pl.$cstart;
  #$cend=$pr.$cend;
  @cx= sort { $main::a->[0]<=>$main::b->[0] } @cx;
  #$cx[0]->[0]=$pl.$cx[0]->[0];
  #$cx[-1]->[1]=$pr.$cx[-1]->[1];
  if ($cstrand eq '-') { #reverse complement features
     ($cstart, $cend)=($cend, $cstart);
     #@cx= sort { $main::b->[0]<=>$main::a->[0] } @cx;
     @cx = map { $_=[$_->[1], $_->[0]] } @cx;
     @cx= reverse(@cx);
     }
 my $s="\t";
 my $cid='';
 my $infrsim;
 if ($cdescr) {
   my $hdescr;
   ($cid, $hdescr)=split(/ /,$cdescr,2);
   $cdescr=$hdescr if $hdescr;
   $cdescr=~s/^gid:\S+\s+//;
   $cdescr=~s/^CDS:\d+\-\d+\s+//;
   $cdescr=~s/UniRef100_\S+\s+//;
   if ($cid=~m/^UPr?\|(\w+)/) {
     my $upid=$1;
     $cid=up2gb($upid);
     if ($cid) {
      $infrsim='AA sequence:INSD:'.$cid;
      $cid.=' ';
      }    
     $cdescr=~s/\([^\)]+\)\s*//g;
     }
    elsif ($cid=~m/[rpvinfmod]+\|([\w\.]+)/) { #refseq entry
     $cid=$1;
     my $refacc=$cid;
     $refacc.='.1' unless $refacc=~m/\.\d+$/;
     $infrsim='RNA sequence:RefSeq:'.$refacc;
     $cid.=' ';
     }
    else {
     $cid='';
     } 
   }
 #print FTBL join($s, '', '', '','locus_tag', $cparent)."\n";
 my $note;
 my $product='hypothetical mRNA';
 my $cdsproduct='hypothetical protein';
 my $gproduct='';
 if ($mcov>50 && $cdescr) {
    my $d=$cdescr; $d=~s/UniRef100_\S+\s+//;
    $d=~s/\s+\{[ \.\,\-\"\'\w]+\}$//;
    # make words start with lowercase:
    $d =~ s/([A-Z][a-z]{5,})\b/\L$1/g;
    $product=$d.' (predicted)';
    $cdescr=~s/\{([ \.\,\-\"\'\w]+)\}$/\($1\)/;
    #$note='supporting evidence includes similarity to '.$cid.$cdescr;
    
    $note='similar to '.$cid.$cdescr;
    $cdsproduct=$product;
    }
   else {
    $gene_name='';
    }
 my $lpid; #local protein id
 if ($gene_name) {
   my $gnum=++$gn{$gene_name};
   if ($gnum>1) {
    $gene_name.='_'.$gnum;
    }
   $gene_name=~s/[ \-_]?predicted$//;
   $lpid=$orgcode.'_'.$gene_name;
   #$gene_name.='_predicted';
   }
  else {
   my $jname=$cparent;
   $jname=~tr/_//d;
   $jname=~s/tjsm/jsm/;
   $lpid=$orgcode.'_jsm'.uc(sprintf('%x',$cstart)).($cstrand eq '+' ? 'f':'r');
   #$gene_name='jigsaw_prediction_'.$jname;   
   $gene_name=$lpid;
   }
 print FTBL join($s, '<'.$cstart, '>'.$cend, 'gene')."\n";
 print FTBL join($s, '', '', '','gene', $gene_name)."\n";
 my $jspredline=join($s,'','','','inference','ab initio prediction:JIGSAW:3.2')."\n";
 print FTBL $jspredline;
 if ($infrsim) {
  print FTBL join($s,'','','','inference','similar to '.$infrsim)."\n";
  }

 #print FTBL join($s, '', '', '','locus_tag', $lpid)."\n";
 my $x=shift(@cx);
 
 print FTBL join($s, '<'.$$x[0], $$x[1], 'mRNA')."\n"; 
 my $endrna= $cx[-1]->[1];
 $cx[-1]->[1]='>'.$endrna;
 foreach my $xc (@cx) {
  print FTBL join($s, $$xc[0], $$xc[1])."\n";
  }
 print FTBL $jspredline;
 #print FTBL join($s, '', '', '','transcript_id', $lpid.'-t')."\n";
 print FTBL join($s, '', '', '','product', $product)."\n";
 #print FTBL join($s, '', '', '','note', $note)."\n";
 print FTBL join($s, $$x[0], $$x[1], 'CDS')."\n";
 $cx[-1]->[1]=$endrna;
 foreach my $xc (@cx) {
  print FTBL join($s, $$xc[0], $$xc[1])."\n";
  }
 print FTBL $jspredline; 
 if ($infrsim) {
  print FTBL join($s,'','','','inference','similar to '.$infrsim)."\n";
  }
 # print FTBL join($s, '', '', '','codon_start', 1)."\n";
 print FTBL join($s, '', '', '','protein_id', 'gnl|NISC-CON|'.$lpid)."\n";
 print FTBL join($s, '', '', '','product', $cdsproduct)."\n";
 print FTBL join($s, '', '', '','note', $note)."\n" if $note;

}
 
sub up2gb {
 my $uid=shift(@_);
 $ds->exec($sth, $uid);
 my $first;
 while (my $r=$ds->fetch($sth)) {
   $first=$$r[0] unless $first;
   }
 ($first)=($first=~m/^([\w\.]+)/) if $first;
 return $first;
}