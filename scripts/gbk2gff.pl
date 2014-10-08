#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q{ Usage:
 gb2gff.pl [-A] [-f <fasta_output>] [-o <out_gff_file>] genbank_file(s)
 
 Options:
  -A  convert ALL features (default is: only convert genes/mRNAs
  -f  also write a fasta version of the base sequence(s) in
      the file <fasta_output>
};
umask 0002;
getopts('Af:o:') || die($usage."\n");
my $allfeats=$Getopt::Std::opt_A;
my $outfile=$Getopt::Std::opt_o;
my $fastaout=$Getopt::Std::opt_f;
my %gnames;
my $fakeGI=121000;
my $vegasrc=0;
if ($fastaout) {
 open(FASTAOUT, '>'.$fastaout)|| die ("Error creating file $fastaout\n");
 }
if ($outfile) {
 open(GFFOUT, '>'.$outfile) || die ("Error creating file $outfile\n");
 select(GFFOUT);
 }
my ($bseq, $fseq, $bdef, $bGI, $bACC, $blen, $freading, $sreading);
my ($gname, $gdescr, $mrna);
while (<>) {
 chomp;
NEXTLINE:
 if (m/^LOCUS/) {
   ($bseq, $blen)=(m/^LOCUS\s+(\S+)\s+(\d+)/);
   die("Error parsing LOCUS tag info: $_\n") unless $blen>1;
   }
 elsif (m/^DEFINITION\s+(.+)/) {
   $bdef=$1;
   }
 elsif (m/^ACCESSION\s+(\S+)/) { 
   $bACC=$1;
   }
 elsif (m/^VERSION/) {
   if (m/chromosome:VEGA/) {
     $bGI = $fakeGI++;
     $vegasrc=1;
     }
    else {
     ($bGI)=(m/\bGI:(\d+)/);
     die("Error parsing GI# from VERSION tag ($_)") unless $bGI>1;
     }
   }
 elsif (m/^FEATURES\s+/) {
   $freading=1;
   ($gname, $gdescr, $mrna)=('','','');
   next;
   }
 elsif (m/^BASE COUNT/) {
   $freading=0;
   next;
   }
 elsif (m/^ORIGIN\s+/) {
   $freading=0;
   $sreading=1;
   next;
   }
 elsif (m/^\/\/$/) {
   $freading=0;
   $sreading=0;
   if ($fastaout) {
     print FASTAOUT &fastafmt(\$fseq, $bseq.' gi|'.$bGI.' '.$bdef);
     }
   $fseq='';
   $bdef='';
   $bseq='';
   }
   
 #-- remaining lines parsing:
 if ($freading) {

   my ($featname, $frange)=(m/^\s+([\w\-]+)\s+(.+)/);
   die("Error parsing feature at: $_\n") unless $featname && $frange;
   ($gname, $gdescr)='' if $featname eq 'gene';

   if ($frange=~m/\(/ && $frange!~m/\)$/) { #wrapped composite intervals line
     do {
        $_=<>;
        last unless $_;
        chomp;
        $frange.=$_;
        } until (m/\)$/);
     } #wrapped composite intervals;
   my %attrs; 
   my $endattr;
   #print STDERR "..reading attributes for $featname ($frange)\n";
   while (<>) { #read all attributes for this feature
     unless (m/^\s+\/([^=]+)=(.+)/) {
        $endattr=1;
        last;
        }
     chomp;
     my ($fattr, $fattr_dta)=&parseAttr($1,$2); #this will read the next line(s) if necessary
     #$fattr_dta =~ tr/"//d; #"
     $attrs{$fattr}=$fattr_dta if ($fattr ne 'translation');
     } #attribute parsing loop
   if ($vegasrc) {
    #print STDERR "############# FEATNAME = $featname \n";
    if ($featname eq 'gene') {
      $gname=$attrs{'locus_tag'};
      $gdescr=$attrs{'note'};
      $mrna='';
      }
   if ($featname eq 'mRNA') {
     my $note=$attrs{'note'};
     ($mrna)=($note=~m/transcript_id=([\-\.\w]+)/);
     #print STDERR "##############>>> MRNA === > $mrna ($note)\n";
     }
    }
   if (!$allfeats && $featname ne 'mRNA' && $featname ne 'CDS' && $featname ne 'gene') {
     goto NEXTLINE if $endattr;
     next;
     }
   &writeFeature($featname, $frange, \%attrs, $gname, $gdescr, $mrna);
   goto NEXTLINE if $endattr;
   next;
   }
   
 if ($fastaout && $sreading && m/^\s+\d+\s+(.+)/) {
   my $s=$1;
   $s=~tr/ \t\n//d;
   $fseq.=$s;
   }
   
}#while input lines

close(FASTAOUT) if $fastaout;

sub parseAttr {
 my ($a, $adta)=@_;
 #print STDERR "[[ parseAttr($a, $adta) ]]\n";

 $adta=~s/^\s+//;$adta=~s/\s+$//; #trim
 #die("Invalid attribute string for $fn\:$a ($adta)\n") unless $adta=~m/^"/;
 if ($adta=~m/^"/) {
   while ($adta!~m/\"$/) {
     $_=<>;
     chomp;
     tr/ \t/  /s; #remove extra spaces
     $adta.=$_;
     }
  }
 $adta=~s/^\"//;$adta=~s/\"$//; #remove quotes
 return ($a, $adta);
}

sub writeFeature {
 my ($fname, $frange, $attrs, $gname, $gdescr, $mrna)=@_;
 #the rest of @_ are attributes
 my $comment='';
 my $gene=$$attrs{'gene'};
 my $id=$gene;
 my $product=$$attrs{'product'};
 if ($fname eq 'mRNA') {
    $id=$$attrs{'transcript_id'} || $gene;
    }
  elsif ($fname eq 'CDS') {
    $id=$$attrs{'protein_id'} || $gene;
    }
 my $strand=($frange=~s/complement\(//) ? '-' : '+';
 $frange=~s/join\(//;
 $frange=~s/\)+$//;
 my @ranges=split(/\,/,$frange);
 my @ex;
 foreach my $r (@ranges) {
  my ($start, $end)=($r=~m/(\d+)\.\.(\d+)/);
  die("Error parsing ranges for $fname ($frange)\n") unless $start>0 && $end>$start;
  push(@ex, [$start, $end, 0]); #start, end, frame;
  }
 @ex= sort { $main::a->[0] <=> $main::b->[0] } @ex;
 my ($min, $max)=($ex[0]->[0], $ex[-1]->[1]);
 #  if ($fname eq 'mRNA' || $fname eq 'CDS') { #isoform testing
 #    my $p= $$attrs{'product'}; #identify isoforms by product.. :(
 #    my $isoform;
 #    if ($p=~m/isoform (\S+)/i) {
 #      $isoform=$1;
 #      $gene.='-'.$isoform;  
 #      }
 #    }
 if ($fname eq 'gene' || $fname eq 'mRNA' || $fname eq 'CDS') {   
   my $gene_name=$gname || $gene;
   $id=$mrna if $mrna;
   $comment="ID=$id;Name=$gene_name";
   $comment.=";geneId=$gname" if $gname;
   $comment.=";descr=\"$gdescr\"" if $gdescr;
   my $numexons = @ex;
   $comment.=";exonCount=$numexons;" if $numexons>0;
   $comment.=";product=\"$product\"" if $product;
   my $gffeat=$fname;
   $gffeat='mRNA' if $fname eq 'CDS';

   print join("\t",$bseq, 'gb2gff', $gffeat, $min, $max, '.',$strand,'.',$comment)."\n" 
        unless $fname eq 'CDS';
   return if $fname eq 'gene';
   $fname='exon' if $fname eq 'mRNA';
   }
 if ($fname eq 'exon' || $fname eq 'CDS') {
   my $n=1; 
   #print "## CDS $gene info: ".join('; ',@_)."\n" if ($fname eq 'CDS');
   my $frame=0; ## actually, this should be: codon_start - 1
   if ($fname eq 'CDS') { # compute frame
     if ($strand eq '-') {  @ex=reverse(@ex);  }
     foreach my $e (@ex) {
        $$e[2]=$frame;
        $frame = ($frame + ($$e[1] - $$e[0] + 1)) % 3;
        }
     if ($strand eq '-') {  @ex=reverse(@ex);  }
     }
   foreach my $e (@ex) {
        my ($start, $stop, $fr)=@$e;
        print join("\t",$bseq, 'gb2gff', $fname, $start, $stop, '.', $strand, $fr,
                   "Parent=$id")."\n";
        #$frame = ($fname eq 'CDS') ? ($frame + ($stop - $start + 1)) % 3 : '.';
        $n++;
        }
    }
  else { # other non-mRNA features
    $comment=join(';',@_);
    foreach my $e (@ex) {
        print join("\t",$bseq, 'gb2gff', $fname, $$e[0], $$e[1], '.',$strand,'.',
                   "Parent=$id")."\n";
        }
    }
}

sub fastafmt { # (seqref, defline, linelen)
 my $seqr=shift(@_);
 my $defline=shift(@_);
 my $linelen=shift(@_) || 60;;
 my $rec;
 $rec=">$defline\n" if $defline;
 my $seqlen=length($$seqr);
 my $pos=0;
 while ($pos<$seqlen) {
   $rec.= substr($$seqr,$pos,$linelen)."\n";
   $pos+=$linelen;
   }
 return $rec;
 }


