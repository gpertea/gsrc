#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;
my $usage = q{ 
 Script for parsing GenBank protein records GenPept File Format.
 Usage:
  genpept2fa.pl [-P][-R] [-f <dbxref.lst>] genbank_format_stream
 
 Use -P option in order to exclude the PREDICTED mRNAs and gene models.
 Use -R option to only write the REVIEWED and VALIDATED entries 
 Use -f option to provide a list of identifiers and output
 only those records matching any of the provided dbxref IDs 
 };

umask 0002;
getopts('PRf:') || die($usage."\n");
my %refseqtype=(
'INFERRED'=> 'inf', # usually not yet supported by experimental evidence
'MODEL' => 'mod', #usually predicted by GNOMON
'PREDICTED' => 'p', # only CDS is predicted in most cases, the mRNA is experimentally confirmed
                    # (it becomes 'px' then)
'PROVISIONAL' => 'pv', #usually based on experimental evidence, but not "reviewed" by NCBI
'VALIDATED' => 'v',  # preliminary reviewed, but did not go through "final review"
'REVIEWED' => 'r'
);

#my $outfile=$Getopt::Std::opt_g;
#my $fastaout=$Getopt::Std::opt_F;
my %typeflt=('mod'=>1,'p'=>1,'inf'=>1,'px'=>1,'pv'=>1, 'v'=>1, 'r'=>1, 'gb'=>1);

my $nopred=$Getopt::Std::opt_P;
my $onlyrev=$Getopt::Std::opt_R;

if ($nopred || $onlyrev) {
 @typeflt{'mod', 'inf', 'p'}=(0,0,0);
 }
 
if ($onlyrev) {
 @typeflt{'pv','px'}=(0,0);
 }

my %gnames;
#if ($fastaout) {
# open(FASTAOUT, '>'.$fastaout)|| die ("Error creating file $fastaout\n");
# }
#if ($outfile) {
# open(GFFOUT, '>'.$outfile) || die ("Error creating file $outfile\n");
# select(GFFOUT);
# }
my ($bseq, $fseq, $descr, $bGI, $bACC, $blen, $org, $comment, $cds, $cdsprod,
    $genesym, $geneID, $freading, $seqreading, $has_dbxref);

my $fxref=$Getopt::Std::opt_f;
my %xr;
if ($fxref) {
 open(XREF, $fxref) || die ("Error opening file $fxref\n");
 while (<XREF>) {
  chomp;
  next unless length($_)>0;
  my ($v)=(m/(\S+)/);
  $xr{$v}=1;
  }
 close(XREF);
}

while (<>) {
 chomp;
NEXTLINE:
 if (m/^LOCUS/) {
   ($bseq, $blen)=(m/^LOCUS\s+(\S+)\s+(\d+)/);
   die("Error parsing LOCUS tag info: $_\n") unless $blen>1;
   }
 elsif (m/^DEFINITION\s+(.+)/) {
   $descr=$1;
   while (<>) {
     chomp;
     if (m/^[A-Z]+\s+/) {
        $descr=~tr/\t \n\r/ /s;
        goto NEXTLINE;
        }
     $descr.=' '.$_;
     }
   $descr=~tr/\t \n\r/ /s;
   next;
   }
 elsif (m/^ACCESSION\s+(\S+)/) {
   $has_dbxref=$fxref ? 0 : 1;
   $bACC=$1;
   next;
   }   
 elsif (m/^VERSION/) {
   ($bACC, $bGI)=(m/([\w\.]+)\s+GI:(\d+)/);
   die("Error parsing GI# from VERSION tag ($_)") unless $bGI>1;
   next;
   }
 elsif (m/^SOURCE\s+(.+)/) {
   $org=$1;
   next;   
   }
 elsif (m/^COMMENT\s+(.+)/) {
   $comment=$1;
   while (<>) {
     chomp;
     if (m/^\s*[A-Z]+\s+/) {
        goto NEXTLINE;
        }
     $comment.=' '.$_;
     }
   next;
   }
 elsif (m/^FEATURES\s+/) {
   $freading=1;
   next;
   }
# elsif (m/^BASE COUNT/) {
#   $freading=0;
#   next;
#   }
 elsif (m/^CONTIG\s+/) {
   $freading=0;
   $seqreading=0;   
   }
 elsif (m/^ORIGIN\s+/) {
   $freading=0;
   $seqreading=1;
   next;
   }
 elsif (m/^\/\/$/) {
   $freading=0;
   $seqreading=0;
   my $defline=$bACC;
   my ($ctype)=($comment=~m/^([A-Z]+)\s+/);
   my $rtype;
   if ($ctype && ($rtype=$refseqtype{$ctype})) {
     $rtype.='x' if $comment=~m/is supported by experimental/;
     #$defline=$rtype.'|'.$bACC;
     }
    else {
     $rtype='gb';
     }
   $defline=$rtype.'|'.$bACC;
   # non-refseq entries will have rtype set to 
   # so they will NOT be shown
   my $noskip=(length($fseq)>9 && $typeflt{$rtype}==1 && $has_dbxref);
   if ($noskip)  {
     $descr=~s/\,\s+mRNA\.$//;
     $descr=~s/ \[[\w\,\. ]+\]\.$//;
     $descr=~s/\.$//;
     $descr=~s/\s*\(\Q$genesym\E\)// if $genesym;
     $descr=~s/^$org\s+//;
     if ($cdsprod) {
        $cdsprod=~tr/\t \n\r/ /s;
        if (!$descr || $descr=~m/hypothetical/i) {
            $descr=$cdsprod;
            }
           else {
            #print STDERR "{$descr}\n{$cdsprod}\n";
            $descr.=' product:'.$cdsprod
              unless (index(lc($descr), lc($cdsprod))>=0 || $cdsprod=~m/hypothetical protein/i);
            }
         }
     $defline.=' gid:'.$genesym if $genesym;
     $defline.=' '.$geneID if $geneID; 
     $defline.=' '.$cds if $cds;
     $defline.=' '.$descr;
     $defline.=' {'.$org.'}';
     $defline =~ tr/ \t/  /s;
     print &fastafmt(\$fseq, $defline);
     }
   $fseq='';
   $descr='';
   $bseq='';
   $genesym='';
   $geneID='';
   $has_dbxref=$fxref ? 0 : 1;
   $comment='';
   $cds='';
   $cdsprod='';
   $comment='';
   next;
   }
   
 #-- remaining lines parsing:
 
 if ($freading) {
   my ($featname, $frange)=(m/^\s+([\w\'\-]+)\s+(.+)/);
   die("Error parsing $bACC feature at: $_\n") unless $featname && $frange;
   if ($frange=~m/\(/ && $frange!~m/\)$/) { #wrapped composite intervals line
     do {
       $_=<>;
       last unless $_;
       chomp;
       $frange.=$_;
       } until (m/\)$/);
     } #unwrapped the composite intervals (multi-exons)
   #if ($featname eq 'CDS') {
   # $cds.='|' if $cds;
   # $frange=~s/\.\./-/g;
   # $cds.='CDS:'.$frange;
   # }
   
   my @attrs; 
   my $endattr;
   #print STDERR "..reading attributes for $featname ($frange)\n";
   while (<>) { #read all attributes for this feature
      unless (m/^\s+\/\w+/) {
        $endattr=1;
        last;
        }
      chomp;
      if (m/^\s+\/([\w\-]+)=(.+)/) { # /attr=value
         my ($a, $av)=($1, $2);
         my ($fattr, $fattr_dta)=&parseAttr($featname, $a,$av); #this will read the next line(s) if necessary
         if ($fattr eq 'gene' && !$genesym) {
           $genesym=$fattr_dta;
           #$gname =~ tr/"//d; #"
           next;
           }
         if ($featname eq 'Protein' && $fattr eq 'product') {  
           $cdsprod=$fattr_dta;
           #    if (!$descr || $descr=~m/hypothetical/i) {
           #            $descr=$fattr_dta;
           #            }
           #    else { $cdsprod=' product:'.$fattr_dta; }
           next;
           }
         if ($has_dbxref==0 && $fattr eq 'db_xref') {
           if ($has_dbxref=exists($xr{$fattr_dta})) {
             $geneID=$fattr_dta;
             }
           next;
           }
         if ($fattr eq 'organism') {
           $org=$fattr_dta;
           #$org=~tr/"//d; #"
           }
         } # /attr=attrvalue parsing
         # store all attributes and their values -- for the current feature only                   
         #elsif ($fattr ne 'translation') {
         # push(@attrs, $fattr.'='.$fattr_dta);
         # }
      } #attribute parsing loop
   #&writeFeature($featname, $frange, $gname, @attrs);
   goto NEXTLINE if $endattr;
   next;
   }# FEATURES section reading
   
 if ($seqreading && m/^\s+\d+\s+(.+)/) {
   my $s=$1;
   $s=~tr/ \t\n//d;
   $fseq.=uc($s);
   }
   
}#while input lines

#close(FASTAOUT) if $fastaout;


sub parseAttr {
 my ($fn, $a, $adta)=@_;
 #print STDERR "[[ parseAttr($fn, $a, $adta) ]]\n";
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
 $adta=~tr/\t \n\r/ /s;
 return ($a, $adta);
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
