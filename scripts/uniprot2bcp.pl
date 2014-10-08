#!/usr/bin/perl
use strict;
use Getopt::Std;
#use LWP::Simple;
use FindBin;use lib $FindBin::Bin;
use dbSession;
my $usage = qq/Usage:
 uniprot2bcp.pl [-t <taxons.tab>] [-F] uniref100_clusters.pfa [<uniprot_data_file_or_stream>]
 
 Parses the UNIPROT data into bcp files for loading into geanno db.

 Also creates a uniprot fasta file (uniprot.fa) with all the
 uniprot entries found in the input data\/stream

 If no database connection is possible, use -t option to provide a file
 with taxon information; <taxons.tab> can be built with 
 pullsql -b common\@IGM -o taxons.tab 
 .. with the query:
  select tax_id, name, c_name from taxon
  
 If -F option is given, ONLY the fasta file will be created. 
/;
umask 0002;
getopts('t:F') || die($usage."\n");
my $fastaOnly=$Getopt::Std::opt_F;
my $ftaxon=$Getopt::Std::opt_t;
my @tax;  # $tax[taxon_id] = [ taxon_name, common_name ] #ncbi taxon ids
my %taxnames;# taxon_name => first tax_id found
my %tmiss; # UniProt taxon_id => uniprot organism name
my %tfixed; # uniprot taxon_id => NCBI taxon_id (through name)
my @tmisslst; # list of missed taxon_ids, in the order they were found
if ($ftaxon) {
 open(TAX, $ftaxon)||die("Error opening $ftaxon\n");
 while (<TAX>) {
  chomp;
  my @t=split(/\t/);
  $tax[int($t[0])]=[$t[1],$t[2]];
  my $tname=lc($t[1]);
  $tname=~s/^[\s\-\.]+//;
  $tname=~s/[\s\-\;\.]+$//;
  $taxnames{$tname}=int($t[0]) unless exists $taxnames{$tname};
  }
 close(TAX);
 print STDERR "Taxonomy data loaded (".scalar(@tax)." entries).\n";
 }
#create the bcp files for loading

my $funiref=shift(@ARGV) || die("$usage\n");

my $ds;
unless ($ftaxon) {
 $ds=dbSession->new('common@IGM');
 }

my %uniref; # uniprot_id => uniref_cluster_id
open(UNIREF, $funiref) || die("Error opening $funiref file!\n");
while (<UNIREF>) {
 chomp;
 s/^>//;
 my @ids=split;
 unless (@ids>1) {
   print STDERR "WARNING: cluster: '$_' has no members!\n";
   next;
   }
 my $clname=shift(@ids);
 foreach my $upid (@ids) {
   $uniref{$upid}=$clname;
   }
 }
close(UNIREF);
print STDERR "UniRef clusters loaded, now parsing the main uniprot KB data..\n";

open(UPFA, '>uniprot.fa') || die ("Error creating file uniprot.fa!\n");
unless ($fastaOnly) {
 open(UNIPROT, '| pbzip2 -p2 -c > uniprot.bcp.bz2') || die ("Error creating file uniprot.bcp!\n");
 open(UNAMES, '| pbzip2 -p2 -c > uniprot_names.bcp.bz2') || die ("Error creating file uniprot_names.bcp!\n");
 open(UXREF,  '| pbzip2 -p2 -c > uniprot_xref.bcp.bz2') || die ("Error creating file uniprot_xref.bcp!\n");
}

#open(FGO, $file) || die ("Error opening $file\n");
my $sth;
unless ($ftaxon) {
 $sth=$ds->prep('select name, c_name from taxon where tax_id=?');
 
 }
my $c_id; #current uniprot ID
my $c_taxon;
my $c_rev=0;
my $c_len;
my $c_priacc;
my $c_descr;
my $c_gene;
my $c_uniref;
my $c_htaxon;
my $seq;
my $inSeq=0;
my $xcode=0;
my ($c_orgname, $c_comname);
while (<>) {
 if (m{^//}) {
   writeUP() if $c_id;
   next;
   }
if ($inSeq) {
  #within SQ record
  chomp;
  tr/\t //d;
  $seq.=$_;
  next;
  }
if (m/^SQ\s+SEQUENCE\s+(\d+)/) {
  my $seqlen=$1;
  die("Error: invalid length for $c_id ($c_len vs $seqlen)!\n") if ($seqlen != $c_len);
  $inSeq=1;
  next;
  }
if (m/^ID\s+(\S+)\s+(\S+)\s+(\d+)/) {
  my ($id, $rev, $len)=($1,$2,$3);
  writeUP() if $c_id;
  ($c_id, $c_len)=($id, $len);
  $c_rev=($rev=~/^Reviewed/)?'1':'0';
  next;
  }# ID
if (m/^AC\s+(\S.+)/) {
  my $a=$1;
  chomp;
  my @accs=split(/;\s*/, $a);
  $c_priacc=$accs[0] unless $c_priacc;
  unless ($fastaOnly) {
    foreach my $acc (@accs) {
     next unless $acc;
     print UXREF join("\t",$c_id, 'UP_ACC', '', $acc, '')."\n";
     }
    }
  next;
  }# AC
if (m/^DE\s+(\S.+)/) {
  my $d=$1;
  $c_descr.=' '.$d;
  next;
  }
if (m/^PE\s+(\d+)\:/) {
  $xcode=$1;
  next;
  }
if (m/^GN\s+(\S.+)/) {
  my $d=$1;
  my ($gname)=($d=~m/\bName\s*=\s*([\w\-\.]+)/);
  if ($gname) {
    $c_gene=$gname unless $c_gene;
    unless ($fastaOnly) {
     print UXREF join("\t", $c_id, 'UP_GENE', '', $gname, '')."\n";
     my ($gsyns)=($d=~m/\bSynonyms\s*=\s*(\S[\w\-\,\.\s]+)\;/);
     if ($gsyns) {
       my @gsyn=split(/\,\s*/,$gsyns);
       foreach my $gs (@gsyn) {
         print UNAMES join("\t",$gname, 'GSYN', $gs)."\n";
         }
       } #synonyms found
     }
    }#a gname found
  next;
  }#GN
if (s/^OS\s+//) {
   unless ($c_orgname) {
      ($c_orgname)=(m/^([\w\. \-]+)/);
       $c_orgname=~s/[\s\-]+$//;
       $c_orgname=~s/[\s\-\.\;]+$//;
      ($c_comname)=(m/\(([\w\. \-]+)/);
   }
  next;
  }
if (m/^OX\s+NCBI_TaxID\s*=\s*(\d+)/) {
 my $t=$1;
 $c_taxon=$t unless $c_taxon;
 print UXREF join("\t", $c_id, 'UP_Taxon', $t, '', '')."\n"
    unless $fastaOnly; #taxon
 next;
 }#OX (Taxon)

if (m/^OH\s+NCBI_TaxID\s*=\s*(\d+)/) {
 my $t=$1;
 $c_htaxon=$t unless $c_htaxon;
 print UXREF join("\t", $c_id, 'UP_Host', $t, '', '')."\n"
  unless $fastaOnly; #host taxon
 next;
 } #OH (Host taxon)
 
if (m/^DR\s+(\S.+)/) { # db xref 
 my $d=$1;
 $d=~s/\-?\.\s*$//;
 $d=~s/\;\s*$//;
 my @xd=split(/\;\s*/,$d);
 @xd=grep { !/^\-/ } @xd;
 my $xdb=shift(@xd);
 if ($xdb eq 'GO') {
   my $goid=shift(@xd);
   my ($gonum)=($goid=~/GO:(\d+)/);
   $gonum=int($gonum);
   die("Error parsing GO ID from DR line: $_\n") unless $gonum>0;
   my $gobranch=shift(@xd);
   my $go_ann_type=shift(@xd);
   print UXREF join("\t", $c_id, $xdb, $gonum, $goid, $go_ann_type.' '.$gobranch)."\n"
     unless $fastaOnly;
   }
 else {
   my $acc=shift(@xd);
   my $data=join(' ', @xd);
   print UXREF join("\t", $c_id, $xdb, '', $acc, $data)."\n" unless $fastaOnly;
   }
 next;
 } #DR 
} #while

writeUP() if $c_id;
unless ($fastaOnly) {
 close(UNIPROT);
 close(UNAMES);
 close(UXREF);
 }
close(UPFA);

if (@tmisslst>0) {
  open(TMISS, '>missing_taxons.tab');
  foreach my $t (@tmisslst) {
    print TMISS "$t\t".$tmiss{$t}."\n";
  }
  close(TMISS);
}

unless ($ftaxon) {
  $ds->logout();
  }
#-- now run bcpin

#my $bcpcmd='bcpin -TI -b geanno@NEOSYBASE';
#system($bcpcmd. ' uniprot.bcp uniprot_names.bcp uniprot_xref') 
#  && die("Error at bcpin!\n");

sub writeUP {
 # c_descr now can contain strange RecName: Full=...; and AltName: Full=...;
 # keep only the RecName Full, remove anything else
 if ($c_descr=~m/Name:\s+Full=\s*([^\;]+)/) { # RecName or SubName
   $c_descr=$1;
   }
  # elsif ($c_descr=~m/SubName:\s+Full=\s*([^\;]+)/) {
  # $c_descr=$1;
  # }
   
 prepFields($c_descr);
 # -- c_descr can probably be insanely long.. truncate it here to the max field length
 $c_uniref=$uniref{$c_id} || $uniref{$c_priacc};
 $c_taxon=int($c_taxon);
 my $upre=($c_rev==1)?'UPr|':'UP|';
 #$upre.='e' if $xcode<3;
 my $defline='';
 $defline.='|gid|'.$c_gene if $c_gene;
 $defline.=' '.$c_uniref if $c_uniref;
 $defline.=' '.$c_descr;
 my ($oname, $coname);
 if ($ftaxon) {
    my $td=$tax[$c_taxon];
    if ($td) {
       ($oname, $coname)=@$td;
       }
     #else {
     #  print STDERR "Warning: no data for taxon $c_taxon".
     #        " was loaded ($c_id|acc|$c_priacc $defline)!\n";
     #  }
    }
  else {
   $ds->exec($sth, $c_taxon);
   while (my $rr=$ds->fetch($sth)) { ($oname, $coname)=@$rr; };
   }
 my $n_taxon=$c_taxon; #NCBI taxon -> default is the same
 unless ($oname) {
   $n_taxon=$tfixed{$c_taxon};
   unless ($n_taxon) { #not fixed before
      unless (exists($tmiss{$c_taxon})) {
       #see if we can retrieve by name
       if ($c_orgname) {
          my $tname=lc($c_orgname);
          $tname=~s/^[\s\-\.]+//;
          $tname=~s/[\s\-\;\.]+$//;
          $n_taxon=$taxnames{$tname};
          if ($n_taxon) {
               #fixed it
               $tfixed{$c_taxon}=$n_taxon;
               }
             else { #can't find NCBI taxon ID
               $tmiss{$c_taxon}=$c_orgname;
               push(@tmisslst, $c_taxon);
               print STDERR "Warning: no NCBI data for taxon id $c_taxon".
                  " ($c_id|acc|$c_priacc $defline) - using '$c_orgname' instead.\n";
               }
          }
          else {
           #unexpected: no organism name
           $tmiss{$c_taxon}='NA';
           push(@tmisslst, $c_taxon);
           print STDERR "Warning: no organism name at all (?!) for taxon $c_taxon\n";
          }
       }
      }
   $oname=$c_orgname;
   } #no NCBI taxon id match
 unless ($fastaOnly) {
  my $wtaxon=$n_taxon || $c_taxon;
  print UNIPROT join("\t",
   $c_id, $wtaxon, $c_rev, $c_len, $c_priacc, $c_descr, $c_gene, $xcode, $c_uniref, $c_htaxon)."\n";
  }
 $defline.=' {'.$oname.'}' if $oname;
 $defline=~tr/;"\t/.' /; #"
 print UPFA '>'.$upre.$c_id.'|acc|'.$c_priacc.$defline."\n";
 print UPFA fastafmt(\$seq);
 $inSeq=0;
 ($c_id, $c_taxon, $c_rev, $c_len, $c_priacc, $c_descr, $c_gene, $xcode, $c_uniref, $c_htaxon)=
 (undef,    undef,    0,        undef,  undef,   undef,   undef,   0,      undef,     undef  );
 ($c_orgname, $c_comname)=(undef,undef);
 $seq=''; 
}

sub fastafmt {
 my $s=shift;
 my $slen=length($$s);
 my @lines=unpack("A60" x (int(($slen-1)/60)+1),$$s);
 return join("\n",@lines)."\n";
 }
 
sub prepFields {
 foreach (@_) {
  s/^\s+//;
  s/\s+$//;
  s/\.$//;
  tr/\t \n/   /s;
  }
 }
