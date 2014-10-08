#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 gff_get_introns.pl [-m <minlen>] <gff_data_stream..>
 Extract intron coordinates from given GFF data stream.
 Output this tabulated fomat:
 
 <chr><strand> <start> <end> <intron_id> <transcript_info>
 
 Options:
  -m minium intron length to extract (default 100000)
/;
umask 0002;
getopts('m:o:') || die($usage."\n");
my $minlen=$Getopt::Std::opt_m || 100000;
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# ---------------------      0         1          2      3        4      5       6       7
my %recs; # recID =>  [ chr.strand, feat_type,  gname, tdescr, fstart, fend, [@exons], [@cds] ]
             #   
my $curtag; #chromosome and strand for current model
my ($gname, $tdescr);
my @exd; #exons for current model
while (<>) {
   next if m/^\s*#/;
   chomp;
   my ($chr, $track, $f, $fstart, $fend, $fscore, $strand, $frame, $lnum)=split(/\t/);
   next unless $fstart>1 && $lnum;
   next if $f eq 'gene'; # Warning: skipping any 'gene' features, unconditionally
   my $gff3_ID;
   my $gff3_Parent;
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
         die("Error: duplicate feature $gff3_ID on $chr\n") if (exists($recs{"$chr|$gff3_ID"}));
         my $recID="$chr|$gff3_ID";
         $curtag=$chr.$strand;
         $recs{$recID} = [$curtag, $f, $gname, $tdescr, $fstart, $fend, [], [] ];
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
   $curtag=$chr.$strand;
   $recID=$chr.'|'.$recID;
   my $ld = $recs{$recID};
   if ($ld) { #existing entry
     my $i=($f eq 'CDS') ? 7 : 6;
     my ($lstart, $lend)=($$ld[4], $$ld[5]);
     $$ld[4]=$fstart if $fstart<$lstart;
     $$ld[5]=$fend if $fend>$lend;
     push(@{$$ld[$i]}, [$fstart, $fend, $fscore]);
     }
    else { # first time seeing this locus/gene
     $recs{$recID} = ($f eq 'CDS') ? 
           [$curtag, $f, $gname, $tdescr, $fstart, $fend, [], [[$fstart, $fend, $fscore]] ] :
           [$curtag, $f, $gname, $tdescr, $fstart, $fend, [[$fstart, $fend, $fscore]], [] ] ;
     }
} #while <>

writeModels();

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************
sub writeModels {
 return if keys(%recs)==0;
 while ( my ($l, $ld)=each(%recs) ) {
   # my $ld=$recs{$l} || die ("Error: locus $l found in list but not in hash!\n");
   my ($chrstrand, $ftype, $gname, $descr, $lstart, $lend, $er, $cr) = @$ld;
   my ($mstart,$mend)=($lstart, $lend);
   my @ex;
   if (@$er<1 && @$cr>0) { 
     @ex = sort { $main::a->[0] <=> $main::b->[0] } @$cr;
     }
    else {
     @ex = sort { $main::a->[0] <=> $main::b->[0] } @$er;
     }
   # -- now look at the introns
   next if @ex<2;
   my @in; #introns to be printed stored here
   #my $icount=0; #intron output counter
   ($mstart, $mend) = ($ex[0]->[0], $ex[-1]->[1]);
   for (my $i=1;$i<@ex;$i++) {
     my ($istart, $iend)=($ex[$i-1]->[1]+1,$ex[$i]->[0]-1);
     next unless ($iend-$istart+1>=$minlen);
     push(@in, [$istart, $iend]);
     #$icount++;
     my $ann="";
     $ann=" gene:$gname" if $gname;
     $ann.=" product:$descr" if $descr;
     print join("\t", $chrstrand, $istart, $iend, $l.".i$i", $ann)."\n";
     }
  } #for each stored transcript
}

