#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 gap_mod.pl [-b[{A|G|C|X}:]out_fasta_repl] [{ -o gapinfo.txt \
    | -I gapinfo.txt [-c SoapScaf2Ctg.txt]] }] fasta_file

Unless -I option is used, for each gap found in the input sequences 
the following tab-delimited columns are shown (or printed in file gapinfo.txt):
 1. sequence name
 2. sequence location of the gap
 3. length of the gap
 4. 10 bases to the left of the gap (prefix)
 5. 10 bases to the right of the gap (sufix)

If -r option is provided, gaps are "filled" with a single base (default T)
replacing all Ns and the modified sequence is written to out_fasta_repl

If the gapinfo file is provided as input (-I option), then a restoration 
of the gaps in the given fasta file is attempted (-b can be used to specify
the base that was used for gap filling before); the SoapScaf2Ctg file (-c 
option) must be produced by ctg2Soap.pl and should be used to provide info
about the location and the orientation of the original contigs (assuming 
that fasta_file contains scaffolds constructed from the original contigs).
/;
umask 0002;
getopts('b:I:c:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $replIn=$Getopt::Std::opt_I;
my $replOut=$Getopt::Std::opt_b;
my $scaf2ctg=$Getopt::Std::opt_c;

my $replBase='T';

if ($replOut=~m/^([ACGX]):(.+)/) {
 $replBase=$1;
 $replOut=$2;
}
### - for gap restore case (-I):
my %gapInfo;    # scaffold_id => [ [ $pos, $len, $prefix, $suffix ], [ ... ] ]
my @gapScaffs; # just an enumeration of scaffolds IDs in gapInfo

my %restoreTargets; # new_scaffold_id => [ [$gapPos, $gapSeqPrint, $Npos, $gaplen], ...]
       # $gapPos=-1 for true new scaffolds, otherwise it's adjusted for prefix/suffix/strand
       # $gapSeqPrint includes prefix and suffix seq
       # $Npos is the gap offset (position) in $gapSeqPrint
if ($replIn) {
  #load gap info
  $replBase=$Getopt::Std::opt_b if $Getopt::Std::opt_b;
  open(INF, $replIn) || die("Error reading gap info file $replIn\n");
  while(<INF>) {
   chomp;
   my ($scaff, $pos, $len, $prefix, $suffix)=split(/\t/);
   if (exists($gapInfo{$scaff})) {
     push(@{$gapInfo{$scaff}}, [$pos, $len, $prefix, $suffix]);
   }
   else {
     $gapInfo{$scaff}=[ [$pos, $len, $prefix, $suffix] ];
     push(@gapScaffs, $scaff);
   }
   #$gapSeq{$scaff.'>'.$pos}=$prefix.($replBase x $len).$suffix;
   #$gapRestore{$scaff.'>'.$pos}=$prefix.('N' x $len).$suffix;
   #$gapPos{$scaff}=$pos-length($prefix);
  }
  close(INF);
  if ($scaf2ctg) { #new scaffold mapping to old contigs/scaffolds
    open(SINFO, $scaf2ctg) || die("Error reading scaf2ctg info file $scaf2ctg\n");
    while(<SINFO>) {
     chomp;
     my @t=split(/\t/);
     my $newscaf=shift(@t);
     my @tgaps; #target gaps in this new scaffold: list of [ $gPos, $gapSeqPrint, $Nofs, $NLen]
     if (@t>1) { # new scaffold
       foreach my $c (@t) {
         my ($oldscaf, $strand, $start, $end)=($c=~m/^(\S+) ([\+\-])\[(\d+)\:(\d+)\]/);
         my @ogaps=getGaps($oldscaf, $strand);
         push(@tgaps, @ogaps) if @ogaps>0;
       } # for each component contig (scaffold)
     $restoreTargets{$newscaf}=[@tgaps] if @tgaps>0;
     }
     else { #straight copy, not a new scaffold, just a new name
       @tgaps=getGaps($t[0], '+', 1);
       $restoreTargets{$newscaf}=[@tgaps] if @tgaps>0;
     }
    }
    close(SINFO);
  }#if $scaf2ctg info given
  else { #same contigs given for gap restore
    foreach my $scaf (@gapScaffs) {
       my @tgaps=getGaps($scaf, '+', 1);
       die("Error: no gaps retrieved for $scaf ?!?\n") unless @tgaps>0;
       $restoreTargets{$scaf}= [ @tgaps ];
    }
  }
}

if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
if ($replOut) {
  open(OUTR, '>'.$replOut) || die("Error creating output fasta file $replOut\n");
}
# --
my $seqname;
my $seqhdr;
my $seq;
if ($replIn) { #restore gap
  while (<>) {
    if (m/^>(\S+)(.*)/) {
      my ($sn,$sh)=($1, $2);
      gapRestore();
      ($seqname, $seqhdr)=($sn,$sh);
    } else {
     tr/\n\r\t //d;
     $seq.=uc($_);
    }
  }
  gapRestore();
}
else { #!generate gap info
  while (<>) {
    if (m/^>(\S+)(.*)/) {
      my ($sn,$sh)=($1, $2);
      procSeq();
      ($seqname, $seqhdr)=($sn,$sh);
    } else {
     tr/\n\r\t //d;
     $seq.=uc($_);
    }
  }
  procSeq();
}

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

if ($replOut) {
 close(OUTR);
}

#---------------------

sub getGaps {
  my ($oldscaf, $strand, $usepos)=@_;
  my @tgaps; #list of [$gapPos, $gapSeqPrint, $Npos, $gaplen]
  my $gd=$gapInfo{$oldscaf};
  if ($gd) { #was it a gapped scaffold previously?
     foreach my $ginfo (@$gd) {
        my ($gpos, $glen, $gprefix, $gsuffix)=@$ginfo;
        my $gapseq=$gprefix.($replBase x $glen).$gsuffix;
        my $posadj=length($gprefix);
        if ($strand eq '-') {
          $gapseq=reverseComplement($gapseq);
          $posadj=length($gsuffix);
        }
        if ($usepos) { $gpos-=$posadj+1; }
                else { $gpos=-1; }
        push(@tgaps, [$gpos, $gapseq, $posadj, $glen] );
     }
  }
  return (@tgaps);
}


sub procSeq {
  return unless $seqname && $seq;
  #print "[DBG]:>SEQ>$seqname\t".length($seq)."\n";
  my @polys; #array of [$poly, $p_start, $p_end] for each base in @PChars (A,C,G,T,N)
  my $cseq=$seq;
  my $p=0;
  while (length($cseq)>=2) {
       if ($cseq=~m/(NN+)/) {
        my ($poly, $p_start, $p_end)=($1, $-[0], $+[0]);
        $p+=$p_start;
        my $pl=$p_end-$p_start;
        my $lstart=$p_start-10;
        my $rstart=$p_end;
        my $rend=$rstart+10;
        my $maxend=length($cseq)-1;
        $rend=$maxend if $rend>$maxend;
        my $rlen=$rend-$rstart;
        $lstart=0 if $lstart<0;
        my $llen=$p_start-$lstart;
        #                0     1   2(prefix)                      3(suffix)
        push(@polys, [ $p+1,  $pl, substr($cseq, $lstart, $llen), substr($cseq, $rstart, $rlen)]);
        $cseq=substr($cseq, $p_end);
        $p+=$pl;
       }
       else { last; }
  } #while gap
  foreach my $pd (@polys) {
    #print "[DBG]:$seqname:poly match=$poly($p_start-$p_end), seqlen=".length($seq)."\n";
    my ($pos, $len, $prefix, $suffix)=@$pd;
    print join("\t",$seqname, $pos, $len, $prefix, $suffix)."\n";
    if ($replOut) {
      substr($seq, $pos-1, $len)=$replBase x $len;
    }
  }
  if ($replOut) {
     print OUTR ">$seqname$seqhdr\n".
       join("\n", unpack("A100"x((length($seq)-1)/100+1),$seq))."\n";
  }
  $seq='';
  $seqname='';
  
}

sub gapRestore {
  return unless $seqname && $seq;
  my $td=$restoreTargets{$seqname};
  my $dbg=1;
  if ($td) {# this scaffold had gaps, or it's build from contigs with gaps
    foreach my $gd (@$td) {
      my ($gPos, $gapSeqPrint, $Npos, $gaplen)=@$gd;
      ##debug only
      if ($seqname eq 'scaffold1') {
       print STDERR "searching ($dbg of ".scalar(@$td)."): $gapSeqPrint\n";
       $dbg++;
      }
      my $printlen=length($gapSeqPrint);
      if ($gPos>=0) {
        die("Error: gap sequence print not found at expected location! ($seqname, $gPos)\n")
          if substr($seq, $gPos, $printlen) ne $gapSeqPrint;
      }
      else { #just search for it in the whole sequence
        $gPos=index($seq, $gapSeqPrint);
        if ($gPos<0) {
           print STDERR "Warning: couldn't restore gap (len $gaplen) in scaffold $seqname\n";
           next;
        }
      }
     my $gapRSeq=$gapSeqPrint;
     substr($gapRSeq, $Npos, $gaplen) = 'N' x $gaplen; # restore gap 
     substr($seq,  $gPos, $printlen) = $gapRSeq; 
    }#foreach former gap
  } # had gaps
  
  print ">$seqname$seqhdr\n".
       join("\n", unpack("A100"x((length($seq)-1)/100+1),$seq))."\n";
  $seq='';
  $seqname='';
  $seqhdr='';
}

sub reverseComplement {
  my $s=reverse($_[0]);
  $s =~ tr/AaCcTtGgUuMmRrWwSsYyKkVvHhDdBb/TtGgAaCcAaKkYyWwSsRrMmBbDdHhVv/;
  return $s;
 }
