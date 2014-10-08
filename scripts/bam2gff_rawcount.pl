#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 bam2gff_rawcount.pl [-o <outfile] <bamfile> <transcripts.gff>
 
 Reports the raw number of read mappings overlapping each transcript 
 found in the given <gff> file (exon-wise)
 
 BAM file MUST have been indexed prior to running this script.
 (using 'samtools index <bamfile>')
 
/;
umask 0002;
getopts('o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $bamfile=shift(@ARGV) || die "$usage\n";
my $f_gtf=shift(@ARGV) || die "$usage\n";
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }


my %gffrecs;
# ---------------  0     1       2       3         4      5       6       7    
#     recID =>  [ chr, strand, fstart, fend,  [@exons], [@cds], geneID, \@attrs ]

my $gffh;
open($gffh, $f_gtf) || die("Error: failed to open $f_gtf. $!\n");
print STDERR "Loading GFF data from $f_gtf..\n";
my @tlist; #list of transcript_IDs, in order they are encountered

loadGff($gffh, \%gffrecs, \@tlist);
# we don't need to sort by location, we'll go through every transcript anyway
#my @tsorted = sort sortByLoc @tlist;


foreach my $tid (@tlist) {
 #for every transcript, print the number of read mappings which overlap at least one exon
 #samtools view will be called for every transcript, so this can be slow
 my $td=$gffrecs{$tid};
 my $ovlcount=bamOvlCount($bamfile, $tid, $td);
 print join("\t",$tid, $$td[1], $$td[0].':'.$$td[2].'-'.$$td[3], $$td[6], $ovlcount)."\n";
 }


#-------------
sub loadGff {
 my ($fh, $recs, $reclist)=@_;
 # ---------------------  0     1       2       3         4      5       6       7    
 #           recID =>  [ chr, strand, fstart, fend,  [@exons], [@cds], geneID, \@attrs ]
 # reclist = just a list of IDs in order they are discovered
 while (<$fh>) {
   next if m/^\s*#/;
   chomp;
   my ($chr, $track, $f, $fstart, $fend, $fscore, $strand, $frame, $attrs)=split(/\t/);
   next unless $fstart>1 && $attrs;
   next if $f eq 'gene' || $f =~ m/locus/; # Warning: skipping any 'gene' or 'locus' features, unconditionally
   my ($tid, $geneid, @lattr);
   my %hattr;
   #if ($attrs=~m/nearest_ref[= ]+"?([^;"]+)/) {
   #   $refid=$1;
   #   }
   #if ($attrs=~m/class_code[\.= ]+"?([\w=])/) {
   #   $classcode=$1;
   #   }
   @lattr=(split(/;\s*/, $attrs));
   my $gff3_ID;
   my $gff3_Parent;
   ($fstart, $fend)=($fend, $fstart) if $fend<$fstart;
   
   ($gff3_ID)=($attrs=~m/\bID=([^;]+)/);
   ($gff3_Parent)=($attrs=~m/\bParent=([^;]+)/);
   
   if ($gff3_ID || $gff3_Parent) {
      # GFF3 format
      @lattr = map { [split(/\s*=\s*/, $_, 2)] } @lattr;
      %hattr = map { lc($_->[0]) => $_->[1]  } @lattr;
      
      $gff3_ID =~ tr/"//d;
      $gff3_Parent =~ tr/"//d;
      $gff3_Parent='' if ($f eq 'mRNA');
      $geneid=$hattr{'gene_id'} || $hattr{'geneid'} || $hattr{'gene_name'};
      $geneid='' if $geneid eq $gff3_ID;
      if ($gff3_ID && !$gff3_Parent) { #top level feature
         $tid=$gff3_ID;
         die("Error: duplicate feature $gff3_ID on $chr\n") # if (exists($$recs{"$chr|$gff3_ID"}));
            if (exists($$recs{$gff3_ID}));
         #my $recID=$gff3_ID;
         #$geneid=$refid if $refid;
         $$recs{$tid} = [$chr, $strand, $fstart, $fend, [], [], $geneid, [@lattr]];
         push(@$reclist, $tid);
         next;
         } # parent/top-level feature
      $tid=$gff3_Parent;
      } #GFF3
     else { #GTF format
      if ($track=~/^jigsaw/ && $attrs=~m/^\d+$/) {
        $tid=$chr.'.jsm.'.$attrs;
        }
       else {
        @lattr = map { [split(/\s+/,$_,2)] } @lattr;
        %hattr = map { lc($_->[0]) => $_->[1]  } @lattr;
        $tid  = $hattr{'transcript_id'} || die("Error: invalid GTF (transcript_id missing: $_)\n");
        $tid =~ tr/"//d; #"
        }
      if ( $geneid = $hattr{'gene_id'} ) {
        $geneid=~tr/"//d; #"
        }
      # gtf with parent 'transcript' feature
      if ($f eq 'transcript') {
           die("Error: duplicate feature $tid on $chr\n") # if (exists($$recs{"$chr|$gff3_ID"}));
            if (exists($$recs{$tid}));
           $$recs{$tid} = [$chr, $strand, $fstart, $fend, [], [], $geneid, [@lattr]];
           push(@$reclist, $tid);
           next;
           }
      } #GTF
   #only 'exon' and 'CDS' sub-features are recognized
   next unless $f eq 'exon' || $f eq 'CDS'; 
   # -----------exon/CDS line here:
   my $ld = $$recs{$tid};
   my $sidx=($f eq 'CDS') ? 5 : 4;
   if ($ld) { #update existing entry
      my ($lstart, $lend)=($$ld[2], $$ld[3]);
      $$ld[2] = $fstart if $fstart < $lstart;
      $$ld[3] = $fend if $fend > $lend;
      }
    else { # first time seeing this transcript
      $$recs{$tid} = [$chr, $strand, $fstart, $fend, [], [], $geneid, [@lattr] ];
      push(@$reclist, $tid);
      }
   push(@{$$ld[$sidx]}, [ $fstart, $fend, $fscore, $frame ]);
 } #while <$fh>
} # loadGff

sub sortByLoc {
 my $da=$gffrecs{$a};
 my $db=$gffrecs{$b};
 if ($$da[0] eq $$db[0]) {
    return ($$da[2]==$$db[2]) ? $$da[3] <=> $$db[3] : $$da[2] <=> $$db[2] ;
    }
  else { return $$da[0] cmp $$db[0] ; }
}

sub bamOvlCount { # get the number of reads overlapping a transcript
 my ($fbam, $tid, $tdata)=@_;
 my ($chr, $cstrand, $cstart, $cend, $exons, $cds, $geneid, $attrs)=@$tdata;
 my $range=$chr.':'.$cstart.'-'.$cend;
 my $segdata=$exons;
 $segdata=$cds unless @$segdata>0;
 my @exons;
 if (@$segdata>0) {
    @exons= sort { $a->[0]<=>$b->[0] } @$segdata;
    }
  else {
    @exons=([$cstart, $cend]);
    }
 
 open(SAMPIPE, "samtools view $fbam $range |") || die ("Error opening samtools pipe ($!)\n");
 my $numovl=0; #number of reads overlapping at least one exon of the current transcript
 while(<SAMPIPE>) {
   my $samline=$_;
   chomp;
   my ($qname, $flags, $gseq, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual, @extra)=
      split(/\t/);
   # if ($strand) {
   #       my $mstrand= (($flags & 0x10)==0) ? '+' : '-';
   #       next if $mstrand ne $strand;
   #       }
   #now extract the CIGAR segments
   my @cigdata=($cigar=~m/(\d+[A-Z,=])/g);
   my ($mstart, $mend);
   my $curpos=$pos;
   $mstart=$pos;
   foreach my $cd (@cigdata) {
     my $code=chop($cd);
     #now $cd has the length of the CIGAR operation
     if ($code eq 'N') { #gap found (possible intron)
        #process previous interval
        if ($mend && checkOverlap($mstart, $mend, \@exons)) {
           $numovl++;
           $mend=0;
           last; #exon overlap found, no need to check the other exons
           }
           
        $curpos+=$cd;
        $mstart=$curpos; #start of the next exon
        $mend=0;
        next;
        }
      if ($code eq 'M' || $code eq 'D') {
        #only advance genomic position for match and delete operations
        $curpos+=$cd;
        $mend=$curpos-1; #advance the end of the exon
        }
     }
   #check the last interval
   if ($mend && checkOverlap($mstart, $mend, \@exons)) {
         #$hasOvl=1;
         $numovl++;
         }
   } # while <SAMPIPE>
  close(SAMPIPE);
  return $numovl; 
}

sub checkOverlap { #check if segment $a-$b overlaps any of the exons in $rx
 my ($a, $b, $rx)=@_;
 return 0 if ($a>$$rx[-1]->[1] || $b<$$rx[0]->[0]); # not overlapping the transcript region at all
 foreach my $x (@$rx) {
   #check for exon overlap
   return 1 if ($a<=$$x[1] && $b>=$$x[0]);
   return 0 if $b<$$x[0]; #@$rx is sorted, so we can give up
   }
}
