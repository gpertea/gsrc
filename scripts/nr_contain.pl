#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 nr_contain.pl [-M] [-o <outfile.cl>] <chr_mapping_intervals>

 Performs containment clustering of genomic intervals\/mappings
 Input is expected to be sorted by chromosome, start, end
 
 Expected input format:
 
 chr[strand] start end mappingID mapping_defline
 
 Options: 

 -M   merge overlapping intervals into single segments
 
/;
# perl -ne \
# '@t=split(/\t/);($c,$a,$b,$s)=m/ ([^\:]+)\:(\d+)\-(\d+)([\-\+])/;print join("\t",$c.$s,$a,$b,$t[0],$t[2])."\n"' \
# hg19_primary_transcripts.fsize | sort -k1,1 -k2,2n -k3,3n > hg19_primary_transcripts.tab
# perl -pe 's/[\-\+]\t/\t/' hg19_primary_transcripts.tab |sort -k1,1 -k2,2n -k3,3n > hg19_pri_rna.nostrand.tab


#After clustering, get the sequence of the container intervals like this:

# grep '^>' containers.cl | perl -ne 'chomp;@t=split(/\s+/,$_,5);$t[0]=~s/^>//; $t[0]=~s/^[^\|]+\|//; \
#   print "$t[1]\t$t[2]-$t[3]\t$t[0]|$t[1]:$t[2]-$t[3] $t[4]\n"' | seqmanip -f - hg19.fa > hg19_cseqs.fa

umask 0002;
getopts('Mo:') || die($usage."\n");
my $merge=$Getopt::Std::opt_M;
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
my %idata; #interval data: idata{id}=>[chr, start, end, annotation]
my $max_end; #max seen right end for current chromosome, still buffered
my $max_end_id; #id of maximum right end entry in buffered region
my $max_end_idx; #idx of maximum right end entry in @reg
my @reg; # data in buffered region starting < $max_end[2]
 # array of [id, chr, start, end, included_ids...]
 #            0   1     2      3   4...
my $chr; 
while (<>) {
 chomp;
 my ($c, $cstart, $cend, $id, $idann)=split(/\t/); 
 if ($chr ne $c) {
    flush_region();
    %idata=();
    $chr=$c;
    }
   else { 
    flush_region() if $cstart>=$max_end; 
    }
 $idann=~s/\w+[\+\-]?\:\d+\-\d+[\+\-]?//;
 $idata{$id}=[$c, $cstart, $cend, $idann];
 unless ($max_end) {
   push(@reg, [$id, $c, $cstart, $cend]);
   ($max_end, $max_end_id, $max_end_idx)=($cend, $id, $#reg);
   next;
   }
 if ($merge) {
   #check for overlaps with other intervals in @reg
   if ($cstart<$max_end-3) {
      push(@{$reg[$max_end_idx]}, $id);
      if ($cend>$max_end) {
         $max_end=$cend;
         $reg[$max_end_idx]->[3]=$max_end; #adjust ending of cluster region
         }
      next;
      }
   }
 else {
   # check for containment IN other intervals of @reg
   if ($cend<=$max_end+3) {
      #merge into max_id
      push(@{$reg[$max_end_idx]}, $id);
      if ($max_end<$cend) {
          $max_end=$cend;
          $reg[$max_end_idx]->[3]=$max_end;
          }
      next; #do not store as independent member of @reg
      }
   # check for containment OF other intervals of @reg, take over
   if ($cstart<=$reg[$max_end_idx]->[1]+3 && $cend>$max_end) {
      $max_end=$cend;
      push(@{$reg[$max_end_idx]}, $reg[$max_end_idx]->[0]);
      $reg[$max_end_idx]->[0]=$id;
      $reg[$max_end_idx]->[3]=$max_end;
      next;
      }
   }
 # no overlap/containment, add to region buffer
 push(@reg, [$id, $c, $cstart, $cend]);
 if ($cend>$max_end) {
   ($max_end, $max_end_id, $max_end_idx)=($cend, $id, $#reg);
   }
}

flush_region();

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

sub flush_region {
 return unless @reg>0;
 foreach my $d (@reg) {
  my @t=@$d;
  my $clid=shift(@t);
  next unless $clid; # could be a contained one, skip it
  print ">$clid ".shift(@t).' '.shift(@t).' '.shift(@t);
  my $ida=$idata{$clid} || die("Error retrieving data for container id $clid!\n");
  print " ".$$ida[3] if $$ida[3];
  print "\n";
  if (@t>0) {
    foreach my $id (@t) {
       my $dv=$idata{$id} || 
          die("Error retrieving data for id $id!\n");
       print join(' ',$id, $$dv[1], $$dv[2], $$d[3])."\n";
       }
    }
  }
 # reset region data
 ($max_end, $max_end_id, $max_end_idx)=(undef,undef,undef);
 @reg=();
}
