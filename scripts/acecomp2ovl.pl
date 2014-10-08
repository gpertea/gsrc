#!/usr/bin/perl
use strict;

#input is the components file from ace2fasta.pl
#output is contig name followed by a list of overlap regions in contig coordinates

#my @ctglst;
#my %ctgdata;
use Getopt::Std;
my $usage = q/Usage:
acecomp2ovl.pl [-o <outfile>] [-F] < ace_file_comp
 Takes the component file produced by ace2fasta.pl 
 and produces a file with overlap regions of interest.
 
Options:
 -o <outfile> : write output to <outfile>
 -F           : consider overlaps only between components
                with names matching a hard coded filter in this script
/;
# comp2ovl.pl [-r <component_masked_seqs.fa.fai] < ace_file_comp
# Optionally can check for lower-case repeats in the components 
# (-r option) and thus exclude regions which are mostly repeats.

umask 0002;
getopts('Fr:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
#my $repfile=$Getopt::Std::opt_r;
my $ovlfilter=$Getopt::Std::opt_F;

#current contig data:
my ($ctg, $ctgn, $ctglen);
my @comps; # list of component contributions
my @allcomps; # unfiltered list of components (always available)
#my @fltcomps; # if $ovlfilter, list of component
my @regs;  #
while (<>) {
 if (m/^>(\S+)\s+(\d+)\s+(\d+)/) {
  my ($c, $cn, $cl)=($1,$2,$3);
  printContig($ctg, $ctgn, $ctglen, \@regs, \@comps, \@allcomps) if $ctg;
  @comps=();
  @allcomps=();
  @regs=();
  ($ctg, $ctgn, $ctglen)=($c, $cn, $cl);
  #push(@ctglst, $ctg);
  #$refctgdata=[$ctgn, $ctglen, [] ];
  #$ctgdata{$ctg}=$refctgdata;
  next;
 }
 if (m/^(\S+)\s+(\d+)\s+([\+\-])\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
  my ($read, $len, $strand, $rstart, $rend, $cstart, $cend)=($1, $2, $3, $4, $5, $6, $7);
  push(@allcomps, [$cstart, $cend, $read, $strand]);
  if ($ovlfilter) {
   #define custom filter here for component names
   next unless ($read=~m/^3ds_anchor_/ || $read=~m/^CONTIG\:\d+/);
  }
  #print STDERR ']'.join(" ",($read, $len, $strand, $rstart, $rend, $cstart, $cend))."\n";
  my @ovls;
  if (@comps>0) {
      #find all overlaps of this read with other reads
      foreach my $cspan (@comps) {
        my ($ol, $or)=getOvl([$cstart, $cend], $cspan);
        if ($or) { #overlap detected
          #print STDERR "\tovl $ol-$or found\n";
          my ($ml, $mr);
          if (@ovls>0) {
            #chances are this overlaps the last ovl found
            ($ml,$mr)=mergeOvl([$ol, $or], $ovls[-1]);
          }
          if ($mr) { 
            #print STDERR "\textend prev ovl $ovls[-1]->[0]-$ovls[-1]->[1] to $ml-$mr\n";
            $ovls[-1]=[$ml, $mr]; 
            }
          else { push(@ovls, [$ol, $or]); }; #new overlap
        } #overlap detected
      } #for each component read
      if (@ovls>0) {
       @ovls=sort { $a->[0] <=> $b->[0] } @ovls;
       my $check=1;
       while ($check) {
        for (my $i=1;$i<@ovls;$i++) {
          my ($ml, $mr)=mergeOvl($ovls[$i-1], $ovls[$i]);
          if ($mr) { splice(@ovls, $i, 1); $ovls[$i-1]=[$ml, $mr]; last; }
        }
        $check=0;
       }
      } #make sure @ovls is sorted and merged
  } #check overlap with existing components
  #print STDERR "ovls = ";printRegs(\@ovls);
  #print STDERR "vs regs = ";printRegs(\@regs);
  push(@comps, [$cstart, $cend, $read, $strand]);
  my $merged=0;
  if (@regs>0) {
    if (@ovls>0) {
      foreach my $reg (@regs) {
        for (my $o=0;$o<@ovls;$o++) {
          my ($ml, $mr)=mergeOvl($reg, $ovls[$o]);
          if ($mr) {
             #print STDERR "\tmerging reg $reg->[0]-$reg->[1] w/ ovl $ovls[$o]->[0]-$ovls[$o]->[1]\n";
             ($reg->[0], $reg->[1])=($ml,$mr);
             $merged=1;
             splice(@ovls, $o, 1);
             last;
          }
        }
      }
      push(@regs, @ovls);
      #
      @regs=sort { $a->[0] <=> $b->[0] } @regs;
      my $check=1;
      while ($check) {
       for (my $i=1;$i<@regs;$i++) {
         my ($ml, $mr)=mergeOvl($regs[$i-1], $regs[$i]);
         if ($mr) { splice(@regs, $i, 1); $regs[$i-1]=[$ml, $mr]; last; }
       }
       $check=0;
      } #made sure @regs is sorted and non-redundant
    } #have overlaps with previous read components
  } #check overlap with existing regions
  else {
   push(@regs, @ovls);
   #$refctgdata->[2]=[@regs];
  }
  #print STDERR "regs = ";printRegs(\@regs);
 }# component(read) info
}

printContig($ctg, $ctgn, $ctglen, \@regs, \@comps, \@allcomps) if $ctg;

if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#==== subroutines ====

sub printContig {
 my ($ctg, $rcount, $clen, $regs, $comps, $allcomps)=@_;
 #print STDERR "<<<<<<>>>>>>>$ctg\t$rcount\t$clen\t";
 print "$ctg\t$rcount\t$clen\t";
 printRegs($regs, \@allcomps);
 if ($comps && $ovlfilter && @$comps > 0) {
   print "\t";
   my @complst;
   foreach my $comp (@$comps) {
    push(@complst, $$comp[2].'['.$$comp[3].']:'.$$comp[0].'-'.$$comp[1]);
   }
   print join(',',@complst);
 }
 print "\n";
}

sub printRegs {
 my ($regs, $comps)=@_;
 my @rs;
 foreach my $r ( @$regs ) {
   my $el=$r->[0].'-'.$r->[1].' ('.int($$r[1]-$$r[0]+1).')';
   if ($comps && @$comps>0) {
     my @thlst; #list of components running through the merge region
     foreach my $cd (@comps) {
       if ($$cd[0] < $$r[0]-6 && $$cd[1]>$$r[1]+6) {
         push(@thlst, $$cd[2].':'.$$cd[0].'-'.$$cd[1]);
       }
     }
     if (@thlst>0) {
       $el.=' ['.join('~',@thlst).']';
     }
   }
   push(@rs, $el);
 }
 print join(',',@rs);
}

sub getOvl { # returns intersection
 my ($rstart, $rend)=@{$_[0]};
 my ($start, $end)=@{$_[1]};
 if ($start<$rstart) {
   return (0,0) if ($rstart>$end);
   return ($rend>$end) ? ($rstart, $end) : ($rstart, $rend);
   }
  else { # $rstart<=$start
   return (0,0) if ($start>$rend);
   return ($rend<$end) ? ($start, $rend) : ($start, $end);
   }
}

sub mergeOvl { # returns union
 my ($rstart, $rend)=@{$_[0]};
 my ($start, $end)=@{$_[1]};
 if ($start<$rstart) {
   return (0,0) if ($rstart>$end);
   return ($rend>$end) ? ($start, $rend) : ($start, $end);
   }
  else { # $rstart<=$start
   return (0,0) if ($start>$rend);
   return ($rend<$end) ? ($rstart, $end) : ($rstart, $rend);
   }
}
