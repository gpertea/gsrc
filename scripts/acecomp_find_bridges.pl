#!/usr/bin/perl
use strict;

#input is the components file from ace2fasta.pl
#output is contig name followed by a list of overlap regions in contig coordinates

#my @ctglst;
#my %ctgdata;
use Getopt::Std;
my $usage = q/Usage:
acecomp_find_bridges.pl [-o <outfile>] [-F] < ace_file_comp
 Takes the component file produced by ace2fasta.pl 
 and annotates it showing overlaps between components of interest
 and signaling bridges between them, if any.
 
 Components of interest are those having names matching a 
 hard coded filter in this script.

Options:
 -o <outfile> : write output to <outfile>
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
#my $ovlfilter=$Getopt::Std::opt_F;

#current contig data:
my ($ctg, $ctgn, $ctglen);
#my @comps; # list of component contributions [start, end, readname, strand]
my @allcomps; # unfiltered list of components (always available)
#my @regs;  #regions of overlapping CoIs (Components of Interest)
my @clusters; #clusters of CoIs within a contig; list of [$start, $end, \@comps, \@regs]
while (<>) {
 if (m/^>(\S+)\s+(\d+)\s+(\d+)/) {
  my ($c, $cn, $cl)=($1,$2,$3);
  printContig($ctg, $ctgn, $ctglen, \@clusters, \@allcomps) if $ctg;
  @clusters=();
  @allcomps=();
  #@regs=();
  ($ctg, $ctgn, $ctglen)=($c, $cn, $cl);
  #push(@ctglst, $ctg);
  #$refctgdata=[$ctgn, $ctglen, [] ];
  #$ctgdata{$ctg}=$refctgdata;
  next;
 }
 if (m/^(\S+)\s+(\d+)\s+([\+\-])\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
  my ($read, $len, $strand, $rstart, $rend, $cstart, $cend)=($1, $2, $3, $4, $5, $6, $7);
  push(@allcomps, [$cstart, $cend, $read, $strand]);
  my $isCoI=($read=~m/^3ds_anchor/ || $read=~m/^CONTIG\:\d+/);
  #print STDERR ']'.join(" ",($read, $len, $strand, $rstart, $rend, $cstart, $cend))."\n";
  next unless $isCoI;
  my @ovlOldClusters; #how many old @clusters are being overlapped by this comp
  if (@clusters>0) {
      #find all overlaps of this read with other reads
      my $clidx=0;
      foreach my $cluster (@clusters) {
        $clidx++;
        last if $$cluster[0]>$cend;
        my ($covl, $covr)=getOvl([$cstart, $cend], $cluster);
        next unless $covr;
        #overlap existing cluster
        push(@ovlOldClusters, $clidx-1);
        my $comps=$cluster->[2];
        my ($newl, $newr)=mergeOvl([$cstart, $cend], $cluster);
        if ($newr>0) {
           #extend cluster to include this comp
           ($$cluster[0], $$cluster[1]) =  ($newl, $newr);
           push(@{$$cluster[2]}, [$cstart, $cend, $read, $strand]) if @ovlOldClusters==1;
           # no need to add the seq info to more than the first cluster because of the merge below
        }
        else {
          die("Error (design bug!): mergeOvl($cstart, $cend | $$cluster[0], $$cluster[1]) failed after getOvl()\n"); 
        }
        my $ovls=$cluster->[3];
        foreach my $cspan (@$comps) {
          my ($ol, $or)=getOvl([$cstart, $cend], $cspan);
          if ($or) { #overlap detected
            #print STDERR "\tovl $ol-$or found\n";
            my ($ml, $mr);
            if (@$ovls>0) {
              #chances are this overlaps the last ovl found
              ($ml,$mr)=mergeOvl([$ol, $or], $$ovls[-1]);
            }
            if ($mr) { 
              #print STDERR "\textend prev ovl $ovls[-1]->[0]-$ovls[-1]->[1] to $ml-$mr\n";
              $$ovls[-1]=[$ml, $mr]; 
              }
            else { push(@$ovls, [$ol, $or]); }; #new overlap
          } #overlap detected
        } #for each component read
      } #for each cluster
      #merge all found overlaps in each cluster
  } #check overlap with existing components
  if (@ovlOldClusters==0) { #new cluster from current comp
      push(@clusters, [$cstart, $cend, [[$cstart, $cend, $read, $strand]], [] ]) ;
  }
  elsif (@ovlOldClusters>1) {
      #merge all clusters overlapping the new comp
      my $mcl=$clusters[shift(@ovlOldClusters)];
      foreach my $clidx (@ovlOldClusters) {
        my ($ml, $mr)=mergeOvl($mcl, $clusters[$clidx]);
        die("Error (design bug): cannot merge clusters sharing $read ($strand) $cstart-$cend\n")
          if ($mr <=0 );
       ($$mcl[0], $$mcl[1])=($ml, $mr);
       push(@{$$mcl[2]}, @{$clusters[$clidx]->[2]});
       push(@{$$mcl[3]}, @{$clusters[$clidx]->[3]}); #append overlaps
      }
      #another pass to just delete the merged clusters
      my $adj=0;
      foreach my $clidx (@ovlOldClusters) {
        splice(@clusters, $clidx-$adj, 1);
        $adj++;
      }
  }
  if (@clusters>1) {
      #make sure clusters are sorted by coordinate
      @clusters = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @clusters;
  }
  #merge overlap regions in each cluster
  foreach my $cl (@clusters) {
    if (@{$cl->[3]}>0) {
     my @ovls=sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$cl->[3]};
     my $check=1;
     while ($check) {
      for (my $i=1;$i<@ovls;$i++) {
        my ($ml, $mr)=mergeOvl($ovls[$i-1], $ovls[$i]);
        if ($mr) { splice(@ovls, $i, 1); $ovls[$i-1]=[$ml, $mr]; last; }
      }
      $check=0;
     }
     $cl->[3] = [@ovls];
   }
  } #making sure @ovls is sorted and merged in each cluster

#   my $merged=0;
#   if (@regs>0) {
#     if (@ovls>0) {
#       foreach my $reg (@regs) {
#         for (my $o=0;$o<@ovls;$o++) {
#           my ($ml, $mr)=mergeOvl($reg, $ovls[$o]);
#           if ($mr) {
#              ($reg->[0], $reg->[1])=($ml,$mr);
#              $merged=1;
#              splice(@ovls, $o, 1);
#              last;
#           }
#         }
#       }
#       push(@regs, @ovls);
#       #
#       @regs=sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @regs;
#       my $check=1;
#       while ($check) {
#        for (my $i=1;$i<@regs;$i++) {
#          my ($ml, $mr)=mergeOvl($regs[$i-1], $regs[$i]);
#          if ($mr) { splice(@regs, $i, 1); $regs[$i-1]=[$ml, $mr]; last; }
#        }
#        $check=0;
#       } #made sure @regs is sorted and non-redundant
#     } #have overlaps with previous read components
#   } #check overlap with existing regions
#   else {
#    push(@regs, @ovls);
#   }
  #print STDERR "regs = ";printRegs(\@regs);
 }# component(read) info
}

printContig($ctg, $ctgn, $ctglen, \@clusters, \@allcomps) if $ctg;

if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#==== subroutines ====

sub printContig {
 my ($ctg, $rcount, $clen, $clusters, $allcomps)=@_;
 #print STDERR "<<<<<<>>>>>>>$ctg\t$rcount\t$clen\t";
 if (@$clusters==0) {
   print STDERR "Warning: mcontig $ctg has no components of interest?!\n";
   return;
   }
 my $bridge=(@$clusters>1)? "\tlink" : '';
 if ($bridge) {
  my $adj=0;
  for (my $i=1;$i<@$clusters;$i++) {
   #TODO: check for adjacent numbers
    foreach my $ra (@{$$clusters[$i-1]->[2]}) {
       my $na;
       if ($$ra[2]=~m/3ds_anchor[ed]*_(\d+)$/) {
         $na=$1;
       }
       else { next; }
       foreach my $rb (@{$$clusters[$i]->[2]}) {
         my $nb;
         if ($$rb[2]=~m/3ds_anchor[ed]*_(\d+)$/) {
           $nb=$1;
           if (abs(int($na)-int($nb))==1) {
             $adj=1;
             last;
           }
         }
         else { next; }
       }
    last if $adj;
    }
  last if $adj;
  }
  $bridge="\tgapfill" if $adj;
 }
 print ">$ctg\t$rcount\t$clen$bridge\n";
 my $i=1;
 foreach my $cl (@$clusters) {
   print "CL$i($$cl[0]-$$cl[1])\t";
   $i++;
   my $regs=$cl->[3];
   my $comps=$cl->[2];
   printRegs($regs, \@allcomps);
   my $ovlfilter=1;
    if ($comps && $ovlfilter && @$comps > 0) {
     print "\t";
     my @complst;
     foreach my $comp (@$comps) {
      push(@complst, $$comp[2].'['.$$comp[3].']:'.$$comp[0].'-'.$$comp[1]);
     }
     print join(',',@complst);
   }
   print "\n";
 } #for each cluster
}

sub printRegs {
 my ($regs, $comps)=@_;
 my @rs;
 foreach my $r ( @$regs ) {
   my $el=$r->[0].'-'.$r->[1].'('.int($$r[1]-$$r[0]+1).')';
   if ($comps && @$comps>0) {
     my @thlst; #list of non-CoI components running through the merge region
     foreach my $cd (@$comps) {
       next if $$cd[2]=~m/^3ds_anchor/ || $$cd[2]=~m/^CONTIG\:\d+/;
       if ($$cd[0] < $$r[0]-6 && $$cd[1]>$$r[1]+6) {
         push(@thlst, $$cd[2].':'.$$cd[0].'-'.$$cd[1]);
       }
     }
     if (@thlst>0) {
       $el.='<'.join(';',@thlst).'>';
     }
   }
   push(@rs, $el);
 }
 print @rs>0 ? join(',',@rs) : 'no_merge';
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
