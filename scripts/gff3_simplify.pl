#!/usr/bin/env perl
#
#converts newer,complex GFF3 file (e.g. from TAIR)
#  to the simpler format expected by Cufflinks/Tophat
#
# Usage: 
#    gff3_simplify.pl < input.gff3 > output.gff3

#http://www.sequenceontology.org/resources/gff3.html
# Unsupported features by Cufflinks/Tophat GFF3 parser:
#  * multiple parents (e.g. exons shared by multiple isoforms)
#  * gene name can only be found in the parent "gene" feature
#  * exon and CDS can be parented directly by a "gene" feature
#  * spliced non-coding transcripts the parent feature is "noncoding_transcript" 
#     with child "exon" features (this is actually supported)

#corrections made:
#  *create dummy RNA IDs for orphan CDS/exons or those owned directly by a gene feature
#  *discard gene features but keep their name in the gene_name attribute of their transcripts
#  *discard any non-transcript or redudant features (i.e. keep only ?RNA, exon and CDS features)
#  *copy shared exons/CDS to their parent transcripts as needed (to avoid multiple parents for a feature)


use strict;
my %genes; # id => name
my %rnas;
# t_id=>[gseq, src, "ftype fstart fend score strand phase", [exons], [cds], name, gene_name]
#          0    1     2                                        3       4     5      6
my @tlist;
my $tdummy; #dummy transcripts counter
while (<>) {
  chomp;
  my ($ctg, $src, $f, $fstart, $fend, $score, $strand, $phase, $attrs)=split(/\t/);
  next unless $attrs;
  my ($id)=($attrs=~m/ID=([^;]+)/);
  my ($name)=($attrs=~m/Name=([^;]+)/);
  if ($f eq 'gene') {
    my $gname=$name || $id;
    $genes{$id}=[$gname, join("\t", $fstart, $fend, $score, $strand, $phase)];
    next;
    }
  my ($parent)=($attrs=~m/Parent=([^;]+)/);
  if ($f=~m/RNA|transcript/) {
    my $gname=$genes{$parent}->[0] if exists($genes{$parent});
    $rnas{$id}=[$ctg, $src, join("\t", $f, $fstart, $fend, $score, $strand, $phase),
                             [], [], $name, $gname];
    push(@tlist, $id);
    next;
    }
  next unless $f eq 'exon' || $f eq 'CDS';
  my @parents=split(/\,/,$parent);
  my $x= ($f eq 'CDS') ? 4 : 3;
  foreach my $p (@parents) {
     my $td=$rnas{$p};
     if ($td) {
         push(@{$$td[$x]}, join("\t",$fstart, $fend, $score, $strand, $phase));
         }
       elsif (exists($genes{$p})) {
            my $grna=$ctg.'|'.$genes{$p}->[0].'|RNA';
            if (!exists($rnas{$grna})) {
                $rnas{$grna}=[$ctg, $src, $genes{$p}->[1], [], [], $genes{$p}->[0]];
                push(@tlist, $grna);
                }
            push(@{$rnas{$grna}->[$x]}, join("\t",$fstart, $fend, $score, $strand, $phase));
            }
     }
  }

foreach my $tid (@tlist) {
 my $td=$rnas{$tid};
 my $tattrs="ID=$tid";
 $tattrs.=";Name=".$$td[6].';gene_name='.$$td[6] if ($$td[6]);
 print join("\t", @$td[0..2], $tattrs)."\n";
 foreach my $xd (@{$$td[3]}) {
    print join("\t", $$td[0], $$td[1], 'exon', $xd, "Parent=$tid")."\n";
    }
 foreach my $xd (@{$$td[4]}) {
    print join("\t", $$td[0], $$td[1], 'CDS', $xd, "Parent=$tid")."\n";
    }
 }

