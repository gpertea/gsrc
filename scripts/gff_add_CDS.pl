#!/usr/bin/perl
use strict;
my $usage=q/
gff_add_CDS.pl <gff_w_CDS.gff> <gff_input..>

Adds the CDS records from a GFF file to another GFF input 
which presumably lacks it (or if the input already has CDS
records, they will be discarded).
/;
my %cdsdata; # CDSParent => [ "CDSseg1_data", "CDSseg2_data", ..]

#loads the CDS data
my $cdsfile=shift(@ARGV);
my $cdscount;
open(CDSGFF, $cdsfile) || die ("Error opening CDS file $cdsfile\n");
while (<CDSGFF>) {
  my @t=split('\t');
  next unless $t[2] eq 'CDS';
  chomp($t[8]);
  my ($p)=($t[8]=~m/\bParent=([^;]+)/);
  die("Error parsing Parent for CDS segment at $_\n") unless $p;
  $t[8]=~s/\bParent=[^;]+;?//;
  push(@{$cdsdata{$p}}, join("\t",@t[2..8]));
  $cdscount++;
}
close CDSGFF;

#print STDERR "..loaded $cdscount CDS segments..\n";

while (<>) {
 my @t=split('\t');
 print $_;
 if ($t[2] eq 'mRNA') {
   # print CDS data here, if found
   my ($id)=($t[8]=~m/\bID=([^;]+)/);
   my $cds=$cdsdata{$id};
   unless ($cds) {
      #search for the special "acc" attribute we may store there for Drupal
      my ($acc)=($t[8]=~m/\bacc=([^;]+)/);
      $cds=$cdsdata{$acc};      
      }
   if ($cds) {
      foreach my $cdseg (@$cds) {
         #$cdseg.="Parent=$id"
         # unless ($cdseg=~s/\t([^\t]+)$/\tParent=$id;$1/);
         $cdseg=~s/\t([^\t]*)$/\tParent=$id;$1/;
         print join("\t",@t[0..1],$cdseg)."\n";
         }
      }
   }
}
