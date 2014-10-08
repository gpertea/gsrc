#!/usr/bin/perl
#this sorts the concatenated deflines after nrdb -d1
#such that "reviewed" deflines come first 'UPr|'
use strict;

#-------- local array holdin uniformative annotation patterns:
my @uninformative=(
'\bunknown\b',
#'unknown protein\b',
'\bhypothetical\b',
'unnamed protein product',
'open reading frame',
'\borf\b',
'\bputative\b',
'\bhomologue\b',
'\bsimilar to',
'^expressed sequence \S+$',
'\bHA\d{4}\b',
'\bDKFZP\S+\b',
'PROTEIN FOR MGC:\d+',
'PROTEIN FOR IMAGE:\d+',
'\bR\d{5}\_\d\b',
'\bPRO\d{4}\b',
# 'KIAA\d+ GENE PRODUCT',
#'KIAA\d+ PROTEIN',
'\bKIAA\d+\b',
'\bHSPC\d+\b',
# HSPC\d+ PROTEIN
#'\bC\d+ORF\d+\b',
'FLJ\d+ PROTEIN',
'\bDJ\d+[A-Z]\d+(\.\d+)*',
'NOVEL PROTEIN',
'CG\d+ PROTEIN',
'CG\d+ GENE PRODUCT',
'^\s*CG\d+\s*$',
'CGI\-\d+ PROTEIN',
'CGI\-\d+',
'CDNA:? FLJ\d+ FIS, CLONE \w+',
'BA\d+[A-Z]\d+[A-Z]?\.\d(\.\d)?',
#'\bRIKEN CDNA .{10} GENE\b',
'\bRIKEN.+?CDNA\b',
'MRNA, COMPLETE CDS, CLONE:\d+(\+\d[A-Z])?\-\d+',
'MRNA, COMPLETE CDS, CLONE:SMAP\d+\-\w+',
'BRAIN CDNA, CLONE MNCB-\d+',
'.{10}RIK PROTEIN',
'^MY\d{3}\s*$',
'MY\d{3} PROTEIN^',
'^probable\b',
'BRAIN MY\d{3}$',
'NPD\d{3} PROTEIN',
'[A-Z]\d{2}[A-Z0-9]+\.\d+ PROTEIN',
'WUGSC:H_\w+\.\w+ PROTEIN',
#'DNA SEGMENT, CHR [0-9XY]+, WAYNE STATE UNIVERSITY \d+, EXPRESSED',
#'DNA SEGMENT, CHR [0-9XY]+, KL MOHLKE \d+',
#'DNA SEGMENT, CHR [0-9XY]+, BAYLOR \d+',
'\bDNA SEGMENT\b',
'PROTEIN HSPC\d+',
#'HYPOTHETICAL [\.\d]+\s*KDA PROTEIN \S+ IN CHROMOSOME \S+',
'EG:[0-9A-Z\.]+ PROTEIN',
'GENOMIC DNA, CHROMOSOME \d+, P1 CLONE:\S+',
'[^,]+, RIKEN FULL-LENGTH ENRICHED LIBRARY, CLONE:.{10}, FULL INSERT SEQUENCE',
'ZK\d+\.\d+ PROTEIN',
'\bEST \w+',
'B2 ELEMENT'
);

$/="\n>";
my $excluded=0;
my $written=0;
while (<>) {
 chomp;
 s/^>//;
 my ($defline, $seq)=(m/^([^\n]+)(.+)/s);
 my $fseq=$seq;
 $seq=~tr/\n\r\t //d;
 if (length($seq)<10) {
    $excluded++;
    next;
    }
 #print STDERR "defline: $defline\nseq: $seq\n";
 #next;
 
 ## temp fix for my silly mistake:
 # $defline=~s/UPr\|/UX\|/sg;
 # $defline=~s/UP\|/UPr\|/sg;
 # $defline=~s/UX\|/UP\|/sg;
 ##
 
 if ($defline=~m/\x01/) {
  my @d=split(/\x01/, $defline);
  @d=sort cmpDeflines @d;
  my $f=shift(@d);
  foreach (@d) { s/ UniRef\w+// };
  $defline=join("\x01",$f,@d);
  }
 print '>'.$defline.$fseq."\n";
 $written++;
}
my $wtotal=$excluded+$written;
print STDERR "Done. $wtotal total entries processed ($excluded excluded for being too short)\n";

#===============================================
# bool isInformative($description) 
# expects only the descripts - not the accession
#===============================================
sub isInformative {
 local $_=$_[0];
 s/^\s+//g;s/\s+$//g;
 return 0 if length($_)<2;
 foreach my $pat (@uninformative) {
   if (m/$pat/i) {
     return 0;
     }
   }
return 1;
}


sub cmpDeflines {
 my $a=$main::a;
 my $b=$main::b;
 my $va=2 if ($a=~m/^UPr\|/ || $a=~m/\bNP_\d+/);
 $va+=isInformative($a);
 my $vb=2 if ($b=~m/^UPr\|/ || $b=~m/\bNP_\d+/);
 $vb+=isInformative($b);
 return ($vb <=> $va);
}
