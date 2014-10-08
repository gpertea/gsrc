#!/usr/bin/perl
use strict;
use Getopt::Std;

my $usage=q{
Usage:

  tran6frames.pl [-O] [-t <tag>] [-R] <seqfile>
  
  -R : the sequence is masked for repeats 
       (with lower case letters) and the translation
       will ignore (skip) such regions
  -O : use the offset information in the defline of each fasta 
       record for generating the base name of the translated
       sequences; the generated ID will have the format:
       <original_seqid>|<original_offset>|<tran_start>_<tran_end>
  -t   add the provided string <tag> to the seq ID generated for each
       record; the new ID will have the format:
       <original_seqid>|<tag>|<tran_start>_<tran_end>
};

getopts('ORt:') || die($usage."\n");

my ($useoffset, $repeat, $nametag)=
  ($Getopt::Std::opt_O, $Getopt::Std::opt_R, $Getopt::Std::opt_t);

if ($useoffset && $nametag) {
  die("$usage\n Error: options -t and -O are exclusive!\n");
  }

my %codons=(
'AAA'=>'K', 'AAC'=>'N', 'AAG'=>'K', 'AAR'=>'K', 'AAT'=>'N',
'AAY'=>'N', 'ACA'=>'T', 'ACB'=>'T', 'ACC'=>'T', 'ACD'=>'T',
'ACG'=>'T', 'ACH'=>'T', 'ACK'=>'T', 'ACM'=>'T', 'ACN'=>'T',
'ACR'=>'T', 'ACS'=>'T', 'ACT'=>'T', 'ACV'=>'T', 'ACW'=>'T',
'ACY'=>'T', 'AGA'=>'R', 'AGC'=>'S', 'AGG'=>'R', 'AGR'=>'R',
'AGT'=>'S', 'AGY'=>'S', 'ATA'=>'I', 'ATC'=>'I', 'ATG'=>'M',
'ATH'=>'I', 'ATM'=>'I', 'ATT'=>'I', 'ATW'=>'I', 'ATY'=>'I',
'CAA'=>'Q', 'CAC'=>'H', 'CAG'=>'Q', 'CAR'=>'Q', 'CAT'=>'H',
'CAY'=>'H', 'CCA'=>'P', 'CCB'=>'P', 'CCC'=>'P', 'CCD'=>'P',
'CCG'=>'P', 'CCH'=>'P', 'CCK'=>'P', 'CCM'=>'P', 'CCN'=>'P',
'CCR'=>'P', 'CCS'=>'P', 'CCT'=>'P', 'CCV'=>'P', 'CCW'=>'P',
'CCY'=>'P', 'CGA'=>'R', 'CGB'=>'R', 'CGC'=>'R', 'CGD'=>'R',
'CGG'=>'R', 'CGH'=>'R', 'CGK'=>'R', 'CGM'=>'R', 'CGN'=>'R',
'CGR'=>'R', 'CGS'=>'R', 'CGT'=>'R', 'CGV'=>'R', 'CGW'=>'R',
'CGY'=>'R', 'CTA'=>'L', 'CTB'=>'L', 'CTC'=>'L', 'CTD'=>'L',
'CTG'=>'L', 'CTH'=>'L', 'CTK'=>'L', 'CTM'=>'L', 'CTN'=>'L',
'CTR'=>'L', 'CTS'=>'L', 'CTT'=>'L', 'CTV'=>'L', 'CTW'=>'L',
'CTY'=>'L', 'GAA'=>'E', 'GAC'=>'D', 'GAG'=>'E', 'GAR'=>'E',
'GAT'=>'D', 'GAY'=>'D', 'GCA'=>'A', 'GCB'=>'A', 'GCC'=>'A',
'GCD'=>'A', 'GCG'=>'A', 'GCH'=>'A', 'GCK'=>'A', 'GCM'=>'A',
'GCN'=>'A', 'GCR'=>'A', 'GCS'=>'A', 'GCT'=>'A', 'GCV'=>'A',
'GCW'=>'A', 'GCY'=>'A', 'GGA'=>'G', 'GGB'=>'G', 'GGC'=>'G',
'GGD'=>'G', 'GGG'=>'G', 'GGH'=>'G', 'GGK'=>'G', 'GGM'=>'G',
'GGN'=>'G', 'GGR'=>'G', 'GGS'=>'G', 'GGT'=>'G', 'GGV'=>'G',
'GGW'=>'G', 'GGY'=>'G', 'GTA'=>'V', 'GTB'=>'V', 'GTC'=>'V',
'GTD'=>'V', 'GTG'=>'V', 'GTH'=>'V', 'GTK'=>'V', 'GTM'=>'V',
'GTN'=>'V', 'GTR'=>'V', 'GTS'=>'V', 'GTT'=>'V', 'GTV'=>'V',
'GTW'=>'V', 'GTY'=>'V', 'MGA'=>'R', 'MGG'=>'R', 'MGR'=>'R',
'NNN'=>'X', 'RAY'=>'B', 'SAR'=>'Z', 'TAA'=>'.', 'TAC'=>'Y',
'TAG'=>'.', 'TAR'=>'.', 'TAT'=>'Y', 'TAY'=>'Y', 'TCA'=>'S',
'TCB'=>'S', 'TCC'=>'S', 'TCD'=>'S', 'TCG'=>'S', 'TCH'=>'S',
'TCK'=>'S', 'TCM'=>'S', 'TCN'=>'S', 'TCR'=>'S', 'TCS'=>'S',
'TCT'=>'S', 'TCV'=>'S', 'TCW'=>'S', 'TCY'=>'S', 'TGA'=>'.',
'TGC'=>'C', 'TGG'=>'W', 'TGT'=>'C', 'TGY'=>'C', 'TRA'=>'.',
'TTA'=>'L', 'TTC'=>'F', 'TTG'=>'L', 'TTR'=>'L', 'TTT'=>'F',
'TTY'=>'F', 'XXX'=>'X', 'YTA'=>'L', 'YTG'=>'L', 'YTR'=>'L'
);

#open(F,$seqfile) || die("Error opening $seqfile!\n");

$/="\n>";
while(<>) {
    chomp;
    s/^>//;
    next unless length($_)>0;
    my ($name, $descr, $seq)=(m/^(\S+)[ \t\x01]*(.*?)\n(.+)/s);
    die "ERROR 21: Wrong FASTA format: .$_." unless $name && $seq;
    $seq =~ tr/ \t\n\r//d;
    my ($ofs)=($descr=~m/^(\d+)/);
    if ($nametag) {
      $name.='|'.$nametag;
      }
    elsif ($useoffset && defined($ofs)) {
      $name.='|'.$ofs;
      }
    process($name,$seq);
}

sub process {
    my ($name,$seq)=@_;

    my $len=length($seq);
    #my $index=0;
    #$index=translate($name,$seq,$len,$index,1);
    translate($name,$seq,$len,1);
    $seq= reverse $seq;
    $seq =~ tr/AaCcTtGgUuMmRrWwSsYyKkVvHhDdBb/TtGgAaCcAaKkYyWwSsRrMmBbDdHhVv/;
    #translate($name,$seq,$len,$index,-1);
    translate($name,$seq,$len);
}

sub translate {
    my ($name,$seq,$len,$fwd)=@_;
    my $i=0;
    my @start;

    $start[0]=0;
    $start[1]=1;
    $start[2]=2;

    my @prot;
    my @aa;
    my (@numaa, @numX);
    while($i+3<$len) {
     for(my $j=0;$j<3;$j++) {
        my $cod=substr($seq,$i+$j,3);
        $aa[$j]=$repeat ? $codons{$cod} : $codons{uc($cod)};
        $aa[$j]='X' unless $aa[$j];
        $numX[$j]++ if $aa[$j] eq 'X';
        if($aa[$j] eq '.') {
         if($i+$j-$start[$j]>=60 && $numaa[$j]-$numX[$j]>6) {
            print ">$name|";
            my ($cstart, $cend)=$fwd ?($start[$j]+1,    $i+$j):
                                      ($len-$start[$j], $len-$i-$j+1);
            print join('_',$cstart, $cend)."\n";
            #$index++;
            #print $prot[$j],"\n";
            &printfa($prot[$j]);
            }
         $start[$j]=$i+$j+3;
         $prot[$j]="";
         $numaa[$j]=0;
         $numX[$j]=0;
         } #stop found
         elsif($aa[$j]) {
          $prot[$j].=$aa[$j];
          $numaa[$j]++;
          } 
       } #for $j
     $i+=3;
     } #while $i
    #process the last aa translation (after the last stop found)
    for(my $j=0;$j<3;$j++) {
      my $lenp=length($prot[$j]);
      if($start[$j]<$len && $lenp>=20 && $numaa[$j]-$numX[$j]>6) {
         print ">$name|";
         #$index++;
         my ($cstart, $cend)= $fwd ? ($start[$j]+1,    $start[$j]+3*$lenp) :
                                     ($len-$start[$j], $len-$start[$j]-3*$lenp+1);
         $cend=1 if $cend<1;
         #if ($cstart==389) {
         #  print STDERR join(", ", "len=$len", "fwd=$fwd","start[j]=$start[$j] (j=$j)", "lenp=$lenp")."\n";
         #  }                               
         print join('_',$cstart, $cend)."\n";
         &printfa($prot[$j]);
         $prot[$j]="";
         }
     }
 #return $index;
}

sub printfa {
 my $seqlen=length($_[0]);
 my $pos=0;
 while ($pos<$seqlen) {
   print substr($_[0],$pos,60)."\n";
   $pos+=60;
   }
 }
