#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 utrpaste.pl [-a 'newID'] [-o <output.gff3>] <jigsaw_out.jgff> <exon_mappings.gff3>
 
 Augment the gene models found in <jigsaw_out.jgff> with
 exon mappings found in <mappings.gff3> (only initial\/terminal exons, 
 if they overlap the CDS)
 
/;
umask 0002;
getopts('a:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $newacc=$Getopt::Std::opt_a;
my ($jgff,$mapgff)=@ARGV;
foreach (@ARGV) {
   die ("Error locating file $_\n") unless -f $_;
   }
open(MGFF, $mapgff) || die ("Error opening mappings file $mapgff\n");
#                        0     1       2    3       4
my @mapdata; #list of [ ID, strand, start, end, [ @exons] ]
             #  where each exon element is [exon_start, exon_end]
my $mi=-1; #current index in @mapdata
my $curid;
while (<MGFF>) {
  my ($chr, $v, $f, $fstart, $fend, $fscore, $strand, $frame, $info)=split(/\t/);
  next unless ($f eq 'exon'); #! ignore anything else but exons
  ($fstart, $fend)=($fend, $fstart) if ($fstart>$fend);
  my ($id)=($info=~m/Parent=([^;]+)/);
  die("Error getting mapping ID at:\n$_\n") unless $id;
  if ($id ne $curid) {
     #start a new record;
     push(@mapdata, [ $id, $strand, 0,0, []]); 
     $mi=$#mapdata;
     $curid=$id;
     }
  push( @{$mapdata[$mi]->[4]} , [$fstart, $fend] );
  }
close(MGFF);

foreach my $md (@mapdata) {
  my @exons=sort { $main::a->[0] <=> $main::b->[0] } @{$md->[4]};
  $md->[2]=$exons[0]->[0];
  $md->[3]=$exons[-1]->[1];
  $md->[4]=[@exons];
  }
  
# sort mappings by start coordinate
@mapdata=sort { $main::a->[2] <=> $main::b->[2] } @mapdata;

if ($outfile) {
  open(FOUT, '>'.$outfile) || die("Error creating file $outfile\n");
  select(FOUT);
  }
#now parse and augment the jgff
open(JGFF, $jgff) || die ("Error opening mappings file $jgff\n");
$curid='';
my @cds=0; #current CDS chain: list of [$start, $end, $frame]
my @curdata;
while (<JGFF>) {
  next if length($_)<4 || m/^# /;
  my ($chr, $track, $f, $fstart, $fend, $fscore, $strand, $phase, $jdata)=split(/\t/);  
  next unless (($f eq 'CDS') || ($track=~m/jigsaw/ && $f=~m/\-exon$/));  
  my ($id)=($jdata=~m/^([^;]+); gene_score=/);
  unless ($id) {
    ($id)=($jdata=~m/Parent=([^;]+)/);    
    unless ($id) { #try GTF too
       ($id)=($jdata=~m/transcript_id "([^"]+)/);
       }
    die("Error getting model ID at:\n$_\n") unless $id;
    }
  if ($id ne $curid) {  
    extendCDS(@curdata, \@cds) if $curid; #checks @mapdata
    @cds=();
    $curid=$id;
    my $wid = $newacc || $curid;
    my ($jnum)=($curid=~m/\.(\d+)$/);
    $jnum=1 unless $jnum;
    @curdata=($wid, $jnum, $chr, $track, $strand);
    }
   push(@cds, [$fstart, $fend, $phase]);   
  }
# last record:
extendCDS(@curdata, \@cds) if @cds>0; #checks @mapdata
close(JGFF);

if ($outfile) { close(FOUT); }

#================ SUBROUTINES ============

sub extendCDS {
 my ($id, $jnum, $chr, $track, $strand, $cdsl)=@_;
 #find mapdata overlapping first and last exons of this model
 my @cds = sort { $main::a->[0] <=> $main::b->[0] } @$cdsl;
 my $cstart=$cds[0]->[0];
 my $cend=$cds[-1]->[1];
 my $mrnaid=$id.'.j'.$jnum;
 
 my $found=0;
 my @exons;

 foreach my $m (@mapdata) {
   my ($mid, $mstrand, $mstart, $mend, $mexons)=@$m;
   #mapping should overlap on the same strand
   #print STDERR "===> comparing $mstart-$mend ($mstrand) w $cstart-$cend ($strand)\n";
   next if ($mstrand ne $strand || $mstart>$cend || $cstart>$mend);
   #mapping ends should cover a larger interval than CDS
   next if (($mstart>$cstart) || ($mend<$cend));
   # at least one of the ends should stretch farther
   next unless (($mstart<$cstart) || ($mend>$cend));
   # -- found it:
   $found=1;
   # augment here, and print the mRNA line too, then exons
   # left end :
   my @lexons;
   my $lovl;
   for (my $i=0;$i<@$mexons;$i++) {
     my ($estart, $eend)=@{$$mexons[$i]};
     #print STDERR "left end cmp: $estart-$eend vs CDS $cstart-$cds[0]->[1]\n";
     if ($estart<=$cstart && $eend>$cstart) {
        # if they really overlapped properly, $estart is <= cds start 
        #print STDERR "                         YES, adding exon $estart-$cds[0]->[1]^^!\n";
        push(@lexons, [$estart, $cds[0]->[1]]);
        $lovl=1;
        last;
        }
     push(@lexons,[$estart, $eend]);
     }
   #right end:
   my @rexons;
   my $rovl;
   for (my $i=@$mexons-1;$i>=0;$i--) {
     my ($estart, $eend)=@{$$mexons[$i]};
     #print STDERR "right end cmp: $estart-$eend vs CDS $cds[-1]->[0]-$cend\n";
     if ($eend>=$cend && $estart<$cend) {
        # if they really overlapped properly, $eend is >= cds end
        #print STDERR "                         YES, adding exon $cds[-1]->[0]-$eend!\n";
        unshift(@rexons, [$cds[-1]->[0], $eend]);
        $rovl=1;
        last;
        }
     unshift(@rexons,[$estart, $eend]);
     }
 
   if ($rovl && $lovl) {
      #special single exon case:
      if ( @cds==1 ) {
         #we have to merge rexon and lexons
         $lexons[-1]->[1]=$rexons[0]->[1];         
         shift(@rexons); #remove the first rexon
         push(@exons, @lexons, @rexons);
         }
       else { #normal multi-exon case
         push(@exons, @lexons);
         foreach my $cd (@cds) {
           push(@exons, [$$cd[0], $$cd[1]]) 
               if ($cd->[0]>$lexons[-1]->[1] && $cd->[1]<$rexons[0]->[0])
           }
         #for (my $i=1;$i<@cds;$i++) {
         #   push(@exons, [$cds[$i]->[0], $cds[$i]->[1]]);
         #   }
         push(@exons, @rexons);
         }
      # print mRNA and exons
      print join("\t", $chr, $track, 'mRNA', $exons[0]->[0], $exons[-1]->[1], '.', $strand, '.',
                    'ID='.$mrnaid)."\n";
      foreach my $e (@exons) {
        print join("\t", $chr, $track, 'exon', $$e[0], $$e[1], '.', $strand, '.',
                    'Parent='.$mrnaid)."\n";
        }
      last;
      }
    else { $found = 0; }
   }
 if (!$found) {
   #print mRNA line and the old CDS as is
   print join("\t", $chr, $track, 'mRNA', $cstart, $cend, '.', $strand, '.',
                    'ID='.$mrnaid)."\n";
   }
 #print CDS segments:
 foreach my $c (@cds) {   
     print join("\t", $chr, $track, 'CDS', $$c[0], $$c[1], '.', $strand, $$c[2],
                    'Parent='.$mrnaid)."\n";
     
     }
}
