#!/usr/bin/perl
use strict;
use Getopt::Std;
#use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 btab2cov.pl [-T] [-p <min_ps>] [-c <min_cov>] [-n <max_hits>] <input_btab>
 
 WARNING: it assumes that the input btatb is sorted by <qry_name>!
 
 It will consider all overlaps of at least <min_ps> percent similarity
 (default: 70) to assess the coverage of the query sequences
 
 Outputs coverage information like this, for any query sequence with 
 an overall percentual coverage of at least <min_cov>% (default 50):
 
 <qry_seq> <strand> <qry_cov%> <top5hits>
 
 ..where:
 <qry_cov%> is the percentual coverage of <qry_seq> length;
 <top5hits> is a comma delimited list of max. 5 subj. sequences with the best 
 overall alignment scores
 
/;
#umask 0002;
getopts('Tn:p:c:') || die($usage."\n");
my $mincov=$Getopt::Std::opt_c || 50;
my $minpsim=$Getopt::Std::opt_p || 70;
my $tophits=$Getopt::Std::opt_n;
my $mgtab=$Getopt::Std::opt_T;
my @qhp; #HSPs on positive strand: list of [ql, qr, sname]
my @qhm; #HSPs on negative strand:list of  [ql, qr, sname]
my @hitsp; # list of [score, sname] for positive strand;
my @hitsm; # list of [score, sname] for negative strand;
my %topm; # names of top n hits on plus strand
my %topp; # names of top n hits on minus strand
my ($curquery, $curlen);
if ($mgtab) { #mgblast TAB format
while (<>) {
 next if m/^\s*#/;
 chomp;
 #       0    1      2       3     4      5       6    7     8      9
 my ($qname, $qlen, $q_5,   $q_3, $sname, $slen, $s_5, $s_3, $pid, $score,
 #    10     11       12          13
    $e_val, $strand, $gapdataq, $gapdatas)=split(/\t/);
 if ($curquery ne $qname) {
   flushQ() if $curquery;
   $curquery=$qname;
   $curlen=$qlen;
   @qhp=();
   @qhm=();
   @hitsp=();
   @hitsm=();
   }
 next unless $pid>=$minpsim;
 if ($strand eq '-') {
   #minus strand:
   push(@hitsm, [$score, $sname]);   
   ($q_5, $q_3)=($q_3,$q_5) if $q_5>$q_3;
   push(@qhm, [$q_5, $q_3, $sname]);
   }
 else {
   #plus strand:
   push(@hitsp, [$score, $sname]);
   ($q_5, $q_3)=($q_3,$q_5) if $q_5>$q_3;
   push(@qhp, [$q_5, $q_3, $sname]);
   }
 } #while readline 
} #mgblast TAB format

else { # BTAB format

while (<>) {
 next if m/^\s*#/;
 chomp;
 #       0      1      2       3     4     5       6    7     8      9
 my ($qname, $date, $qlen, $method, $db, $sname, $q_5, $q_3, $s_5, $s_3,
 #    10    11     12      13      14    15        16       17       18
    $pid, $psim, $score, $fofsb, $fofse, $sdescr, $frame, $strand, $slen, 
 #    19    20      21
    $e_val, $scov, $hsps)=split(/\t/);
 if ($curquery ne $qname) {
   flushQ() if $curquery;
   $curquery=$qname;
   $curlen=$qlen;
   @qhp=();
   @qhm=();
   @hitsp=();
   @hitsm=();
   }
 my @hsp=split(/\~/,$hsps);
 if ($strand eq '-' || lc($strand) eq 'minus') {
   #minus strand:
   push(@hitsm, [$score, $sname]); 
   foreach my $h (@hsp) {
     my ($q5,$q3, $h5, $h3, $p)=($h=~/(\d+)\-(\d+)\:(\d+)\-(\d+)\|([\d\.]+)/);
     next unless $p>=$minpsim;
     die("Error parsing segment data $h for btab line:\n$_\n") unless $p>1;
     ($q5, $q3)=($q3,$q5) if $q5>$q3;
     push(@qhm, [$q5, $q3, $sname]);
     }
   }
 else {
   #plus strand:
   push(@hitsp, [$score, $sname]);
   foreach my $h (@hsp) {
     my ($q5,$q3, $h5, $h3, $p)=($h=~/(\d+)\-(\d+)\:(\d+)\-(\d+)\|([\d\.]+)/);
     next unless $p>=$minpsim;
     die("Error parsing segment data $h for btab line:\n$_\n") unless $p>1;
     #($q5, $q3)=($q3,$q5) if $q5>$q3;
     push(@qhp, [$q5, $q3, $sname]);
     }
   }
 }
} #BTAB case

flushQ() if $curquery;

sub flushQ {
 #compute coverage per strand
 my ($cp, $cm)=(0,0);
 my $rhits; #points to either @hitsp or @hitsm
 undef %topp;
 undef %topm;
 if (@hitsp) {
   @hitsp = sort { $main::b->[0]<=>$main::a->[0] } @hitsp; #sort hits by score
   if ($tophits) {
     foreach my $hd (@hitsp) {
        $topp{$hd->[1]}=1;
        last if keys(%topp)>=$tophits;
        }
     }
   }
 if (@hitsm) {
   @hitsm = sort { $main::b->[0]<=>$main::a->[0] } @hitsm;
    if ($tophits) {
      foreach my $hd (@hitsm) {
        $topm{$hd->[1]}=1;
        last if keys(%topm)>=$tophits;
        }
    }
   }
 if (@qhp) {
   if ($tophits) { @qhp = grep { exists($topp{$_->[2]}) } @qhp };
   #compute coverage on positive strand
   $cp=getCov($curlen, \@qhp);
   }
 if (@qhm) {
  if ($tophits) { @qhm = grep { exists($topm{$_->[2]}) } @qhm };
   #compute coverage on negative strand
   $cm=getCov($curlen, \@qhm);
   }
 if ($cm>=$mincov) {
   my @h;
   foreach my $p (@hitsm) { if ($tophits) { 
                                 push(@h, $$p[1]) if exists($topm{$$p[1]});
                             } else { push(@h, $$p[1]); }
                            last if @h>=5; 
                            }
   print join("\t",$curquery, '-', $cm, join(',',@h))."\n";
   }
 if ($cp>=$mincov) {
   my @h;
   foreach my $p (@hitsp) { 
                      if ($tophits) {
                        push(@h, $$p[1]) if exists($topp{$$p[1]});
                        }
                       else { 
                        push(@h, $$p[1]); 
                        }
                      last if @h>=5; 
                      }
   print join("\t",$curquery, '+', $cp, join(',',@h))."\n";
   }
}


sub getCov {
 my ($len, $r)=@_;
 my @a=sort { $main::a->[0]<=>$main::b->[0] } @$r;
 my $wasovl=1;
 my @m = map { [$_->[0], $_->[1]] } @a;
 WASOVL: while ($wasovl) {
  $wasovl=0;
  for (my $i=0;$i<@m-1;$i++) {
     for (my $j=$i+1;$j<@m;$j++) {
      next if $i==$j;
      my ($l1, $r1)=($m[$i]->[0], $m[$i]->[1]);
      my ($l2, $r2)=($m[$j]->[0], $m[$j]->[1]);
      #we know that l2>=l1 for sure
      if ($l2<=$r1) {
        #intersection detected;
        $wasovl=1;
        $m[$i]->[1] = ($r1>$r2)?$r1:$r2; #max right coord
        splice(@m, $j);
        next WASOVL;
        } #intersection
     } #j
   } #i
 }#while wasovl
 my $cov=0;
 map { $cov+=$_->[1]-$_->[0]+1 } @m;
 return sprintf('%d',($cov*100.00)/$len);
}
