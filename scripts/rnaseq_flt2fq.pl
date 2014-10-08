#!/usr/bin/perl
use strict;
use Getopt::Std;

my $usage = q/Usage:
  rnaseq_flt2fq.pl [-Q] [-D] [-o <out.fq>] [-p <name_prefix>] <flt_file>

  Use -D to apply a minimal dust filter and not output any reads
  that are more than 50% low-complexity

/;
umask 0002;
getopts('QDo:p:') || die($usage."\n");
my $prefix=$Getopt::Std::opt_p || 'MIR0';
my $dust=$Getopt::Std::opt_D;
my $isfq=$Getopt::Std::opt_Q;
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --

my $minlen=16;
my $rc=0; #counter
my $totalreads=0;
my $discarded=0;
my %h; # $h{read}=[count, 'quals'];
# also computes qualities by averaging the values
while (<>) {
 chomp;
 next unless length($_)>0;
 my ($seq, $qv);
 if ($isfq) {
    die("formatting Error! (\@ missing)\n") unless m/^\@/;
    $seq=<>;
    chomp($seq);
    $_=<>;
    die("formatting Error! (+ missing)\n") unless m/^\+/;
    $qv=<>;
    chomp($qv);
    }
  else { #simple 2-column format: sequence, quality values
   ($seq, $qv)=split;
    }
 $totalreads++;
 die("Error: quality value length doesn't match sequence length!\n")
   unless (length($seq)==length($qv));
 my $plen=length($seq);
 if ($seq=~s/^n+//i) {
   $qv=substr($qv,$plen-length($seq));
   }
 if (length($seq)<16) {
   $discarded++;
   next;
   }
 $rc++;
 storeRead(uc($seq), $qv);
 #printfq($seq, $qv);
}


# sprintf('MIRD1%08d',$rc),
 #print join("\n",'@'.$id, uc($seq), '+',$quals)."\n";

print STDERR "$discarded reads discarded out of initial $totalreads reads\n";
my $newcount=keys(%h);
print STDERR "After collapsing duplicates, the new total is: $newcount reads\n";
my $counter=0;
my $dustremoved=0;
while (my ($seq, $sd)=each(%h)) {
 my $id=sprintf($prefix.'%08d',$counter+1);
 if ($dust) {
   my $mseq=$seq;
   $mseq =~ s/((\w+?)\2{5,})/'N' x length $1/oeg;
   my $masked=($mseq=~tr/N//); #bases masked
   if ($masked>(length($seq)>>1)) { #>50% masked
     $dustremoved++;
     next;
     }
   }
 print join("\n", '@'.$id.'_x'.$sd->[0], $seq, '+',$sd->[1])."\n";
 $counter++;
}
my $dustmsg='';
$dustmsg=" (dust filter discarded $dustremoved reads)" if ($dustremoved>0);
print STDERR "$counter fastq records written.$dustmsg\n";

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#============ Subroutines:

sub storeRead {
 my ($seq, $qv)=@_;
 my @q=unpack('C*', $qv); #get ascii values
 my $slen=length($qv);
 my $sd=$h{$seq};
 if ($sd) { # same read found before
    $sd->[0]++;
    my @oldq=unpack('C*',$sd->[1]); #old qualities
    my $i=0;
    my $quals='';
    foreach my $q0 (@oldq) {
       $quals.= chr( int(($q[$i]-31+$q0)/2) ); #average qualities
       $i++;
       }
    $sd->[1]=$quals; #update qualities
    }
  else {
    my $quals='';
    foreach my $c (@q) {
      $quals .= chr($c-31);
      }
    $h{$seq}=[1, $quals];
    }
}
