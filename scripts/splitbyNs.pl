#!/usr/bin/perl
use Getopt::Std;
use strict;
my $usage=q{
 splitbyNs.pl [-n <minNs>] [-M] [-G] [-o <output_name>] <input_fasta>
 
 Looks for regions of at least <minNs> Ns (default 10) and splits 
 the input sequence there.
Options:
 -M  create separate files for every split sequence
 -G  show gap info during splitting
};
getopts('MGo:n:') || die($usage."\n");
my $outputfile=$Getopt::Std::opt_o;
my $multiple=$Getopt::Std::opt_M;
my $showGaps=$Getopt::Std::opt_G;
my $minNs=$Getopt::Std::opt_n || 10;
my $Ns='N' x $minNs;

my $outfile;
my $fastafile = shift @ARGV || die($usage."\n");

if ($outputfile) {
  if (!$multiple) {
     $outfile=$outputfile;
     open(OUTF, '>'.$outputfile) || 
      die("Error creating output file $outputfile\n");
     select(OUTF);
     }
}
else {
  $outputfile='stdout' if $multiple;
}

my $n=1;
$/="\n>";

open(F,$fastafile);

while(<F>) {
    s/^>//;
    chomp;
    next unless $_;
    my ($name, $descr, $seq)=(m/^(\S+)[ \t\x01]*(.*?)\n(.+)/s);
    die "ERROR 21: Wrong FASTA format: .$_." unless $name && $seq;
    $seq=~tr/\t \n\r//d;
    #$seq=lc($seq);
    $seq=~tr/n/N/;
    $n=procScaff($name,\$seq, $n);
}

close(F);

if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

sub procScaff {
	my ($name, $seq, $n)=@_;
	my $laststart=0;
	my $count=1;
	my $pos;
	while(($pos=index($$seq,$Ns,$laststart))>0) {
	  my $newseq=substr($$seq,$laststart,$pos-$laststart);
	  writeCtg($name, \$newseq, $laststart, $n, $count);
	  $n++;
	  $count++;
	  $laststart=skipns($seq,$pos);
	  if ($showGaps) {
	    printf STDERR 'gap:%d-%d'."\t".'%d'."\n", $pos+1, $laststart, $laststart-$pos;
	  }
	}
	#write last line;
	my $newseq=substr($$seq,$laststart);
	writeCtg($name, \$newseq, $laststart, $n, $count);
	return $n;
}

sub writeCtg {
  my ($name, $rnewseq, $laststart, $n, $count)=@_;
  my $slen=length($$rnewseq);
  my ($ctgstart, $ctgend)=($laststart+1, $laststart+$slen);
  my @seqlines=unpack("A72" x (int(($slen-1)/72)+1),$$rnewseq);
  if ($multiple) {
      my $newout=$outputfile.".".sprintf('%06d',$n);
      open(O,">$newout");
      print O ">$name $ctgstart-$ctgend\n";
      print O join("\n",@seqlines)."\n";
      #"$$rnewseq\n";
      close(O);
  }
  else {
      print ">$name".'_'.sprintf('%04d',$count)." $ctgstart-$ctgend\n";
      # $$rnewseq\n";
      print join("\n",@seqlines)."\n";
  }
}

sub fastafmt {
 my $s=shift;
 my $slen=length($$s);
 my @lines=unpack("A72" x (int(($slen-1)/72)+1),$$s);
 return join("\n",@lines)."\n";
 }


sub skipns {
   my($seq, $pos)=@_;
   while(substr($$seq,$pos,1) eq 'N') { $pos++;}
   return($pos);
}
