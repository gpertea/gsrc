#!/usr/bin/perl
use strict;
use Getopt::Std;
use File::Basename;
use Cwd qw(abs_path cwd);
umask 0002;

my $output;
my $myfasta;
my $evidence_file;
my $treedir;
my $idir_list;
my $addPrefix;
my $useLinComb;
my $noOverWrite;
my $penStr;
my $inLen;

my $usage = q~
 Wrapper script for running jigsaw (not the training part).

    run_jigsaw.pl [options]

 Options:
 -e <evidence_fie>   - evidence descriptor file (required)
 -l <dirlist_file>   - file with full paths to the genomic sequence and data
                       directories, one per line (required)
 -o <output_file>    - name of output file for gene predictions (required)
 -f <fasta_suffix>   - the filename suffix of the genomic sequence file
                       (default 'fa')
 -d <training_dir>   - directory with training data (required, unless -L)
 -n <intron_penalty> - specify maximum intron length and penalty for 
                       exceeding the length
 -g <grid_engine>    - grid engine to use: 'smp', 'sge' or 'condor';
                       (default is: run locally, do not use the grid)
 -c <numCPUs>        - max grid nodes/CPUs to use (default 20; requires -g)
 -m <e-mail>         - e-mail to notify when all jobs are finished
 -i <min_intron_len> - specify minimum intron length
 -C                  - do not clobber existing output file (default is to 
                       overwrite)
 -L                  - run the linear combiner instead

 Miscellaneous Options:
       -h, -help        - print this help message
       -V               - obtain program version
~;

my $cmdline="$0 ".join(' ',@ARGV);
getopts('LCm:c:g:f:e:d:l:p:n:o:') || die($usage."\n");

($evidence_file, $output, $treedir, $idir_list, $myfasta) = 
($Getopt::Std::opt_e, $Getopt::Std::opt_o, $Getopt::Std::opt_d, $Getopt::Std::opt_l, 
$Getopt::Std::opt_f);
$myfasta='fa' unless $myfasta;
$useLinComb=$Getopt::Std::opt_L;
$penStr=$Getopt::Std::opt_n;
$inLen=$Getopt::Std::opt_i;
$noOverWrite=$Getopt::Std::opt_C;
my $mailnotify=$Getopt::Std::opt_m;
my $gridengine=$Getopt::Std::opt_g;
my $maxCPUs=$Getopt::Std::opt_c || 20;


$addPrefix = 1;

# if ( defined $addPrefix ) {
#   $addPrefix = 0;
# } else {
#   $addPrefix = 1;
# }

die("$usage\n Error: not all the required parameters where given!\n") 
  unless $output && -f $idir_list && -f $evidence_file;

die("$usage\nError: training directory not given!\n") if (!-d $treedir && !$useLinComb);

print STDERR "#Command line:\n$cmdline\n";

my $penArg=$penStr ? "-n \"$penStr\"" : '';

my $linArg = $useLinComb ? '-l' : '';

my $ilArg="";

#if ( defined $penStr ) {
   #$ilArg="-i $inLen";
#}

my $evidence;
open(FILE,$evidence_file) || die("Cannot open $evidence_file: $!");
while(my $line = <FILE>) {
  chomp($line);
  ## ignore the curation line, not used when actually running combiner
  $evidence.="$line:" unless ($line =~ /\scuration\s*/ );
}
close(FILE);

my @testOn;
open(FILE,$idir_list) || die "unable to open [$idir_list]\n";
while(my $line = <FILE>) {
  chomp($line);
  push(@testOn,$line) if $line;
}
close(FILE);
my $combid = basename($evidence_file);

my $CMDQFile='.runJIGSAW.'.$combid.'.'.$$.'.gridcmds';

my @cmdqueue;
my $treearg = $treedir ? '-d '.abs_path($treedir) : '';
for(my $cnt = 0; $cnt < @testOn; $cnt++) {
  my $dir = $testOn[$cnt];
  my $dname = basename($dir);
  my $prefix;
  if( $addPrefix ) {
    $prefix = "$dir/$dname.";
  } else {
    $prefix = "$dir/";
  } 
  my $combd = "${prefix}${combid}.run";
  writeEvidenceFile($combd, $evidence, $prefix);
  my $combout = "${prefix}$output";
  #my $cmd=;
  push(@cmdqueue, "jigsaw -f ${prefix}$myfasta $treearg $linArg -e $combd -m $combout $penArg");
     
     #print "[$cmd]\n";
     #if (!defined($noOverWrite) || ! -e $combout || -z $combout )  {
     #  runCmd($cmd);
     #  }
}

close(RUNQ);
$mailnotify = $mailnotify ? "-m '$mailnotify'": '';
if ($gridengine) {
 open(RUNQ, '>'.$CMDQFile) || die("Error creating $CMDQFile!\n");
 foreach my $l (@cmdqueue) {
  print RUNQ $l."\n";
  }
 close(RUNQ);
 &runCmd("gridx -g $gridengine -p $maxCPUs -S -f $CMDQFile $mailnotify -O logs_jigs_run -q");
 }
else { # run locally:
 foreach my $cmd (@cmdqueue) {
  &runCmd($cmd);
  }
 }

exit(0);

sub writeEvidenceFile {
  my ($tfname,$evlst,$prefix) = @_;
  my @evlst = split(/:/,$evlst);
  my $evdata = $tfname.'.run';
  my $status = 0;
  open(EVDATA, '>'.$tfname) || die ("Cannot create evidence file: $tfname\n");
  for(my $iter = 0; $iter <= $#evlst; $iter++) {
       print EVDATA "$prefix$evlst[$iter]\n";
     }
  close(EVDATA);
}

sub runCmd {
 my $cmd = $_[0];
 print STDERR "#>running: $cmd\n";
 system($cmd);
}
