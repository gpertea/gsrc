#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 rnaseq_cufflinks.pl <libprefix> 
 Prepares and runs tophat+cufflinks on a pair of fastq files
 e.g. 
 rnaseq_cufflinks.pl RIVA 
 
/;

my %libs=('BJAB'=>['GCB',250],
'Dogum'=>['GCB', 250],
'DOHH2'=>['GCB', 250],
'Gumbus'=>['GCB',250],
'HT'=>['GCB',250],
'Karpas422'=>['GCB',250],
'Ly1'=>['GCB',150],
'Ly18'=>['GCB',250],
'Ly19'=>['GCB',250],
'Ly2'=>['GCB',250],
'Ly4'=>['GCB',250],
'Ly7'=>['GCB',250],
'Ly8'=>['GCB',150],
'Mieu'=>['GCB',250],
'Pfieffer'=>['GCB',250],
'ULA'=>['GCB',250],
# ABC:
'HBL1'=>['ABC',150],
'TMD8'=>['ABC',150],
'Ly10'=>['ABC',150],
'Ly3'=>['ABC',250],
'RIVA'=>['ABC',250],
'SUDHL2'=>['ABC',150],
'U2932'=>['ABC',150]
);

my $tcmd='tophat --solexa1.3-quals --coverage-search --microexon-search --closure-search -p 8 ';
#my $tophat1=$tcmd.' -r 36 -o RIVA_tophat_untrimmed hg19 RIVA_1.fq RIVA_2.fq > & log_tophat_untrimmed.log &';
#my $tophat2=$tcmd.' -r 80 -o RIVA_cut84_tophat hg19 RIVA_1.unmapped_cut84.fq RIVA_2.unmapped_cut84.fq > & log_tophat_cut84.log &';
my $gzdir='/fs/szattic-asmg6/NCI_rnaseq/Cell_line_seq';
my $workdir='/fs/szannotation/human_rnaseq/NCI';
umask 0002;
chdir($workdir) || die("Error switching to working directory: $workdir\n");

getopts('o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --

my $fbase=$ARGV[0];
my $libd=$libs{$fbase};
die("$usage Error: unrecognized sample ($fbase)!\n") unless $libd;
my ($gz1, $gz2)=($fbase.'_1_sequence.txt.gz',$fbase.'_2_sequence.txt.gz');
foreach my $gz ($gz1, $gz2) {
 die("Error: file $gzdir/$gz not found!\n") unless -f "$gzdir/$gz";
}
my ($fq1, $fq2)=($fbase.'_1.fq', $fbase.'_2.fq');
my ($fqu1, $fqu2)=($fbase.'_unmapped_1.fq', $fbase.'_unmapped_2.fq');
my ($fqustats1, $fqustats2)=($fbase.'_unmapped_1.qstats', $fbase.'_unmapped_2.qstats');
my ($fqcut1, $fqcut2)=($fbase.'_trimmed_1.fq', $fbase.'_trimmed_2.fq');
if (-f $fq1) {
  mlog("$fq1 already exists, skipping decompression.");
  }
 else {
  mlog("decompressing $gz1 ..");
  system("gzip -cd $gzdir/$gz1 > $fq1");
  }
unless (-f "$fq1.cidx") {
  mlog("indexing $fq1..");
  system("cdbfasta -Q $fq1");
  system("cdbyank -l $fq1.cidx | sort -u > $fq1.lst");
  }
# fq2
if (-f $fq2) {
  mlog("$fq2 already exists, skipping decompression.");
  }
 else {
  mlog("decompressing $gz2 ..");
  system("gzip -cd $gzdir/$gz2 > $fq2");
  }
unless (-f "$fq2.cidx") {
  mlog("indexing $fq2..");
  system("cdbfasta -Q $fq2");
  system("cdbyank -l $fq2.cidx | sort -u > $fq2.lst");
  }
unless (-s "$fbase.lst" ) {
  system("sort -o $fbase.lst -m $fq1.lst $fq2.lst");
  }
my $rlen=int(`head -2 $fq1 | tail -1 | wc -L`);
my $frgsize=$libd->[1]; #
my $innerdist=$frgsize-$rlen*2;
mlog("read length: $rlen, fragment len: $frgsize (inner dist=$innerdist)");
#die("Just for now, exit\n");
my $tophat_dir=$fbase.'_tophat_untrimmed';
if (-d $tophat_dir) {
  mlog("dir $tophat_dir exists, skipping first tophat run.");
  }
 else {
  mkdir($tophat_dir) || die("Error creating directory $tophat_dir ($!)\n");
  my $tophatcmd="$tcmd -r $innerdist -o $tophat_dir hg19 $fq1 $fq2 >& $tophat_dir/tophat.log";
  mlog('Starting tophat #1 run on the untrimmed data:', ' '.$tophatcmd);
  system($tophatcmd);
  mlog('tophat run #1 done.');
  }
unless (-s "$tophat_dir/mapped_reads.lst") {
  mlog("generating the list of mapped reads ($tophat_dir/mapped_reads.lst");
  sysrun("samtools view $tophat_dir/accepted_hits.bam | sam2reads.pl | sort -u > $tophat_dir/mapped_reads.lst");
  }
# - generate the fq file with the unmapped reads
unless (-f $fqu1) {
 mlog("extracting unmapped reads from $fq1"); 
 sysrun("comm -13 $tophat_dir/mapped_reads.lst $fq1.lst | cdbyank $fq1.cidx > $fqu1");
}
unless (-f $fqu2) {
 mlog("extracting unmapped reads from $fq2"); 
 sysrun("comm -13 $tophat_dir/mapped_reads.lst $fq2.lst | cdbyank $fq2.cidx > $fqu2");
}
#analyse the qv distribution for $fqu1 and $fqu2, find the cutoff point
unless (-f $fqustats1) {
 mlog("getting quality stats for $fqu1"); 
 sysrun("fastx_quality_stats -i $fqu1 -o $fqustats1");
}

unless (-f $fqustats2) {
 mlog("getting quality stats for $fqu2"); 
 sysrun("fastx_quality_stats -i $fqu2 -o $fqustats2");
}

my ($q1cut, $q2cut);
open(FQ, $fqustats1);
while (<FQ>) {
 my @t=split();
 next unless $t[0]>5;
 $q1cut=$t[0];
 last if $t[5]<10
 }
close(FQ);

open(FQ, $fqustats2);
while (<FQ>) {
 my @t=split();
 next unless $t[0]>5;
 $q2cut=$t[0];
 last if $t[5]<10
 }
close(FQ);
  
my $cutlen= ($q1cut<$q2cut ? $q1cut : $q2cut);
die("Error getting cut length ($q1cut, $q2cut)!\n") if $cutlen<20;
unless ( -f $fqcut1 ) {
 mlog("trimming $fqu1 to length $cutlen (into $fqcut1)"); 
 sysrun("fastx_trimmer -f 1 -l $cutlen -i $fqu1 -o $fqcut1");
 }
unless ( -f $fqcut2 ) {
 mlog("trimming $fqu2 to length $cutlen (into $fqcut2)"); 
 sysrun("fastx_trimmer -f 1 -l $cutlen -i $fqu2 -o $fqcut2");
 }

my $tophat_dir2=$fbase.'_tophat_cut';
my $cutinnerdist=$frgsize-$cutlen*2;
mlog("cut read length: $cutlen, fragment len: $frgsize (inner dist=$cutinnerdist)");

if (-d $tophat_dir2) {
  mlog("dir $tophat_dir2 exists, skipping 2nd tophat run.");
  }
 else {
  mkdir($tophat_dir2) || die("Error creating directory $tophat_dir ($!)\n");
  my $tophatcmd="$tcmd -r $cutinnerdist -o $tophat_dir2 hg19 $fqcut1 $fqcut2 >& $tophat_dir2/tophat.log";
  mlog('Starting tophat #2 run on the trimmed data:', ' '.$tophatcmd);
  system($tophatcmd);
  mlog("tophat run #2 done (exit code $?)");
  }
#join resulting bam files for cufflinks
chdir($tophat_dir) || die("Error at chdir into $tophat_dir!\n");
unless (-f 'merged_hits.bam') {
 mlog("merging mappings into $tophat_dir/merged_hits.bam");
 sysrun("samtools merge merged_hits.bam accepted_hits.bam ".
        "../$tophat_dir2/accepted_hits.bam");
}
unless (-d 'cufflinks_out') {
 my $idist=int(($innerdist+$cutinnerdist)/2);
 #my $cufflinks_cmd="cufflinks.new -m $idist -s 34 -p 8 -o cufflinks_out ".
 #          ' merged_hits.bam >& log_cufflinks.log';
 my $cufflinks_cmd="cufflinks -m $frgsize -s 60 -p 8 -o cufflinks_out ".
           ' merged_hits.bam >& log_cufflinks.log';
 mlog("running cufflinks command:\n$cufflinks_cmd");
 system($cufflinks_cmd);
 mlog("cufflinks done (exit code $?)");
}
# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }
print STDERR "done.\n";
#************ Subroutines **************

sub mlog {
foreach my $m (@_) {
  print STDERR $m."\n";
  }
}

sub sysrun {
 my ($cmd)=@_; 
 my $errmsg = `($cmd) 2>&1`;
 my $exitcode=$?;
 if ($exitcode || ($errmsg=~/Error|Segmentation|Fail|Invalid|Cannot|Unable/si)) {
  print STDERR "!Error at:\n$cmd\n";
  print STDERR "    Exit code: $exitcode, message:\n$errmsg\n";
  #unlink(@todel);
  exit(1);
  }
}
