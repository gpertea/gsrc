#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 mem_usage.pl <pid>
/;
umask 0002;
getopts('o:') || die($usage."\n");
my $pid=$ARGV[0] || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my %cache;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my $smaps=getMemUsage($pid);

sub getMemUsage {
  my $pid=shift;
  my $name='/proc/'.$pid.'/smaps';
  open(my $f, $name) || die("Cannot open $name: $!\n");
  my $l;
  my %smaps;
  my ($vma_start, $vma_end, $perms, $file_off, $dev_major, $dev_minor,
      $inode, $file_name, $is_deleted);
  while (defined($l=<$f>)) {
    if ($l=~/([\da-f]+)-([\da-f]+)\s        # range
             ([r\-])([w\-])([x\-])([sp])\s  # access mode
             ([\da-f]+)\s                   # page offset in file
             ([\da-f]+):([\da-f]+)\s        # device
             (\d+)\s*                       # inode
             (.*?)                          # file name
              (\s\(deleted\))?$
            /xi ) {
       $vma_start=hex $1;
       $vma_end=hex $2;
       unless (exists($cache{$vma_start."\0".$vma_end})) {
          $cache{$vma_start."\0".$vma_end}=1;
          $perms=$3.$4.$5.$6;
          ($file_off,$dev_major,$dev_minor)=(hex($7), hex($8), hex($9));
          $inode=$10;
          $file_name=$11;
          $is_deleted=defined($12);
          }
       }
      elsif ($l=~/^(\w+):\s*(\d+) kB$/) {
        my $m=lc $1;
        $m=~s/\s/_/g;
        $smaps{$m}+=$2;
        }
    } # each line in /proc/<pid>/smaps
  close($f);
  return \%smaps;
}

foreach my $k (keys(%$smaps)) {
  printf('%16s:%9s'."\n",$k,fmtInt($smaps{$k}));
  }
# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************
sub fmtInt {
 (my $num = shift) =~ s/\G(\d{1,3})(?=(?:\d\d\d)+(?:\.|$))/$1,/g;
 return $num;
}
