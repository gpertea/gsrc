#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 batchcmds.pl [<shellcmds_file>] -f <datafile> [-t <delim>]
 
 <shellcmds_file>  is a file containing an shell command(s) to be executed 
               for each line in the <datafile>; if missing, the command is 
               taken interactively from the standard input (use Ctrl+D to
               end the text input). Columns taken from <datafile> should 
               be specified by placeholders :0 , :1 , .. :9
               global authentication file, see note below)
    -f         <datafile> should contain lines with space delimited fields;
               each field will replace the corresponding placeholders 
               specified in the shellcmds_file (:0, :1 .. :9)
 /;
umask 0002;
my $qfile;
$qfile=shift if (substr($ARGV[0],0,1) ne '-');

getopts('t:f:h') || die($usage."\n");
#my $outfile=$Getopt::Std::opt_o;
if ($Getopt::Std::opt_h) { print $usage; exit;}

die($usage."Error: shell command file not found!\n")
               unless (!$qfile || -e $qfile);

my $datafile=$Getopt::Std::opt_f;
die($usage."Error: data file $datafile not found.") unless (-e $datafile);

my $delim=$Getopt::Std::opt_t;
my $query='';
if ($qfile) { #qfile given
    local $/=undef;#one sip
    open(INFILE, '<'.$qfile);
    $query=<INFILE>;
    close(INFILE);
    }
 else {
   print STDERR ">Enter shell command(s) to be executed for each line of\n".
                " $datafile (use :1 .. :9 column placeholders; press Ctrl+D to end):\n";
   local $/="\n";
   $query.=$_ while (<STDIN>); 
   }
die("No placeholders found in the given commands. Aborting..\n") 
    unless ($query =~ m/\:\d+/);

$query=~s/\\[\n\r]+/ /sg;
$query=~s/[\n\r]+/\x01/sg;
my @cmds=split(/\x01/, $query);
open(XLSTFILE, $datafile) || die ("Error opening data file $datafile !\n");
while (<XLSTFILE>) {
 chomp;
 #s/'/''/g;
 next unless $_;
 my @pdata;
 if ($delim) {
    @pdata=split(/$delim/);
    }
   else { 
    @pdata=split(/\s+/);
    }
  #print STDERR ">dataline: $pdata[0] $pdata[1]\n";
  foreach my $c (@cmds) {
     next if $c=~/^\s*$/s;
     my $cmd=$c; # copy value
     $cmd=~s/\:(\d+)\b/$pdata[$1]/sg;
          #the parameters are numbered from 0, not from 1 !
          #so the first field from the file will be :0
     # execute
     #print STDERR "  .. exec: $cmd\n";
     my $xcode=system($cmd);
     #sleep(1);
     print STDERR "Warning: error detected ($?) for command $cmd!\n" if $xcode;
     }
 } #for each param line
close(XLSTFILE);

#================ SUBROUTINES ============


