#!/usr/bin/perl
use strict;
use Getopt::Std;
use LWP::Simple;
use FindBin;use lib $FindBin::Bin;
#use dbSession;
my $url='http://www.expasy.org/cgi-bin/lists?dbxref.txt';
my $usage = qq/Usage:
 dbxref2bcp.pl [-f <input_dbxref.txt>]
 
 Parses the UniProt's dbxref.txt info preparing it for loading into geanno db, 
 table xrefdbs. The required data file is either given explicitely or 
 downloaded from the embedded url:
 $url
/;

umask 0002;
getopts('f:') || die($usage."\n");
my $file=$Getopt::Std::opt_f;
unless ($file) {
 print STDERR "downloading file..\n";
 $file='dbxref_uniprot.txt';
 unlink($file);
 die("Error: file $file already exists (couldn't remove!)\n") 
   if -f $file;
 die "Error downloading with getstore()!\n" 
   if is_error(getstore($url, $file));
 die "Error: file $file is non-existent or zero size after retrieval!\n"
   unless -s $file;
 print STDERR "Download OK.\n";
}

my @knowntags=('Abbrev', 'Name', 'Cat', 'LinkTp', 'Ref',       'Server', 'Db_URL', 'Note');
#               xrefdb    name    cat    linktp    linktp_info  server    db_url    note
my %known;
@known{@knowntags}=(1) x scalar(@knowntags);

open(INF, $file) || die("Error opening file $file!\n");

my ($curtag, $curval);
open(OUTF, '>xrefdbs.bcp') || die("Error creating file xrefdbs.bcp!\n");

my %val;
while (<INF>) {
 if (m/^$/) {
   #end record
   if ($curtag) {
     storeCurVal();
     putValues();
     ($curtag,$curval)=(undef,undef);
     }
   next;
   }
 if (m/^Abbrev:\s*(.+)/) {
   ($curtag, $curval)=('Abbrev',$1);
   putValues();
   next;
   }
 next unless $curtag;
 chomp;
 if (m/^(\S+)\s*:(.+)/) { #new tag
   my ($t, $v)=($1,$2);
   die("Error: unrecognized tag $t at '$_'!\n") unless exists $known{$t} || $t eq 'AC';   
   storeCurVal() unless $t eq 'AC';
   $curtag=$t;
   $curval=$v;
   if ($curtag eq 'LinkTp') {
     if ($v=~m/^\s*(\w[\w ]+)\;\s+(\S.+)/) {
       $curval=$1;
       $val{$curtag}=$curval;
       $curtag='linfo';
       $curval=$2;
       }
     }
   }
 elsif (m/^\s+(.+)/) {
   my $vx=$1;
   $vx=~s/^\s+//;$vx=~s/\s+$//;
   $curval.=' '.$vx;
   }
 chomp;
}

storeCurVal() if ($curtag);
putValues();

close(INF);
close(OUTF);

# map { print STDERR $_."\t".$known{$_}."\n" } @knowntags;

sub putValues {
 return unless keys(%val);
 my @vals = map { $val{$_} } @knowntags;
 print OUTF join("\t",@vals)."\n";
 undef(%val);
}

sub storeCurVal {
 $curval=~s/^\s+//; $curval=~s/\s+$//;
 $val{$curtag}=$curval;
 $known{$curtag}=length($curval) if $known{$curtag}<length($curval);
 }

