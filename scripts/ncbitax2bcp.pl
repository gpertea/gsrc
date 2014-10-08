#!/usr/bin/perl
use strict;
#use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q{Usage:

Download ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
And unpack it in the current directory

Then run $0 which will create
the bcp files for loading into the database.
};
umask 0002;
my %merged; #insert "merged" entries too, even though their data is duplicate
#getopts('o:') || die($usage."\n");
#my $outfile=$Getopt::Std::opt_o;
die("$usage\n") if $ARGV[0]=~m/^\-[\-help]*$/;
die("$usage Error: cannot find nodes.dmp!\n") unless -f 'nodes.dmp';
my $cmd='cut -f5 nodes.dmp | sort -u |';
open(TAXC, '>taxon_cat.bcp') || die "Error creating taxon_cat.bcp file!\n";

my %ranks; #rank_name => rank_id
open(RANK, $cmd) || die "Error opening $cmd pipe!\n";
my ($rank_id, $rank_name);
while (<RANK>) {
 chomp;
 s/^\s+//;s/\s+$//;
 next unless $_;
 $rank_id++; 
 $rank_name=$_;
 $ranks{$rank_name}=$rank_id;
 print TAXC join("\t", 'R', $rank_id, '', $rank_name, '')."\n";
 }
close(RANK);

open(MRG, 'merged.dmp') || die("Error opening merged.dmp!\n");
while (<MRG>) {
 chomp;
 next unless $_;
 s/\t\|$//;
 my @t=split(/\t\|\t/);
 push(@{$merged{$t[1]}},$t[0]);
}
close(MRG);

open(DIV, 'division.dmp') || die("Error opening division.dmp!\n");
while (<DIV>) {
 chomp;
 next unless $_;
 s/\t\|$//;
 my @t=split(/\t\|\t/);
 print TAXC join("\t", 'D', @t[0..3])."\n";
}

close(DIV);
close(TAXC);
open(NAMES, 'names.dmp') || die("Error opening names.dmp!\n");
open(TAXN, '>taxon_names.bcp') || die("Error creating taxon_names.bcp!\n");
my %names; # tax_id => [sci_name, com_name]
while(<NAMES>) {
 chomp;
 next unless $_;
 s/\t\|$//;
 my @t=split(/\t\|\t/);
 print TAXN join("\t", @t)."\n";
 if ($t[3]=~/\bcommon name/) {
   my $d=$names{$t[0]};
   $names{$t[0]}= $d ? [$$d[0], $t[1]] : ['',$t[1]];   
   }
  elsif ($t[3]=~/\bscientific name/) {
   my $d=$names{$t[0]};
   $names{$t[0]}= $d ? [$t[1], $$d[1]] : [$t[1]];
   } 
 }
close(NAMES);
close(TAXN);
open(TAX, '>taxon.bcp') || die("Error creating taxon.bcp!\n");
open(NODES, 'nodes.dmp') || die("Error opening names.dmp!\n");
while(<NODES>) {
 chomp;
 next unless $_;
 s/\t\|$//;
 my @t=split(/\t\|\t/);
 my $rank_id=$ranks{$t[2]};
 die("Error: invalid rank $t[2] from nodes.dmp line:\n$_\n")
   unless defined($rank_id);
 my $d=$names{$t[0]};
 die("Error: names not found for taxon $t[0] at nodes.dmp line:\n$_\n")
   unless $d && $$d[0];
 my $sci_name=$$d[0];
 my $com_name=$$d[1] || '';
 print TAX join("\t",$t[0],$t[1],$rank_id, $t[4], $t[5],
                       $sci_name, $com_name)."\n";
 if (my $md=$merged{$t[0]}) {
   foreach my $tid (@$md) {
     print TAX join("\t",$tid,$t[1],$rank_id, $t[4], $t[5],
                       $sci_name, $com_name)."\n";
     }
   }
}
close(NODES);
close(TAX);
