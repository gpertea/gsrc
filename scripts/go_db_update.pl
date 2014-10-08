#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q{Usage:
 go_db_update.pl [-f <path/to/go_daily-termdb-data.gz>]

 Prepares data from go_daily-termdb-data.gz in "bcp" tab delimited
 format for the following tables:

 go_term
 go_term2term
 go_synonym
 go_xdb
 go_dbxref
 go_xref
 go_rel_prop
 go_rel_comp
};
umask 0002;
getopts('f:o:') || die($usage."\n");
my $infile=$Getopt::Std::opt_f || 'go_daily-termdb-data.gz';
die "$usage Error: cannot locate input file $infile\n" unless -f $infile;
# --
# - plan of action: parse term and term_definition data in memory so we can link them
# -  parse the other tables as they are (no changes)
# --
my %go_tables=('term'=>'go_term',  # special case - schema changed
               'term_definition'=>'definition', # special case: schema included in go_term
               'term2term'=>'go_term2term', 
               'term_synonym'=>'go_synonym',
               'relation_properties'=>'go_rel_prop',
               'relation_composition'=>'go_rel_comp',
               'db'=>'go_xdb',
               'dbxref'=>'go_dbxref',
               'term_synonym'=>'go_synonym',
               'term_dbxref'=>'go_xref'
               );
my %go_written; # table_name => written_flag
my @go_terms; # list of [id, name, term_type, acc, is_obsolete, is_root, is_relation]
my %go_defs; # id=>definition
open(GZ, "gzip -cd $infile |") || die ("Error starting gzip -cd $infile pipe !\n");

while (<GZ>) {
 if (m/^\s*INSERT\s+INTO\s+\`?([\-\.\w]+)/i) {
   my $table=$1;
   my $mytable=$go_tables{$table};
   next unless $mytable;
   chomp;
   s/^\s*INSERT\s+INTO\s+\`?[\-\.\w]+\`?\s+VALUES\s*\(//i;
   my @r=split(/\),\(/);
   $r[-1]=~tr/\n\r//d;
   $r[-1]=~s/\)$//;
   my $bcpwrite=($table ne 'term' && $table ne 'term_definition');
   if ($bcpwrite) {
     if (exists($go_written{$mytable})) {
        open(BCP, ">>$mytable.bcp") || die "Error writing to $mytable.bcp\n";
        }
      else {
        open(BCP, ">$mytable.bcp") || die "Error creating $mytable.bcp\n";
        $go_written{$mytable}=1;
        }
      }
   foreach my $rd (@r) {
     # for each row
     $rd=~s/\\'/\x81/g;
     $rd=~s/\\"/\x82/g;
     $rd=~tr/\\//d; #there is no reason to protect anything else
     $rd=~tr/\n\t\r / /s; #WARNING: tabs are clearly not allowed within field data!
     $rd=~s/\,\s+'/,'/g;
     $rd=~s/'\s+\,/',/g;
      my @f=();
      push(@f, defined($1) ? $1:$3) 
          while ($rd=~m/'([^'\\]*(\\.[^'\\]*)*)'|([^,]+)/g);
      foreach my $s (@f) {
        $s=~s/^['"]//;
        $s=~s/['"]$//;
        $s=~tr/\x81\x82/'"/; #'
        }
      if ($table eq 'term') {
        push(@go_terms,[@f]);
        next;
        }
      if ($table eq 'term_definition') {
        $go_defs{$f[0]}=$f[1];
        next;
        }
      if ($bcpwrite) {
        print BCP join("\t",@f)."\n";
        }
     } #for each row data
    close(BCP) if ($bcpwrite);
   } #if INSERT line
 } #while GZ lines
close(GZ);
# write go_term table
open(BCP, '>go_term.bcp') || die("Error creating file go_term.bcp\n");
foreach my $td (@go_terms) {
  my $tdef=$go_defs{$$td[0]} || '';
  print BCP join("\t",@$td,$tdef)."\n";
  }
close(BCP);
# --
#************ Subroutines **************

