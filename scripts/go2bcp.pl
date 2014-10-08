#!/usr/bin/perl
use strict;
use Getopt::Std;
use LWP::Simple;
use FindBin;use lib $FindBin::Bin;
use dbSession;
my $url=
  'http://www.geneontology.org/ontology/gene_ontology_edit.obo';
my $usage = qq/Usage:
 go2bcp.pl [-f <input_obo_file>]
 
 Parses and loads the GO data in geanno db.
 The GO data are either downloaded from the embedded url:
 $url
 ..or, if -f option is given, from a local .obo file downloaded previously.
/;
umask 0002;
getopts('f:') || die($usage."\n");
my $file=$Getopt::Std::opt_f;
unless ($file) {
 print STDERR "downloading file..\n";
 $file='gene_ontology_data.obo';
 unlink($file);
 die("Error: file $file already exists (couldn't remove!)\n") 
   if -f $file;
 die "Error downloading with getstore()!\n" 
   if is_error(getstore($url, $file));
 die "Error: file $file is non-existent or zero size after retrieval!\n"
   unless -s $file;
 print STDERR "Download OK.\n";
}

my %relcode=( 'part_of' => 'P',
               'is_a' => 'I',
               'alt_id'=> 'A',
               'upd_ann'=>'U'
              );
my %namespace=(
'biological_process'=>'P',
'cellular_component'=>'C',
'molecular_function'=>'F'
 );              
#create the bcp files for loading

open(GOTERM, '>go_term.bcp') || die ("Error creating file go_term.bcp!\n");
open(GOLINK, '>go_link.bcp') || die ("Error creating file go_link.bcp!\n");
open(GOSYN,  '>go_synonyms.bcp') || die ("Error creating file go_synonyms.bcp!\n");
open(GOXREF, '>go_xref.bcp') || die("Error creating file go_xref.bcp!\n");
open(GOSET,  '>go_subset.bcp') || die("Error creating file go_subset.bcp!\n");

#my $ds=dbSession->new('geanno@NEOSYBASE');

open(FGO, $file) || die ("Error opening $file\n");

my $c_id; #current go term being parsed
my $c_type; # 'P', 'F', or 'C'
my $c_name;
my $c_def;
my $c_obsolete=0;
my $c_comment;
while (<FGO>) {
 if (m/^\s*\[Term\]/) {
   goTerm() if $c_id;
   my $id;
   do { $_=<FGO>;
       ($id)=(/^id:\s*GO\:(\d+)/);
     } until $id || !$_;
   die("Couldn't parse current term ID!\n") unless $id;
   $c_id=int($id);
   ($c_name, $c_type, $c_def, $c_comment)=(undef,undef,undef,undef);
   $c_obsolete=0;
   next;
   }
 next unless $c_id;
 next if m/^\s*$/;
 chomp;
 s/\s+$//;tr/\t/ /s;
 if (/^name:\s*(.+)/) {
   $c_name=$1;
   }
 elsif (/^namespace:\s*(\S+)/) {
   $c_type=$namespace{lc($1)};
   }
 elsif (/^def:\s*(.+)/) {
   $c_def=$1;
   $c_def =~ tr/"//d; #"#--
   }
 elsif (/^is_obsolete:\s*true/) {
   $c_obsolete=1;
   }
 elsif (/^relationship:\s*(\S+)\s*GO:(\d+)/) {
  my ($rel, $id2)=($1, int($2));
  $rel=$relcode{lc($rel)} || 
     die("Error: unrecognized relationship code '$rel' for term $c_id ($_)\n");
  print GOLINK join("\t",$c_id, $rel, $id2)."\n";
  }
 elsif (/^alt_id:\s*GO:(\d+)/) {
  my $id2=int($1);
  print GOLINK join("\t",$c_id, $relcode{'alt_id'}, $id2)."\n";
  }
 elsif (/^is_a:\s*GO:(\d+)/) {
  my $id2=int($1);
  print GOLINK join("\t",$c_id, $relcode{'is_a'}, $id2)."\n";
  }
 elsif (/^subset:\s*(\S+)/) {
  my $subset=$1;
  print GOSET join("\t",$c_id, $subset)."\n";
  }
 elsif (/^synonym:/) {
  my ($s,$t)=(m/\"(.+?)\"\s*(\w+)/);
  print STDERR "WARNING: parsing synonym & its type at $c_id, $_\n" unless $t;
  print GOSYN join("\t",$c_id, $s, $t)."\n" if $t;
  }
 elsif (m/^xref:\s*(\S+):(.+)$/) {
  my ($xdb, $xref)=($1,$2);
  $xref=~s/^\s+//;$xref=~s/\s+$//;
  print GOXREF join("\t",$c_id, $xdb, $xref)."\n";
  }
 elsif (/^comment:\s*(.+)/) {
  my $comment=$1;$comment=~tr/\t\n/  /s;
  $c_comment=$comment;
  my @upds=($comment=~m/to\s*update\s*annotation.+?GO:(\d+)/g);
  foreach my $u_id (@upds) {
    print GOLINK join("\t",$c_id, 'U', int($u_id))."\n";
    }
  }
 
} #while
goTerm() if $c_id;
close(FGO);

close(GOTERM);
close(GOLINK);
close(GOSYN);
close(GOXREF);
close(GOSET);


#-- now run bcpin

my $bcpcmd='bcpin -TI -b geanno@NEOSYBASE';
system($bcpcmd. ' go_term.bcp go_link.bcp go_synonyms.bcp go_xref.bcp go_subset.bcp') 
  && die("Error at bcpin!\n");

sub goTerm {
 print GOTERM join("\t",$c_id, $c_type, $c_name, $c_def, $c_obsolete, $c_comment)."\n";
}
