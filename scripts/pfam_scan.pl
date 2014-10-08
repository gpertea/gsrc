#!/usr/bin/perl -w

# Copyright (c) 2002-2005 Genome Research Ltd
# Distributed under the same terms as the Pfam database
# See ftp://ftp.sanger.ac.uk/pub/databases/Pfam/COPYRIGHT
# for more details.
#
# Written by Sam Griffiths-Jones borrowing heavily from some
# old code by Ewan Birney.  Contact pfam@sanger.ac.uk for
# help.
#
# version 0.5 (16/3/2005)

=head1 NAME

pfam_scan.pl - search protein fasta sequences against the Pfam
library of HMMs.

=head1 VERSION

This is version 0.5 of pfam_scan.pl.  See the history section for
recent changes.

Behaviour of recent versions is a significantly different from 0.1.
From version 0.5, overlapping matches to families within the same clan
are removed, keeping the best scoring hit.  This behaviour can be
overridden with the --overlap option.  From version 0.2, we can use
BLAST to preprocess the input sequences with the --fast option, so we
only have to search a subset of sequences against a subset of HMMs
using hmmpfam.  For puritanical reasons we don't do this by default!
Read the notes about this below.

This version has been tested with Perl 5.6.1, Pfam 10.0 (through
13.0), Bioperl 1.2 and HMMER 2.3.1.  It should work with any versions
higher than these.

=head1 REQUIREMENTS

 - this script
 - Perl 5.6 or higher (and maybe lower)
 - The Pfam database (downloadable from
   ftp://ftp.sanger.ac.uk/pub/databases/Pfam/)
 - HMMER software (from http://hmmer.wustl.edu/)
 - NCBI BLAST binaries (from http://www.ncbi.nlm.nih.gov/Ftp/)
 - Bioperl (from http://bio.perl.org/)

The Bioperl modules directory must be in your perl library path, and
the HMMER and BLAST binaries must be in your executable path.

You also need to be able to read and write to /tmp on your machine.

Some of these requirements are easily circumvented, but this script
should at least give you a start.

=head1 HOW TO INSTALL PFAM LOCALLY

1. Get the Pfam database from
   ftp://ftp.sanger.ac.uk/pub/databases/Pfam/.  In particular you need
   the files Pfam-A.fasta, Pfam_ls, Pfam_fs, and Pfam-A.seed.

2. Unzip them if necessary
    $ gunzip Pfam*.gz

3. Grab and install HMMER, NCBI BLAST and Bioperl, and make sure your
   paths etc are set up properly.

4. Index Pfam-A.fasta for BLAST searches
    $ formatdb -i Pfam-A.fasta -p T

5. Index the Pfam_ls and Pfam_fs libraries for HMM fetching
    $ hmmindex Pfam_ls
    $ hmmindex Pfam_fs


=head1 SEARCHING PFAM

This script is really just a wrapper around hmmpfam.

Run pfam_scan.pl -h to get a list of options.  Probably the only thing
to worry about is supplying the -d option with the location of your
downloaded Pfam database.  Or you can set the PFAMDB environment
variable to point to the right place and things should work without
-d.  And you should decide whether or not to use --fast.

A few things to note: 

--fast uses BLAST as a preprocessor to reduce the amount of compute we
have to do with hmmpfam.  This is known to reduce sensitivity in the
case of a very small number of families (those whose length is
exceptionally short, like the XYPPX repeat).  If you're annotating
genomes then you *probably* don't care too much about these families.
Omiting this option may give you a small added sensitivity, but with a
rough 10 fold time cost.  If you want to exactly replicate the Pfam
web site results or distributed data, you probably shouldn't use this.

Overlapping above-threshold hits to families within the same clan are
removed -- only the best scoring hit is kept.  You can override this
behaviour with the --overlap option.

Pfam provides two sets of models, called ls and fs models, for whole
domain and fragment searches.  This wrapper basically returns all hits
to the ls models, and then adds to these all non-overlapping hits to
the fragment models.  This mimics the behaviour of Pfam web site
searches.  You can choose to search only one set of models with the
--mode option.

Unless you want to grub around in the noise you should probably use
the default thresholds - these are hand curated for every family by
the Pfam team, such that we believe false positives will not score
above these levels.  The consequence is that some families may miss
members.

You may want to adjust the threshold used for the preprocessing BLAST
search (default evalue 10).  Raising this to 50 will slow everything
down a bit but may gain you a little sensitivity.  Lowering the evalue
cutoff will speed things up but with an inevitable sensitivity cost.

It is important that each sequence in the fasta file has a unique
identifier.  Note that the fasta header format should be:

>identifier  <optional description>

so the identifier should not contain whitespace.

The format of the output is:

<seq id> <seq start> <seq end> <hmm acc> <hmm start> <hmm end> <bit score> <evalue> <hmm name>

hmmpfam returns scores for sequence and domain matches seperately.
For simplicity, the single line for each domain format returned here
reports domain scores.

=head1 BUGS

Many options are not rigorously tested.  Error messages are
uninformative.  The documentation is inadequate.  You may find it
useful.  You may not.

=head1 HISTORY

Version     Main changes
-------     ------------

0.5         Removes overlapping above-threshold hits to families
            within the same clan. --overlap overrides.

0.4         Work-around for hmmpfam bug/feature that reports hits
            above domain threshold even if the sequence doesn't 
            score above the sequence threshold.

0.3         Fix minor bugs to be compatable with HMM versions in
            Pfam 13.

0.2         --fast option to use BLAST preprocessing for significant
            speed-up.

0.1         First effort, simply wraps up hmmpfam without doing
            anything clever.

=head1 CONTACT

This script is copyright (c) Genome Research Ltd 2002-2005.  Please
contact pfam@sanger.ac.uk for help.

=cut

use strict;
use lib '/fs/sz-user-supported/common/lib/bioperl-1.5.0';
use Getopt::Long;
use IO::File;
use Bio::SeqIO;
use Bio::SearchIO;


sub show_help {
    print STDERR <<EOF;

$0: search a fasta file against Pfam

Usage: $0 <options> fasta_file
    Options
        -h             : show this help
	-d <dir>       : directory location of Pfam flatfiles
	-o <file>      : output file, otherwise send to STDOUT
        --fast         : use BLAST preprocessing

    Expert options
	--overlap      : show overlapping hits within clan member families
	--mode <ls|fs> : search only ls or fs models (default both)
	-e <n>         : specify hmmpfam evalue cutoff, else use Pfam cutoffs
	-t <n>         : specify hmmpfam bits cutoff, else use Pfam cutoffs
	-be <n>        : specify blast evalue cutoff (default 1)

    Output format is:
        <seq id> <seq start> <seq end> <hmm acc> <hmm start> <hmm end> <bit score> <evalue> <hmm name>

EOF
}

my( $help,
    $ecut,
    $bcut,
    $outfile,
    $pfamdir,
    $mode,
    $fast,
    $noclean,
    $verbose,
    $overlap,
    );

# defaults
my $blastecut = 1;

if( $ENV{ 'PFAMDB' } ) {
    $pfamdir = $ENV{ 'PFAMDB' };
}
else {
    $pfamdir = '.';
}

&GetOptions( "mode=s"     => \$mode,
	     "fast"       => \$fast,
	     "o=s"        => \$outfile,
	     "e=s"        => \$ecut,
	     "t=s"        => \$bcut,
	     "h"          => \$help,
	     "be=s"       => \$blastecut,
	     "d=s"        => \$pfamdir,
	     "noclean"    => \$noclean,
	     "v"          => \$verbose,
	     "overlap"    => \$overlap,
	     );

my $fafile = shift;

if( $help ) {
    &show_help();
    exit(1);
}
if( !-s $fafile ) {
    &show_help();
    print STDERR "FATAL: can't find fasta file [$fafile]\n\n";
    exit(1);
}   
if( !-d $pfamdir ) {
    &show_help();
    print STDERR "FATAL: database location [$pfamdir] is not a directory\n\n";
    exit(1);
}    
if( !-s "$pfamdir/Pfam_ls" and ( $mode ne "fs" ) ) {
    &show_help();
    print STDERR "FATAL: database location [$pfamdir] does not contain Pfam_ls\n\n";
    exit(1);
}    
if( !-s "$pfamdir/Pfam_fs" and ( $mode ne "ls" ) ) {
    &show_help();
    print STDERR "FATAL: database location [$pfamdir] does not contain Pfam_fs\n\n";
    exit(1);
}    
if( $fast ) {
    if( !-s "$pfamdir/Pfam-A.fasta" ) {
	&show_help();
	print STDERR "FATAL: database location [$pfamdir] does not contain Pfam-A.fasta.\nYou cannot use --fast without this file.\n\n";
	exit(1);
    }    
    if( !-s "$pfamdir/Pfam-A.fasta.phr" ) {
	&show_help();
	print STDERR "FATAL: you need to index Pfam-A.fasta using formatdb to use --fast\n\n";
	exit(1);
    }    
    if( !-s "$pfamdir/Pfam_ls.ssi" or !-s "$pfamdir/Pfam_fs.ssi" ) {
	&show_help();
	print STDERR "FATAL: you need to index Pfam_ls and Pfam_fs using hmmindex to use --fast\n\n";
	exit(1);
    }	
}

if( !$overlap and !-s "$pfamdir/Pfam-C" ) {
    print STDERR "WARNING: database directory [$pfamdir] does not contain Pfam-C.\nYou may get overlapping hits to families in the same clan\nwithout this file.\n\n";
}



my $maxidlength = 1;
open( FA, $fafile ) or die "FATAL: can't open fasta file [$fafile]\n";
while(<FA>) {
    if( /^\>(\S+)/ ) {
	$maxidlength = length( $1 ) if( $maxidlength < length( $1 ) );
    }
}
close FA;

my $options = "";
if( $ecut ) {
    $options .= "-E $ecut";
}
elsif( $bcut ) {
    $options .= "-T $bcut";
}
else {
    $options .= "--cut_ga";
}

# map Pfam accessions to ids
# expensive, but we have to read Pfam-A.seed or Pfam-A.fasta
my %accmap;
if( -s "$pfamdir/Pfam-A.seed" ) {
    my $id;
    open( SEED, "$pfamdir/Pfam-A.seed" ) or die "FATAL: can't open Pfam-A.seed file\n";
    while(<SEED>) {
	if( /^\#=GF ID\s+(\S+)/ ) {
	    $id = $1;
	}
	if( /^\#=GF AC\s+(PF\d+\.?\d*)/ ) {
	    $accmap{$id} = $1;
	}
    }
}
elsif( -s "$pfamdir/Pfam-A.fasta" ) {
    open( FASTA, "$pfamdir/Pfam-A.fasta" ) or die;
    while(<FASTA>) {
	if( /^\>\S+\s+(PF\d+\.?\d*)\;(\S+)\;/ ) {
	    $accmap{$2} = $1;
	}
    }
    close FASTA;
}
else {
    die "FATAL: can't get Pfam id/accession mapping.  I need either\nPfam-A.seed or Pfam-A.fasta.\n";
}
$options .= " -Z ".scalar( keys %accmap );


my( %clanmap, $clanacc );
# read the clan mappings in
if( !$overlap ) {
    open( CLAN, "$pfamdir/Pfam-C" ) or die;
    while(<CLAN>) {
	if( /^AC\s+(\S+)/ ) {
	    $clanacc = $1;
	}
	if( /^MB\s+(\S+)\;/ ) {
	    $clanmap{$1} = $clanacc;
	}
    }
    close CLAN;
}


# The organisation here is not intuitive.  This is the first way I
# thought of reusing code nicely for both the full hmmpfam search
# and the mini database searches created from blast results

my( $hmm_seq_map, @hmmlist, %seq );

if( $fast ) {
    # read sequences in
    my $seqin = Bio::SeqIO -> new( '-file'   => $fafile,
				   '-format' => 'Fasta' );
    while( my $seq = $seqin->next_seq() ) {
	if( exists $seq{ $seq->id } ) {
	    die "FATAL: We're going to struggle here - you have two sequences\nwith the same name in your fasta file\n";
	}
	$seq{ $seq->id } = $seq;
    }

    # then run blast searches
    my $blastdb  = "$pfamdir/Pfam-A.fasta";
    my $blastcmd = "blastall -p blastp -i $fafile -d $blastdb -e $blastecut -F F -b 100000 -v 100000";
    print STDERR "running [$blastcmd > tmp_$$.blast] ....\n" if( $verbose );
    system "$blastcmd > tmp_$$.blast" and die "FATAL: failed to run blastall - is it in your path?\n";
    print STDERR "parsing blast output ....\n" if( $verbose );
    $hmm_seq_map = parse_blast( "tmp_$$.blast" );
    unless( $noclean ) {
	unlink( "tmp_$$.blast" ) or die "FATAL: failed to remove temporary files - can you write to tmp?\n";
    }
    @hmmlist = keys %{ $hmm_seq_map };
}

unless( $fast ) {
    @hmmlist = (1);  # just so we get through the following while loop once
}

my $allresults = HMMResults->new();

while( my $hmmacc = shift @hmmlist ) {
    my $seqfile;
    if( $fast ) {
	$seqfile = "tmp_$$.fa";

	open( SOUT, ">$seqfile" ) or die "FATAL: failed to create temporary files - can you write to /tmp?\n";
	my $seqout = Bio::SeqIO -> new( '-fh'     => \*SOUT,
					'-format' => 'Fasta' );
    
	foreach my $id ( keys %{ $hmm_seq_map->{$hmmacc} } ) {
	    if( not exists $seq{ $id } ) {
		warn "can't find [$id] in your sequence file\n";
	    }
	    else {
		$seqout -> write_seq( $seq{ $id } );
		print STDERR "searching [$id] against [$hmmacc]\n" if( $verbose );
	    }
	}
	close SOUT;
    }
    else {
	$seqfile = $fafile;
    }

    foreach my $m ( "ls", "fs" ) {
	my $hmmfile;
	if( $mode ) {
	    next unless( $m eq $mode );
	}
	if( $fast ) {
	    $hmmfile = "tmp_$$.hmm_$m";
	    system "hmmfetch $pfamdir/Pfam_$m $hmmacc > $hmmfile" and die "FATAL: failed to fetch [$hmmacc] from $pfamdir/Pfam_$m\n";
	}
	else {
	    $hmmfile = "$pfamdir/Pfam_$m";
	}

	my $fh = new IO::File;
	my $results = new HMMResults;

	$fh -> open( "hmmpfam --compat -A 0 $options $hmmfile $seqfile |" );
	$results -> parse_hmmpfam( $fh );
	$fh -> close;
	$allresults -> merge_results( $results );
	if( -e "tmp_$$.hmm_$m" and not $noclean ) {
	    unlink "tmp_$$.hmm_$m" or die "FATAL: failed to remove temporary files - can you write to /tmp?\n";
	}
    }
    if( -e "tmp_$$.fa" and not $noclean ) {
	unlink "tmp_$$.fa" or die "FATAL: failed to remove temporary files - can you write to /tmp?\n";
    }
}

if( !$overlap ) {
    $allresults = $allresults->remove_overlaps( \%clanmap );
}


if( $outfile ) {
    open( O, ">$outfile" ) or die "FATAL: can't write to your output file [$outfile]\n";
    $allresults -> write_ascii_out( \*O );
    close O;
}
else {
    $allresults -> write_ascii_out( \*STDOUT );
}

exit(0);



###############
# SUBROUTINES #
###############

sub parse_blast {
    my $file = shift;
    my %hmm_seq_map;
    my %already;

    my $searchin = Bio::SearchIO -> new( '-file'   => $file,
					 '-format' => 'Blast' );

    while( my $result = $searchin->next_result() ) {
        while( my $hit = $result -> next_hit() ) {
            my( $acc, $id ) = $hit->description =~ /(PF\d+\.?\d*)\;(\S+)\;/;
#	    print STDERR "seq [",$result->query_name(),"] hmm [$acc]\n" if( $verbose );
	    $hmm_seq_map{$acc} -> { $result->query_name() } = 1;
	}
    }

    return \%hmm_seq_map;
}


############
# PACKAGES #
############

package HMMResults;

sub new {
    my $ref = shift;
    my $class = ref($ref) || $ref;
    my $self = {
	'domain'  => [],
	'seq'     => {},
	'hmmname' => undef };
    bless $self, $class;
    return $self;
}

sub hmmname {
    my $self  = shift;
    my $value = shift;
    $self->{'hmmname'} = $value if( defined $value );
    return $self->{'hmmname'};
}

sub addHMMUnit {
    my $self = shift;
    my $unit = shift;
    my $name = $unit->seqname();

    if( !exists $self->{'seq'}->{$name} ) {
	warn "Adding a domain of $name but with no HMMSequence. Will be kept in domain array but not added to a HMMSequence";
    } else {
	$self->{'seq'}->{$name}->addHMMUnit($unit);
    }
    push( @{$self->{'domain'}}, $unit );
}

sub eachHMMUnit {
    my $self = shift;
    return @{$self->{'domain'}};
}

sub addHMMSequence {
    my $self = shift;
    my $seq  = shift;
    my $name = $seq->name();
    if( exists $self->{'seq'}->{$name} ) {
	warn "You already have $name in HMMResults. Replacing by a new entry!";
    }
    $self->{'seq'}->{$name} = $seq;
}

sub eachHMMSequence {
    my $self = shift;
    my (@array,$name);

    foreach $name ( keys %{$self->{'seq'}} ) {
	push( @array, $self->{'seq'}->{$name} );
    }

    return @array;
}

sub getHMMSequence {
    my $self = shift;
    my $name = shift;
    return $self->{'seq'}->{$name};
}

sub parse_hmmpfam {
    my $self = shift;
    my $file = shift;

  QUERY: {
      my( %hash, $seqname );
      while(<$file>) {
	  /^Scores for sequence/ && last;
	  if( /^Query:\s+(\S+)/ ) {
	      $seqname = $1;
	      my $seq = HMMSequence->new();
	      $seq->name($seqname);
	      $self->addHMMSequence($seq);
	  }
      }
	  
      while(<$file>) {
	  /^Parsed for domains/ && last;
	  if( my( $id, $sc, $ev, $nd ) = /^\s*(\S+).+?\s(\S+)\s+(\S+)\s+(\d+)\s*$/ ) {
	      $hash{$id} = $sc;
	  }
      }
	  
      while(<$file>) {
	  /^\/\// && redo QUERY;
	  
	  if( my( $id, $sqfrom, $sqto, $hmmf, $hmmt, $sc, $ev ) = 
	      /(\S+)\s+\S+\s+(\d+)\s+(\d+).+?(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)\s*$/ ) {
		  
	      if( !exists($hash{$id}) ) {
		  # hmmpfam seems to return hits above the domain score that aren't
		  # above the sequence score.  We don't want these, so skip if we
		  # haven't seen a sequence score.
		  next;
#		  warn("HMMResults parsing error in hmmpfam for $id - can't find sequence score");
	      }
	      
	      my $unit = HMMUnit->new();
	      $unit ->   seqname($seqname);
	      $unit ->   hmmname($id);
	      $unit ->    hmmacc($accmap{$id});
	      $unit -> start_seq($sqfrom);
	      $unit ->   end_seq($sqto);
	      $unit -> start_hmm($hmmf);
	      $unit ->   end_hmm($hmmt);
	      $unit ->      bits($sc);
	      $unit ->    evalue($ev);
	      $unit->seqbits($hash{$id});
	      
	      # this should find it's own sequence!
	      $self->addHMMUnit($unit);
	  }
      }
  }
    return $self;
}


sub get_unit_nse {
    my $self = shift;
    my $seqname = shift;
    my $domname = shift;
    my $start = shift;
    my $end = shift;
    my($seq,$unit);
    
    $seq = $self->getHMMSequence($seqname);
    if( !defined $seq ) {
	warn("Could not get sequence name $seqname - so can't get its unit");
	return undef;
    }
    foreach $unit ( $seq->eachHMMUnit() ) {
	if( $unit->hmmname() eq $domname && $unit->start_seq() == $start &&  $unit->end_seq() == $end ) {
	    return $unit;
	}
    }

    return undef;
}

sub write_ascii_out {
    my $self = shift;
    my $fh = shift;

    if( !defined $fh) {
	$fh = \*STDOUT;
    }

    my @units = sort { $a->seqname cmp $b->seqname ||
		       $a->start_seq() <=> $b->start_seq } $self->eachHMMUnit();

    foreach my $unit ( @units ) {
	print $fh sprintf( "%-".$maxidlength."s  %4d %4d %-11s %4d %4d %7s  %8s  %s\n",
			   $unit->seqname(),
			   $unit->start_seq(),
			   $unit->end_seq(),
			   $unit->hmmacc,
			   $unit->start_hmm,
			   $unit->end_hmm,
			   $unit->bits,
			   $unit->evalue,
			   $unit->hmmname );
    }	    
}


sub remove_overlaps {
    # self is unchanged
    # so need to grab the returned results object
    my $self = shift;
    # if clanmap is specified then remove overlaps between
    # families within a clan
    my $clanmap = shift;

    my $new = HMMResults->new();

  UNIT:
    foreach my $unit ( sort { $b->bits <=> $a->bits } $self->eachHMMUnit() ) {
	if( not $new->getHMMSequence( $unit->seqname ) ) {
	    my $seq = HMMSequence->new();
	    $seq->name( $unit->seqname );
	    $seq->bits( $unit->seqbits );
	    $new->addHMMSequence( $seq );
	}
	else {
	    # check overlaps
	    my ($acc1) = $accmap{ $unit->hmmname } =~ /(\S+)\.\d+/;   # grrrrr
	    foreach my $u ( $new->getHMMSequence( $unit->seqname )->eachHMMUnit ) {
		my ($acc2) = $accmap{ $u->hmmname } =~ /(\S+)\.\d+/;

		if( $unit->hmmname eq $u->hmmname ) {
		    next UNIT if( $unit->overlap( $u ) );
		}

		if( ( $clanmap and 
		      exists $clanmap->{ $acc1 } and 
		      exists $clanmap->{ $acc2 } and 
		      $clanmap->{ $acc1 } eq $clanmap->{ $acc2 } ) ) {
		    next UNIT if( $unit->overlap( $u ) );
		}
	    }
	}

	$new->addHMMUnit( $unit );
    }

    return $new;
}

    
sub merge_results {
    # merge two HMMResults objects, preferring those from self over res in case of overlap
    my $self = shift;
    my $res  = shift;

  UNIT:
    foreach my $unit ( $res->eachHMMUnit() ) {
	my $seqname = $unit->seqname();
	my $seqbits = $unit->seqbits();
	my $oldseq  = $self->getHMMSequence($seqname);
	if( not $oldseq ) {
	    my $seq = HMMSequence->new();
	    $seq->name($seqname);
	    $seq->bits($seqbits);
	    $self->addHMMSequence($seq);
	}
	else {
	    foreach my $oldunit ( $oldseq->eachHMMUnit() ) {
		if( ( $unit->hmmname() eq $oldunit->hmmname() ) and $unit->overlap( $oldunit ) ) {
		    next UNIT;
		}
	    }
	}
	$self->addHMMUnit($unit);
    }
    return $self;
}


######################

package HMMSequence;

sub new {
    my $ref = shift;
    my $class = ref($ref) || $ref;
    my $self = {
	'name'   => undef,
	'desc'   => undef,
	'bits'   => undef,
	'evalue' => undef,
	'domain' => [] };
    bless $self, $class;
    return $self;
}

sub name {
    my $self = shift;
    my $name = shift;

    if( defined $name ) {
        $self->{'name'} = $name;
    }
    return $self->{'name'};
}

sub desc {
    my $self = shift;
    my $desc = shift;

    if( defined $desc ) {
        $self->{'desc'} = $desc;
    }
    return $self->{'desc'};
}

sub bits {
    my $self = shift;
    my $bits = shift;

    if( defined $bits ) {
        $self->{'bits'} = $bits;
    }
    return $self->{'bits'};
}

sub evalue {
    my $self   = shift;
    my $evalue = shift;

    if( defined $evalue ) {
        $self->{'evalue'} = $evalue;
    }
    return $self->{'evalue'};
}

sub addHMMUnit {
    my $self = shift;
    my $unit = shift;

    push(@{$self->{'domain'}},$unit); 
}

sub eachHMMUnit {
    my $self = shift;
    
    return @{$self->{'domain'}};
}


###################

package HMMUnit;

sub new {
    my $ref = shift;
    my $class = ref($ref) || $ref;
    my $self = {
	seqname    => undef,
	seq_range  => new Range,
	hmmname    => undef,
	hmmacc     => undef,
	hmm_range  => new Range,
	bits       => undef,
	evalue     => undef,
	prob       => undef,
	seqbits    => undef,
	alignlines => [],
	mode       => undef };

    bless $self, $class;
    return $self;
}

sub seqname {
    my $self  = shift;
    my $value = shift;
    $self->{'seqname'} = $value if( defined $value );
    return $self->{'seqname'};
}

sub hmmname {
    my $self  = shift;
    my $value = shift;
    $self->{'hmmname'} = $value if( defined $value );
    return $self->{'hmmname'};
}

sub hmmacc {
    my $self  = shift;
    my $value = shift;
    $self->{'hmmacc'} = $value if( defined $value );
    return $self->{'hmmacc'};
}

sub bits {
    my $self  = shift;
    my $value = shift;
    $self->{'bits'} = $value if( defined $value );
    return $self->{'bits'};
}

sub evalue {
    my $self  = shift;
    my $value = shift;
    $self->{'evalue'} = $value if( defined $value );
    return $self->{'evalue'};
}

sub seqbits {
    my $self  = shift;
    my $value = shift;
    $self->{'seqbits'} = $value if( defined $value );
    return $self->{'seqbits'};
}

sub mode {
    my $self  = shift;
    my $value = shift;
    $self->{'mode'} = $value if( defined $value );
    return $self->{'mode'};
}

sub add_alignment_line {
    my $self = shift;
    my $line = shift;
    push(@{$self->{'alignlines'}},$line);
}

sub each_alignment_line {
    my $self = shift;
    return @{$self->{'alignlines'}};
}

sub get_nse {
    my $self = shift;
    my $sep1 = shift;
    my $sep2 = shift;

    if( !defined $sep2 ) {
	$sep2 = "-";
    }
    if( !defined $sep1 ) {
	$sep1 = "/";
    }

    return sprintf("%s%s%d%s%d",$self->seqname,$sep1,$self->start_seq,$sep2,$self->end_seq);
}


sub start_seq {
    my $self = shift;
    my $start = shift;

    if( !defined $start ) {
	return $self->{'seq_range'}->start();
    }
    $self->{'seq_range'}->start($start);
    return $start;
}

sub end_seq {
    my $self = shift;
    my $end = shift;

    if( !defined $end ) {
	return $self->{'seq_range'}->end();
    }
    $self->{'seq_range'}->end($end);
    return $end;

}


sub start_hmm {
    my $self = shift;
    my $start = shift;

    if( !defined $start ) {
	return $self->{'hmm_range'}->start();
    }
    $self->{'hmm_range'}->start($start);
    return $start;
}

sub end_hmm {
    my $self = shift;
    my $end = shift;

    if( !defined $end ) {
	return $self->{'hmm_range'}->end();
    }
    $self->{'hmm_range'}->end($end);
    return $end;

}


sub overlap {
    # does $self overlap with $unit?
    my $self = shift;
    my $unit = shift;
    my( $u1, $u2 ) = sort { $a->start_seq <=> $b->start_seq } ( $self, $unit );

    if( $u2->start_seq <= $u1->end_seq ) {
	return 1;
    }

    return 0;
}


#################

package Range;

sub new {
    my $ref = shift;
    my $class = ref($ref) || $ref;

    my $self = {
	start => undef,
	end   => undef
	};

    bless $self, $class;

    if(@_ == 2) {
	$self->start(shift);
	$self->end(shift);
    } elsif (@_) {
	print STDERR "Usage: new Range()\n\tnew Range(start, stop)\n";
    }

    return $self;
}

sub start {
    my $self  = shift;
    my $value = shift;
    $self->{'start'} = $value if( defined $value );
    return $self->{'start'};
}

sub end {
    my $self  = shift;
    my $value = shift;
    $self->{'end'} = $value if( defined $value );
    return $self->{'end'};
}


