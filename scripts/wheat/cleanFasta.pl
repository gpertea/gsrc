#!/usr/bin/perl -w

use strict;

use File::Basename;
use AMOS::AmosFoundation;
use AMOS::AmosLib;
use AMOS::ParseFasta;

my $base = new AMOS::AmosFoundation();    
die("ERROR: problems creating an AMOS::Foundation object\n") unless($base);

my $VERSION = '$Revision: 1.1 $ ';

my $PRG = fileparse($0);

my $HELP = qq~
Program that ...
	
	USAGE: 
		$PRG
	
	INPUT:
		
		options:
		  ...
		
	OUTPUT:
	
	EXAMPLE:
		input:
		  ...
		
		output:
		  ...
	
	RETURN CODES:   
		0 - on success, 1 - on failure.
		
~;

$base->setVersion($VERSION);
$base->setHelpText($HELP);

###############################################################################
#             
# Main program       
#
###############################################################################

MAIN:         
{
	my %opt;		# stores command line option values	
	$opt{len}=100;
	$opt{min}=0;

	my $result = $base->getOptions(
		"l=i"	=> \$opt{len},
		"m=i"	=> \$opt{min},
		"acc"	=> \$opt{acc}
	);
	$result    or $base->bail("Command line parsing failed");
	$opt{len}  or $base->bail("Command line -len argument error");
	
	####################		
		
	my $fh = new AMOS::ParseFasta(\*STDIN);
        $fh or $base->bail("Fasta format expected in $fh: $!\n");
	
	while (my ($head, $seq) = $fh->getRecord())
    	{       
		next if($opt{min} and $opt{min}>length($seq));

		$head =~ /^(\S+)/ or $base->bail("Sequence head format error\n");

		if($opt{acc} and $head=~/^gi\|\d+\|\w+\|(\w+)\S*\s*(.*)/) { print ">$1 $2\n"         }
		else                                            { print ">",$head,"\n"; }
		
		my @seq=split //,$seq;
		my $len=scalar(@seq);
		my $i=0;
		
		if($opt{len})
		{
			$seq="";
			while($i<$len)
			{
				my $j=($i+$opt{len}-1<$len-1)?$i+$opt{len}-1:$len-1;
				
				$seq.= join "",@seq[$i..$j];
				$seq.="\n";
				
				$i+=$opt{len};
			}
		}	
		
		print uc($seq);
	}
		        					
	exit 0;
}

