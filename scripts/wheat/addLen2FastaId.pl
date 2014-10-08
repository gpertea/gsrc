#!/usr/bin/perl -w

use strict;

use File::Basename;
use AMOS::AmosFoundation;
use AMOS::AmosLib;
use AMOS::ParseFasta;

my $base = new AMOS::AmosFoundation();    
die("ERROR: problems creating an AMOS::Foundation object\n") unless($base);

my $VERSION = '$Revision: 1.2 $ ';

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

	my $result = $base->getOptions(
 		);
	$result or $base->bail("Command line parsing failed");
	
	####################		
		
        my $fh = new AMOS::ParseFasta(\*STDIN);	
        $fh or $base->bail("Fasta format expected in $fh: $!\n");
	
	while (my ($head, $seq) = $fh->getRecord())
    	{       
		$head =~ /^(\S+)(.*)/ or $base->bail("Sequence head format error\n");

		print ">",$1," ",length($seq),$2,"\n";
		print $seq,"\n";
	}
		        					
	exit 0;
}

