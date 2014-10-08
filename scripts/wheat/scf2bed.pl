#!/usr/bin/perl -w

use strict;

use File::Basename;
use AMOS::AmosFoundation;
use AMOS::AmosLib;

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
	my (%opt);
	my $result = $base->getOptions(	);
	$result or $base->bail("Command line parsing failed");

	########################################################################

	my $scaf;
        while(<>)
        {
		my @F=split;
		if(/>(\w+)/) 
		{ 
			$scaf=$1 ;
		}
		else
		{
			print join "\t",($scaf,$F[1],$F[1]+$F[3],$F[0],$F[3],$F[2]);
			print "\n";
		}		
        }
     					
	exit 0;
}
