#!/usr/bin/perl  -w

use strict;

use scheduler;
use File::Basename;
use AMOS::AmosFoundation;
use AMOS::AmosLib;

my $base = new AMOS::AmosFoundation();    
die("ERROR: problems creating an AMOS::Foundation object\n") unless($base);

my $VERSION = '$Revision: 1.1 $ ';

my $PRG = fileparse($0);

my $HELP = qq~
Program that launches concurrent jobs 
	
	USAGE: 
		$PRG -concurrency <n> command_file
	
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
	my $concurrency=8;
	my $result = $base->getOptions(
		"concurrency=i"	=>	\$concurrency
	);
	$result or $base->bail("Command line parsing failed");

	while(<>)
	{
		chomp;

		next unless($_);
		next if(/^#/);

                &scheduler::schedulerSubmit($_);
        }

        &scheduler::schedulerSetNumberOfProcesses($concurrency);
        &scheduler::schedulerFinish();
			        					
	exit 0;
}


