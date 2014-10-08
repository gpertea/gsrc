#!/usr/bin/perl -w

use strict;

use File::Basename;
use AMOS::AmosFoundation;
use AMOS::AmosLib;

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
	my %max;
	my %line;
	my ($i,$j)=(0,1);
	my $abs;

	my $result = $base->getOptions(
		"i=i"	=>	\$i,
		"j=i"	=>	\$j,
		"abs"	=>	\$abs
	);
	$result or $base->bail("Command line parsing failed");
			
	while(<>)
	{
		chomp;
		next if(/^$/ or /^#/);

		my @f=split;		
		my ($key,$val) =($f[$i],$f[$j]);
		
		if(!defined($max{$key}) or $max{$key}<$val)
		{
 			 $val=abs($val) if($abs);
			 $max{$key}=$val ;
			 $line{$key}=$_;
		}
	}		        	
	
	foreach my $key (sort {$max{$b}<=>$max{$a}} keys %max) 
	{ 
		print $line{$key},"\n"; 
	}
					
	exit 0;
}

