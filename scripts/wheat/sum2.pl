#!/usr/bin/perl -w
 
use strict;
use warnings;

# TIGR Modules
use TIGR::Foundation;

my $tigr_tf = new TIGR::Foundation;
my $PRG = $tigr_tf->getProgramInfo('name');
my $REV="1.0";
my @DEPENDS=("TIGR::Foundation");

# help info
my $HELPTEXT = qq~
Program that ...

Usage: $PRG file [options]
	  
  INPUT:   
  
  file:	
  
  options:

	-h|help		- Print this help and exit;
	-V|version	- Print the version and exit;
	-depend		- Print the program and database dependency list;
	-debug <level>	- Set the debug <level> (0, non-debug by default); 
 
  OUTPUT:  
~;

my $MOREHELP = qq~
Return Codes:   0 - on success, 1 - on failure.
~;


###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	# define variables
	my %options;
	my %sum;
	my %count;
	my %max;
	my @keys;
	my ($i,$j)=(0,1);
	my $min;
		
	# Configure TIGR Foundation
	$tigr_tf->setHelpInfo($HELPTEXT.$MOREHELP);
        $tigr_tf->setUsageInfo($HELPTEXT);
        $tigr_tf->setVersionInfo($REV);
        $tigr_tf->addDependInfo(@DEPENDS);
	
	# validate input parameters
	my $result = $tigr_tf->TIGR_GetOptions(
		"i=i"	=>	\$i,
		"j=i"	=>	\$j,
		"min=f"	=>	\$min
	);
	$tigr_tf->printUsageInfoAndExit() if (!$result);	

	while(<>)
	{
		my @f=split;
		next unless(@f);

		my $key=$f[$i];
		my $val=$f[$j];
		
		if(!defined($sum{$key}))
		{
			push @keys,$key;
			$sum{$key}=0;
		}
		$sum{$key}+=$val;
		$count{$key}++;
		$max{$key}=$val if(!defined($max{$key}) or $max{$key}<$val);
	}

	#print join "\t",("key","sum","count","mean"); print "\n";
	
	foreach my $key ( @keys) 
	{		
		#print join "\t",($key,$sum{$key},$count{$key},$sum{$key}/$count{$key});  print "\n";
		#print join "\t",($key,$count{$key},$max{$key},$sum{$key});  print "\n";

		next if(defined($min) and $sum{$key}<$min);

		print join "\t",($key,$sum{$key});  print "\n";
	}
			
	exit 0;
}
