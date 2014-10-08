#!/usr/bin/perl -w
 
use strict;
use warnings;

#use lib "/nfshomes/dpuiu/JCVI/common" ;

# TIGR Modules
use TIGR::Foundation;

my $tigr_tf = new TIGR::Foundation;
my $PRG = $tigr_tf->getProgramInfo('name');
my $REV="1.0";
my @DEPENDS=("TIGR::Foundation");

# help info
my $HELPTEXT = qq~
Program that counts the number of elements in a certain column.
Default to the first column.

Usage: $PRG file [options]
	  
  INPUT:   
  	
  file:			- TAB delimited format	
  
  options:

  	-h|help		- Print this help and exit;
	-V|version	- Print the version and exit;
	-depend		- Print the program and database dependency list;
	-debug <level>	- Set the debug <level> (0, non-debug by default); 

	-i|c int 	- column number starting from 0
	-min int	- display only instances that occur at least min times
	-Max int        - display only instances that occur at most Max times
	-p real		- show percentages (default 0)
	-t string	- display text
	-u		- test for uniqueness
 
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
	my %options;

	my %count;
	my $count;
	$options{percentage}=0;
	
	# Configure TIGR Foundation
	$tigr_tf->setHelpInfo($HELPTEXT.$MOREHELP);
        $tigr_tf->setUsageInfo($HELPTEXT);
        $tigr_tf->setVersionInfo($REV);
        $tigr_tf->addDependInfo(@DEPENDS);
	
	# validate input parameters
	my $result = $tigr_tf->TIGR_GetOptions(
		"percentage=s"	=>	\$options{percentage},		
		"i|col=s"		=> 	\$options{col},
		"min=s"		=> 	\$options{min},
		"Max=s"	   	=>      \$options{max},
	);
	$tigr_tf->printUsageInfoAndExit() if (!$result);

	######################################################
	
	while(<>)
	{
		chomp;
		next if(/^$/) ;                     
		my $key=$_;

		if(defined($options{col}))
		{
			my @f=split;
			$key=""; 
			$key=$f[$options{col}] if(defined($f[$options{col}]));
		}

		$count{$key}++;
		$count++;
	}

	##########################################################

	foreach my $key (sort {$count{$b}<=>$count{$a}} keys %count)
	{
		next if (defined($options{min}) and $count{$key}<$options{min});
		next if (defined($options{max}) and $count{$key}>$options{max});

                print $key,"\t", $count{$key};
                printf("\t%.4f",$count{$key}/$count) if($options{percentage});
                print "\n";
	}

	exit 0;
}
