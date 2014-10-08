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
Program that sorts numerically a TAB file according to 2 columns

Usage: $PRG file [options]
	  
  INPUT:   
	TAB delimited file; The sorting columns (start from 0) must have numeric values  
  
  options:

	-i <n>		- First  column to sort after (Default 0)
  	-j <n>		- Second column to sort after (Default 1)
	
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

sub mySort
{
	my ($a1,$b1);
	
	$a1=$1 if($a=~/^(-?\d+)$/);
	$b1=$1 if($b=~/^(-?\d+)$/);
	
	if($a1 and $b1) { return $a1<=>$b1 ; }	
	else            { return $a cmp $b ; }
}

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	my ($i,$j,$absi,$absj)=(0,1);
	my $noPrefix;
	my %h;
	
	# Configure TIGR Foundation
	$tigr_tf->setHelpInfo($HELPTEXT.$MOREHELP);
        $tigr_tf->setUsageInfo($HELPTEXT);
        $tigr_tf->setVersionInfo($REV);
        $tigr_tf->addDependInfo(@DEPENDS);
	
	# validate input parameters
	my $result = $tigr_tf->TIGR_GetOptions(
		"i=i"		=>	\$i,
		"j=i"		=>	\$j,
		"noPrefix"	=>	\$noPrefix
	);
	$tigr_tf->printUsageInfoAndExit() if (!$result);


	##########################################################
	
	#parse input file
	while(<>)
	{		
		if(/^$/)    { next }		
		else
		{		
			my @f=split;
			$f[0]=~s/^#//;

			if($noPrefix)
			{
				$f[$i]=$1 if($f[$i]=~/.+\.(.+)/);
				$f[$j]=$1 if($f[$j]=~/.+\.(.+)/);
			}
	
			die unless(defined($f[$i]));
			die unless(defined($f[$j]));
	
			$h{$f[$i]}{$f[$j]}=$_;
		}
	}

	##########################################################

	foreach my $ki (sort mySort keys %h)
	{
		foreach my $kj (sort mySort keys %{$h{$ki}})
		{
			print $h{$ki}{$kj};
		}
	}
	
	exit 0;
}
