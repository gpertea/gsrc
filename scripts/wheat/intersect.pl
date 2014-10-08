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
Program that computes the difference of 2 data sets

Usage: $PRG file1 file2 [options]

	INPUT:   
		2 tab delimited files
	
	  options:
	  	-i number:  column number (starting from 0; default 0) 
				in the 1st file
		-j number:  column number (starting from 0; default 0) 
			in the 2nd file
			
	OUTPUT:  
		all the records in file1 that have a correspondent in file2
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
	my $i=0;
	my $i2;
	my $j=0;	 
	my %h;
	my $filter;
	my $ignoreVersion;
	my $s='\s+';

	# Configure TIGR Foundation
	$tigr_tf->setHelpInfo($HELPTEXT.$MOREHELP);
        $tigr_tf->setUsageInfo($HELPTEXT);
        $tigr_tf->setVersionInfo($REV);
        $tigr_tf->addDependInfo(@DEPENDS);
	
	# validate input parameters
	my $result = $tigr_tf->TIGR_GetOptions(
		"i=i"	=> \$i,
		"i2=i"	=> \$i2,
		"j=i"	=> \$j,
		"f=s"	=> \$filter,
		"s=s"	=>\$s,
		"ignoreVersion"	=> \$ignoreVersion
		
	);
	$tigr_tf->printUsageInfoAndExit() if (!$result);
	
	#########################################
	
	if(defined($filter)) { open(FILTER,$filter) or $tigr_tf->bail("Cannot open input file".$!) ;	}
	else                 { open(FILTER,$ARGV[1]) or $tigr_tf->bail("Cannot open input file".$!) ; }
	
	while(<FILTER>)
	{
    		my @f=split /$s/;
		next unless(defined($f[$j])); 
		$f[$j]=$1 if($ignoreVersion and $f[$j]=~/(.+)\.\d+$/);		
    		$h{$f[$j]}=1;
	}
	close(FILTER);
	last unless(%h);

	#########################################
	
	if(defined($filter)) {}
	else                 { open(STDIN,$ARGV[0]) or $tigr_tf->bail("Cannot open input file".$!) ;	}
	
	while(<STDIN>)
	{
		chomp;

    		my @f=split /$s/;
		next if(!defined($f[$i]));
		next if(defined($i2) and !defined($f[$i2]));	
		$f[$i]=$1 if($ignoreVersion and $f[$i]=~/(.+)\.\d+$/);
    		if(defined($h{$f[$i]})) { print $_,"\n" }
	}
	
	exit 0;
}
