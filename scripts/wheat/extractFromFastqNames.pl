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
Program that extracts a set of sequences from a FASTQ file

Usage: $PRG -f filter_file < fastq_file 

        options:
	
	-f <n>		- Filter file: contains the ids of the sequences to filter, one/line
	-i <n>          - Column number    (Default 0)
	-suffix string  - Add this suffix to all ids
        -prefix string  - Add this prefix to all ids
	-negate		- Negate

	-h|help		- Print this help and exit;
	-V|version	- Print the version and exit;
	-depend		- Print the program and database dependency list;
	-debug <level>	- Set the debug <level> (0, non-debug by default); 
 
	OUTPUT:  
		Fasta sequences at the console
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
	my %opt; $opt{i}=0;
	my %h;
	my $id;
	
	# Configure TIGR Foundation
	$tigr_tf->setHelpInfo($HELPTEXT.$MOREHELP);
        $tigr_tf->setUsageInfo($HELPTEXT);
        $tigr_tf->setVersionInfo($REV);
        $tigr_tf->addDependInfo(@DEPENDS);
	
	# validate input parameters
	my $result = $tigr_tf->TIGR_GetOptions(
		"f=s"   =>      \$opt{f},
		"i=i"	=>	\$opt{i},
		"preffix=s"	=> \$opt{prefix},
		"suffix=s"	=> \$opt{suffix},
		"negate"	=> \$opt{negate}
	);
	$tigr_tf->printUsageInfoAndExit() if (!$result);
	
	#########################################################
	# read the name file
	defined($opt{f}) or $tigr_tf->bail("Input file missing") ;
	open(IN,$opt{f}) or $tigr_tf->bail("Cannot open input file ".$!) ;
	while(<IN>)
	{
		next if/^#/;

		my @f=split;
		next unless(@f);
		
		$id=$f[$opt{i}];
		$id=$opt{prefix}.$id if(defined($opt{prefix}));
		$id=$id.$opt{suffix} if(defined($opt{suffix}));
		$h{$id}=1;
	}
	close(IN);

	#########################################################
	# read the Fastq file
	while(<>)
	{	    		
		if($.%4==1)
		{ 
			/^@(\S+)/ or die "ERROR: $_";
			$id=$1 ;
		}

		print $_ if($h{$id} xor $opt{negate});
	}
                      
	exit 0;
}
