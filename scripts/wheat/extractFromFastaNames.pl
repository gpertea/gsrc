#!/usr/bin/perl -w
 
use strict;
use warnings;
$|=1;

# TIGR Modules
use TIGR::Foundation;

my $tigr_tf = new TIGR::Foundation;
my $PRG = $tigr_tf->getProgramInfo('name');
my $REV="1.0";
my @DEPENDS=("TIGR::Foundation");

# help info
my $HELPTEXT = qq~
Program that extracts a set of sequences from a multi-Fasta file

Usage: $PRG fasta_file seq_name_file [options]

	INPUT:   
		the multi-Fasta file
		file containing the sequence names in the first column
		
        options:
	
	-c <n>		- Column number    (Default 0)
	-n <n>		- Occurance nucmer (Default undef)
	-w              - word match

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
	my %opt;
	my %h;
	$opt{i}=0;
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
		"prefix=s"      => \$opt{prefix},
                "suffix=s"      => \$opt{suffix},
                "negate"        => \$opt{negate},
		"word"		=> \$opt{word}
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
                $id.=$opt{suffix}    if(defined($opt{suffix}));
                $h{$id}=1;


	}
	close(IN);

	#########################################################
	# read the Fasta file
	while(<>)
	{	
    		
		chomp;
		if($opt{word} and />(\w+)(.*)/ or />(\S+)(.*)/)
		{
			#/>gi\|\d+\|\w+\|(\S+)\|/ or die "ERROR:$_";
			$id=$1 ;
		}

		print $_,"\n" if($h{$id} xor $opt{negate});
	}
                      
	exit 0;
}
