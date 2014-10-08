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
	-suffix		- add /[12] suffixes

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
				
	# Configure TIGR Foundation
	$tigr_tf->setHelpInfo($HELPTEXT.$MOREHELP);
        $tigr_tf->setUsageInfo($HELPTEXT);
        $tigr_tf->setVersionInfo($REV);
        $tigr_tf->addDependInfo(@DEPENDS);
	
	# validate input parameters
	my $result = $tigr_tf->TIGR_GetOptions(
		"trim"	=>  \$options{trim},
		"suffix" => \$options{suffix},
		"minLen=i" => \$options{minLen},
		"uc"	    => \$options{uc},
		"454"	=> \$options{454},
		"RF"	=> \$options{RF}

	);
	$tigr_tf->printUsageInfoAndExit() if (!$result);	
	my ($qry,$seq,$qlt,$dir);

	####################################
	
	while(<>)
	{
		# 0                   1   2          3        4    5    6  7  8  9                                     10                                    11...		
		# SRR022906.1.1       0   NC_000913  4579753  255  36M  *  0  0  TGGTCGGTGTTTTTGTTCTCTTCGCTGTCTGCCACT  IIIIIIIDIIIIII1IIIIIIIIIII%III?I&8II  ...
	
		next if /^\@/;

		my @F=split;
		next if(defined($qry) and ($qry eq $F[0]) and ($dir eq $F[1]&0x40));

		$qry=$F[0];
		$seq=$F[9];
		$qlt=$F[10];
		$dir=$F[1]&0x40;

		#if($options{trim})
		if(1)
		{
			if($F[5]=~/^(\d+)S/)
			{
				my $trim=$1;
				$seq=substr($seq,$trim,length);
			}

			if($F[5]=~/M(\d+)S$/)
        	        {
                	        my $trim=$1;
                        	$seq=substr($seq,0,length($seq)-$trim);
        	        }
		}

		next if($options{minLen} and $options{minLen}>length($seq));

		#test if reversed alignment
		if($F[1] & 0x10)
		{        
			$seq=~tr/acgtACGT/tgcaTGCA/;
                       	$seq=scalar reverse($seq);
		}			

		$seq=uc($seq) if($options{uc});

		# all seqs become FR
		if($options{454} and $options{RF})
                {
                        $seq=~tr/ACGT/TGCA/;
                        $seq=scalar reverse($seq);
               	}

		print ">$qry";
		if($options{suffix})
		{
			if($F[1]&0x40) { print "/1"; }
			else           { print "/2"; }
		}
		elsif($options{454})
		{
			$qry=~/^(.........)/;
			my $lib=$1;

			if($F[1]&0x40) { print "/1 template=$qry library=$lib dir=F"; }
                        else           { print "/2 template=$qry library=$lib dir=R"; }
		}
		
		print "\n$seq\n";
	}	

	exit 0;
}
