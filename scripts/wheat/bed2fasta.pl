#!/usr/bin/perl -w
 
use strict;
use warnings;

# TIGR Modules
use TIGR::Foundation;
use TIGR::FASTAreader;
use TIGR::FASTArecord;

my $tigr_tf = new TIGR::Foundation;
my $PRG = $tigr_tf->getProgramInfo('name');
my $REV="1.0";
my @DEPENDS=("TIGR::Foundation","TIGR::AsmLib");

my $GAP_LEN=20;

# help info
my $HELPTEXT = qq~
Program that generates a file.scaffold.fasta out of a file.scaffold & file.fasta 

Usage: $PRG file.fasta [options]
	  
  INPUT:  
	file.fasta
	file.scaff  
  
  options:
	-s		- file.scaff (required)

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
	my %options;
        my ($contigId,%contigSeq);
	
	# Configure TIGR Foundation
	$tigr_tf->setHelpInfo($HELPTEXT.$MOREHELP);
        $tigr_tf->setUsageInfo($HELPTEXT);
        $tigr_tf->setVersionInfo($REV);
        $tigr_tf->addDependInfo(@DEPENDS);
	
	# validate input parameters
	my $result = $tigr_tf->TIGR_GetOptions(
		"bed=s"	=>      \$options{bed},
	);
	$tigr_tf->printUsageInfoAndExit() if (!$result);	
	
	########################################################################

	while(<>)
	{
		chomp;

		if(/^>(\S+)/)
		{
			$contigId=$1;
			$contigSeq{$contigId}="";
		}
		else
		{
			$contigSeq{$contigId}.=$_;
		}
	}



	########################################################################
	
	my @P=();
	my $scafSeq="";
	open(IN,$options{bed}) or $tigr_tf->bail("Cannot open or read file $options{bed}");
	while(<IN>)
	{
                #scaffold1	0	239     jcf7180062029523        239     -
                #scaffold1	416     706     deg7180045116455        290     -

		my @F=split;

		if(@P and $P[0] ne $F[0])
		{
			print ">$P[0]\n$scafSeq\n";
			$scafSeq="";
			@P=();
		}

		my $gapLen=0;
		$gapLen=$F[1]-$P[2] if(@P);
		$gapLen=20 if($gapLen<0);

		$scafSeq.='N'x$gapLen if($gapLen);
		
		my $contigSeq=$contigSeq{$F[3]};
		if($F[5] eq "-")
		{
                      	$contigSeq=reverse($contigSeq);
                       	$contigSeq=~tr/ACGTacgt/TGCAtgca/;
		}
		$scafSeq.=$contigSeq;

		@P=@F;
	}
	
	if(@P)
	{
		print ">$P[0]\n$scafSeq\n";
	}
			
	exit 0;
}
