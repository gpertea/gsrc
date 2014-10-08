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
	my %opt;	
	my @colMax;
	#$opt{sep}='\t';

	# Configure TIGR Foundation
	$tigr_tf->setHelpInfo($HELPTEXT.$MOREHELP);
        $tigr_tf->setUsageInfo($HELPTEXT);
        $tigr_tf->setVersionInfo($REV);
        $tigr_tf->addDependInfo(@DEPENDS);
	
	# validate input parameters
	my $result = $tigr_tf->TIGR_GetOptions(
		"i=s"	=>	\$opt{i},
		"sep=s"	=>	\$opt{sep},
		"skip=i" =>     \$opt{skip},
		"o"	=>	\$opt{o},
		"integers" =>   \$opt{integers},
	);
	$tigr_tf->printUsageInfoAndExit() if (!$result);	

	###########################################################################
	
	my @lines=<>;
	if($opt{skip})
	{
		foreach (0..$opt{skip}-1)
		{
			print shift @lines;
		}
	}
	
	foreach (@lines)
	{		
		next if(/^$/);	
		next if(/^\//);
		next if(/^#/);

		s/^\s+//;
		
		my @f;
		if(defined($opt{sep})) { @f=split /$opt{sep}/; }
		else                   { @f=split /\s+/; }
		
		my $n=scalar(@f);	
		$n=$opt{i} if(defined($opt{i}));	
				
		foreach my $i (0..$n-1)
		{
			$f[$i]=int($f[$i]+0.5) if($opt{integers} and $f[$i]=~/(.*\d+)\.(\d+)$/);
			
			my $colLen=length($f[$i]);
			$colMax[$i]=$colLen if(!defined($colMax[$i]) or $colLen>$colMax[$i]);	
		}	
	}

	$colMax[0]+=2;	
	
	###########################################################################
	
	foreach (@lines)
	{
		if(/^#/) { print $_; next; }
		chomp;
		s/^\s+//;

		my @f;		
                if(defined($opt{sep})) { @f=split /$opt{sep}/; }
                else                   { @f=split /\s+/; }

		my $n=scalar(@f);		
		$n=$opt{i} if(defined($opt{i}));
		
		print "  " if($opt{o});
		
		foreach my $i (0..$n-1)
		{
			$f[$i]=int($f[$i]+0.5) if($opt{integers} and $f[$i]=~/(.*\d+)\.(\d+)$/);
			printf("%-".$colMax[$i]."s  ",$f[$i]);			
		}				

		if(defined($opt{i}))
		{
			print join " ",@f[$opt{i} .. @f-1];
		}			
	
		print "\n";
	}
		
	###########################################################################
			
	exit 0;
}
