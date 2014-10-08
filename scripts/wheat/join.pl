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
Program that joins two files based on the first column;
The files do not have to have the same number of lines 

Usage: $PRG file1 file2
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
	my (%h0,%h1);
	my ($cols0,$cols1);
	my ($missing0,$missing1);

	my %options;
	my $i=0;
	my $j=0;
	my ($i1,$j1,$jh);

	$options{empty}="";
	
	# Configure TIGR Foundation
	$tigr_tf->setHelpInfo($HELPTEXT.$MOREHELP);
        $tigr_tf->setUsageInfo($HELPTEXT);
        $tigr_tf->setVersionInfo($REV);
        $tigr_tf->addDependInfo(@DEPENDS);
	
	# validate input parameters
	my $result = $tigr_tf->TIGR_GetOptions(	
		"i=i"	=>	\$i,
		"i1=i"	=>	\$i1,
		"j=i"	=>	\$j,
		"j1=i"	=>	\$j1,
		"jh=s"	=>	\$jh,
		"all"	=>	\$options{all},
		"0"	=>	\$options{0},
		"1"     =>      \$options{1},
		"empty=s" =>	\$options{empty},
		"f=s"	=> 	\$options{filter}
	);
	$tigr_tf->printUsageInfoAndExit() if (!$result);
	
	###########################################################################

	if(defined($options{filter})) { open(FILTER,$options{filter}) or $tigr_tf->bail("Cannot open input file".$!) ;	}
	else                          { open(FILTER,$ARGV[1]) or $tigr_tf->bail("Cannot open input file".$!) ; }
	  
	if($jh)
	{
		my @f=split /\s+/,$jh;
		my $name=splice(@f,$j,1);
                $name.=".".splice(@f,$j1-1,1) if(defined($j1));

                $h1{$name}=join "\t",@f;
        }

	while(<FILTER>)
	{
    		chomp;
		next if(/^$/);
		next if(/^#/ and $.>1);
		
    		my @f=split;
    		die "ERROR: $_" if(@f<$j);
	
		$cols1=scalar(@f);

		my $name=splice(@f,$j,1);
		$name.=".".splice(@f,$j1-1,1) if(defined($j1));
		
    		$h1{$name}=join "\t",@f;
	}
	close(FILTER);

	$missing1=($options{empty}."\t") x ($cols1-1);
	
	###########################################################################
	
	if(defined($options{filter})) {}
	else                          { open(STDIN,$ARGV[0]) or $tigr_tf->bail("Cannot open input file".$!) ;	}
	
	# read 1st file
	while(<STDIN>)
	{
    		chomp;
		if(/^$/ or /^#/)
		{
			print $_,"\n";
			next;
		}
		
    		my @f=split;		
		die "ERROR: $_" if(@f<$i);

		$cols0=scalar(@f);

		my $name=$f[$i];
		$name.=".".splice(@f,$i1,1) if(defined($i1));
		
		if(defined($h1{$name}))             { print $_,"\t",$h1{$name},"\n"; }	
		elsif($options{all} or $options{0}) { print $_,"\t",$missing1,"\n"   }
	
		$h0{$name}=1;		
	}

	$missing0=($options{empty}."\t") x ($cols0-1);

	###########################################################################

	if($options{all} or $options{1})
	{
		foreach my $name (keys %h1)
		{
			if(!defined($h0{$name})) { print $name,"\t",$missing0,"\t",$h1{$name},"\n"; }
		}
	}
	
	exit 0;
}
