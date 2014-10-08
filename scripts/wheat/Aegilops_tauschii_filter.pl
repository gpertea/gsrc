#!/usr/bin/perl -w
 
use strict;
use warnings;

# TIGR Modules
use TIGR::Foundation;

my $tigr_tf = new TIGR::Foundation;
my $PRG = $tigr_tf->getProgramInfo('name');
my $VERSION="1.0";
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

# Configure TIGR Foundation
$tigr_tf->setHelpInfo($HELPTEXT.$MOREHELP);
$tigr_tf->setUsageInfo($HELPTEXT);
$tigr_tf->setVersionInfo($VERSION);
$tigr_tf->addDependInfo(@DEPENDS);

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	# define variables
	my %opt;
	$opt{min}=127;

	# validate input parameters
	my $result = $tigr_tf->TIGR_GetOptions(
		"filter=s"	=> \$opt{filter},
		"min=i"		=> \$opt{min},
		"max=i"		=> \$opt{max},
		"k=i"		=> \$opt{k}
	);
	$tigr_tf->printUsageInfoAndExit() if (!$result);	
	$tigr_tf->printUsageInfoAndExit() if (!$opt{filter});
	$tigr_tf->printUsageInfoAndExit() if (!$opt{max});
	$tigr_tf->printUsageInfoAndExit() if (!$opt{k});

	###############################################

	my %scf;
	my %pos;
	my %end;

        open(FILTER,$opt{filter}) or die $!;
        while(<FILTER>)
        {
                my @f=split;
		if($f[2]>=$opt{min} and $f[2]<=$opt{max})
		{
	                $scf{$f[0]}=$f[1];
			$pos{$f[0]}=$f[2];
			$end{$f[0]}=$f[-1];
		}
        }
        close(FILTER);

	###############################################

	my ($id1,$id2,$seq1,$seq2,$qlt1,$qlt2);

        open(IN1,$ARGV[0])     or die $!;
        open(IN2,$ARGV[1])     or die $!;
        open(OUT1,">$ARGV[2]") or die $!;
        open(OUT2,">$ARGV[3]") or die $!;

        while(<IN1>)
        {
                $id1=$_;   
		$seq1=<IN1>;
		<IN1>;
		$qlt1=<IN1>;

                $id2=<IN2>;              
                $seq2=<IN2>;
                <IN2>;
                $qlt2=<IN2>;

		my $subseq1=substr($seq1,0,$opt{k});
		my $subseq2=substr($seq2,0,$opt{k});

		if($scf{$subseq1} and $scf{$subseq2} and $scf{$subseq1} ne $scf{$subseq2}) 
		{
			if($pos{$subseq1}+$pos{$subseq2}<=$opt{max}*1.2)
			{
				$id1=~/^@(\S+)/; $id1="@".$1." ".$scf{$subseq1}.".".$end{$subseq1}." ".$pos{$subseq1}."\n";
				$id2=~/^@(\S+)/; $id2="@".$1." ".$scf{$subseq2}.".".$end{$subseq2}." ".$pos{$subseq2}."\n";
				print OUT1 $id1,$seq1,"+\n",$qlt1;
				print OUT2 $id2,$seq2,"+\n",$qlt2;
			}
		}
	}
	close(IN1);
	close(IN2);
        close(OUT1);
        close(OUT2);

	exit 0;
}
