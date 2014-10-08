#!/usr/bin/perl -w

use strict;

use File::Basename;
use AMOS::AmosFoundation;
use AMOS::AmosLib;

my $base = new AMOS::AmosFoundation();    
die("ERROR: problems creating an AMOS::Foundation object\n") unless($base);

my $VERSION = '$Revision: 1.2 $ ';

my $PRG = fileparse($0);

my $HELP = qq~
Program that ...
	
	USAGE: 
		$PRG
	
	INPUT:
		
		options:
		  ...
		
	OUTPUT:
	
	EXAMPLE:
		input:
		  ...
		
		output:
		  ...
	
	RETURN CODES:   
		0 - on success, 1 - on failure.
		
SAM format
1 QNAME String  : Query template NAME
2 FLAG Int [0,2^16-1] bitwise FLAG
3 RNAME String  : Reference sequence NAME
4 POS Int [0,2^29-1] 1-based leftmost mapping POSition
5 MAPQ Int [0,2^8-1] MAPping Quality
6 CIGAR String \*|([0-9]+[MIDNSHPX=])+ CIGAR string
7 RNEXT String : Ref. name of the mate/next segment
8 PNEXT Int [0,2^29-1] Position of the mate/next segment
9 TLEN Int [-2^29+1,2^29-1] observed Template LENgth
10 SEQ String : segment SEQuence
11 QUAL String : ASCII of Phred-scaled base QUALity+33

FLAG
0x1 template having multiple segments in sequencing
0x2 each segment properly aligned according to the aligner
0x4 segment unmapped; Bit 0x4 is the only reliable place to tell whether the segment is unmapped. 
0x8 next segment in the template unmapped
0x10 SEQ being reverse complemented
0x20 SEQ of the next segment in the template being reversed
0x40 the 1st segment in the template
0x80 the last segment in the template
0x100 secondary alignment
0x200 not passing quality controls
0x400 PCR or optical duplicate

~;

$base->setVersion($VERSION);
$base->setHelpText($HELP);
$|=1;

###############################################################################
#             
# Main program       
#
###############################################################################

MAIN:         
{
	my %options;
	$options{min_align}=0;
	$options{min_score}=0;

	my $result = $base->getOptions(	
		"nh"	        => \$options{no_header} , 
		"min_align=i"	=> \$options{min_align},
		"min_score=i"   => \$options{min_score},
		"max_mismatch=i"	 => \$options{max_mismatch},
		);
	$result or $base->bail("Command line parsing failed");
			   
	###########################################

	while(<>)
	{
		if(/^@/) 
		{
			print $_ unless($options{no_header});
		}
		else
		{
			next if(/^\[/);		
        	        my @F=split /\t/;	
			
			my $skip=0;

			$skip=1 if($F[1] & 0x4)   ;   #no alignment
			$skip=1 if($F[1] & 0x100) ;   #secondary alignment
	
			my $align=length($F[9]);
			$align-=$1 if($F[5]=~/^(\d+)S/);
			$align-=$1 if($F[5]=~/(\d+)S$/);

			$skip=1 if($align<$options{min_align});
                        $skip=1 if($F[4]<$options{min_score});
			$skip=1 if(defined($options{max_mismatch}) and /NM;i;(\d+)/ and $1>$options{max_mismatch});
			
			print $_ unless($skip);
		}
	}
     					
	exit 0;
}
