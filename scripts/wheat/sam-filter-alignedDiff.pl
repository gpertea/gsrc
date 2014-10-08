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

###############################################################################
#             
# Main program       
#
###############################################################################

MAIN:         
{
	my %options;
	
	my (@P,$P,%len);
	my $result = $base->getOptions(
		"min_score=i"	=> \$options{min_score},
		"prefix|s"	=>  \$options{prefix}
	);
	$result or $base->bail("Command line parsing failed");
			   
	###########################################

	while(<>)
	{
		if (/^\@/) 
		{	
			print $_;
		}	
                else 
		{		
        	        my @F=split /\t/;
	
			next if(! $F[1] & 0x1 || $F[1] & 0x2 || $F[1] & 0x4 || $F[1] & 0x100 || (@P and $P[1] & 0x4)) ;
			next if(defined($options{min_score}) and $F[4]<$options{min_score}) ; 
		        if(@P and $P[0] eq $F[0] and $P[1] & 0x40 and $F[1] & 0x80 and ($P[2] ne $F[2]))
			
			{

				if(!$options{prefix})
				{
					print $P,$_ ;
				}
				else
				{
					$P[2]=~/^(.+?)\./; my $pprefix=$1;
					$F[2]=~/^(.+?)\./; my $fprefix=$1;
				
					if($pprefix eq $fprefix)
					{
						 print $P,$_ ;
					}
				}
			}

			$P=$_;
			@P=@F;
			
		}
	}
     					
	exit 0;
}

