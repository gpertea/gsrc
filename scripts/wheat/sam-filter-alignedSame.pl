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
                "min_dist=i"            => \$options{min_dist},
                "max_dist=i"            => \$options{max_dist},		
		"dist=i"		=> \$options{dist},
		"RF"			=> \$options{RF},
		"FR"                    => \$options{FR},
		"FF"                    => \$options{FF},
		"RR"                    => \$options{RR},
	);
	$result or $base->bail("Command line parsing failed");
			   
	if(defined($options{dist}) and  !defined($options{min_dist}) and  !defined($options{max_dist}) )	
	{
		$options{min_dist}=0.7*$options{dist};
		$options{max_dist}=1.3*$options{dist};
	}

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
			@P=@F;

                        next if($F[1] & 0x4)   ;   #not mapped
			next if($F[1] & 0x8)   ;   #next segment unmapped
                        next if($F[1] & 0x100) ;   #skip secondary alignments
			next if($F[6] ne "=");

			next if(defined($options{min_dist}) and abs($F[8])<$options{min_dist});
			next	if(defined($options{max_dist}) and abs($F[8])>$options{max_dist});


		        if(defined($options{RF})) # mp
        	        {
				next unless($F[1] & 0x10 and !($F[1]&0x20) and $F[3]<=$F[7]  or !($F[1] & 0x10) and ($F[1]&0x20) and $F[3]>=$F[7]); 
                	}
                        elsif(defined($options{FR})) # pe
                        {
                                next unless($F[1] & 0x10 and !($F[1]&0x20) and $F[3]>=$F[7]  or !($F[1] & 0x10) and ($F[1]&0x20) and $F[3]<=$F[7]); # ???
                        }
			elsif(defined($options{FF})) # pe
                        {
                                next unless(!($F[1]&0x10) and !($F[1]&0x20)); 
                        }
                       elsif(defined($options{RR})) # pe
                        {
                                next unless(($F[1]&0x10) and ($F[1]&0x20));
                        }

			
			print $_ ;
		}
	}
     					
	exit 0;
}

