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
	my %count;

	my (@P,$P,%len);
	my $result = $base->getOptions(
		"min_dist=i"		=> \$options{min_dist},
		"max_dist=i"            => \$options{max_dist},
		"min_refLen=i"  	=> \$options{min_refLen},
		"min_score=i"           => \$options{min_score}, 
		"RF"			=> \$options{RF}
	);
	$result or $base->bail("Command line parsing failed");

	###########################################

	#issue: short ctg
	#SRR124243.1549598	65	contig00001	2424	39	89M	contig00007	217	0	ATTTTCCTGAATTTGAAAATCTTAGTACTCAAAACGTGCGACCAATAGTGATCGTATTGAACTACGGTCACATGTATTTAGGAAAGATG	GGGGGGGGGGGGGG=GGGFGGGFGGFFFAGED?EEE?EEACC?>C=C@>CD=AEE5;::<AA5?B-5C:C>,B;>@4@@445<7<287+	NM:i:0	AS:i:89	XS:i:0
	#SRR124243.1549598	129	contig00007	217	40	1H90M	contig00001	2424	0	AGCTTCAAGCCTTTCGGAATAGTGCAAACTGCACTGCACATAGCTCCGTGCAGTCTACACTATTCCTCAAGGCTTAAAGCTAAGCAACGT	D=DAD?DDB5==>B=5A:C=EBCAEEAEEB=-CCC5?:A5?AAABDADBA=ABB:E?A-A>.:=;C>A?C-AC:5@;B>>C=ECC+7@%%	NM:i:0	AS:i:90	XS:i:0

	while(<>)
	{
		my @F=split /\t/;

		if (/^\@/) 
		{	
			$len{$1}=$2 if(/\@SQ\s+SN:(\S+)\s+LN:(\d+)/);
			#print $_;
			next;
		}	
			
		next if($F[1] & 0x4);		#no alignment
		next if(defined($options{min_refLen}) and $len{$F[2]}<$options{min_refLen});
                next if(defined($options{min_score}) and $F[4]<$options{min_score});

		my ($dist,$dist5,$dist3);		

		$dist5=$F[3];
                $dist3=$len{$F[2]}-$F[3]-length($F[9]);
		$dist3-=$1 if($F[5]=~/^(\d+)S/);
		$dist3-=$1 if($F[5]=~/(\d+)S$/);

		if(defined($options{RF})) # mp
		{
			if($F[1]&0x10) { $dist=$dist3 }
			else           { $dist=$dist5 }
		}
		else
		{		
	                if($F[1]&0x10) { $dist=$dist5 }
                       	else           { $dist=$dist3 }

		}

		next if(defined($options{max_dist}) and $dist>$options{max_dist});
		#next if(defined($options{min_dist}) and $dist5<$options{min_dist} and $dist3<$options{min_dist}); ???
		next if(defined($options{min_dist}) and $dist<$options{min_dist});

		chomp;
		print $_,"\tCO:Z:$dist\n";
	}
     					
	exit 0;
}

