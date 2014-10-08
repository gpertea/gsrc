#!/usr/bin/perl -w

use strict;

use AMOS::AmosFoundation;
use AMOS::AmosLib;
use AMOS::ParseFasta;
use Digest::MD5 qw(md5 md5_hex md5_base64);

my $base = new AMOS::AmosFoundation();    
die("ERROR: problems creating an AMOS::Foundation object\n") unless($base);

my $VERSION = '$Revision: 1.2 $ ';

my $HELP = q~
	Program that splits format scaffolds into contigs. 
	Input/Output data in FASTA format.
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
	$options{flank}=0;

	my $result = $base->getOptions(
		"bed"	=> \$options{bed},
		"flank=i" => \$options{flank}
	);

	####################		
		
        my $fh = new AMOS::ParseFasta(\*STDIN);
        $fh or $base->bail("Fasta format expected in $fh: $!\n");
	
	while (my ($head, $seq) = $fh->getRecord())
    	{       
		my @head=split /\s+/,$head;
		$head[0]=$1 if($head[0]=~/^\S*\|(\w+)\.\d+/);

		my $id=1;
		my $offset=0;
		
		while($seq=~/([ACGT]+)(N+)([ACGT]+.*)/i)
		{
			my ($five,$gap,$three)=($1,$2,$3);
			
			$offset+=length($five);
			my $len=length($gap);

			if($options{bed}) { print "$head[0]\t",$offset-$options{flank},"\t",$offset+$len+$options{flank},"\t$head[0].gap.$id\n" } 
			else              
			{ 
				if($options{flank})
				{
					$gap=substr($five,length($five)-$options{flank}).$gap.substr($three,0,$options{flank});
				}
				print ">$head[0].gap.$id $len\n$gap\n";                          
			}
	
			$offset+=length($gap);
			$seq=$three;			
			$id++;	
		}
	}
		        					
	exit 0;
}

