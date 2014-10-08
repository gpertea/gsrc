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
Program that 

Usage: $PRG file [options]
	  
  INPUT:   
  
  file:	
  
  options:

  	-min		- min
	-max		- max
	-i		- col	
	-j		- col [optional]
	-t		- description
	-g 		- genome size (for N50 computation)
	-Esize          - show estimated size
	-nh		- no header
	-all		- account for all

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
$tigr_tf->setVersionInfo($REV);
$tigr_tf->addDependInfo(@DEPENDS);

###############################################################################
#
# Main program
#
###############################################################################

MAIN:
{
	my $i=0;
	my $j;
	
	my $sum=0;
	my $sum2=0;
	my $sump;
	
	my $n50;	
	my @list;
	my %count;
	my @keys;
	
	my ($min,$max,$mean,$mode);
	my ($q1,$q2,$q3,$qMin,$qMax)=(0,0,0);
	#my $q90=0;
	
	my $text="";
	my $nh=0;
	my $zero;		
		
	my $count=0;

	my $genomeSize;
	my $IQR;
	my $len;
	my $Esize;
	my $all;
	my $nl;	

	# validate input parameters
	my $result = $tigr_tf->TIGR_GetOptions(
		"min=s"		=>	\$min,
		"max=s"		=>	\$max,
		"i=i"		=> 	\$i,
		"j=i"		=> 	\$j,
		"text=s"	=>	\$text,
		"nheader"	=>	\$nh,
		"nline"		=>	\$nl,
		"genomeSize=i"	=>	\$genomeSize,
		"IQR"		=>	\$IQR,	
		"Esize"		=>      \$Esize,
		"all=i"		=>	\$all
	);
	$tigr_tf->printUsageInfoAndExit() if (!$result);
	
	$text=$ARGV[0] if(!$text and defined($ARGV[0]));

	# parse input
	while(<>)
	{
		next if(/^$/ or /^#/ or /^@/);
		chomp;
		my @f=split;
		next unless(@f);
		
		$tigr_tf->bail("Undefined column $i in $_") if(!defined($f[$i]));
		$tigr_tf->bail("Undefined column $j in $_") if(defined($j) and !defined($f[$j]));
		
		my $key=$f[$i];				
		my $val=1;
		$val=$f[$j] if(defined($j));		

		next unless(is_integer($key) or is_float($key));
		next unless(is_integer($val) or is_float($val));
		
		next if(defined($min) and $key<$min);
		next if(defined($max) and $key>$max);
		
		$count{$key}+=$val;
		$count+=$val;		
		$sum+=$key*$val;
	}
	
	# get mean	
	exit 0 unless($count) ;

	if(defined($all) and $count<$all)
	{
		$count{0}+=$all-$count;
		$count=$all;
	}
		
	$mean=$sum/$count;
	
	# get all keys	
	@keys=sort {$a <=> $b} keys %count;
	$min=$keys[0];
	$max=$keys[-1];
					
	# get q1
	$sump=0;
	foreach my $key (@keys)
	{									
		#if($sump<=($count/10))  { $q90=$key; }
		if($sump<=($count/4))   { $q1=$key; }	
		if($sump<=($count/2))   { $q2=$key; }				
		if($sump<=($count*3/4)) { $q3=$key; }	
		else                    { last;     }
			
		$sump+=$count{$key};

		$mode=$key if(!defined($mode) or $count{$mode}<$count{$key});
	}
		
	if($IQR)
	{				
		$qMin=$q1-1.5*($q3-$q1);
		$qMax=$q3+1.5*($q3-$q1);
	}
	
	# get n50 	
	$genomeSize=$sum unless($genomeSize);
	if($genomeSize)
	{
		$genomeSize=$sum unless($genomeSize);

		$sump=0;
		foreach my $key (reverse @keys)
		{
			$n50=$key;
			$sump+=$key*$count{$key};
			last if($sump>$genomeSize/2);
		}

		$n50=0 if($sump<$genomeSize/2);
	}

	if($genomeSize and $Esize)
	{
		$sump=0;
                foreach my $key (@keys)
               	{
                        $sump+=$key*$key*$count{$key};
               	}

		$Esize=($sump/($genomeSize));
	}

	unless($nh)
	{
		printf ("%-20s ",".")        if(defined($text));			
		printf ("%-10s ","elem");
		
		printf ("%-10s ","<=$zero")  if(defined($zero));			 		
		printf ("%-10s ",">$zero")   if(defined($zero)); 
        	
		printf ("%-6s ",  "min");
                #printf ("%-6s ",  "q90");
		printf ("%-6s ",  "q1");
		printf ("%-6s ",  "q2");
		printf ("%-6s ",  "q3");
		printf ("%-10s ",  "max");
		printf ("%-10s ",  "mean");  
                #printf ("%-10s ",  "mode")   if(defined($mode));              
		printf ("%-15s",   "sum");
               printf ("%-10s ",  "n50")    if(defined($n50));
		
		printf ("%-10s ",  "qMin")   if(defined($IQR));		 
		printf ("%-10s ",  "qMax")   if(defined($IQR));		 						
		printf ("%-8s ",  "Esize")   if(defined($Esize));

		print "\n";
	}
	
	printf ("%-20s ",$text)              if(defined($text));	
	printf ("%-10s ",$count);
	
	printf ("%-6s ",   $min);
        #printf ("%-6s ",   $q90);
	printf ("%-6s ",   $q1);
	printf ("%-6s ",   $q2);
	printf ("%-6s ",   $q3);
	printf ("%-10s ",  $max);
	#printf ("%-10.2f ",$mean);
    printf ("%-10d ",$mean);
        #printf ("%-10s ",  $mode)            if(defined($mode));
	printf ("%-15s",   $sum);
        printf ("%-10s ",  $n50)             if(defined($n50));
	printf ("%-10.2f ",  $qMin)          if(defined($IQR));		 
	printf ("%-10.2f ",  $qMax)          if(defined($IQR));		 
        #printf ("%-8.6f ",   $Esize)        if(defined($Esize));
	printf ("%-8i ",   $Esize)           if(defined($Esize));
	
	print "\n" unless(defined($nl));
		
	exit 0;
}

####################################################

sub is_integer {
   defined $_[0] && $_[0] =~ /^[+-]?\d+$/;
}

sub is_float {
   defined $_[0] && $_[0] =~ /^[+-]?\d+(\.\d+)?$/;
}
