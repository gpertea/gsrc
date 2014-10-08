#!/usr/bin/perl -w
                                                                                                                                                             
use strict;
use warnings;

my $name;

while(<>)
{
	chomp;

	my @f=split;
	next unless(defined($f[0]));

	if(/^>(.+)/)	
	{ 
		print "\n" if(defined($name));
		
		$name=$1;
		print $1,"\t";
	}
	else	
	{ 
		print $_;
	} 
}

print "\n";
