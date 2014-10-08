#!/usr/bin/perl -w
                                                                                                                                                             
use strict;
use warnings;

while(<>)
{
	my @f=split;

	my $n=scalar(@f);
	next unless($n);

	print ">";
	print join " ",@f[0..$n-2];
	print "\n";

	print $f[$n-1];
	print "\n";
}

