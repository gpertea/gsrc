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

Usage: $PRG fasta_file [options]
	  
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


###############################################################################
#
# Main program
#
###############################################################################


MAIN:
{
	my %options;
	my (%seq,%len,$id,@ids,%info);

	# Configure TIGR Foundation
	$tigr_tf->setHelpInfo($HELPTEXT.$MOREHELP);
        $tigr_tf->setUsageInfo($HELPTEXT);
        $tigr_tf->setVersionInfo($REV);
        $tigr_tf->addDependInfo(@DEPENDS);
	
	# validate input parameters
	my $result = $tigr_tf->TIGR_GetOptions(
		"rev"	=>	\$options{rev},
		"newId"	=>	\$options{newId},
		"prefix=s"	=>      \$options{prefix},
		"suffix=s"	=>	\$options{suffix},
	);
	$tigr_tf->printUsageInfoAndExit() if (!$result);	
	
	##############################

        while (<>)
	{
		if(/>(\S+)(.*)/)
		{
			$id=$1;
			$seq{$id}="";
			$info{$id}=$2;
		}
		else
		{
			$seq{$id}.=$_;
			chomp;

			$len{$id}+=length($_);
		}
	}

	##############################

	@ids=sort {$len{$a}<=>$len{$b}} keys %len;
	@ids=reverse(@ids) if($options{rev});

	my $count=1;
	foreach $id (@ids)
	{
		print ">";
		if(defined($options{prefix})) { print $options{prefix},"." }
		if($options{newId})  { print "$count"; $count++; }
		else                 { print "$id"; }
		if(defined($options{suffix})) { print ".",$options{suffix} }

		print $info{$id} if(defined($info{$id}));
		
		print "\n",$seq{$id},"\n";
	}

	exit 0;
}

