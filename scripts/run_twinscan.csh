#!/bin/tcsh -f
if (x$1 == 'x') then
  echo "Usage: run_twinscan.csh fasta_seq"
  exit 1
endif
set plat=`uname -p`
set infile=$1
set outfile=$infile.tsgff
set tsdir='fs/sz-user-supported/Linux-'$plat
#set tsdir=fs/szannotation/$plat/bin
if (-d "/sz$tsdir") then
 set tsbin="/sz$tsdir/twinscan"
else
 set tsdir="/$tsdir/twinscan"
endif

# - hard coded parameter file
set hmm=$tsdir/parameters/human_iscan-9993-genes-09-13-2004.zhmm

$tsdir/bin/iscan $hmm $infile > $outfile
