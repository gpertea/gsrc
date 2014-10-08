#!/bin/tcsh -f
if ($1"x" == "x") then
 echo "Usage: run_glimmerhmm.csh <genomic_fasta> "
 echo ""
 echo "If a full path to the <genomic_fasta> file is given"
 echo "it will cd to that directory first"
 exit
 endif

set fa=$1
#-- change these here if you need different parameters
set fs=/szfs
set ghmmbin=$fs/szgenefinding/devel/GlimmerHMM/bin/glimmerhmm_linux
set cktrain=$fs/szgenefinding/devel/GlimmerHMM/trained_dir/chicken/train09_09
#--
set fname=$fa:t
set fb=$fname:r
set fadir=$fa:h
echo "Switching to $fadir"
if ($fadir =~ /*) cd $fadir
if ( -s $fname ) then
 $ghmmbin $fname $cktrain -f -g -o $fb.ghmm.gff3
else
 echo "Error: file $fadir/$fname does not exist!\n"
endif
