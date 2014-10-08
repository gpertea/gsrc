#!/bin/tcsh -f


set cdb=refseq.mammalian_np.CDS.fa.cidx

if ($1"x" == "x") then
 echo "Usage: prep_sim4.csh fasta_files.."
 echo ""
 echo "Must be run from the sim4search subdirectory"
 echo "and expects a gblat_*.gff and/or gmap_*.gff3 file with gblat scanning results"
 echo "and $cdb in the current directory to pull RefSeqs for sim4cc"
 echo ""
 echo "Example: prep_sim4.csh ../ga_??.fa"
 exit 1
 endif

# prepare the name translation table

set f=$1;

if (! -f $cdb ) then
 echo "Error: missing file in current directory: $cdb"
 exit 1
endif

set cdir=`pwd`
set pdir=$cdir:h
set gblat_gff=`ls gblat*.gf{f,f3}`

if (! -f $gblat_gff ) then
 echo "Error: no single gblat*.gf{f,f3} file found in current directory"
 echo "make sure you're in <organism>/sim4search and the gblat results are there" 
 exit 1
endif
set infiles=$gblat_gff

set gmap_gff=`ls gmap*.gf{f,f3}`
if ( -f $gmap_gff ) then
 #echo "Error: no single gmap*.gf{f,f3} file found in current directory"
 #echo "make sure you're in <organism>/sim4search and the gblat/gmap results are there" 
 #exit 1
 set infiles="$infiles $gmap_gff"
endif

set fp=$f:h
set ftx="$fp/seqn2fname.lst"
set dlst=""
set flst=""
/bin/rm -f $ftx
foreach f ( $* )
 set fn=$f:t
 set bd=$fn:r
 set sn=`head -1 $f | cut -f1 -d ' ' | cut -b2-`
 #if (! -d $bd) mkdir $bd
 echo "$sn\t$bd" >> $ftx
 #/bin/rm -f $bd/$bd.refseq.gmap.gff3
 #/bin/rm -f $bd/$bd.unipr.pmap.gff3
 set flst = "$flst $fn"
end
echo "$fp/seqn2fname.lst created."
cat $infiles | trcol.pl -G $ftx | sort -u > anchored.lst
set simbat=sim4_grid.cmds
/bin/rm -f $simbat
gawk '{ print "qsim4cc -p75 -c75 -o wrk_sim4cc.gff3 -q \047"$2"\047"" '$cdir/$cdb' '$pdir/'"$1".fa"}' \
 anchored.lst > $simbat

echo "$simbat file created, use gridx with it, e.g.:"
echo "gridx -q -N -O grdsim4logs -f sim4_grid.cmds -p40 -m gpertea"
