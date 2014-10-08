#!/bin/tcsh -f

#set pmap = 'unipr.blat.gff3'
set pmap = 'pexo_exon.gff3'
#set gmap = 'refseq_exon.gff3'
set gmap = 'sim4.gff3'
set ptrack='p2g_exo'
set mtrack='sim4cc'
set protdb = '~/ann/protdb/unipr_mammals.fa.cidx'

if ($2'x' == 'x') then
  echo "Usage: rebuild_ann.csh [<step>] <jigsaw_suffix> <dirs..>"
  echo "  Must be run in the '.../encode/<species>/jigsaw' directory "
  echo "  where the <dirs..> are"
  echo "  <jigsaw_suffix>.jgff is/was used as jigsaw output suffix"
  echo " Example: "
  echo "   rebuild_ann.csh 1 wpmap ra_0?"
  echo ""
  echo '  WARNING: make sure variables $pmap, $ptrack, $gmap and $mtrack are correct!'
  echo "       (currently they are: $pmap (track: $ptrack) and $gmap (track: $mtrack))"
  echo "<step> is an optional numeric argument denoting the stage at which"
  echo " reprocessing should start (default is 2):" 
  echo "   0: rerun jigsaw with <dirs> and evidence_<jigsaw_suffix> "
  echo "   1: combine the 3 jigsaw output file for each <dirs..>"
  echo "   2: (default) annotate jigsaw results from <dir>.recon_<jigsaw_suffix>.jgff.gff3"
  echo "   3: build the submit & gview files from <dir>.ann.recon_<jigsaw_suffix>.jgff.gff3"
  exit 1
  endif

set jstart=$1
set jstageinfo="annotation of predicted gene models"
set stest=`echo "$jstart" | egrep -c '^[0-4]$'`
if ( ${%jstart} > 1 || $stest < 1 ) then
  set jstart=2
else
  shift
  if ($jstart == 0) set jstageinfo="running jigsaw 3 times"
  if ($jstart == 1) set jstageinfo="combine the 3 jigsaw runs"
  if ($jstart == 3) set jstageinfo="generate submit and gview files"
endif
set jsuf=$1
set jsufbase=$jsuf
shift

set gvdir='/mylocal/geo/httpd/html/data/enc'
if ( ! -d $gvdir ) then
echo "Error: $gvdir not found!"
exit 1
endif

set startdir=$PWD  #this doesn't resolve simbolic links
set basedir=`pwd` 
set updir=$basedir:h # should be ~/ann/encode/spdir (remove ./jigsaw)
set spdir=$updir:t #species directory name

set tgtdir=$gvdir/$spdir
if ( -d $tgtdir"_$jsufbase" ) then
   set tgtdir=$tgtdir"_$jsufbase"
 endif
if ( ! -d $tgtdir) then
 mkdir -p $tgtdir
 endif

set map=$updir/fa2encode.map
if ( ! -f $map ) then
  echo "Error: $map file not found!"
  exit 1
endif


echo "# Starting from stage $jstart ($jstageinfo)"
if ($jstart < 2 ) then
  #locate jigsaw's evidence file
  set evf=`ls {ev,weight}*_$jsuf.jigsaw exp*/{exp,weight}*_$jsuf.jigsaw`
  if ( $%evf > 0 ) then
    set evf=$basedir/$evf[1]
    echo "Using jigsaw evidence file: $evf"
  else
    echo "Error: no evidence file found $evf"
    exit 1
  endif
  #build the dirlist
endif
set jsuf="$jsuf.jgff"
#echo "Jigsaw file suffix: $jsuf"

foreach d ( $argv ) # { <------- for each directory
 echo "====================> $d <==================="
 if ( ! -d $d ) then
   echo "ERROR: directory $d not found!"
   exit 2
 endif

set j_gff = $d.recon_$jsuf.gff3
set j_anngff = $d.ann.recon_$jsuf.gff3
 
if ($jstart < 2) then
  #running linear jigsaw here.. change this if you want something else
  #prepare dirlist
  set jdirf=dirlst-$d.jigsaw
  echo "$basedir/$d" > $jdirf
  if ($jstart < 1) then
   echo "$basedir/${d}f" >> $jdirf
   echo "$basedir/${d}r" >> $jdirf
   set jcmd="run_jigsaw.pl -e $evf -l $jdirf -L -o $jsuf"
   echo "..running: $jcmd"
   $jcmd
  endif
  #-- now recombine
  set df=$d'f'
  set dr=$d'r'
  recon_jigsaw.pl -P $d/$d.fa -o $j_gff $d/$d.$jsuf $df/$df.$jsuf $dr/$dr.$jsuf
  #validate the file one more time:
  set badcds=$d.recon_$jsuf.bad_JIGSAW.gtf
  gffilter -g $d/$d.fa -b $badcds $j_gff
  if (-s $badcds) then 
     echo "**** ERROR found at CDS validation for $j_gff! (see $badcds)"
   else
     /bin/rm -f $badcds
   endif
endif # jstart is 0 or 1

#stage 2: annotation
if ($jstart < 3) then 
  if (-s $j_gff) then
   /bin/rm -f $d.[pg]map.i{it,fa}
   #gff2iit -o $d.pmap -t $ptrack $d/$d.$pmap
    gtf2gff -t $ptrack $d/$d.$pmap  > $d.iit_pmap.gff3
    gtf2gff -t 'pmap'  $d/$d.pmap_exon.gff3 >> $d.iit_pmap.gff3
    gff2iit -o $d.pmap $d.iit_pmap.gff3
   #gff2iit -o $d.gmap -t $mtrack $d/$d.$gmap
   gtf2gff -t $mtrack $d/$d.$gmap | perl -pe 's/\.mrna(\d+)/\.m$1/' > $d.iit_gmap.gff3
   gtf2gff -t 'gmap'  $d/$d.gmap_exon.gff3  >> $d.iit_gmap.gff3
   gff2iit -o $d.gmap -t $mtrack $d.iit_gmap.gff3

   #echo "gffann.pl -o $j_anngff -t $d -m $d.gmap.iit -p $d.pmap.iit -P $protdb $j_gff" 
   #gffann.pl -o $j_anngff -t $d -m $d.gmap.iit -p $d.pmap.iit -P $protdb $j_gff
   #echo "gffann.pl -o $j_anngff -t $d -m $d.gmap.iit -p $d.pmap.iit -P $protdb $j_gff"
   gffann.pl -o $j_anngff -t $d -m $d.gmap.iit -p $d.pmap.iit -P $protdb $j_gff
   #  -- temporary: throw in all the mapping evidence
   gtf2gff -t $mtrack $d/$d.$gmap  > $d.top5gmap.gff3
   gtf2gff -t $ptrack $d/$d.$pmap  > $d.top5pmap.gff3
   gtf2gff -t 'gmap'  $d/$d.gmap_exon.gff3  >> $d.top5gmap.gff3

   gtf2gff -t 'pmap'  $d/$d.pmap_exon.gff3 >> $d.top5pmap.gff3
   # gtf2gff -t 'blat'  $d/$d.pblat.gff3 >> $d.top5pmap.gff3
  else #no predictions; just show whatever mappings we got
   cp $j_gff $j_anngff
   gtf2gff -t $mtrack $d/$d.$gmap | perl -pe 's/\.mrna(\d+)/\.m$1/' > $d.top5gmap.gff3
   gtf2gff -t $ptrack $d/$d.$pmap  > $d.top5pmap.gff3
   gtf2gff -t 'gmap'  $d/$d.gmap_exon.gff3  >> $d.top5gmap.gff3
   gtf2gff -t 'pmap'  $d/$d.pmap_exon.gff3 >> $d.top5pmap.gff3
   # gtf2gff -t 'blat'  $d/$d.pblat.gff3 >> $d.top5pmap.gff3
  endif
endif # jstart was 0,1 or 2

# stage 3: prepare submit and gview files
set fo=$d:s/_//
if ( ! -d submit ) then
   mkdir submit
   endif
echo "jgff2tbl -f $d/$d.fa -A -o submit/$fo $j_anngff"
jgff2tbl -f $d/$d.fa -A -o submit/$fo $j_anngff
# -- now prepare and copy the gview files
if ( ! -f $gvdir/tracks.def ) then
 echo "Error: file $gvdir/tracks.def not found!"
 exit 3 
endif
  set m=`grep $d < $map | tr "\t" '~'`
  set m = $m:s/~/ /
  set ar=($m)
  set fn=$ar[1]
  if ($fn != $d) then
    echo "Error: grep for $d in $map didn't work properly!"
    exit 1
  endif
  set en=$ar[2]
  set fsub=$d:s/_//
  #echo "> preparing $en ($fn) for gview in $gvdir/$spdir/"
  set tgf=../$en.gv.gff3
  set bkm=../$en.bkm
  /bin/rm -f $tgf
  # set gnum=`cat $j_anngff | grep -c mRNA`
  # set gl=`cat $j_anngff | perl -ne 'print $1."," if m/GeneId=([\w\-]+)/'`
  # echo "$en\t$fn\t$gnum\t$gl" >> $tgtdir/$sm
  if ( -f submit/$fsub.tbl.gff3 ) then
   ~/ann/bin/fltgff4gv -o $tgf -b $bkm submit/$fsub.tbl.gff3
  else
   ~/ann/bin/fltgff4gv -t jigsaw -o $tgf -b $bkm $j_anngff
  endif
  #~/ann/bin/fltgff4gv -t gmap < $fn/$fn.refseq_exon.gff3 >> $tgf
  #~/ann/bin/fltgff4gv -t gmap < $fn/$fn.sim4exon.gff >> $tgf
  ~/ann/bin/fltgff4gv < $fn.top5gmap.gff3 >> $tgf
  ~/ann/bin/fltgff4gv < $fn.top5pmap.gff3 >> $tgf
  sed 's/GlimmerHMM/glimmerhmm/' < $fn/$fn.glimmerhmm_pdom.gff >> $tgf
  cat $fn/$fn.twinscan.gff >> $tgf
  echo "running: gff2gview -d $tgtdir -o $en -T $gvdir/tracks.def $fn/$fn.fa $tgf"
  gff2gview -d $tgtdir -o $en -T $gvdir/tracks.def $fn/$fn.fa $tgf
  if ( -f $bkm ) then
   cp $bkm $tgtdir
  endif
  #ls -l $tgtdir/$en.*
  # prepare the sequin files for reviewing
  #mkdir sequin
  ln -s ../submit/$fo.{tbl,fsa} sequin/
  tbl2asn -p sequin -t template.sequin
end # } for each directory

#-- rebuild the gene counts

cd $tgtdir
set gtab='gcount.tab'
/bin/rm -f $gtab
foreach bkm ( *.bkm )
 set gnum=`wc -l < $bkm`
 set gl=`cut -f2 -d '|' $bkm | cut -f1 | perl -ne 'chomp;print "$_," unless m/\d+_jsm\w+$/'`
 set en=$bkm:r
 set fn=`grep $en < $map | cut -f1`
 echo "$en\t$fn\t$gnum\t$gl" >> $gtab
end #for each gene bookmark file
cd $basedir
#echo "done."
exit 0

