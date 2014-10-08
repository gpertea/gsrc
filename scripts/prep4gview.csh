#!/bin/tcsh -f
set gvdir='/mylocal/geo/httpd/html/data/enc'
set map='fa2encode.map'
if ($1'x' == 'x') then
 echo "Must be run from the encode base directory"
 echo ""
 echo "Usage: prep4gview.csh <jigsaw_output_suffix> <organism_dirs>.. "
 echo ""
 echo " e.g.: prep4gview.csh gblatp.jgff rabbit"
 echo "Please check if output directory is set properly"
 echo "(currently $gvdir)"
 #echo "Expected suffix for annotated jigsaw results: .$jsr"
 echo ""
 echo "Requires map file $map to be present under each <dir>"
 echo "Uses common track definition file: $gvdir/tracks.def"
 exit 1
endif

set jr=$1
shift

set jsr="ann.recon_$jr.gff3"
if ( ! -f $gvdir/tracks.def ) then
 echo "Error: file $gvdir/tracks.def not found!"
 exit 3 
endif
set sm='gcount.tab'
foreach d ( $* )
 if ( ! -d $d ) then
  echo "Error: no directory $d found!"
  exit 2
 endif
 cd $d
 echo "..processing $d"
 if ( ! -d $gvdir/$d) then
   mkdir -p $gvdir/$d
   endif
 /bin/rm -f $gvdir/$d/$sm
 foreach m ( `cat $map | tr "\t" '~'` )
  set m = $m:s/~/ /
  set ar=($m)
  set fn=$ar[1]
  set en=$ar[2]
  set fsub=$fn:s/_//
  echo "> processing $en ($fn)..."
  set tgf=$en.gv.gff3
  set bkm=$en.bkm
  /bin/rm -f $tgf
  set gnum=`cat jigsaw/$fn.$jsr | grep -c mRNA`
  set gl=`cat jigsaw/$fn.$jsr | perl -ne 'print $1."," if m/GeneId=([\w\-]+)/'`
  echo "$en\t$fn\t$gnum\t$gl" >> $gvdir/$d/$sm
  if ( -f jigsaw/submit/$fsub.tbl.gff3 ) then
   ~/ann/bin/fltgff4gv -o $tgf -b $bkm jigsaw/submit/$fsub.tbl.gff3
  else
   ~/ann/bin/fltgff4gv -t jigsaw -o $tgf -b $bkm jigsaw/$fn.$jsr
  endif
  #~/ann/bin/fltgff4gv -t gmap < jigsaw/$fn/$fn.refseq_exon.gff3 >> $tgf
  #~/ann/bin/fltgff4gv -t gmap < jigsaw/$fn/$fn.sim4exon.gff >> $tgf
  ~/ann/bin/fltgff4gv < jigsaw/$fn.top5gmap.gff3 >> $tgf
  ~/ann/bin/fltgff4gv < jigsaw/$fn.top5pmap.gff3 >> $tgf
  sed 's/GlimmerHMM/glimmerhmm/' < jigsaw/$fn/$fn.glimmerhmm_pdom.gff >> $tgf
  cat jigsaw/$fn/$fn.twinscan.gff >> $tgf
  echo "running: gff2gview -d $gvdir/$d -o $en -T $gvdir/tracks.def jigsaw/$fn/$fn.fa $tgf"
  gff2gview -d $gvdir/$d -o $en -T $gvdir/tracks.def jigsaw/$fn/$fn.fa $tgf
  if ( -f $bkm ) then
   cp $bkm $gvdir/$d
   endif
  end
  
 cd ..
 
end
