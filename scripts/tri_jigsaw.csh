#!/bin/tcsh -f

#set gmap = 'refseq_exon.gff3'
set gmap = 'sim4.gff3'
#set pmap = 'unipr.blat.gff3'
set pmap = 'pexo_exon.gff3'

set mtrack='sim4cc'
set ptrack='p2g_exo'

if ($2'x' == 'x') then
  echo "Usage: tri_jigsaw.csh <dirlist_file> <jigsaw_output_suffix>"
  echo "  <dirlist_file> is the one corresponding to the 'both strands' run"
  echo "  <jigsaw_output_suffix> is the one that was used when jigsaw was run"
  echo "         (e.g. refseq.jgff) "
  echo '  WARNING: make sure variables $pmap and $gmap are correct!'
  echo "          (currently they are: $pmap and $gmap)"
  exit 1
  endif

set dirlist=$1
set jsuf=$2
echo "..running:\ntri_jigsaw.csh $*"
set protdb = '~/ann/protdb/unipr_mammals.fa.cidx'
set predcount=0
foreach d ( `cat $dirlist` )
 set bname=$d:t
 set dbase=$d:h
 #echo "dirbase=$dbase (bname=$bname)"
 set fname=$bname'f'
 set rname=$bname'r'
 cd $dbase
 
 set j_gff3 = $bname.recon_$jsuf.gff3
 set j_anngff3 = $bname.ann.recon_$jsuf.gff3
 echo "recon_jigsaw.pl -P $bname/$bname.fa -o $j_gff3 $bname/$bname.$jsuf $fname/$fname.$jsuf $rname/$rname.$jsuf"
 recon_jigsaw.pl -P $bname/$bname.fa -o $j_gff3 $bname/$bname.$jsuf $fname/$fname.$jsuf $rname/$rname.$jsuf
 #ls -al $bname/$bname.pmap.fltOK.gtf
 set pcount=`cat $bname/$bname.$jsuf | grep -v '^#' |  cut -f9 | cut -f1 -d ';' | sort -u | wc -l`
 @ predcount = $predcount + $pcount
 set d=$bname
 if (-s $j_gff3) then
   #--validate the resulting gff
   # set badcds=$bname.recon_$jsuf.badCDS.gtf
   # gffilter -g $bname/$bname.fa -b $badcds $j_gff3
   # if (-s $badcds) then 
   #   echo "**** ERROR found at CDS validation for $j_gff3! (see $badcds)"
   #   else
   #   /bin/rm -f $badcds
   #   endif
   /bin/rm -f $d.[pg]map.i{it,fa}
   /bin/rm -f $d.pmap.iit $d.pmap.ifa
    #gff2iit -o $d.pmap -t $ptrack $d/$d.$pmap
    gtf2gff -t $ptrack $d/$d.$pmap  > $d.iit_pmap.gff3
    gtf2gff -t 'pmap'  $d/$d.pmap_exon.gff3 >> $d.iit_pmap.gff3
    gff2iit -o $d.pmap $d.iit_pmap.gff3
   /bin/rm -f $d.gmap.iit $d.gmap.ifa
    #gff2iit -o $d.gmap -t $mtrack $d/$d.$gmap
   gtf2gff -t $mtrack $d/$d.$gmap  | perl -pe 's/\.mrna(\d+)/\.m$1/' > $d.iit_gmap.gff3
   gtf2gff -t 'gmap'  $d/$d.gmap_exon.gff3  >> $d.iit_gmap.gff3
   gff2iit -o $d.gmap -t $mtrack $d.iit_gmap.gff3
    
   #echo "gffann.pl -o $j_anngff3 -t $d -m $d.gmap.iit -p $d.pmap.iit -P $protdb $j_gff3"
   gffann.pl -o $j_anngff3 -t $d -m $d.gmap.iit -p $d.pmap.iit -P $protdb $j_gff3
   #  -- temporary: throw in all the mapping evidence
   gtf2gff -t $mtrack $d/$d.$gmap  > $d.top5gmap.gff3
   gtf2gff -t $ptrack $d/$d.$pmap  > $d.top5pmap.gff3
   gtf2gff -t 'gmap'  $d/$d.gmap_exon.gff3  >> $d.top5gmap.gff3

   gtf2gff -t 'pmap'  $d/$d.pmap_exon.gff3 >> $d.top5pmap.gff3
   gtf2gff -t 'blat'  $d/$d.pblat.gff3 >> $d.top5pmap.gff3
   
  else #no predictions.. just show whatever evidence is there
   cp $j_gff3 $j_anngff3
   gtf2gff -t $mtrack $d/$d.$gmap  | perl -pe 's/\.mrna(\d+)/\.m$1/' > $d.top5gmap.gff3
   gtf2gff -t $ptrack $d/$d.$pmap  > $d.top5pmap.gff3
   gtf2gff -t 'gmap'  $d/$d.gmap_exon.gff3  >> $d.top5gmap.gff3
   gtf2gff -t 'pmap'  $d/$d.pmap_exon.gff3 >> $d.top5pmap.gff3
   gtf2gff -t 'blat'  $d/$d.pblat.gff3 >> $d.top5pmap.gff3
  endif
 
 cd ..
end

echo "Total: $predcount gene models predicted."
