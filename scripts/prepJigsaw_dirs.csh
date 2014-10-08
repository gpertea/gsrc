#!/bin/tcsh -f
if ($1"x" == "x") then
 echo "Usage: prepJigsaw_dirs.csh fasta_files.."
 echo ""
 echo "Must be run from the jigsaw subdirectory"
 echo "and expects ../glimmerhmm/, ../mrnasearch and ../protsearch/ to"
 echo " contain the appropriate evidence files"
 echo ""
 echo "Example: prepJigsaw ../ga_??.fa"
 exit
 endif

# prepare the name translation table


set f=$1
set fp=$f:h
set ftx="$fp/seqn2fname.lst"
set dlst=""
set flst=""
/bin/rm -f $ftx
foreach f ( $* )
 set fn=$f:t
 set bd=$fn:r
 set sn=`head -1 $f | cut -f1 -d ' ' | cut -b2-`
 if (! -d $bd) mkdir $bd
 echo "$sn\t$bd" >> $ftx
 /bin/rm -f $bd/$bd.pmap.gff3
 /bin/rm -f $bd/$bd.pblat.gff3
 /bin/rm -f $bd/$bd.pexo.gff3
 /bin/rm -f $bd/$bd.gmap.gff3
 /bin/rm -f $bd/$bd.sim4.gff3
 set flst = "$flst $fn"
end

#set 
set pmap_res=`/bin/ls ../protsearch/*pmap*.gff3`
if ( $#pmap_res > 0 ) then
  set pmap_res=$pmap_res[1]
  echo "Runnining gff_split for pmap results ($pmap_res).."
  gff_split -m $ftx -s pmap.gff3 $pmap_res
endif

#set 
set pexo_res=`/bin/ls ../protsearch/exo*.gff3`
if ( $#pexo_res > 0 ) then
  set pexo_res=$pexo_res[1]
  echo "Runnining gff_split on exonerate p2g results ($pexo_res).."
  gff_split -m $ftx -s pexo.gff3 $pexo_res
endif


set pblt_res=`/bin/ls ../protsearch/*blat*.gff3`
if ( $#pblt_res > 0 ) then
  set pblt_res=$pblt_res[1]
  echo "Runnining gff_split on blat p2g results ($pblt_res).."
  gff_split -m $ftx -s pblat.gff3 $pblt_res
endif


set gmap_res=`/bin/ls ../sim4search/*gmap*.gff3`
if ( $#gmap_res > 0 ) then
  set gmap_res=$gmap_res[1]
  echo "Runnining gff_split on gmap results ($gmap_res).."
  gff_split -m $ftx -s gmap.gff3 $gmap_res
endif

set sim4_res=`/bin/ls ../sim4search/sim4*.gff3`
if ( $#sim4_res > 0 ) then
  set sim4_res=$sim4_res[1]
  echo "Runnining gff_split on sim4cc results ($sim4_res).."
  gff_split -m $ftx -s sim4.gff3 $sim4_res
endif


foreach fn ( ${flst} )
 set d=$fn:r
 if ( ! -d $d ) continue
 echo "processing dir: $d"
 cp $fp/$fn $d/$d.fa
 cp ../glimmerhmm/$d.glimmerhmm*.gff $d/
 cp ../glimmerhmm/$d.twinscan.gff $d/
 
 cd $d
 #gmap
 if ( $gmap_res"x" != "x" ) then
  if ( -s $d.gmap.gff3 ) then
   gffilter -g $d.fa -a $d.gmap_CDS.gtf $d.gmap.gff3
   perl -ne '@t=split("\t");print $_ if $t[2] eq "mRNA" || $t[2] eq "exon"' \
       < $d.gmap.gff3 > $d.gmap_exon.gff3
   else
    touch $d.gmap_CDS.gtf $d.gmap_exon.gff3
    echo "Warning: gmap file $d.seqfseq.gmap.gff3 not found."
   endif
 endif
 #pmap
 if ( $pmap_res"x" != "x" ) then
   if ( -s $d.pmap.gff3 ) then
    gffilter -g $d.fa -a $d.pmap_CDS.gtf $d.pmap.gff3
    perl -ne '@t=split("\t");print $_ if $t[2] eq "mRNA" || $t[2] eq "exon"' \
      < $d.pmap.gff3 > $d.pmap_exon.gff3
    else
     touch $d.pmap_CDS.gtf $d.pmap_exon.gff3
     echo "Warning: pmap file $d.pmap.gff3 not found."
    endif
 endif 
 #exonerate p2g
 if ( $pexo_res"x" != "x" ) then
   if ( -s $d.pexo.gff3 ) then
    gffilter -g $d.fa -a $d.pexo_CDS.gtf $d.pexo.gff3
    perl -ne '@t=split("\t");print $_ if $t[2] eq "mRNA" || $t[2] eq "exon"' \
      < $d.pexo.gff3 > $d.pexo_exon.gff3
    else
     touch $d.pexo_CDS.gtf $d.pexo_exon.gff3
     echo "Warning: exonerate file $d.pexo.gff3 not found."
    endif
 endif 
 #blat p2g
 if ( $pblt_res"x" != "x" ) then
   if ( -s $d.pblat.gff3 ) then
    gffilter -g $d.fa -a $d.pblat_CDS.gtf $d.pblat.gff3
    else
     touch $d.pblat_CDS.gtf
     echo "Warning: gblat uniprot file $d.pblat.gff3 not found."
    endif
 endif 
 if ( $sim4_res"x" != "x" ) then
   if ( -s $d.sim4.gff3 ) then
    gffilter -g $d.fa -a $d.sim4_CDS.gtf $d.sim4.gff3
    else
     touch $d.sim4_CDS.gtf
     echo "Warning: sim4cc file $d.sim4.gff3 not found."
    endif
 endif 
 cd ..
end

echo "-= Done =-"
exit
