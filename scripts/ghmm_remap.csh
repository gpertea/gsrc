#!/bin/tcsh
# if ($1'x' == 'x') then
#  echo "Usage: ghmm_remap.csh <ghmm..."
#  echo 
#  exit 1
# endif


foreach f ( *.isoFa.[1-9] *.isoFa.[1-9][0-9] )
 set ghpd = $f:s/.isoFa./.isoFa.ghmm.pdom./
 set gh = $f:s/.isoFa./.isoFa.ghmm./
 set ghopd = $f:r
 set ghopd = $ghopd:s/.isoFa/.glimmerhmm_pdom.gff/
 set gho = $f:r
 set gho = $gho:s/.isoFa/.glimmerhmm.gff/
 set ofs = `head -1 $f | cut -f2 -d ' '`
 mapgffback.pl $ghpd $ofs >> $ghopd
 mapgffback.pl $gh $ofs >> $gho
 #exit
end
