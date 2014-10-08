#!/bin/tcsh
foreach f ( *.isoFa.[1-9] *.isoFa.[1-9][0-9] *.isoFa.[1-9][0-9][0-9] )
 #set fnum = $f:e
 set ts = $f'.twsgtf'
 set tsgff = $ts:r
 set tsgff = $tsgff'.tws.gff'
 set tso = $f:s/.isoFa/.twinscan.gff/
 set tso = $tso:r
 if (-s $ts) then
  echo "processing $f [=> $tsgff, => $tso ]"
 
  twinscan2gff.pl $ts > $tsgff
  set ofs = `head -1 $f | cut -f2 -d ' '`
  echo "             ofs = $ofs"
  mapgffback.pl $tsgff $ofs >> $tso
  endif
end
