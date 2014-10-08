#!/bin/tcsh -f
if ($2'x' == 'x') then
  echo "Usage: prep_strands.csh <dirlist_file> <evidence_file>"
  exit 1
  endif
set dirlist=$1
set evf=$2
set dlf='f_'$dirlist;
set dlr='r_'$dirlist;
rm -rf $dlf $dlr
foreach d ( `cat $dirlist` )
 set bname=$d:t
 set dbase=$d:h
 echo " processing $bname .."
 set df=$d'f'
 set dr=$d'r'
 echo $df >> $dlf
 echo $dr >> $dlr
 if (! -d $df) mkdir $df
 if (! -d $dr) mkdir $dr
 cp $d/$bname.fa $df/$bname'f.fa'
 cp $d/$bname.fa $dr/$bname'r.fa'
 foreach f ( `cut -f1 -d ' ' $evf` )
   set src=$d/$bname.$f
   #ls -l $f
   set nf=$bname'f.'$f
   set nr=$bname'r.'$f
   #echo "$nf   $nr"
   gawk ' {if ($7=="+") print $0 }' < $src | sed "s/^$bname/${bname}f/" > $df/$nf
   gawk ' {if ($7=="-") print $0 }' < $src | sed "s/^$bname/${bname}r/" > $dr/$nr
 end
end
