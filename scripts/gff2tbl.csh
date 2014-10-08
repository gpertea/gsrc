#!/bin/tcsh -f

if ($1'x' == 'x') then
 echo "Usage: gff2tbl <jigsaw_output_suffix> <dirs>.. "
 echo "   Example:"
 echo " gff2tbl gblatp.jgff ba_?? "
 echo ""
 exit 1
endif
if ( ! -d submit ) then
 mkdir submit
endif
set jres=$1
shift
set fm="../fa2encode.map"
if ( ! -f $fm ) then
 foreach fb ( $* )
   set fa=../$fb.fa
   set eid=`head -1 $fa | perl -ne 'print $1 if (/ENCODE region (\w+)/)'`
   set acc=`head -1 $fa | perl -ne 'print $1 if (/^>(\w+)/)'`
   echo "$fb\t$eid" >> $fm
 end
endif

foreach d ( $* )
 if ( ! -d $d ) then
  echo "Warning: skipping $d, it's not a directory"
  continue
 endif
 set f=$d.ann.recon_$jres.gff3
 echo "processing file: $f"
 set fo=$d:s/_//
 echo "running: jgff2tbl -f $d/$d.fa -A -o submit/$fo $f"
 jgff2tbl -f $d/$d.fa -A -o submit/$fo $f
end
