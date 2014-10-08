#!/bin/tcsh -f
if ("x"$1 == "x") then
 echo "Usage: cuffcompare_stats.txt <transcripts.tmap>"
 exit 1
endif
set t=$1
#echo "Significantly matching transfrags distribution:"
set tdist=`grep -P '\t[=cj]\t' $t.tmap | cut -f3 | sort | uniq -c`
#echo "$tdist[1]\t$tdist[2]"
#echo "$tdist[3]\t$tdist[4]"
#echo "$tdist[5]\t$tdist[6]"

#echo "Perfectly matched gene transcripts:"
echo "transfrags\tcode\tgenes"
set g1=`cut -f1,3 $t.tmap | grep -P '\t=$' | sort -u | wc -l`
set g2=`cut -f1,3 $t.tmap | grep -P '\tc$' | sort -u | wc -l`
set g3=`cut -f1,3 $t.tmap | grep -P '\tj$' | sort -u | wc -l`
echo "$tdist[1]    \t$tdist[2]\t$g1"
echo "$tdist[3]    \t$tdist[4]\t$g2"
echo "$tdist[5]    \t$tdist[6]\t$g3"

