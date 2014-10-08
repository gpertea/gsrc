#!/bin/bash -eux

#63.0 or 127.0
C=$1

/data1/igm3/sw/bin/infoseq -noheading -only -name -len -pgc -desc -filter asm-1.contig | sort -nk2  | sed 's|_| |g' | awk '{print $1,$2,$3,$7}' > asm-1.contig.info 

cat asm-1.scaf | grep "^>" -A 1 | \
 perl -ane 'if(/>(\S+)/) { $id=$1 } elsif(/UP (.+)/) { @F=split /\s/,$1 ; foreach (@F) { /(\d+):(\d+)/ ;print "$1 $id\n" if($2>1) }} ;'  >  asm-1.scaf.up
cat asm-1.scaf | grep "^>" -B 1 | \
 perl -ane 'if(/>(\S+)/) { $id=$1 } elsif(/DOWN (.+)#UP/ or /DOWN (.+)/) { @F=split /\s/,$1 ; foreach (@F) { /(\d+):(\d+)/ ;print "$1 $id\n" if($2>1)}} ;'  >   asm-1.scaf.down

cat asm-1.scaf.up asm-1.scaf.down  | count.pl -m 2        | awk '{print $1,$3,$2}' > asm-1.scaf.rep
cat asm-1.scaf.{up,down}           | count.pl -i 0 -min 3 | difference.pl - asm-1.scaf.rep | \
 awk '{print $1,$2}' >> asm-1.scaf.rep
cat asm-1.contig.info | awk '{print $1,$2,$3,$4}' | grep "$C$" | difference.pl - asm-1.scaf.rep | \
 difference.pl - asm-1.scaf  | awk '{print $1, 0}' >> asm-1.scaf.rep

join.pl asm-1.contig.info asm-1.scaf.rep | sort -nk1 -r  > asm-1.scaf.rep.tmp ; mv asm-1.scaf.rep.tmp asm-1.scaf.rep

#join.pl  asm-1.contig.info asm-1.scaf.rep  | sort -nk1 -r | sed 's|^|C|' > asm-1.scaf.rep.tmp ; mv asm-1.scaf.rep.tmp asm-1.scaf.rep
#cat asm-1.scaf.rep  | perl -ane 'print $_ if(@F==6);' | tee asm-1.scaf.rep2
