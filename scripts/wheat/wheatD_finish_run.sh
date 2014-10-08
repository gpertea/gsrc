#!/bin/bash -eu
#must be run in the main run directory, with ./asm/<pool#>/ subdirectories
test -f pools_count_mea_std.info

#R=${PWD##*/}
export R=`echo ${PWD##*/} | sed 's|run||' | sed 's|_| |' | awk '{print $1}'`
N=`cat pools_count_mea_std.info | grep -v ^# | wc -l`
RN="$R:$N"

export D=`date "+%m%d%Y"`

###########################################

grep -v '^#' pools_count_mea_std.info   | getSummary.pl -i 1 -t mates | awk '{print $9}' | sed 's|sum|#mates|' > run_summary.mates.txt
grep -v '^#' pools_count_mea_std.info   | getSummary.pl -i 2 -t meaIns |  awk '{print $8}' | sed 's|mean|meaIns|' > run_summary.meaIns.txt


cat asm/*/asm.fa   > asm.fa
cat asm/*/asm.info > asm.info
scf2gap.pl < asm.fa | /data1/igm3/sw/bin/infoseq  -noheading -only -name -len -filter > asm.gap.len
cat asm.info | getSummary.pl -i 1 -min 2000 -t $RN  | sed 's|^\.|run:#pools|' | sed 's|:| |' | pretty.pl  > run_summary.txt

grep -m 1 ^scafSeq asm/????/summary.txt  | perl -pe 's|asm/(\w+)/summary.txt:scafSeq|$1|' > pools_summary.txt
join.pl pools_summary.txt pools_count_mea_std.info  > pools_summary.txt.tmp; mv pools_summary.txt.tmp pools_summary.txt
egrep -m 1 '^[1-9]' asm/????/asm-1.scaf.rep | perl -ane '$F[0]=~s/^asm\/(\w+)\/.+/$1/; print $F[0]." ".$F[1]."\n"' | join.pl pools_summary.txt - -all -empty 0 > pools_summary.txt.tmp ; mv pools_summary.txt.tmp pools_summary.txt
grep -c "^>" asm/*/BES.seq | perl -pe 'tr/:/ /;s/^asm\/(\w+)\/\S+/$1/' | join.pl pools_summary.txt - -all -empty 0 > pools_summary.txt.tmp ; mv pools_summary.txt.tmp pools_summary.txt
wc -l asm/*/BES.coords | fgrep coords | perl -ane '$F[1]=~s/^asm\/(\w+)\/.+/$1/; print $F[1]." ".$F[0]."\n"' | join.pl pools_summary.txt - -all -empty 0 > pools_summary.txt.tmp ; mv pools_summary.txt.tmp pools_summary.txt
cat pools_summary.txt | perl -ane 'chomp ; print $_," 0 0\n";' > pools_summary.txt.tmp ; mv pools_summary.txt.tmp pools_summary.txt
cat asm.gap.len | sed 's|\.| |' | count.pl -i 0 |  join.pl pools_summary.txt - -all -empty 0 > pools_summary.txt.tmp ; mv pools_summary.txt.tmp pools_summary.txt
cat asm.gap.len | sed 's|\.| |' | awk '{print $1,$3}' |  max2.pl |  join.pl pools_summary.txt - -all -empty 0 > pools_summary.txt.tmp ; mv pools_summary.txt.tmp pools_summary.txt
cat asm.gap.len | sed 's|\.| |' | awk '{print $1,$3}' |  sum2.pl |  join.pl pools_summary.txt - -all -empty 0 > pools_summary.txt.tmp ; mv pools_summary.txt.tmp pools_summary.txt
ls -l --time-style="+%m%d%Y" asm/????/asm-3.scafSeq | perl -ane 'print "$ENV{R} $F[5] v1 $ENV{D}\n";'  | paste pools_summary.txt - | pretty.pl > pools_summary.txt.tmp ; mv pools_summary.txt.tmp pools_summary.txt
cat asm.info     | getSummary.pl -i 1 -min 2000 -t $RN  | sed 's|^\.|run:#pools|' | sed 's|:| |' | paste - run_summary.mates.txt run_summary.meaIns.txt  | pretty.pl  > run_summary.txt

###

cat asm/*/asm.wgs.fa  > asm.wgs.fa
cat asm/*/asm.wgs.info > asm.wgs.info
scf2gap.pl < asm.wgs.fa | /data1/igm3/sw/bin/infoseq  -noheading -only -name -len	   -filter > asm.wgs.gap.len
cat asm.wgs.info | getSummary.pl -i 1 -min 2000 -t $RN  | sed 's|^\.|run:#pools|' | sed 's|:| |' | pretty.pl  > run_summary.wgs.txt

grep -m 1 ^scafSeq.wgs asm/????/summary.txt  | perl -pe 's|asm/(\w+)/summary.txt:scafSeq\.wgs|$1|' > pools_summary.wgs.txt
join.pl pools_summary.wgs.txt pools_count_mea_std.info  > pools_summary.wgs.txt.tmp; mv pools_summary.wgs.txt.tmp pools_summary.wgs.txt
egrep -m 1 ^[1-9] asm/????/asm-1.scaf.rep | perl -ane '$F[0]=~s/^asm\/(\w+)\/.+/$1/; print $F[0]." ".$F[1]."\n"' | join.pl pools_summary.wgs.txt - -all -empty 0 > pools_summary.wgs.txt.tmp ; mv pools_summary.wgs.txt.tmp pools_summary.wgs.txt
grep -c "^>" asm/*/BES.seq | perl -pe 'tr/:/ /;s/^asm\/(\w+)\/\S+/$1/' | join.pl pools_summary.wgs.txt - -all -empty 0 > pools_summary.wgs.txt.tmp ; mv pools_summary.wgs.txt.tmp pools_summary.wgs.txt
wc -l asm/*/BES.wgs.coords | fgrep coords | perl -ane '$F[1]=~s/^asm\/(\w+)\/.+/$1/; print $F[1]." ".$F[0]."\n"' | join.pl pools_summary.wgs.txt - -all -empty 0 > pools_summary.wgs.txt.tmp ; mv pools_summary.wgs.txt.tmp pools_summary.wgs.txt
cat pools_summary.wgs.txt | perl -ane 'chomp ; print $_," 6\n";' > pools_summary.wgs.txt.tmp ; mv pools_summary.wgs.txt.tmp pools_summary.wgs.txt
grep -c ^@SRR asm/*/SRR*_1.fastq  | perl -pe 's|^asm/(\w+)/|$1 |' | sed 's|:| |' | sum2.pl -i 0 -j -1 | join.pl pools_summary.wgs.txt - -all -empty 0 > pools_summary.wgs.txt.tmp ; mv pools_summary.wgs.txt.tmp pools_summary.wgs.txt
cat asm.wgs.gap.len | sed 's|\.| |' | count.pl -i 0 |  join.pl pools_summary.wgs.txt - -all -empty 0 > pools_summary.wgs.txt.tmp ; mv pools_summary.wgs.txt.tmp pools_summary.wgs.txt
cat asm.wgs.gap.len | sed 's|\.| |' | awk '{print $1,$3}' |  max2.pl |  join.pl pools_summary.wgs.txt - -all -empty 0 > pools_summary.wgs.txt.tmp ; mv pools_summary.wgs.txt.tmp pools_summary.wgs.txt
cat asm.wgs.gap.len | sed 's|\.| |' | awk '{print $1,$3}' |  sum2.pl |  join.pl pools_summary.wgs.txt - -all -empty 0 > pools_summary.wgs.txt.tmp ; mv pools_summary.wgs.txt.tmp pools_summary.wgs.txt
ls -l --time-style="+%m%d%Y" asm/????/asm-6.scafSeq | perl -ane 'print "$ENV{R} $F[5] v1 $ENV{D}\n";'  | paste pools_summary.wgs.txt - | pretty.pl > pools_summary.wgs.txt.tmp ; mv pools_summary.wgs.txt.tmp pools_summary.wgs.txt
cat asm.wgs.info | getSummary.pl -i 1 -min 2000 -t $RN  | sed 's|^\.|run:#pools|' | sed 's|:| |' | paste - run_summary.mates.txt run_summary.meaIns.txt  | pretty.pl  > run_summary.wgs.txt

#rm run_summary.mates.txt  run_summary.meaIns.txt 
