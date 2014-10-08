#!/bin/bash -eux
#lower cvg => maxCvgReported=63.0
N="$1"
if [[ -z "$N" ]]; then
  echo "Error: pool# or 'fix' required as parameter 1" >&2
  exit 1
fi
fixing=""
nsplit=""
BESonly=""
if [[ "$N" == "fix" ]]; then
  #just to fix a too-long scaffold, run this in the pool directory
  N=${PWD##*/}
  fixing="yes"
else
  if [[ "$N" == "split" || "$N" == "nsplit" ]]; then
    N=${PWD##*/}
    fixing="yes"
    nsplit="yes"
  fi
fi
if [[ "$N" == "BES" ]]; then
 N=${PWD##*/}
 BESonly="yes"
fi

if [[ -z "$fixing" ]]; then
  nsplit="yes" #automatically split large gaps
fi


P=12
K=127
k=63
L=155
export D=5
d=1
I=98.5

###########################
if [[ -z "$BESonly" ]]; then

if [ ! -s asm-3.scafSeq.info ] ; then
SOAPdenovo2-127mer pregraph -s SOAPdenovo.config -K $K -p $P -R -o asm-1 -d $d
SOAPdenovo2-127mer contig -R -g asm-1 -D $D
/data1/igm3/sw/bin/infoseq -noheading -only -name -len -pgc -desc -filter asm-1.contig | sort -nk2  | sed 's|_| |g' | \
  awk '{print $1,$2,$3,$7,$9}'  > asm-1.contig.info
SOAPdenovo2-127mer map -s SOAPdenovo.config -g asm-1 -p $P -k $k -f
SOAPdenovo2-127mer scaff -g asm-1 -p $P -F  -w -u
echo ">" >> asm-1.scaf
~/bin/wheat/SOAPdenovo2repeats.sh ${k}.0

GapCloser -a asm-1.scafSeq -b SOAPdenovo.config -o asm-2.scafSeq -l $L  -t $P
/data1/igm3/sw/bin/infoseq -noheading -only -name -len -pgc -desc -filter asm-2.scafSeq | sort -nk2 -r > asm-2.scafSeq.info

#filterFastaLen.pl -Max 1999 < asm-2.scafSeq | cleanFasta.pl > asm-2.scafSeq2K
seqmanip -l100 -x1999 -UP < asm-2.scafSeq > asm-2.scafSeq2K

nucmer -c $k -l $k -maxmatch -nosimplify asm-2.scafSeq2K asm-2.scafSeq -p asm-2.scafSeq 
show-coords -l -r -d -o -I $I asm-2.scafSeq.delta  | grep CONTAINED | count.pl -i -3  > asm-2.scafSeq.delete
show-coords -l -r -d -o -I $I asm-2.scafSeq.delta  | grep IDENTITY | perl -ane 'print $_ if($F[-3] lt $F[-2]);' | count.pl -i -3 >> asm-2.scafSeq.delete
show-coords -l -r -d    -I $I asm-2.scafSeq.delta  | perl -ane 'print if( $F[11]<=255 and $F[11]<$F[12] and $F[6]>=127 );' | count.pl -i -2 >> asm-2.scafSeq.delete
cat asm-2.scafSeq.info  | perl -ane 'print $_ if($F[1]<2000 and $F[-1]<$ENV{D});' >> asm-2.scafSeq.delete
intersect.pl asm-2.scafSeq.info  asm-2.scafSeq.delete  > asm-2.scafSeq.delete.tmp
mv asm-2.scafSeq.delete.tmp asm-2.scafSeq.delete    

extractFromFastaNames.pl -f asm-2.scafSeq.delete < asm-2.scafSeq -n > asm-3.scafSeq
/data1/igm3/sw/bin/infoseq -noheading -only -name -len -pgc -desc -filter asm-3.scafSeq | sort -nk2 -r > asm-3.scafSeq.info

fi

########################
#if [ ! -s summary.txt ]; then
####
if [ -n "$fixing" ] ; then
 test -f asm-3.scafSeq.ctgs.fa
 test -f asm-4.updated.edge
 if [[ -z "$nsplit" ]]; then
   if [ ! -f asm-4.updated.edge.orig ]; then
     mv asm-4.updated.edge asm-4.updated.edge.orig
   fi
   perl -pe '($n)=(/length (\d+)/);if ($n) {$n-=127;s/length \d+/length $n/}' \
    asm-4.updated.edge.orig > asm-4.updated.edge
 else
   if [ -f asm-4.updated.edge.orig ]; then
     cp asm-4.updated.edge.orig asm-4.updated.edge
   fi
 fi
 /bin/rm -f asm-4.scaf*
 set +e
else
 /bin/rm -rf asm-4*
 #prepare -K $K -g asm-4 -c asm-3.scafSeq
 #splitbyNs.pl -n3 -o asm-3.scafSeq.ctgs.fa asm-3.scafSeq
 #finalFusion -D -K $K -g asm-4 -c asm-3.scafSeq.ctgs.fa
 gap_mod.pl -b asm-3.scafSeq.N2T.fa -o asm-3.scafSeq.gapinfo asm-3.scafSeq
 finalFusion -D -K $K -g asm-4 -c asm-3.scafSeq.N2T.fa
 set +e
 SOAPdenovo2-127mer map -s SOAPdenovo.config -g asm-4 -p $P -k $k -f
fi
SOAPdenovo2-127mer scaff -g asm-4 -p $P -F -u -w
set -e
if [ -s asm-4.scaf ]; then
  #restore gaps
  perl -p -i -e 's/^\n$//' asm-4.contig
  ctg2Soap.pl -o asm-4.scaff2ctg.tab asm-3.scafSeq.N2T.fa asm-4
  gap_mod.pl -I asm-3.scafSeq.gapinfo -c asm-4.scaff2ctg.tab asm-4.scafSeq > asm-4.scafSeq.fa
  maxgap=$(show_poly_runs.pl -p60 -G asm-4.scafSeq.fa | sort -k4nr | head -1 | cut -f4)
  if [[ -n "$fixing" && -z "$nsplit" && $maxgap -gt 600 ]]; then
    echo "Error: there is still a gap larger than 600 ($maxgap), try: $0 split'\n" >&2
    exit 1
  fi
  if [[ -n "$nsplit" && $maxgap -gt 600 ]]; then
    echo "$N error: gap too large detected ($maxgap), splitting by 600Ns" >&2
    mv asm-4.scafSeq.fa asm-4.scafSeq.to_splitbyN
    splitbyNs.pl -n600 -o asm-4.scafSeq.fa asm-4.scafSeq.to_splitbyN
  fi
  /data1/igm3/sw/bin/infoseq -noheading -only -name -len -pgc -desc -filter asm-4.scafSeq.fa | sort -nk2 > asm-4.scafSeq.info
  #touch asm-4.scafSeq.flt
  #cat asm-4.scafSeq.fa | fasta2tab.pl  | grep ^s | tab2fasta.pl  | cleanFasta.pl -l 100 >> asm-4.scafSeq.flt
  seqmanip -l100 -q100 -UP asm-4.scafSeq.fa > asm-4.scafSeq.flt
  #cat asm-4.scafSeq.info | grep ^C | sed 's|C||' | extractFromFastaNames.pl -f - asm-4.contig | sed 's|>|>C|' >> asm-4.scafSeq.flt
  
  GapCloser -a asm-4.scafSeq.flt -b SOAPdenovo.config -l $L  -t $P -o asm-5.scafSeq
  /data1/igm3/sw/bin/infoseq -noheading -only -name -len -pgc -desc -filter asm-5.scafSeq | sort -nk2 > asm-5.scafSeq.info
else
  ln -s asm-3.scafSeq asm-5.scafSeq
fi

#cat asm-5.scafSeq | awk '{print $1}' |  sortFastaLen.pl -rev -new -prefix $N | addLen2FastaId.pl | cleanFasta.pl > asm.fa
cat asm-5.scafSeq | awk '{print $1}' |  sortFastaLen.pl -rev -new -prefix $N | seqmanip -l100 -UPD > asm.fa
cat asm.fa  | fasta2tab.pl > asm.tab

/data1/igm3/sw/bin/infoseq -noheading -only -name -len -pgc -filter asm.fa > asm.info
cat asm.info | getSummary.pl -i 1 -min 2000  -t scafSeq > summary.txt
##########################

fi # not BES only

if [ -s /scratch0/igm1/dpuiu/Aegilops_tauschii/Data/BES/fastq/*/$N.seq ]; then
  if  [ ! -f BES.seq ] ; then
   ln -s /scratch0/igm1/dpuiu/Aegilops_tauschii/Data/BES/fastq/*/$N.seq BES.seq
  fi
  if [ -f BES.seq ] ; then
    /data1/igm3/sw/bin/infoseq -noheading -only -name -len -pgc -filter BES.seq > BES.info
    nucmer BES.seq asm.fa -p BES
    /data1/igm3/sw/bin/show-coords -c -l -H -q  -d  BES.delta |  max2.pl -i -2 -j 6 | sort -k20   >  BES.coords
    cat BES.coords  | sort2col.pl -i -1 -j 3 | perl -ane 'print "$F[-1] $F[-2]\n";' | uniq | \
     perl -ane ' if($.==1 or $F[0] ne $P[0]) { print "\n$F[0]"  } print " $F[1]" ; @P=@F; END { print "\n"}' >  BES.tab
    cat asm.fa | awk '{print $1}' | addLen2FastaId.pl | fasta2tab.pl  > asm.tab
    join.pl -all asm.tab BES.tab  | perl -ane '$n=scalar(@F); @F=@F[0,1,3..$n-1,2] if($n>3); print join " ",@F; print "\n";'| \
     tab2fasta.pl | cleanFasta.pl > asm.fa
  fi

  /data1/igm3/sw/bin/infoseq -noheading -only -name -len -pgc -filter -desc asm.fa  | perl -ane '$F[3]=""; print join "\t",@F; print "\n";' > asm.info
fi
exit 0
###########################

#to fix the beginning/end of scf truncation
#cat asm-4.scaf | ~/bin/scf2bed.pl | ~/bin/Pinus_taeda/bed2fasta.pl -bed - asm-4.contig  > asm-4.scafSeq.tmp
cat asm-4.scaf | scf2bed.pl | bed2fasta.pl -bed - asm-4.contig  > asm-4.scafSeq.fix
cat asm-4.scafSeq.fa | grep "^>C" | sed 's|>C||' | extractFromFastaNames.pl -f - asm-4.contig | sed 's|>|>C|' >> asm-4.scafSeq.fix
#mv asm-4.scafSeq.tmp asm-4.scafSeq

