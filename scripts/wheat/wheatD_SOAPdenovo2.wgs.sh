#!/bin/bash -eux
#lower cvg => maxCvgReported=63.0
#echo v.1.03
fixing=""
N="$1"
if [[ -z "$N" ]]; then
  echo "Error: pool#, 'fix' or 'redo' required as parameter 1" >&2
  exit 1
fi
bfix=""
nsplit=""
BESonly=""
if [[ "$N" == "fix" ]]; then
  #just to fix a too-long scaffold, run this in the pool directory
  N=${PWD##*/}
  fixing="yes"
  bfix="-b 1.3"
else
 if [[ "$N" == "efix" || "$N" == "fix_edge" || "$N" == "fixedge" ]]; then
  N=${PWD##*/}
  fixing="yes"
  bfix=""
 else
  if [[ "$N" == "split" || "$N" == "nsplit" ]]; then
    N=${PWD##*/}
    fixing="yes"
    nsplit="yes"
  fi
 fi
fi

if [[ "$N" == "redo" || "$N" == "again" || "$N" == "rerun" ]]; then
  N=${PWD##*/}
  /bin/rm -f asm.*mer.tab
fi

if [[ "$N" == "BES" ]]; then
 N=${PWD##*/}
 BESonly="yes"
fi


if [[ -z "$fixing" ]]; then
  nsplit="yes"
fi
P=16
L=155
export K=75
export M=8600

test -f asm.fa

###########################
if [[ -z "$BESonly" ]]; then


if [ -n "$fixing" ] ; then
 test -f asm.ctgs.fa
 test -f asm-6.updated.edge
 if [ -z $bfix ] ;  then 
   if [ ! -f asm-6.updated.edge.orig ]; then
     mv asm-6.updated.edge asm-6.updated.edge.orig
   fi
   edgedec=$(head -2 asm-6.updated.edge.orig | tail -1 | perl -lne '($v)=(m/length\s+(\d+)/);print $v<50?$v:$v-40')
   perl -pe '($n)=(/length (\d+)/);if ($n) {$n-='$edgedec';s/length \d+/length $n/}' \
   asm-6.updated.edge.orig > asm-6.updated.edge
 else # nsplit or bfix
   if [ -f asm-6.updated.edge.orig ]; then
     cp asm-6.updated.edge.orig asm-6.updated.edge
   fi
 fi
 /bin/rm -f asm-6.scaf*
 set +e
else
 if [ ! -s asm.75mer.tab ] ; then
 cat asm.fa | fasta2tab.pl |  perl -ane \
  '$n=length($F[-1]); foreach(0..$n-$ENV{K}) { print substr($F[-1],$_,$ENV{K}),"\t",$F[0],"\t",$_,"\t",$_+$ENV{K},"\t5\n" }' \
  | grep -v N  > asm.${K}mer.tab
 cat asm.fa | /data1/igm3/sw/bin/revseq -filter | fasta2tab.pl | \
  perl -ane '$n=length($F[-1]); foreach(0..$n-$ENV{K}) { print substr($F[-1],$_,$ENV{K}),"\t",$F[0],"\t",$_,"\t",$_+$ENV{K},"\t3\n" }' \
  | grep -v N  >> asm.${K}mer.tab
 cat asm.${K}mer.tab |  count.pl -i 0 -M 1 |  intersect.pl asm.${K}mer.tab -  \
  | perl -ane 'print $_ if($F[2]<$ENV{M});' > asm.${K}mer.tab.uniq ; mv asm.${K}mer.tab.uniq asm.${K}mer.tab
 cat /home/gpertea/bin/wheat/filter${K}mer.all.sh  | scheduler.pl
 fi

 /bin/rm -rf asm-6* asm.wgs* BES.wgs.coords
 ##prepare -K $K -g asm-6 -c asm.fa
 #splitbyNs.pl -n3 -o asm.ctgs.fa asm.fa
 gap_mod.pl -b asm.N2T.fa -o asm.gapinfo asm.fa
 #finalFusion -D -K $K -g asm-6 -c asm.ctgs.fa
 finalFusion -D -K $K -g asm-6 -c asm.N2T.fa
 set +e
 SOAPdenovo2-127mer map -s SOAPdenovo.wgs.config -g asm-6 -p $P -k $K -f
fi

# SOAPdenovo2-127mer scaff -g asm-6 -p $P -F -u -w -b 1.3 # try lower -b value
SOAPdenovo2-127mer scaff -g asm-6 -p $P -F -u -w $bfix
set -e
################33
fi # no BESonly
if [ -s asm-6.scaf ]
then
  ###########################
  if [[ -z "$BESonly" ]]; then

    #restore gaps
    perl -p -i -e 's/^\n$//' asm-6.contig
    ctg2Soap.pl -o asm-6.scaff2ctg.tab asm.N2T.fa asm-6
    gap_mod.pl -I asm.gapinfo -c asm-6.scaff2ctg.tab asm-6.scafSeq > asm-6.scafSeq.fa

    maxgap=$(show_poly_runs.pl -p60 -G asm-6.scafSeq.fa | sort -k4nr | head -1 | cut -f4)
    if [[ -n "$fixing" && -z "$nsplit" && $maxgap -gt 9000 ]]; then
      echo "Error: there is still a gap larger than 9000 ($maxgap), try: $0 split'\n" >&2
      exit 1
    fi
    if [[ -n "$nsplit" && $maxgap -gt 9000 ]]; then
      echo "$N error: gap too large detected ($maxgap), splitting by 9000Ns" >&2
      mv asm-6.scafSeq.fa asm-6.scafSeq.to_splitbyN
      splitbyNs.pl -n9000 -o asm-6.scafSeq.fa asm-6.scafSeq.to_splitbyN
    fi
    /data1/igm3/sw/bin/infoseq -noheading -only -name -len -pgc -desc -filter asm-6.scafSeq.fa | sort -nk2 > asm-6.scafSeq.info
    #touch asm-6.scafSeq.flt
    #cat asm-6.scafSeq.fa | fasta2tab.pl  | grep ^s | tab2fasta.pl  | cleanFasta.pl -l 100 >> asm-6.scafSeq.flt
    #cat asm-6.scafSeq.info | grep ^C | sed 's|C||' | extractFromFastaNames.pl -f - asm-6.contig | sed 's|>|>C|' >> asm-6.scafSeq.flt
    seqmanip -l100 -q100 -UP asm-6.scafSeq.fa > asm-6.scafSeq.flt

    GapCloser -a asm-6.scafSeq.flt -b SOAPdenovo.config -l $L  -t $P -o asm-7.scafSeq

    cat asm-7.scafSeq | awk '{print $1}' | sortFastaLen.pl -rev -new -prefix $N > asm.wgs.fa
    /data1/igm3/sw/bin/infoseq -noheading -only -name -len -pgc -desc -filter asm.wgs.fa > asm.wgs.info
  fi ### no BESonly
  if [ -s /scratch0/igm1/dpuiu/Aegilops_tauschii/Data/BES/fastq/*/$N.seq ]; then
   if [ -f BES.seq ]; then
     nucmer BES.seq asm.wgs.fa -p BES.wgs
     /data1/igm3/sw/bin/show-coords -c -l -H -q  -d  BES.wgs.delta |  max2.pl -i -2 -j 6 | sort -k20   >  BES.wgs.coords
     cat BES.wgs.coords  | sort2col.pl -i -1 -j 3 | perl -ane 'print "$F[-1] $F[-2]\n";' | uniq  \
      | perl -ane ' if($.==1 or $F[0] ne $P[0]) { print "\n$F[0]"  } print " $F[1]" ; @P=@F; END { print "\n"}' >  BES.wgs.tab
     cat asm.wgs.fa | awk '{print $1}' | addLen2FastaId.pl | fasta2tab.pl  > asm.wgs.tab
     join.pl -all asm.wgs.tab BES.wgs.tab  | perl -ane '$n=scalar(@F); @F=@F[0,1,3..$n-1,2] if($n>3); print join " ",@F; print "\n";'| tab2fasta.pl | cleanFasta.pl > asm.wgs.fa
     /data1/igm3/sw/bin/infoseq -noheading -only -name -len -pgc -filter -desc asm.wgs.fa  | \
        perl -ane '$F[3]=""; print join "\t",@F; print "\n";' > asm.wgs.info
   fi
  fi
  if [[ "$BESonly" ]]; then
    exit 0
  fi
else
  ln -s asm.fa asm.wgs.fa
  ln -s	asm.info asm.wgs.info
  if [ -f BES.coords ]; then
    ln -s BES.coords BES.wgs.coords
  fi
fi

###################################
if [ -n "$fixing" ] ; then
 head -2  summary.txt > summary.txt.tmp
 mv summary.txt.tmp summary.txt
fi
cat asm.wgs.info | getSummary.pl -i 1 -min 2000  -t scafSeq.wgs -nh >> summary.txt



