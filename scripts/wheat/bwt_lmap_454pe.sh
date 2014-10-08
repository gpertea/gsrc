#!/bin/sh
pelib=$1
pedir=/scratch0/igm1/dpuiu/Aegilops_tauschii/Data/454_data_for_3DS.part/PEsff
#cd /home/gpertea/work/wheat/3DS/pool42/asm
bowtie2 --reorder --mm --local --no-unal -L17 --fr -X 6500 -p8 -x 454AllContigs \
  -1 $pedir/$pelib/${pelib}_1.fq -2 $pedir/$pelib/${pelib}_2.fq \
  | sam_flt -p94 -v8 -l42 > ${pelib}_PE.sam 2>bwt_${pelib}.log
#pid1=$!

#wait $pid1

#echo "bowtie2 finished for $pelib."
