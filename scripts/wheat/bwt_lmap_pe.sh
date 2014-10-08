#!/bin/sh
pelib=$1
pedir=/scratch0/igm1/dpuiu/Aegilops_tauschii/Data/454_data_for_3DS.part/PEsff

bowtie2 --reorder --mm --local --no-unal -L16 --fr -p12 -x 454AllContigs \
  -1 $pedir/$pelib/${pelib}_1.fq -2 $pedir/$pelib/${pelib}_2.fq \
  -S ${pelib}_PE.sam &> bwt_${pelib}.log  &
pid1=$!

wait $pid1

echo "bowtie2 finished for $pelib."
