#!/bin/sh
pelib=$1
pedir=/scratch0/igm1/dpuiu/Aegilops_tauschii/Data/454_data_for_3DS.part/PEsff

bowtie2 --reorder --mm --local -L16 --no-unal -x 454AllContigs -p7 -U $pedir/$pelib/${pelib}_1.fq -S ${pelib}_1.sam &
pid1=$!
bowtie2 --reorder --mm --local -L16 --no-unal -x 454AllContigs -p7 -U $pedir/$pelib/${pelib}_2.fq -S ${pelib}_2.sam &
pid2=$!
bowtie2 --reorder --mm --local -L16 --no-unal -x 454AllContigs -p9 -U $pedir/$pelib/${pelib}_u.fq -S ${pelib}_u.sam &
pid3=$!

wait $pid1 $pid2 $pid3

echo "bowtie2 finished for $pelib."
