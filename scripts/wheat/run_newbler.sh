#!/bin/bash
pnum=$1
inf=$2
odir=$3
if [[ -z $odir ]]; then
   echo -e "Usage:\n newbler_run.sh <num_cpus> <input_file> <out_dir> [other_runAssembly_options]"
   echo -e " e.g.:\n newbler_run.sh 8 reads.sff asm_urt -urt"
   exit 1
fi
shift 3
runAssembly -o $odir -cpu $pnum -tr -vt /scratch0/igm1/dpuiu/Aegilops_tauschii/Data/RefSeq/pCC1BAC.fasta \
 -vs /scratch0/igm1/dpuiu/Aegilops_tauschii/Data/RefSeq/EcoliPuc18.fasta -rip -ml 42 -mi 94 -ace "$@" $inf &
echo "runAssembly of $inf started, output will be in $odir"
