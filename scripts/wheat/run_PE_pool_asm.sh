#!/bin/bash -eux
#this MUST be run in the main RUN directory
#with subdirs fastq, filterMC, filterM
R=${PWD##*/}
# parameters: $1=pool#, $2=avg_insert_size
#nonsparse
mkdir -p asm/$1
cd asm/$1

echo "max_rd_len=301" > SOAPdenovo.config
echo "" >> SOAPdenovo.config
echo "[LIB]" >> SOAPdenovo.config
echo "avg_ins=$2" >> SOAPdenovo.config
echo "reverse_seq=0" >> SOAPdenovo.config
echo "asm_flags=3" >> SOAPdenovo.config
echo "map_len=63" >> SOAPdenovo.config
echo "rank=1" >> SOAPdenovo.config
echo "q1=frag_1.fastq" >> SOAPdenovo.config
echo "q2=frag_2.fastq" >> SOAPdenovo.config

#ln -s ~dpuiu/bin/Aegilops_tauschii/SOAPdenovo2.pregraph.sh SOAPdenovo2.sh
set +e
ln -s ../../filterMC/${1}_1.fastq frag_1.fastq
ln -s ../../filterMC/${1}_2.fastq frag_2.fastq
set -e
/home/gpertea/bin/wheat/wheatD_SOAPdenovo2.pregraph.sh $1 >> SD2.pregraph.sh.log 2>&1
cd -
