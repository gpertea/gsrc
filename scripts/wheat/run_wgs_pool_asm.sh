#!/bin/bash -eux
#must be in the top run directory for assembly
cd "asm/$1"

cat SOAPdenovo.config  | sed 's|map_len=63|map_len=75|' > SOAPdenovo.wgs.config
cat /home/gpertea/bin/wheat/SOAPdenovo.wgs.config >> SOAPdenovo.wgs.config
wheatD_SOAPdenovo2.wgs.sh "$1" > SD2.wgs.sh.log 2>&1

cd -
