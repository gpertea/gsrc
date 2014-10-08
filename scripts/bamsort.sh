#!/bin/sh
if [[ -z $1 || "$1" == "-h" || "$1" == "--help" ]]; then
   echo -e "Usage:\nbamsort.sh <base_bam_name> [<numcpus>]"
   echo "Default <numcpus> is 4"
   exit 1
   fi
pmem="7G"
bname="$1"
numcpus="$8"
if [[ -z $numcpus ]]; then
 numcpus=4
fi
samtools-mt sort -@ $numcpus -m $pmem $bname.bam $bname.locsort

