#!/bin/bash
if [[ -z $1 || "$1" == "-h" || "$1" == "--help" ]]; then
   echo -e "Usage:\nbam2fastq.sh file.bam [<numcpus>]"
   echo "Default <numcpus> is 8"
   exit 1
   fi
bam="$1"
numcpus="$2"
if [[ -z $numcpus ]]; then
 numcpus=8
fi
fname=${bam##*/}
#echo "bam   = $bam"
#echo "fname = $fname"
bname=${fname%.*}
echo " Sorting bam into $bname.srt.bam .."

samtools-mt sort -n -l9 -@$numcpus -m 9G $bam $bname.srt

echo " Writing fastq files $bname.[12].fq.gz .."

bam2fastx -o $bname.fq.gz -PAN $bname.srt.bam

