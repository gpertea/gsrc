#!/bin/bash -x

QRY=$1
QR1=${QRY}_1.fastq
QR2=${QRY}_2.fastq
BASEDIR=/home/dpuiu/Aegilops_tauschii/Data/RefSeq
ADAPTER=$BASEDIR/adapters.34bp.fasta
CONTAMINANT=$BASEDIR/bwa-indeces/contaminants.fasta

INDIR=fastq
OUTDIR0=sample
OUTDIR1=filterM
OUTDIR2=filterMC

#P=6
P=16
T=0.001
Q=20
L=32
K=31

test -f ${INDIR}/${QR1}
test -f ${INDIR}/${QR2}
test -f ${ADAPTER}

################################################################################################

#mkdir -p ${OUTDIR0}
#~/bin/FASTQ/sampleFastq12.pl -c 10000 ${INDIR}/${QR1} ${INDIR}/${QR2} ${OUTDIR0}/${QR1} ${OUTDIR0}/${QR2}
#fastx_quality_stats -Q 33 -i ${OUTDIR0}/${QR1} > ${OUTDIR0}/${QRY}_1.stats
#fastx_quality_stats -Q 33 -i ${OUTDIR0}/${QR2} > ${OUTDIR0}/${QRY}_2.stats

mkdir -p ${OUTDIR1}
fastq-mcf -t ${T} -q ${Q} ${ADAPTER} ${INDIR}/${QR1} ${INDIR}/${QR2} -o ${OUTDIR1}/${QR1} -o ${OUTDIR1}/${QR2} -l ${K}
set +x
echo "+ mkdir -p $OUTDIR2" >&2
mkdir -p ${OUTDIR2}
echo "+ bwa mem -k $K -t $P -S $CONTAMINANT $OUTDIR1/$QR1 $OUTDIR1/$QR2 | sam-filter.pl -nh > $OUTDIR2/$QRY.sam" >&2
bwa mem -k $K -t $P -S ${CONTAMINANT} ${OUTDIR1}/${QR1} ${OUTDIR1}/${QR2} \
| sam-filter.pl -nh > ${OUTDIR2}/${QRY}.sam
set -x
extractFromFastqNames.pl -f ${OUTDIR2}/${QRY}.sam ${OUTDIR1}/${QR1} -n > ${OUTDIR2}/${QR1}
extractFromFastqNames.pl -f ${OUTDIR2}/${QRY}.sam ${OUTDIR1}/${QR2} -n > ${OUTDIR2}/${QR2}

