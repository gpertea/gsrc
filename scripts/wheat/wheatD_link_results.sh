#!/bin/bash -eu
#must be run in the main run directory, with ./asm/<pool#>/ subdirectories
test -f pools_count_mea_std.info
test -f pools_mea_std.info

RUNDIR=$(pwd -P)

#Daniela's data dir:
DATADIR=/scratch0/igm1/dpuiu/Aegilops_tauschii/Data
#Daniela's assembly dir:
ASMDIR=/scratch0/igm1/dpuiu/Aegilops_tauschii/Assembly

RUN=${PWD##*/}
export R=`echo ${PWD##*/} | sed 's|run||' | sed 's|_| |' | awk '{print $1}'`
N=`cat pools_count_mea_std.info | grep -v ^# | wc -l`
RN="$R:$N"

export D=`date "+%m%d%Y"`

if [ -d $DATADIR/$RUN ]; then
  echo "Oops: Run directory $RUN already exists in $DATADIR"
  echo "Sort this out first."
  exit 1
fi

if [ -d $ASMDIR/$RUN ]; then
  echo "Oops: Run directory $RUN already exists in $ASMDIR"
  echo "Sort this out first."
  exit 1
fi
declare -a pools=( )
cd asm
for asm in ????/asm.fa; do
 pool=${asm%/*}
 pools=("${pools[@]-}" $pool)
done

mkdir $DATADIR/$RUN
mkdir $ASMDIR/$RUN

cd $DATADIR/$RUN
ln -s $RUNDIR/{original,fastq,filterM,filterMC,pools*.info} .

cd $ASMDIR/$RUN
ln -s $RUNDIR/{pools_count_mea_std.info,*summary*.txt,asm*.fa,asm*.info,asm*.len} .

for pool in ${pools[@]}; do
  mkdir $pool
  cd $pool
  ln -s $RUNDIR/asm/$pool SOAPdenovo2
  cd ..
done


echo "Done here."
exit 0

