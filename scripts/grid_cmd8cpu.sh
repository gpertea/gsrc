#!/bin/sh
#PBS -N gridCmds
#PBS -l walltime=22:00:00
#PBS -j oe
#PBS -m n
#PBS -l nodes=1:ppn=8,pvmem=4gb
if [[ -z $PBS_O_WORKDIR ]]; then
   #submit mode
   if [[ -z $1 || "$1" == "-h" || "$1" == "--help" ]]; then
     echo -e "Usage:\ngrid_cmd8cpu.sh <cmdfile> [<work_path>]"
     exit 1
   fi
   cmdfile="$1"
   wrkpath="$2"
   echo "Submitting commands in $cmdfile as an array job.."
   lnum=$(wc -l < $cmdfile)
   export grid_cmds_file=$(readlink -m "$cmdfile")
   varlist="grid_cmds_file"
   if [[ -z $wrkpath ]]; then
     wrkpath="$PWD"
   fi
   export grid_wrk_dir="$wrkpath"
   #cp $cmdfile $wrkpath     
   varlist="$varlist,grid_wrk_dir"
   echo -e "PBS command:\nqsub -t 1-$lnum -M $USER@jhu.edu-v $varlist -q batch $0"
   cd ~
   qsub -t 1-$lnum -M $USER@jhu.edu -v $varlist -q batch $0
   exit
fi
if [[ $grid_wrk_dir ]]; then
 cd $grid_wrk_dir
else
 cd $PBS_O_WORKDIR
fi
#env
#echo "my var=$grid_cmds_file (array id = $PBS_ARRAYID)"
if [[ "$grid_cmds_file" && $PBS_ARRAYID ]]; then
 #on the job
 l=$PBS_ARRAYID
 cmdline=$(head -$l $grid_cmds_file | tail -1)
 echo -e ">> cmd$l : $cmdline" 1>&2
 eval $cmdline 1>&2
 echo -e "<< cmd$l End." 1>&2
 exit
fi




