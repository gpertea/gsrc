#!/bin/sh
#PBS -N gridCmds
#PBS -l walltime=18:00:00
#PBS -j oe
#PBS -m n
#PBS -l nodes=1:ppn=1,pvmem=8gb

### -- Usage advice --
### Make a copy of this script in order to customize 
### the PBS grid options:
###  -N: modify the name for your grid job (instead of "gridCmds")
###  -l walltime=   : run time requested for each job in the array (18h)
### The grid process is single threaded (ppn=1) and uses up to 8GB RAM;
### If different RAM requirements or CPUs per job are desired use 
### a PBS directive like this: -l nodes=1:ppn=8,pvmem=4gb
### (which means 8cpus per node with 4GB for each CPU (pvmem) so in fact 
###  32GB RAM are being requested)

if [[ -z $PBS_O_WORKDIR ]]; then
   #submit mode
   if [[ -z $1 || "$1" == "-h" || "$1" == "--help" ]]; then
     echo -e "Usage:\n$0 <cmdfile>"
     echo "<cmdfile> has the commands to be run on the grid, one command per line"
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
   script=$(readlink -f $0)
   echo -e "PBS command:\nqsub -t 1-$lnum -M $USER@jhu.edu -v $varlist -q batch $script"
   cd ~
   qsub -t 1-$lnum -M $USER@jhu.edu -v $varlist -q batch $script
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
 #cmdline=$(head -$l $grid_cmds_file | tail -1)
 cmdline=$(awk 'NR=='$l'{print;exit}' $grid_cmds_file)
 echo -e ">> Running cmd-$l : $cmdline" 1>&2
 eval $cmdline 1>&2
 ret_code=$?
 if [[ $ret_code != 0 ]]; then
  printf "Warning: nonzero exit status [%d] for cmd: $cmdline" $ret_code
  echo -e "<< cmd-$l End Error." 1>&2
  exit $ret_code
 fi
 echo -e "<< cmd-$l End OK." 1>&2
 exit
fi




