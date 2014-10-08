#!/bin/sh
bdir=$1
if [[ ! $bdir ]]; then
 echo "Error: no parameters given!"
 exit 3
fi
if [[ ! -d $bdir ]]; then
 mkdir -p $bdir
 if [[ $? -gt 0 ]]; then
   echo "Failed to create directory $bdir" 1>&2
   exit 3
 fi
fi
shift
bn=$1
shift
cleanup_old=$1
shift
ferr=$bdir/$bn.stderr
cd $bdir 2>>$ferr
echo "parameters: $bdir $bn $cleanup_old $@" >> $ferr
for f in "$@" ; do
 # these files are bigger and should have slower turnover -- keep only the last 5 archives
 if [[ $f = $bn* ]]; then
   fb="${f/$bn/}"  # created by current cron job
   c=0
   for fn in $( ls -1t *$fb ) ; do
    (( c++ ))
    if [[ $fn != $f  && $c -gt 5 ]]; then
      #echo "$fn  (..to delete)"
      /bin/rm -f $fn
    #else
    #  echo "$fn <-- will NOT be deleted!"
    fi
   done
 else
   #not a current file, but from a previous backup
   if [[ -f $f ]]; then 
      touch $f
   else # file not found!
    host=$(hostname -s)
    echo "Error on $host: file $f not found among current backups!" >> $ferr
   fi
 fi
done
numfiles=$(ls -1 | wc -l)

if [[ $numfiles -lt 6 ]]; then
 #don't cleanup when there are no previous backups
 cleanup_old=0
fi

if [[ $cleanup_old -gt 0 ]]; then
  #echo "Cleaning up backup files older than 30 days.." 2>>$ferr
  nice -n19 find . -maxdepth 1 -mtime +30 -exec /bin/rm -rf {} \;
fi
cat $ferr
