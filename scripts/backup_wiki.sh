#!/bin/sh
##-- setup
PATH=/bin:/usr/bin:/sbin:/usr/sbin
bn=bck_$(date +%y%m%d_%H-%M)
host=$(hostname -s)
dbuser=$1
dbpass=$2
##!IMPORTANT: do not place '/' at the end of these paths!
#wikidir=/home/html/genomics/igmwiki
wikidir=/data1/igm3/www/ccb.jhu.edu/html/wiki

#rdirs=( "geom.us.to:/mylocal/geo/backups/igmwiki" "igm3.pha.jhu.edu:/home/gpertea/backups/wiki" )
rdirs=( "geom.us.to:/mylocal/geo/backups/ccbwiki" )

ferr=~/backups/wiki/$bn.stderr

cd $wikidir

##- check if LocalSettings.php has changed
files_changed=''
if ! cmp -s ./LocalSettings.php ~/backups/wiki/ccbwiki/LocalSettings.php; then
 cp -p ./LocalSettings.php ~/backups/wiki/ccbwiki/
 files_changed=1
fi

##- make the whole wiki read-only:
nice -n 19 perl -i -pe 's/^[#\s]+(\$wgReadOnly[ =])/$1/' LocalSettings.php
#clean up pChart4mw cache files older than 20 days
cd uploads/pChart4mw_cache
cache_size=$(du -s | cut -f1)
if [[ $cache_size -gt 20000000 ]]; then
  echo "Cleaning up cache.." > $ferr
  ( nice -n19 find . -maxdepth 1 -mtime +20 -exec /bin/rm -rf {} \; ) 2>>$ferr
fi
##### --- see if wiki files have been updated (uploads, extensions, configuration etc.)
cd $wikidir
ff=$bn.ccbwiki_files.tar.bz2
fsync=~/backups/wiki/$bn.rsync.info
(
nice -n19 rsync -Wai --delete --exclude=uploads/pChart4mw_cache/'*' \
 --exclude=LocalSettings.php ./ ~/backups/wiki/ccbwiki > $fsync
) 2>>$ferr
if grep -q -E '^(>f|c)' $fsync 2>/dev/null; then
  #we have changes, build the tar file
  files_changed=1
#else
# rsync -aicI ./LocalSettings.php  ~/backups/wiki/ccbwiki/LocalSettings.php > $fsync
# if grep -q '^>f' $fsync 2>/dev/null; then
#    files_changed=1
# fi
fi
#( nice -n 19 tar cfj ~/backups/wiki/$bfiles uploads \
# images extensions LocalSettings.php ) 2>>$ferr

dbf="$bn.ccbwiki_db.sql.gz"
echo "Backing up database with mysqldump ($dbf).." >> $ferr
( nice -n 19 mysqldump -u "$dbuser" -p"$dbpass" --default-character-set=binary \
   --single-transaction igmwiki -c | nice -n 19 gzip -9 > ~/backups/wiki/$dbf
) 2>>$ferr

### -- xml dump:
xmlf="$bn.ccbwiki_xml.bz2"
echo "Backing up pages as xml ($xmlf).." >> $ferr
( nice -n19 php -d error_reporting=E_ERROR maintenance/dumpBackup.php --current | \
  nice -n19 bzip2 -9 > ~/backups/wiki/$xmlf 
) 2>>$ferr

##- restore wiki to read-write:
nice -n 19 perl -i -pe 's/^\s*(\$wgReadOnly[ =])/## $1/' LocalSettings.php

##- prepare and send the backups to remote backup locations
cd ~/backups/wiki
prevff=$(ls -1t *.ccbwiki_files.tar.bz2 2>/dev/null | head -1)
if [[ -z $prevff ]]; then
 files_changed=1
fi
if [[ $files_changed ]]; then
  echo "Backing up files ($ff).." >> $ferr
  #echo "Files changed, building $ff" >> $ferr
  nice -n 19 tar cfj $ff ccbwiki 2>> $ferr
  c=0
  for fn in $( ls -1t *.ccbwiki_files.tar.bz2 2>/dev/null ) ; do
    (( c++ ))
    if [[ $c -gt 5 ]]; then
      /bin/rm -f $fn
    fi
  done
else
  if [[ $prevff ]]; then
    touch $prevff #make it current
    ff=$prevff
    prevff=''
  fi
fi
# clean up files older than 30 days ONLY IF these backup files are larger
dbfs=$(stat -c'%s' $dbf)
oldbf=$(ls -1t *.ccbwiki_db.sql.gz 2>/dev/null | tail -1)
oldbfs=$(stat -c'%s' $oldbf)

cleanup_old=1
(( oldbfs -= 50000 ))
#echo "Adjusted oldb file size: $oldbfs"

errtext=$(grep -E -i 'error|fail|cannot|space|couldn|fault' $ferr 2>/dev/null)
#echo "errtext='$errtext'"

if [[ $dbfs -lt 1024 || $dbfs -lt $oldbfs || $errtext ]]; then
  cleanup_old=0
  echo "Warning: issues encountered during last backup ($dbf, size $dbfs), see above." >> $ferr
  oldmailer=$(mail -r 2>&1 | grep -F 'invalid option')
  if [[ $oldmailer ]]; then
     mail -s 'ccbwiki backup issue' 'gpertea@jhu.edu' -- -f 'backuper@'$host.jhu.edu < $ferr
  else
     mail -s 'ccbwiki backup issue'  -r 'backuper@'$host.jhu.edu 'gpertea@jhu.edu' < $ferr
  fi
fi

if [[ $cleanup_old -gt 0 ]]; then
  #echo "Cleaning up older backups older than 30 days"
  if [[ $prevff ]]; then 
     /bin/rm -f $prevff
     fi
  nice -n19 find . -maxdepth 1 -name 'bck_*.ccbwiki_*.*' -mtime +30 -exec /bin/rm -rf {} \;
fi
ferr=$bn.ssh.stderr
echo "" > $ferr
scperr=''
for rdest in "${rdirs[@]}" ; do
 scp -q $dbf $rdest/ 2>>$ferr
 if [[ $? != 0 ]]; then scperr=1; fi
 scp -q $xmlf $rdest/ 2>>$ferr
 if [[ $? != 0 ]]; then scperr=1; fi
 if [[ $files_changed ]]; then
   scp -q $ff $rdest/ 2>>$ferr
   if [[ $? != 0 ]]; then scperr=1; fi
 fi
 
 oifs=$IFS
 IFS=':'
 arr=($rdest)
 rhost=${arr[0]}
 rdir=${arr[1]}
 #echo "Host is '$rhost', dir is '$rdir'"
 rbasedir=${rdir%/*} #the remote "backups" directory
 rtdir=${rdir##*/} #target subdirectory under the "backups" directory
 # notify the remote host
 #echo "running: ssh $rhost \"nohup $rbasedir/backup_received.sh '$rdir' '$bn' '$cleanup_old' '$ff' > /dev/null 2>&1\"" >> $ferr
 echo "Notifying host $rhost .." >> $ferr
 ( ssh $rhost "nohup $rbasedir/backup_received.sh '$rdir' '$bn' '$cleanup_old' \
 '$ff'" 2>&1
  ) 1>>$ferr
  
 IFS=$oifs
done
err=$(grep -E -i 'error|fail|cannot|couldn|fault|timeout|discon' $ferr 2>/dev/null)
if [[ $err || $scperr ]]; then
  oldmailer=$(mail -r 2>&1 | grep -F 'invalid option')
  if [[ $oldmailer ]]; then
     mail -s 'ccbwiki backup remote issue' 'gpertea@jhu.edu' -- -f 'backuper@'$host.jhu.edu < $ferr
  else
     mail -s 'ccbwiki backup remote issue'  -r 'backuper@'$host.jhu.edu 'gpertea@jhu.edu' < $ferr
  fi
fi

