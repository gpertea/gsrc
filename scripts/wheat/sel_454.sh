#!/bin/sh
cdir=$(pwd)
for f in PE1sff PE2sff PE3sff PE4sff PE5sff PE6sff PE7sff PE8sff plate1sff plate2sff plate3sff ; do
 #ls -al ${f}_1.sam ${f}_2.sam ${f}_u.sam 
 #cat ${f}_1.flt.sam ${f}_2.flt.sam ${f}_u.flt.sam | grep -v '^@' | cut -f1 | cut -f1 -d '/' | sort -u > $f.sel.lst
 cat ${f}_1.flt.sam ${f}_2.flt.sam | grep -v '^@' | cut -f1 | cut -f1 -d '/' | sort -u > $f.selpe.lst
 cd ../../PEreads/$f/
 sfffile -i $cdir/$f.selpe.lst -o $cdir/PEonly/$f.selpe.sff *.sff
 cd $cdir
done

#sfffile -o allPEsel.sff *.sel.sff
