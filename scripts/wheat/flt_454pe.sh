#!/bin/sh
for f in PE1sff PE2sff PE3sff PE4sff PE5sff PE6sff PE7sff PE8sff plate1sff plate2sff plate3sff ; do
 #ls -al ${f}_1.sam ${f}_2.sam ${f}_u.sam 
 sam_flt.pl -p93 -l36 -v8 -e 6000 < ${f}_1.sam > ${f}_1.flt.sam &
 sam_flt.pl -p93 -l36 -v8 -e 6000 < ${f}_2.sam > ${f}_2.flt.sam &
 sam_flt.pl -p93 -l36 -v8 -e 6000 < ${f}_u.sam > ${f}_u.flt.sam
done
