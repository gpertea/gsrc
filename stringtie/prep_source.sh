#!/bin/sh
ver=$(fgrep '#define VERSION ' stringtie.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
pack=stringtie-$ver
echo "preparing $pack.tar.gz"
echo "-------------------"
/bin/rm -rf $pack
/bin/rm -f $pack.tar.gz
mkdir $pack
mkdir $pack/gclib
#sed 's|../gclib|./gclib|' Makefile > $pack/Makefile
cd samtools-0.1.18
make clean
cd ..
libdir=stringtie-$ver/gclib/
cp -p LICENSE README Makefile stringtie.cpp rlink.cpp rlink.h $pack/
cp -pr samtools-0.1.18 $pack/
/bin/rm -rf $pack/samtools-0.1.18/.svn
cp -p ./gclib/{GVec,GList,GHash}.hh $libdir
cp -p ./gclib/GBitVec.h $libdir
cp -p ./gclib/{GArgs,GBam,GBase,gdna,GStr,gff,codons,GFaSeqGet,proc_mem,GThreads}.{h,cpp} $libdir
tar cvfz $pack.tar.gz $pack
ls -l $pack.tar.gz
#scp $pack.tar.gz igm3:~/src/
