#!/bin/sh
ver=$(fgrep '#define VERSION ' fqtrim.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
pack=fqtrim-$ver
echo "preparing $pack"
echo "-------------------"
/bin/rm -rf $pack
/bin/rm -f $pack.tar.gz
mkdir $pack
mkdir $pack/gclib
sed 's|../gclib|./gclib|' Makefile > $pack/Makefile
libdir=fqtrim-$ver/gclib/
cp -p LICENSE README fqtrim.cpp fqtrim-$ver/
cp -p ../gclib/{GVec,GList,GHash}.hh $libdir
cp -p ../gclib/{GAlnExtend,GArgs,GBase,gdna,GStr,GThreads}.{h,cpp} $libdir
tar cvfz $pack.tar.gz $pack
scp $pack.tar.gz igm3:~/src/
