#!/bin/sh
ver=$(fgrep '#define VERSION ' fqtrim.cpp)
ver=${ver#*\"}
ver=${ver%%\"*}
src=fqtrim-$ver
echo "preparing $src"
echo "-------------------"
/bin/rm -rf $src
/bin/rm -f $src.tar.gz
mkdir $src
mkdir $src/gclib
sed 's|../gclib|./gclib|' Makefile > $src/Makefile
libdir=fqtrim-$ver/gclib/
cp -p LICENSE README fqtrim.cpp fqtrim-$ver/
cp -p ../gclib/{GVec,GList,GHash}.hh $libdir
cp -p ../gclib/{GAlnExtend,GArgs,GBase,gdna,GStr,GThreads}.{h,cpp} $libdir
tar cvfz $src.tar.gz $src
#scp $src.tar.gz igm3:~/src/
linpack=fqtrim-$ver.Linux_x86_64
/bin/rm -rf $linpack
/bin/rm -rf $linpack.tar.gz
mkdir $linpack
cd $src
if [[ $(uname -m) = "x86_64" ]]; then
 echo "Linking statically on x86_64 (only works for gcc 4.5+).."
 export LDFLAGS="-static-libgcc -static-libstdc++"
fi
if [[ $(uname) = "Darwin" ]]; then
 export CFLAGS="-mmacosx-version-min=10.6"
fi

make release
mv fqtrim README LICENSE ../$linpack/
cd ..
pwd
echo "tar cvfz $linpack.tar.gz $linpack"
tar cvfz $linpack.tar.gz $linpack
echo "On igmX machines you can also update the web files:"
echo "cp $linpack.tar.gz $src.tar.gz  ~/html/software/fqtrim/dl/"
echo "perl -i -pe 's/fqtrim-[0-9]+\.\d+\./fqtrim-$ver./g' ~/html/software/fqtrim/index.shtml"
