#!/bin/sh
lab="$1"
lib="$2"
bwtIdx="$3"
threads="$4"

if [ -z "$lab" ]; then
 echo -e "Usage\n wgs_map_bowtie_PEflt.sh {CSHL|BGI} <library> <bowtie_index> <num_threads>"
 exit 1
fi
if [[ "$lab" != "BGI" && "$lab" != "CSHL" ]]; then
 echo "Error: BGI or CSHL must be the 1st parameter"
 echo -e "Usage\n wgs_map_bowtie.sh {CSHL|BGI} <library> <bowtie_index> <num_threads>"
 exit 1
fi
WGS=WGS_$lab
src=/scratch0/igm1/gpertea/wheat/3DS

if [[ $lib && $lib == "list" ]]; then
 head -1 $src/$WGS/info.txt
 grep -v '^#' $src/$WGS/info.txt | sort -k6,6 -nr | more
 exit 1
fi

if [ -z "$threads" ]; then
 echo -e "Usage\n wgs_map_bowtie.sh {CSHL|BGI} <library> <bowtie_index> <num_threads>"
 exit 1
fi

cwd=$(pwd -P)
outdir=$cwd/bwt_wgsmap_${lab}_$lib

#bwtIdx=$cwd/$bwtIdx
bwtIdx=$(readlink -f $bwtIdx)

if [ ! -f $bwtIdx.1.bt2 ]; then
  echo "Error: cannot find bowtie index $bwtIdx"
  exit 2
fi

mgap=$(fgrep -w $lib $src/$WGS/info.txt | cut -f5)
rlen=$(fgrep -w $lib $src/$WGS/info.txt | cut -f4)
if [ ! -d "$outdir" ]; then
 mkdir -p $outdir
fi

fgrep -w $lib $src/$WGS/info.txt > $outdir/lib_info.txt

#mgap=$(awk "BEGIN {printf \"%d\", $mgap * 0.4 + $mgap - $rlen * 2}")
fraglen=$(awk "BEGIN {printf \"%d\", $mgap * 0.4 + $mgap}")
#cseg=$(awk "BEGIN {printf \"%d\", $rlen * 0.8}")
pairedtype=""
if (( $fraglen > 2000 )); then
 pairedtype="--rf"
fi

echo "fraglen = $fraglen"
echo "paired type = $pairedtype"
#echo "mgap=$mgap" >> $outdir/lib_info.txt
#echo "cseg=$cseg" >> $outdir/lib_info.txt

machine=$(uname -n)
cd $outdir
bowtie_options="--reorder --local --mm -L 17 --no-unal \
 -X $fraglen $pairedtype -x $bwtIdx -p $threads "
echo "Running: bowtie2 $bowtie_options \
-1 $src/$WGS/$lib"'_1.fastq' -2 $src/$WGS/$lib'_2.fastq' "| sam_flt -p94 -v8 -l42 > bwt_$lib.sam" > bowtie.cmd
bowtie2 $bowtie_options \
-1 $src/$WGS/$lib'_1.fastq' -2 $src/$WGS/$lib'_2.fastq' | sam_flt -p94 -v8 -l42 > bwt_$lib.sam 2>bowtie.log
# & 
#pid1=$!
#wait $pid1
#echo "Bowtie finished in $outdir on $machine" | mailx -s "Bowtie $lab $lib on $bwtIdx finished on $machine in $outdir" 'gpertea@jhu.edu'
