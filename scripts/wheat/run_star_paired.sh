#!/bin/sh
lab="$1"
lib="$2"
starIdx="$3"
threads="$4"
#starIdx=$cwd/star_merged

if [ -z "$lab" ]; then
 echo -e "Usage\n run_star.sh {CSHL|BGI} <library> <star_genome> <num_threads>"
 exit 1
fi
if [[ "$lab" != "BGI" && "$lab" != "CSHL" ]]; then
 echo "Error: BGI or CSHL must be the 1st parameter"
 echo -e "Usage\n run_star.sh {CSHL|BGI} <library> <star_genome> <num_threads>"
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
 echo -e "Usage\n run_star.sh {CSHL|BGI} <library> <star_genome> <num_threads>"
 exit 1
fi

cwd=$(pwd -P)
outdir=$cwd/star_paired.${lab}_$lib

#starIdx=$cwd/$starIdx
starIdx=$(readlink -f $starIdx)

if [ ! -d $starIdx ]; then
  echo "Error: cannot find dir $starIdx"
  exit 2
fi

mgap=$(fgrep -w $lib $src/$WGS/info.txt | cut -f5)
rlen=$(fgrep -w $lib $src/$WGS/info.txt | cut -f4)
if [ ! -d "$outdir" ]; then
 mkdir -p $outdir
fi

fgrep -w $lib $src/$WGS/info.txt > $outdir/lib_info.txt

#mgap=$(echo "scale = 0; $mgap * 0.2 + $mgap - $rlen * 2" | bc)
mgap=$(awk "BEGIN {printf \"%d\", $mgap * 0.4 + $mgap - $rlen * 2}")
cseg=$(awk "BEGIN {printf \"%d\", $rlen * 0.8}")

echo "mgap=$mgap" >> $outdir/lib_info.txt
echo "cseg=$cseg" >> $outdir/lib_info.txt
#threads=$(awk "BEGIN {printf \"%d\", $threads / 2}")
#if [ ! -L star_$lib ]; then
# ln -s $d .
#fi

#STAR --genomeDir $starIdx --genomeLoad LoadAndKeep --runThreadN $threads \
# --alignIntronMax 1 --alignMatesGapMax $mgap --outFileNamePrefix out_${lib}. \
# --outSAMprimaryFlag AllBestScore --outSAMreadID Number --outSAMorder PairedKeepInputOrder \
# --outFilterScoreMinOverLread 0.85 --outFilterMatchNminOverLread 0.92 --chimSegmentMin $cseg \
# --outFilterMismatchNoverLmax 0.07 --outFilterMismatchNmax 5 --outFilterMultimapScoreRange 6 \
# --readFilesIn $src/$WGS/$lib"_1.fastq" $src/$WGS/$lib"_2.fastq" &> star_$lib.log &
machine=$(uname -n)
cd $outdir
star_options="--genomeDir $starIdx --genomeLoad LoadAndKeep --runThreadN $threads \
 --alignIntronMax 1 --alignIntronMin 2 --alignMatesGapMax $mgap --outFileNamePrefix alnPair. \
 --outSAMprimaryFlag AllBestScore --outSAMreadID Number --outSAMorder PairedKeepInputOrder \
 --outSAMattributes All --outFilterScoreMinOverLread 0.8 --outFilterMatchNminOverLread 0.8 \
 --outFilterMismatchNoverLmax 0.05 --outFilterMismatchNmax 50 \
 --outFilterMultimapNmax 40 --outFilterMultimapScoreRange 2 "

STAR $star_options \
 --readFilesIn $src/$WGS/$lib"_1.fastq" $src/$WGS/$lib"_2.fastq" &> star.log &
pid1=$!
wait $pid1
echo "STAR finished in $outdir on $machine" | mailx -s "STAR $lab $lib finished on $machine" 'gpertea@jhu.edu'
