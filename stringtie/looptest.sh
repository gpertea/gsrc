#!/bin/sh
r=0
c=0
while (( $r==0 )); do
(( c++ ))
echo "run #$c.."
set -vx
./stringtie chr20.polya.bam -o chr20.polya.out -p5 >& chr20.polya.log
r=$?
set +vx
echo " >> result $r."
done
