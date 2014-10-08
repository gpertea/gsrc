#!/bin/sh
rm cuffcompare.hprof*
export HEAP_PROFILE_ALLOCATION_INTERVAL=50000000
export HEAPPROFILE=$PWD/cuffcompare.hprof
# enabling the heap cheker will disable heap profiling
# export HEAPCHECK=normal
exec ./cuffcompare "$@"
