#!/bin/sh
for f in PE1sff PE2sff PE3sff PE4sff PE5sff PE6sff PE7sff PE8sff plate1sff plate2sff plate3sff ; do
echo "Running bowtie2 for $f"
./bwt_lmap.sh $f &> ${f}_bowtie.log
done
