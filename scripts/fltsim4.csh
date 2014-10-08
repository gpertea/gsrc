#!/bin/tcsh -f
if ($1"x" == "x") then
 echo "Usage: fltsim4.csh fasta_files.."
 echo ""
 echo "Must be run from the <organism>/jigsaw directory"
 echo "Expects a *.sim4cds.gff to exist in each of the jigsaw/NN_??/ "
 echo " subdirectories (generated previously by gff_split)"
 echo ""
 echo "Example: fltsim4cc.csh ../ga_??.fa"
 endif

foreach f ( $* )
 set fn=$f:t
 set d=$fn:r
 #echo "fn=$fn, d=$d"
 echo "processing dir: $d"
 cd $d
 if ( -s $d.sim4cds.gff ) then
  gffilter -g $d.fa -a $d.sim4_CDS.gtf $d.sim4cds.gff
  sed 's/\tCDS\t/\texon\t/' < $d.sim4cds.gff > $d.sim4exon.gff
  #perl -ne '@t=split("\t");print $_ if $t[2] eq "mRNA" || $t[2] eq "exon"' \
  #    < $d.sim4cds.gff > $d.refseq_exon.gff3
  else
   touch $d.sim4_CDS.gtf $d.sim4exon.gff
   echo "Warning: sim4cc results not found for $d"
 endif
 cd ..
end

