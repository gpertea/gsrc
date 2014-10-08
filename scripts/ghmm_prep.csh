#!/bin/tcsh

# training dir for glimmerhmm:
set traindir='/fs/szgenefinding/devel/GlimmerHMM/trained_dir/human/train_08_05'
set ghmm_cmds='glimmerhmm_cmds.bat'
if ($1'x' == 'x') then
 echo "Usage: ghmm_prep.csh <fasta_file(s)>..."
 exit 1
endif

cat > run_twinscan.csh << 'EOT'
#/bin/tcsh -f
set plat=`uname -p`
set tsdir='fs/sz-user-supported/Linux-'$plat
if (-d "/sz$tsdir") then
 set tsdir="/sz$tsdir/twinscan"
else
 set tsdir="/$tsdir/twinscan"
endif

set hmm=$tsdir/parameters/human_iscan-9993-genes-09-13-2004.zhmm
foreach inf ( *.isoFa.[1-9] *.isoFa.[1-9][0-9] )
$tsdir/bin/iscan $hmm $inf > $inf.twsgtf
end
'EOT'

##split each of the input sequence files (one seq per .fa file)
 set fprots=prots_to_search.pfa
 # /bin/rm -rf $fprots
 /bin/rm -rf $ghmm_cmds

foreach fa ( $* )
 echo "processing $fa related files.." 
 set nf=$fa:s/.fa/.split.fa/
 splitbyNs.pl $fa $nf 1
 set totnum=1
 set bname=$fa:r
  foreach sf ( $nf.* )
    set isof=$sf:s/.split.fa/.split.iso/
    set ofa=$sf:s/.split.fa/.isoFa/
    set ofa=$ofa:r
    isoscan ~/ann/bin/isoscan_human.hmm $sf 200000 > $isof
 
    
    iso2fasta.pl $sf $isof $ofa $totnum
    
    set inum=`grep -c -v '^#' $isof`
    @ totnum+=$inum
  end
 #-- generate protein domains file
 foreach ifa ( $bname.isoFa.[1-9] $bname.isoFa.[1-9][0-9] )
  set ifnum=$ifa:e; # the number - extension
  set ftag=$bname.$ifnum
  tran6frames.pl -t $ftag $ifa >> $fprots
  set idom = $ifa:s/.isoFa./.pdom./
  set goutpd = $ifa:s/.isoFa./.isoFa.ghmm.pdom./
  set gout = $ifa:s/.isoFa./.isoFa.ghmm./
  echo "glimmerhmm $ifa $traindir -f -p $idom -g -o $goutpd" >> $ghmm_cmds
  echo "glimmerhmm $ifa $traindir -f -g -o $gout" >> $ghmm_cmds
 end
end

echo "All protein translations are in file: $fprots"
echo "Use this file for pfam searches, then run:"
echo "cat gridx-*/wrk_*/*.pfam | map_pfam_domains.pl"
echo "..then run the commands in: $ghmm_cmds"
