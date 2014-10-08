#/bin/tcsh -f
set plat=`uname -p`
set tsdir='fs/sz-user-supported/Linux-'$plat
if (-d "/sz$tsdir") then
 set tsdir="/sz$tsdir/twinscan"
else
 set tsdir="/$tsdir/twinscan"
endif

set hmm=$tsdir/parameters/human_iscan-9993-genes-09-13-2004.zhmm
foreach inf ( *.isoFa.* )
$tsdir/bin/iscan $hmm $inf > $inf.twsgtf
end
