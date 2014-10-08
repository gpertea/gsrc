#!/bin/tcsh -f
set fq=$1
set fesi=$2
set fb=$fesi:r
set fb=$fb:r
set srvlog=srv_$fb.log
set exout=$fb.exout

#goto TEST_LOOP

set port=$GRID_WORKER
@ port += 3804
set sv=`ps auxwww | grep exonerate-server | grep "port $port"`

if ( $#sv > 2 ) then
 echo "Existing server on this machine on port $port, aborting.."
 exit 241
endif

# echo "Launching server.."

exonerate-server --verbosity 0 --maxconnections 1 --proteinwordlen 4 --port $port $fesi >& $srvlog &
set bpid=$!
#echo "Background process ID= $bpid"
while ( `grep listening $srvlog`"x" == "x" )
 sleep 2
end

pexo_slice.pl -S localhost:$port -o $exout $fq

fltexonerate -p 70 -c 50 < $exout > $exout.p70c50.gff3
echo '<Done.' >> /dev/stderr
sleep 2

kill $bpid

exit
#------------------------------


#TEST_LOOP:
#echo "server $fesi started, running exonerate search.."
/bin/rm -f $exout
set params = '-V 0 --model p2g --percent 20 -n 11 --showalignment no --showtargetgff yes --showvulgar no'
set params = "$params --maxintron 400000 --seedrepeat 2 --geneseed 60 -x 50 --proteinwordlen 4 -M 256 "
set params = "$params --ryo "'%ti\t%qi\tmap_end\t%r\t%ql|%qab-%qae\t%ps\t%g\t%qd\n'

set maxc = 2571
@ c = 0 
while ( $c < $maxc ) 
 @ c += 1 
 set excmd = "exonerate $params"
 set excmd = "$excmd --querychunkid $c --querychunktotal $maxc $fq localhost:3804"
 #echo ">run $c of $fesi" >> /dev/stderr
 # exonerate $params --querychunkid $c --querychunktotal $maxc $fq localhost:3804 >> $exout 
 $excmd >> $exout
 if ($status != 0) then
   echo ">Error status for $fesi at --querychunkid $c --querychunktotal $maxc $fq localhost:3804" >> /dev/stderr
 endif 
end

echo ">all $maxc runs completed for $fesi ; Now running fltexonerate on $exout" >> /dev/stderr
#exonerate -V 0 --model p2g --percent 20 -n 11 --showalignment no \
# --maxintron 400000 --seedrepeat 2 --geneseed 60 -x 50 --proteinwordlen 4 -M 256 \
# --ryo '%ti\t%qi\tmap_end\t%r\t%ql|%qab-%qae\t%ps\t%g\t%qd\n' \
# --showtargetgff yes --showvulgar no $fq localhost:3804 > $exout

fltexonerate -p 75 -c 70 < $exout > $exout.p75c70.gff3
echo '<Done.' >> /dev/stderr
sleep 2

kill $bpid
