#!/bin/tcsh -f

exonerate $1 localhost:3043 --model p2g -V 0 --showalignment no --showvulgar no -n 5 --showtargetgff yes \
--seedrepeat 1 --geneseed 250 |& fltexonerate

#\t%ps
#exonerate $1 localhost:3043 --model p2g -V 0 --showalignment no --showvulgar no -n 11 --seedrepeat 1 --geneseed 180 \
#--showcigar yes --showtargetgff yes --percent 70 

#--ryo '%ti\t%qi\tmap_end\t%r\t%et/%ql\t%ps,%pi\t%g\t%em\n'

#--ryo '%ti\t%qi\tmap_end\t%r\t%ql|%qab-%qae\t%ps\t%g\t%qd\n'

