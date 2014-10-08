#!/bin/tcsh -f

if ( $1"x" == "x" ) then
echo "Usage: start_exonerate_server.csh <gseqdb.fa> [T]"
echo " Use any second parameter to start a translated server."
exit 1
endif

set fa=$1
set tran="n"
if ( $2"x" == "x" ) set tran="t"

set fesd=$fa.esd
set fesi=$fa.esi
if ($tran == "t") set fesi = $fa.trans.esi
if ( ! -f $fesd ) then
 fasta2esd --alphabet dna --softmask no $fa $fesd
endif

if ( ! -f $fesi ) then
 if ($tran == "t") then
   esd2esi --translate yes --proteinwordlen 4 --saturatethreshold 80 $fesd $fesi
 else
   esd2esi $fesd $fesi
 endif
endif

exonerate-server --proteinwordlen 4 --port 3804 $fesi
