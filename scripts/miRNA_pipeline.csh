#!/bin/tcsh -f

if ('x'$1 == 'x') then
 echo "Usage: miRNA_pipeline.csh [-f <prefix>] [-q] [-K] input.fq"
 echo ""
 echo "Process a fastq file with miRNA reads"
 echo "\nUse -f option only if the input file is not filtered/dusted"
 echo "Use -q option if the quality values are already phred33"
 echo " (otherwise they are assumed to be phred64)"
 exit 1
 endif
echo "# Command line was:\n$0 $argv:q"
set mdir=/fs/szannotation/microRNA
set unfiltered=0
set q33=0
set skipbowtie=0
OPT_LOOP:
if ('x'$1 == 'x-f') then
  set unfiltered=1
  shift
  set rprefix=$1
  shift
endif
if ('x'$1 == 'x-q') then
 set q33=1
 shift
 goto OPT_LOOP
endif
if ('x'$1 == 'x-K') then
 set skipbowtie=1
 shift
 goto OPT_LOOP
endif

set f=$1
if ('x'$f == 'x') then
  echo "Usage: miRNA_pipeline.csh [-f <prefix>] input.fq"
  echo ""
  exit 1
endif
if ( ! -f $f ) then
 echo "Error: file $f not found!"
 exit 1
endif

set fname=$f:r
set fa=$fname.fa
#set fname=$fname:r

if ($unfiltered == 1 || $skipbowtie == 1 ) then
 set fname = $fname.flt
 set ff=$fname.fq
 if ($skipbowtie == 1) then
   set f=$ff
   set fa=$fname.fa
   goto REFILTER
   endif
 if (-f $ff) then
   echo "Error: file $ff exists. Remove the file and run this again."
   exit 1
   endif
  echo "Running fqtrim.." 
  #echo "fqtrim -5 CGACAGGTTCAGAGTTCTACAGTCCGACGATC -3 TCGTATGCCGTCTTCTGCTTG -l16 -n $rprefix -DCQ -r $fname.trashed.tab -o $ff $f"
  # MUST REMOVE the -Q option here if the input fq is already phred33 
  set topt='DCQ'
  if ($q33 == 1) set topt='DC'
  fqtrim -5 CGACAGGTTCAGAGTTCTACAGTCCGACGATC -3 TCGTATGCCGTCTTCTGCTTG -l16 -n $rprefix -$topt -r $fname.trashed.tab -o $ff $f
  set f=$ff
  echo "        fqtrim done."
endif

set fa=$fname.fa
# convert to fasta
fq2fa $f > $fa
echo "..converted to fasta file $fa"
# scanning known pre-miRNAs
echo "running mgblast vs known microRNA precursors.."
mgblast -d $mdir/miRNA_hsa_pre.fa -i $fa -W12 -D5 -X4 -y8 -FF -p93 | \
    sort -k11,11g -k10,10nr -k9,9nr > $fname.vs_pre.mtab
mgbl_miRNA_cov.pl -P -B -f $fname.pre.mcov.tab -o $fname.pre.mcov $fname.vs_pre.mtab
# scanning known mature miRNAs
echo "running mgblast vs known mature microRNAs.."
mgblast -d $mdir/miRNA_hsa_mature_and_star.fa -i $fa -W12 -D5 -X4 -y8 -FF -p93 | \
   sort -k11,11g -k10,10nr -k9,9nr -k12,12r > $fname.vs_mature.mtab
mgbl_miRNA_cov.pl -U -B -o $fname.mature.mcov $fname.vs_mature.mtab
sort -k3,3nr -o $fname.mature.mcov $fname.mature.mcov
echo ">> $fname.mature.mcov << has the known miRNA expression estimates"
cdbfasta -Q $f
cut -f1 $fname.pre.mcov.tab >$fname.pre.mcov.lst
cdbyank -l $f.cidx | fgrep -v -w -f $fname.pre.mcov.lst | \
    cdbyank $f.cidx > $fname.novel.fq
#echo "Starting bowtie search: "
#echo "bowtie -l16 -v3 -a --best --strata -m30 -p2 --mm hg19 $fname.novel.fq $fname.novel.map"
#bowtie -l16 -v3 -a --best --strata -m30 -p2 --mm hg19 $fname.novel.fq $fname.novel.map
echo "bowtie -l16 -v3 -a --best --strata -m30 -p4 --mm hg19 $fname.novel.fq $fname.novel.map"
bowtie -l16 -v3 -a --best --strata -m30 -p4 --mm hg19 $fname.novel.fq $fname.novel.map
#apply a filter to allow only mappings with mismatches at the 3' end, and sort these entries acoordingly:
REFILTER:
echo " .. filter and sort bowtie mappings.."
rnaseq_mapflt.pl $fname.novel.map | sort -k3,3 -k4,4n > $fname.novel_sorted.map
set totalflt=`cdbyank -l $f.cidx | perl -e '$n=0; while (<>) { ($m)=(m/_x(\d+)$/); $n+=$m } print $n'`
echo "\nTotal reads left after noise/dust/adapter filtering: $totalflt\n(  $fname.mature.mcov:$totalflt  )\n"
echo ">> $fname.novel_sorted.map << has filtered and sorted mappings for novel reads"
#echo "TODO: these mappings must be filtered by all human annotation coordinates"
#echo "    (to check if they overlap any other known RNA features)"
annbygff -r ../human_annotation.non-miRNA.gff3 -o $fname.novel_annRNA.mapann -n $fname.novel_nonann.map $fname.novel_sorted.map

echo "$fname.novel_nonann.map <= has mappings that do not match any other RNA annotations"
echo "..looking for novel miRNAs.."
set fnovel=$fname:r
set fnovel=$fnovel"_novel"
mirfind_excise.pl -o $fnovel /fs/szdata/genomes/human_hg19/all_chrs.fa.cidx $fname.novel_nonann.map
mirfind_score.pl $fnovel.xc.fa $fnovel.xc.xmap > $fnovel.mirfind
# ---- compare them after
# miRNA_cmp.pl -o miRNA_1vs2.cmptab -d diff_exp.cmptab \
#   s1.flt.mature.mcov:6618663 s2.flt.mature.mcov:7297401 >& diff_exp.cmptab.log 
##-- and then: 
# sort -k6,6g -k7,7nr -o diff_exp.cmptab diff_exp.cmptab
# -- to compare novel miRNA prediction sets:
# mirfind_common.pl -t 5478035:6100493 s1_novel.new.mirfind s2_novel.new.mirfind > s1_s2.novel.common.tab
# cut -f1 s1_s2.novel.common.tab | perl -ne 'chomp;@t=split(/\|/);($a,$b)=($t[2]=~m/(\d+)_(\d+)/);print join("\t",$_,$t[0], $t[1], $a, $b)."\n"' > s1_s2.novel.hairpins.tab
# annbygff -N -g ~/ann/Human_gene_data/human_genes_descr.txt -r ~/ann/Human_gene_data/human_mRNAs_only.gff3 -T -o s1_s2.novel_pred.ovlmRNA.tab -n s1_s2.novel_pred.anntab s1_s2.novel.hairpins.tab
#cut -f1 s1_s2.novel_pred.ovlmRNA.tab | fgrep -v -f- s1_s2.novel.common.tab | sort -k1,1 > s1_s2.novel.common.fltsort.tab
# cut -f1,6 s1_s2.novel_pred.anntab > s1_s2.novel.loc2ann
#add_descr_by_join.pl s1_s2.novel.loc2ann s1_s2.novel.common.fltsort.tab | tr -d '"' > s1_s2.novel.tab
