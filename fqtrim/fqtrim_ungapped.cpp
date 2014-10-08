#include "GArgs.h"
#include "GStr.h"
#include "GHash.hh"
#include "GList.hh"
#include <ctype.h>

#define USAGE "Usage:\n\
fqtrim [-5 <5adapter>] [-3 <3adapter>] [-a <min_matchlen>] [-p {64|33}] [-q <minq> [-t <trim_max>]]\\\n\
   [-n <rename_prefix>] [-o <outsuffix>] [-z <zcmd>] [-r <discarded.lst>]\\\n\
   [-l <minlen>] [-C] [-D] [-Q] <input.fq>[,<input_mates.fq>\n\
 \n\
 Trim low quality bases at the 3' end, optionally trim adapter sequence, filter\n\
 for low complexity and collapse duplicate reads\n\
 If read pairs should be trimmed and kept together (i.e. without discarding\n\
 one read in a pair), the two file names should be given delimited by a comma\n\
 or a colon character\n\
\n\
Options:\n\
-n  rename all the reads using the <prefix> followed by a read counter;\n\
    if -C option was given, the suffix \"_x<N>\" is appended, with <N> being\n\
    the read duplication count\n\
-o  write the trimmed/filtered reads to file(s) named <input>.<outsuffix>\n\
    which will be created in the current (working) directory\n\
-5  trim the given adapter or primer sequence at the 5' end of each read\n\
    (e.g. -5 CGACAGGTTCAGAGTTCTACAGTCCGACGATC)\n\
-3  trim the given adapter sequence at the 3' end of each read\n\
    (e.g. -3 TCGTATGCCGTCTTCTGCTTG)\n\
-a  minimum bases to match to adaptor sequence (default 5)\n\
-q  trim bases with quality value lower than <minq> (starting at the 3' end)\n\
-t  for -q option, maximum trimming at the 3' end is limited to <trim_max>\n\
-m  maximum percentage of Ns allowed in a read after trimming (default 7)\n\
-l  minimum \"clean\" length after trimming that a read must have\n\
    in order to pass the filter (default: 16)\n\
-r  write a simple \"trash report\" file listing the discarded reads with a\n\
    one letter code for the reason why they were trashed\n\
-D  apply a low-complexity (dust) filter and discard any read that has over\n\
    50% of its length detected as low complexity\n\
-C  collapse duplicate reads and add a prefix with their count to the read\n\
    name (e.g. for microRNAs)\n\
-p  input is phred64/phred33 (use -p64 or -p33)\n\
-Q  convert quality values to the other Phred qv type\n\
-V  verbose processing\n\
 "
 
//-z  for -o option, the output stream(s) will be first piped into the given\n
//   <zcmd> command, which must output to stdout (e.g. -z 'bzip2 -9 -c')\n

 
// example 3' adapter for miRNAs: TCGTATGCCGTCTTCTGCTTG

//For pair ends sequencing:
//3' : ACACTCTTTCCCTACACGACGCTCTTCCGATCT
//5' : GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG
//FILE* f_out=NULL; //stdout if not provided
//FILE* f_out2=NULL; //for paired reads
//FILE* f_in=NULL; //input fastq (stdin if not provided)
//FILE* f_in2=NULL; //for paired reads
FILE* freport=NULL;

bool debug=false;
bool verbose=false;
bool doCollapse=false;
bool doDust=false;
bool fastaOutput=false;
bool trashReport=false;
//bool rawFormat=false;
int min_read_len=16;
double max_perc_N=7.0;
int dust_cutoff=16;
bool isfasta=false;
bool convert_phred=false;
GStr outsuffix; // -o 
GStr adapter3;
GStr adapter5;
GStr prefix;
GStr zcmd;
int num_trimmed5=0;
int num_trimmed3=0;
int num_trimmedN=0;
int num_trimmedQ=0;
int min_trimmed5=INT_MAX;
int min_trimmed3=INT_MAX;

int qvtrim_qmin=0;
int qvtrim_max=0;  //for -q, do not trim at 3'-end more than this number of bases
int qv_phredtype=0; // could be 64 or 33 (0 means undetermined yet)
int qv_cvtadd=0; //could be -31 or +31

int a3len=0;
int a5len=0;
// adaptor matching metrics -- for extendMatch() function
const int a_m_score=2; //match score
const int a_mis_score=-3; //mismatch
const int a_dropoff_score=7;
int a_min_score=10; //an exact match of 5 bases at the proper ends WILL be trimmed
const int a_min_chain_score=15; //for gapped alignments

class CSegChain;

class CSegPair {
  public:
   GSeg a;
   GSeg b; //the adapter segment
   int score;
   int flags;
   CSegChain* chain;
   CSegPair(int astart=0, int aend=0, int bstart=0, int bend=0, int mscore=0):a(astart,aend),b(bstart, bend) {
      score=mscore;
      if (score==0) score=a.len()*a_m_score;
      flags=0;
      chain=NULL;
      }
   int len() { return  a.len(); }
   bool operator==(CSegPair& d){
      //return (a.start==d.a.start && a.end==d.a.end && b.start==d.b.start && b.end==d.b.end);
      //make equal even segments that are included into one another:
      return (d.a.start>=a.start && d.a.end<=a.end && d.b.start>=b.start && d.b.end<=b.end);
      }
   bool operator>(CSegPair& d){ //ordering based on b (adaptor) start coord and score
     if (b.start==d.b.start) {
        if (score==d.score) {
           //just try to be consistent:
           if (b.end==d.b.end) {
             return (a.start==d.a.start)?(a.end<d.a.end):(a.start<d.a.start);
             }
           return (b.end>d.b.end);
           }
         else return (score<d.score);
        }
     return (b.start>d.b.start);
     }
   bool operator<(CSegPair& d){ //ordering based on b (adaptor) coord
     /*if (b.start==d.b.start && b.end==d.b.end) {
          return (a.start==d.a.start)?(a.end<d.a.end):(a.start<d.a.start);
          }
     return (b.start==d.b.start)?(b.end<d.b.end):(b.start<d.b.start);*/
     if (b.start==d.b.start) {
        if (score==d.score) {
           //just try to be consistent:
           if (b.end==d.b.end) {
             return (a.start==d.a.start)?(a.end>d.a.end):(a.start>d.a.start);
             }
           return (b.end<d.b.end);
           }
         else return (score>d.score);
        }
     return (b.start<d.b.start);
     }
 };

int cmpSegEnds(pointer sa, pointer sb) { //sort by adaptor seg ends AND score
 CSegPair& x = *(CSegPair *)sa;
 CSegPair& y = *(CSegPair *)sb;
 /*
 if (x.b.end==y.b.end) {
     if (x.b.start==y.b.start) {
         if (x.a.end==y.a.end) {
            if (x.a.start==y.a.start) return 0;
            return ((x.a.start>y.a.start) ? -1 : 1);
            }
          else {
            return ((x.a.end>y.a.end) ? -1 : 1);
            }
          }
      else {
       return ((x.b.start>y.b.start) ? -1 : 1);
       }
     }
    else {
     return ((x.b.end>y.b.end) ? -1 : 1);
     }
 */
  if (x.b.end==y.b.end) {
     if (x.score==y.score) {
     if (x.b.start==y.b.start) {
         if (x.a.end==y.a.end) {
            if (x.a.start==y.a.start) return 0;
            return ((x.a.start<y.a.start) ? -1 : 1);
            }
          else {
            return ((x.a.end<y.a.end) ? -1 : 1);
            }
          }
      else {
       return ((x.b.start<y.b.start) ? -1 : 1);
       }
      } else return ((x.score>y.score) ? -1 : 1);
     }
    else {
     return ((x.b.end>y.b.end) ? -1 : 1);
     }

}

class CSegChain:public GList<CSegPair> {
 public:
   uint astart;
   uint aend;
   uint bstart;
   uint bend;
   int score;
   bool endSort;
  CSegChain(bool aln5=false):GList<CSegPair>(true,true,true) {//sorted, free elements, unique
   //as SegPairs are inserted, they will be sorted by a.start coordinate
   score=0;
   astart=MAX_UINT;
   aend=0;
   bstart=MAX_UINT;
   bend=0;
   endSort=aln5;
   if (aln5) { setSorted(cmpSegEnds); }
   } 
 bool operator==(CSegChain& d) {
   //return (score==d.score);
    return (astart==d.astart && aend==d.aend && bstart==d.bstart && bend==d.bend);
   }
 bool operator>(CSegChain& d) { // order based on b (adaptor) coordinate
   //return (score<d.score);
   if (bstart==d.bstart && bend==d.bend) {
          return (astart==d.astart)?(aend>d.aend):(astart>d.astart);
          }
     return (bstart==d.bstart)?(bend>d.bend):(bstart>d.bstart);
   }
 bool operator<(CSegChain& d) {
   //return (score>d.score);
   if (bstart==d.bstart && bend==d.bend) {
          return (astart==d.astart)?(aend<d.aend):(astart<d.astart);
          }
     return (bstart==d.bstart)?(bend<d.bend):(bstart<d.bstart);
   }
 void addSegPair(CSegPair* segp) {
   if (AddIfNew(segp)!=segp) return;
   score+=segp->score;
   if (astart>segp->a.start) astart=segp->a.start;
   if (aend<segp->a.end) aend=segp->a.end;
   if (bstart>segp->b.start) bstart=segp->b.start;
   if (bend<segp->b.end) bend=segp->b.end;
   }
 //for building actual chains:
 bool extendChain(CSegPair* segp) { //segp expected to be "Greater Than" current chain
   int bgap=0;
   int agap=0;
   //if (endSort) {
   if (bstart>segp->b.start) {
      bgap = (int)(bstart-segp->b.end);
      if (abs(bgap)>2) return false;
      agap = (int)(astart-segp->a.end);
      if (abs(agap)>2) return false;
      }
     else {
      bgap = (int) (segp->b.start-bend);
      if (abs(bgap)>2) return false;
      agap = (int)(segp->a.start-aend);
      if (abs(agap)>2) return false;
      }
   if (agap*bgap<0) return false;
   addSegPair(segp);
   score-=abs(agap)+abs(bgap);
   return true;
   }
 };

// element in dhash:
class FqDupRec {
 public:
   int count; //how many of these reads are the same
   int len; //length of qv
   char* firstname; //optional, only if we want to keep the original read names
   char* qv;
   FqDupRec(GStr* qstr=NULL, const char* rname=NULL) {
     len=0;
     qv=NULL;
     firstname=NULL;
     count=0;
     if (qstr!=NULL) {
       qv=Gstrdup(qstr->chars());
       len=qstr->length();
       count++;
       }
     if (rname!=NULL) firstname=Gstrdup(rname);
     }
   ~FqDupRec() {
     GFREE(qv);
     GFREE(firstname);
     }
   void add(GStr& d) { //collapse another record into this one
     if (d.length()!=len)
       GError("Error at FqDupRec::add(): cannot collapse reads with different length!\n");
     count++;
     for (int i=0;i<len;i++)
       qv[i]=(qv[i]+d[i])/2;
     }
 };

void openfw(FILE* &f, GArgs& args, char opt) {
  GStr s=args.getOpt(opt);
  if (!s.is_empty()) {
      if (s=='-') f=stdout;
      else {
       f=fopen(s.chars(),"w");
       if (f==NULL) GError("Error creating file: %s\n", s.chars());
       }
     }
}

#define FWCLOSE(fh) if (fh!=NULL && fh!=stdout) fclose(fh)
#define FRCLOSE(fh) if (fh!=NULL && fh!=stdin) fclose(fh)

GHash<FqDupRec> dhash; //hash to keep track of duplicates

void setupFiles(FILE*& f_in, FILE*& f_in2, FILE*& f_out, FILE*& f_out2, 
                       GStr& s, GStr& infname, GStr& infname2); 
// uses outsuffix to generate output file names and open file handles as needed
 
void writeRead(FILE* f_out, GStr& rname, GStr& rinfo, GStr& rseq, GStr& rqv, int& outcounter);
void trash_report(char trashcode, GStr& rname, FILE* freport);

bool getFastxRec(GLineReader& fq, GStr& rseq, GStr& rqv, 
          GStr& rname, GStr& rinfo, GStr& infname);

char process_read(GStr& rname, GStr& rseq, GStr& rqv, int &l5, int &l3); 
//returns 0 if the read was untouched, 1 if it was trimmed and a trash code if it was trashed

bool ntrim(GStr& rseq, int &l5, int &l3); //returns true if any trimming occured
bool qtrim(GStr& qvs, int &l5, int &l3); //return true if any trimming occured
int dust(GStr& seq);
bool trim_adapter3(GStr& seq, int &l5, int &l3); //returns true if any trimming occured
bool trim_adapter5(GStr& seq, int &l5, int &l3); //returns true if any trimming occured

void convertPhred(char* q, int len);
void convertPhred(GStr& q);

int main(int argc, char * const argv[]) {
  GArgs args(argc, argv, "YQDCVl:d:3:5:m:n:r:p:q:t:o:z:a:");
  int e;
  if ((e=args.isError())>0) {
      GMessage("%s\nInvalid argument: %s\n", USAGE, argv[e]);
      exit(224);
      }
  debug=(args.getOpt('Y')!=NULL);
  verbose=(args.getOpt('V')!=NULL);
  convert_phred=(args.getOpt('Q')!=NULL);
  doCollapse=(args.getOpt('C')!=NULL);
  doDust=(args.getOpt('D')!=NULL);
  /*
  rawFormat=(args.getOpt('R')!=NULL);
  if (rawFormat) {
    GError("Sorry, raw qseq format parsing is not implemented yet!\n");
    }
  */
  prefix=args.getOpt('n');
  GStr s=args.getOpt('l');
  if (!s.is_empty()) 
     min_read_len=s.asInt();
  s=args.getOpt('m');
  if (!s.is_empty()) 
     max_perc_N=s.asDouble();
  s=args.getOpt('d');
  if (!s.is_empty()) {
     dust_cutoff=s.asInt();
     doDust=true;
     }
  s=args.getOpt('q');
  if (!s.is_empty()) {
     qvtrim_qmin=s.asInt();
     }
  s=args.getOpt('t');
  if (!s.is_empty()) {
     qvtrim_max=s.asInt();
     }
  s=args.getOpt('p');
  if (!s.is_empty()) {
     int v=s.asInt();
     if (v==33) {
        qv_phredtype=33;
        qv_cvtadd=31;
        }
      else if (v==64) {
        qv_phredtype=64;
        qv_cvtadd=-31;
        }
       else 
         GMessage("%s\nInvalid value for -p option (can only be 64 or 33)!\n",USAGE);
     }
  if (args.getOpt('3')!=NULL) {
    adapter3=args.getOpt('3');
    adapter3.upper();
    a3len=adapter3.length();
    }
  if (args.getOpt('5')!=NULL) {
    adapter5=args.getOpt('5');
    adapter5.upper();
    a5len=adapter5.length();
    }
  s=args.getOpt('a');
  if (!s.is_empty()) {
     int a_minmatch=s.asInt();
     a_min_score=a_minmatch<<1;
     }

  if (args.getOpt('o')!=NULL) outsuffix=args.getOpt('o');
  trashReport=  (args.getOpt('r')!=NULL);
  int fcount=args.startNonOpt();
  if (fcount==0) {
    GMessage(USAGE);
    exit(224);
    }
   if (fcount>1 && doCollapse) {
    GError("%s Sorry, the -C option only works with a single input.\n", USAGE);
    }
  //openfw(f_out, args, 'o');
  //if (f_out==NULL) f_out=stdout;
  if (trashReport)
    openfw(freport, args, 'r');
  char* infile=NULL;
  while ((infile=args.nextNonOpt())!=NULL) {
    int incounter=0; //counter for input reads
    int outcounter=0; //counter for output reads
    int trash_s=0; //too short from the get go
    int trash_Q=0;
    int trash_N=0;
    int trash_D=0;
    int trash_A3=0;
    int trash_A5=0;
    s=infile;
    GStr infname;
    GStr infname2;
    FILE* f_in=NULL;
    FILE* f_in2=NULL;
    FILE* f_out=NULL;
    FILE* f_out2=NULL;
    bool paired_reads=false;
    setupFiles(f_in, f_in2, f_out, f_out2, s, infname, infname2);
    GLineReader fq(f_in);
    GLineReader* fq2=NULL;
    if (f_in2!=NULL) {
       fq2=new GLineReader(f_in2);
       paired_reads=true;
       }
    GStr rseq, rqv, seqid, seqinfo;
    GStr rseq2, rqv2, seqid2, seqinfo2;
    while (getFastxRec(fq, rseq, rqv, seqid, seqinfo, infname)) {
       incounter++;
       int a5=0, a3=0, b5=0, b3=0;
       char tcode=0, tcode2=0;
       tcode=process_read(seqid, rseq, rqv, a5, a3);
       //if (!doCollapse) {
         if (fq2!=NULL) {
            getFastxRec(*fq2, rseq2, rqv2, seqid2, seqinfo2, infname2);
            if (seqid.substr(0,seqid.length()-1)!=seqid2.substr(0,seqid2.length()-1)) {
               GError("Error: no paired match for read %s vs %s (%s,%s)\n",
                  seqid.chars(), seqid2.chars(), infname.chars(), infname2.chars());
               }
            tcode2=process_read(seqid2, rseq2, rqv2, b5, b3);
            //decide what to do with this pair and print rseq2 if the pair makes it
            if (tcode>1 && tcode2<=1) {
               //"untrash" rseq
               if (a3-a5+1<min_read_len) {
                   a5=1;
                   if (a3<min_read_len) { a3= GMIN(rseq.length()-1, min_read_len+1); }
                   }
               tcode=1;
               }
             else if (tcode<=1 && tcode2>1) {
               //"untrash" rseq2
               if (b3-b5+1<min_read_len) {
                   b5=1;
                   if (b3<min_read_len) { b3= GMIN((rseq2.length()-1),(min_read_len+1)); }
                   }
               tcode2=1;
               }
            if (tcode<=1) { //trimmed or left intact -- write it!
               if (tcode>0) {
                 rseq2=rseq2.substr(b5,b3-b5+1);
                 if (!rqv2.is_empty()) rqv2=rqv2.substr(b5,b3-b5+1);
                 }
               int nocounter=0;
               writeRead(f_out2, seqid2, seqinfo2, rseq2, rqv2, nocounter);
               }
            } //paired read
       // }
       if (tcode>1) { //trashed
          if (tcode=='s') trash_s++;
            else if (tcode=='Q') trash_Q++;
              else if (tcode=='N') trash_N++;
               else if (tcode=='D') trash_D++;
                else if (tcode=='3') trash_A3++;
                 else if (tcode=='5') trash_A5++;
          if (trashReport) trash_report(tcode, seqid, freport);
          }
         else if (!doCollapse) { //write it
          if (tcode>0) {
            rseq=rseq.substr(a5,a3-a5+1);
            if (!rqv.is_empty()) rqv=rqv.substr(a5,a3-a5+1);
            }
          writeRead(f_out, seqid, seqinfo, rseq, rqv, outcounter);
          }
       } //for each fastq record
    delete fq2;
    FRCLOSE(f_in);
    FRCLOSE(f_in2);
    if (doCollapse) {
       outcounter=0;
       int maxdup_count=1;
       char* maxdup_seq=NULL;
       dhash.startIterate();
       FqDupRec* qd=NULL;
       char* seq=NULL;
       while ((qd=dhash.NextData(seq))!=NULL) {
         GStr rseq(seq);
         //do the dusting here
         if (doDust) {
            int dustbases=dust(rseq);
            if (dustbases>(rseq.length()>>1)) {
               if (trashReport && qd->firstname!=NULL) {
                 fprintf(freport, "%s_x%d\tD\n",qd->firstname, qd->count);
                 }
               trash_D+=qd->count;
               continue;
               }
            }
         outcounter++;
         if (qd->count>maxdup_count) { 
            maxdup_count=qd->count;
            maxdup_seq=seq;
            }
         if (isfasta) {
           if (prefix.is_empty()) {
             fprintf(f_out, ">%s_x%d\n%s\n", qd->firstname, qd->count, 
                           rseq.chars());
             }
           else { //use custom read name
             fprintf(f_out, ">%s%08d_x%d\n%s\n", prefix.chars(), outcounter,
                        qd->count, rseq.chars());
             }
           }
         else { //fastq format
          if (convert_phred) convertPhred(qd->qv, qd->len);
          if (prefix.is_empty()) {
            fprintf(f_out, "@%s_x%d\n%s\n+\n%s\n", qd->firstname, qd->count, 
                           rseq.chars(), qd->qv);
            }
          else { //use custom read name
            fprintf(f_out, "@%s%08d_x%d\n%s\n+\n%s\n", prefix.chars(), outcounter,
                        qd->count, rseq.chars(), qd->qv);
            }
           }
         }//for each element of dhash
       if (maxdup_count>1) {
         GMessage("Maximum read multiplicity: x %d (read: %s)\n",maxdup_count, maxdup_seq);
         }
       } //collapse entries
    if (verbose) {
       if (paired_reads) {
           GMessage(">Input files : %s , %s\n", infname.chars(), infname2.chars());
           GMessage("Number of input pairs :%9d\n", incounter);
           GMessage("         Output pairs :%9d\n", outcounter);
           }
         else {
           GMessage(">Input file : %s\n", infname.chars());
           GMessage("Number of input reads :%9d\n", incounter);
           GMessage("         Output reads :%9d\n", outcounter);
           }
       GMessage("------------------------------------\n");
       if (num_trimmed5)
          GMessage("           5' trimmed :%9d  (min. trim: %d)\n", num_trimmed5, min_trimmed5);
       if (num_trimmed3)
          GMessage("           3' trimmed :%9d  (min. trim: %d)\n", num_trimmed3, min_trimmed3);
       GMessage("------------------------------------\n");
       if (trash_s>0)
         GMessage("     Trashed by length:%9d\n", trash_s);
       if (trash_N>0)
         GMessage("         Trashed by N%%:%9d\n", trash_N);
       if (trash_Q>0)
         GMessage("Trashed by low quality:%9d\n", trash_Q);
       if (trash_A5>0)
         GMessage(" Trashed by 5' adapter:%9d\n", trash_A5);
       if (trash_A3>0)
         GMessage(" Trashed by 3' adapter:%9d\n", trash_A3);
       }
    if (trashReport) {
          FWCLOSE(freport);
          }
    FWCLOSE(f_out);
    FWCLOSE(f_out2);
   } //while each input file
 
 //getc(stdin);
}

class NData {
 public:
   int NPos[1024]; //there should be no reads longer than 1K ?
   int NCount;
   int end5;
   int end3;
   int seqlen;
   double perc_N; //percentage of Ns in end5..end3 range only!
   const char* seq;
   bool valid;
   NData() {
    NCount=0;
    end5=0;
    end3=0;
    seq=NULL;
    perc_N=0;
    valid=true;
    }
   void init(GStr& rseq) {
     NCount=0;
     perc_N=0;
     seqlen=rseq.length();
     seq=rseq.chars();
     end5=1;
     end3=seqlen;
     for (int i=0;i<seqlen;i++)
        if (seq[i]=='N') {// if (!ichrInStr(rseq[i], "ACGT")
           NPos[NCount]=i;
           NCount++;
           }
     perc_N=(NCount*100.0)/seqlen;
     valid=true;
     }
  void N_calc() { //only in the region after trimming
     int n=0;
     for (int i=end3-1;i<end5;i++) {
       if (seq[i]=='N') n++;
       }
     perc_N=(n*100.0)/(end5-end3+1);
     }
 };
 
static NData feat;
int perc_lenN=12; // incremental distance from ends, in percentage of
          // sequence length, where N-trimming is done (default:12 %) (autolimited to 20)
          
void N_analyze(int l5, int l3, int p5, int p3) {
/* assumes feat was filled properly */
 int old_dif, t5,t3,v;
 if (l3<l5+2 || p5>p3 ) {
   feat.end5=l5+1;
   feat.end3=l3+1;
   return; 
   }

 t5=feat.NPos[p5]-l5;
 t3=l3-feat.NPos[p3];
 old_dif=p3-p5;
 v=(int)((((double)(l3-l5))*perc_lenN)/100);
 if (v>20) v=20; /* enforce N-search limit for very long reads */ 
 if (t5 < v ) {
   l5=feat.NPos[p5]+1;
   p5++;
   }
 if (t3 < v) {
   l3=feat.NPos[p3]-1;
   p3--;
   }
 /* restNs=p3-p5; number of Ns in the new CLR */
 if (p3-p5==old_dif) { /* no change, return */
           /* don't trim if this may shorten the read too much */
           //if (l5-l3<min_read_len) return;
           feat.end5=l5+1;
           feat.end3=l3+1;
           return;
           }
    else 
      N_analyze(l5,l3, p5,p3);
}


bool qtrim(GStr& qvs, int &l5, int &l3) {
if (qvtrim_qmin==0 || qvs.is_empty()) return false;
if (qv_phredtype==0) {
  //try to guess the Phred type
  int vmin=256, vmax=0;
  for (int i=0;i<qvs.length();i++) {
     if (vmin>qvs[i]) vmin=qvs[i];
     if (vmax<qvs[i]) vmax=qvs[i];
     }
  if (vmin<64) { qv_phredtype=33; qv_cvtadd=31; }
  if (vmax>95) { qv_phredtype=64; qv_cvtadd=-31; }
  if (qv_phredtype==0) {
    GError("Error: couldn't determine Phred type, please use the -p33 or -p64 !\n");
    }
  if (verbose)
    GMessage("Input reads have Phred-%d quality values.\n", (qv_phredtype==33 ? 33 : 64));
  } //guessing Phred type
for (l3=qvs.length()-1;l3>2;l3--) {
  if (qvs[l3]-qv_phredtype>=qvtrim_qmin && qvs[l3-1]-qv_phredtype>=qvtrim_qmin) break;
  }
//just in case, check also the 5' the end (?)
for (l5=0;l5<qvs.length()-3;l5++) {
  if (qvs[l5]-qv_phredtype>=qvtrim_qmin && qvs[l5+1]-qv_phredtype>=qvtrim_qmin) break;
  }
if (qvtrim_max>0) {
  if (qvs.length()-1-l3>qvtrim_max) l3=qvs.length()-1-qvtrim_max;
  if (l5>qvtrim_max) l5=qvtrim_max;
  }
return (l5>0 || l3<qvs.length()-1);
}

bool ntrim(GStr& rseq, int &l5, int &l3) {
 //count Ns in the sequence, trim N-rich ends
 feat.init(rseq);
 l5=feat.end5-1;
 l3=feat.end3-1;
 N_analyze(feat.end5-1, feat.end3-1, 0, feat.NCount-1); 
 if (l5==feat.end5-1 && l3==feat.end3-1) {
    if (feat.perc_N>max_perc_N) {
           feat.valid=false;
           l3=l5+1;
           return true;
           }
      else {
       return false; //nothing changed
       }
    }
 l5=feat.end5-1;
 l3=feat.end3-1;
 if (l3-l5+1<min_read_len) {
   feat.valid=false;
   return true;
   }
 feat.N_calc();
 
 if (feat.perc_N>max_perc_N) {
      feat.valid=false;
      l3=l5+1;
      return true;
      }
 return true;
 }

//--------------- dust functions ----------------
class DNADuster {
 public:
  int dustword; 
  int dustwindow; 
  int dustwindow2; 
  int dustcutoff;
  int mv, iv, jv;
  int counts[32*32*32];
  int iis[32*32*32];
  DNADuster(int cutoff=16, int winsize=32, int wordsize=3) {
    dustword=wordsize;
    dustwindow=winsize;
    dustwindow2 = (winsize>>1);
    dustcutoff=cutoff;
    mv=0;
    iv=0;
    jv=0;
    }
  void setWindowSize(int value) {
    dustwindow = value;
    dustwindow2 = (dustwindow >> 1);
    }
  void setWordSize(int value) {
    dustword=value;
    }
void wo1(int len, const char* s, int ivv) {
  int i, ii, j, v, t, n, n1, sum;
  int js, nis;
  n = 32 * 32 * 32;
  n1 = n - 1;
  nis = 0;
  i = 0;
  ii = 0;
  sum = 0;
  v = 0;
  for (j=0; j < len; j++, s++) {
        ii <<= 5;
        if (*s<=32) {
           i=0;
           continue;
           }
        ii |= *s - 'A'; //assume uppercase!
        ii &= n1;
        i++;
        if (i >= dustword) {
              for (js=0; js < nis && iis[js] != ii; js++) ;
              if (js == nis) {
                    iis[nis] = ii;
                    counts[ii] = 0;
                    nis++;
              }
              if ((t = counts[ii]) > 0) {
                    sum += t;
                    v = 10 * sum / j;
                    if (mv < v) {
                          mv = v;
                          iv = ivv;
                          jv = j;
                    }
              }
              counts[ii]++;
        }
  }
}

int wo(int len, const char* s, int* beg, int* end) {
      int i, l1;
      l1 = len - dustword + 1;
      if (l1 < 0) {
            *beg = 0;
            *end = len - 1;
            return 0;
            }
      mv = 0;
      iv = 0;
      jv = 0;
      for (i=0; i < l1; i++) {
            wo1(len-i, s+i, i);
            }
      *beg = iv;
      *end = iv + jv;
      return mv;
 }

void dust(const char* seq, char* seqmsk, int seqlen, int cutoff=0) { //, maskFunc maskfn) {
  int i, j, l, a, b, v;
  if (cutoff==0) cutoff=dustcutoff;
  a=0;b=0;
  //GMessage("Dust cutoff=%d\n", cutoff);
  for (i=0; i < seqlen; i += dustwindow2) {
        l = (seqlen > i+dustwindow) ? dustwindow : seqlen-i;
        v = wo(l, seq+i, &a, &b);
        if (v > cutoff) {
           //for (j = a; j <= b && j < dustwindow2; j++) {
           for (j = a; j <= b; j++) {
                    seqmsk[i+j]='N';//could be made lowercase instead
                    }
           }
         }
//return first; 
 }
};

static DNADuster duster;

int dust(GStr& rseq) {
 char* seq=Gstrdup(rseq.chars());
 duster.dust(rseq.chars(), seq, rseq.length(), dust_cutoff);
 //check the number of Ns:
 int ncount=0;
 for (int i=0;i<rseq.length();i++) {
   if (seq[i]=='N') ncount++;
   }
 GFREE(seq);
 return ncount;
 }


// ------------------ adapter matching - simple k-mer seed & extend, no indels for now
//when a k-mer match is found, simply try to extend the alignment using a drop-off scheme
//check minimum score and 
//for 3' adapter trimming: 
//     require that the right end of the alignment for either the adaptor OR the read must be 
//     < 3 distance from its right end
// for 5' adapter trimming: 
//     require that the left end of the alignment for either the adaptor OR the read must 
//     be at coordinate < 3 from start

bool extendMatch(const char* a, int alen, int ai, 
                 const char* b, int blen, int bi, int mlen, int& l5, int& l3, CSegChain& segs, bool end5=false) {
 //so the alignment starts at ai in a, bi in b, with a perfect match of length mlen
#ifdef DEBUG
 GStr dbg(b);
#endif 
 //if (debug) {
 //  GMessage(">> in %s\n\textending hit: %s at position %d\n", a, (dbg.substr(bi, mlen)).chars(), ai);
 //  }
 int a_l=ai; //alignment coordinates on a
 int a_r=ai+mlen-1;
 int b_l=bi; //alignment coordinates on b
 int b_r=bi+mlen-1; 
 int ai_maxscore=ai;
 int bi_maxscore=bi;
 int score=mlen*a_m_score;
 int maxscore=score;
 int mism5score=a_mis_score; 
 if (end5 && ai<(alen>>1)) mism5score-=2; // increase penalty for mismatches at 5' end
 //try to extend to the left first, if possible
 while (ai>0 && bi>0) {
   ai--;
   bi--;
   score+= (a[ai]==b[bi])? a_m_score : mism5score;
   if (score>maxscore) {
       ai_maxscore=ai;
       bi_maxscore=bi;
       maxscore=score;
       }
     else if (maxscore-score>a_dropoff_score) break;
   }
 a_l=ai_maxscore;
 b_l=bi_maxscore;
 //if (debug) GMessage("  after l-extend: %*s%s\t\t(score=%d)\n",a_l," ",dbg.substr(b_l,b_r-b_l+1).chars(),maxscore);
 //now extend to the right
 ai_maxscore=a_r;
 bi_maxscore=b_r;
 ai=a_r;
 bi=b_r;
 score=maxscore;
 //sometimes there are extra AAAAs at the end of the read, ignore those
 if (strcmp(&a[alen-4],"AAAA")==0) {
    alen-=3;
    while (a[alen-1]=='A' && alen>ai) alen--;
    }
 while (ai<alen-1 && bi<blen-1) {
   ai++;
   bi++;
   //score+= (a[ai]==b[bi])? a_m_score : a_mis_score;
   if (a[ai]==b[bi]) { //match
      score+=a_m_score;
      if (ai>=alen-2) {
           score+=a_m_score-(alen-ai-1);
           }
      }
    else { //mismatch
      score+=a_mis_score;
      }  
   if (score>maxscore) {
       ai_maxscore=ai;
       bi_maxscore=bi;
       maxscore=score;
       }
     else if (maxscore-score>a_dropoff_score) break;
   }
  a_r=ai_maxscore;
  b_r=bi_maxscore;
  int a_ovh3=alen-a_r-1;
  int b_ovh3=blen-b_r-1;
  int mmovh3=(a_ovh3<b_ovh3)? a_ovh3 : b_ovh3;
  int mmovh5=(a_l<b_l)? a_l : b_l;
  //if (debug) GMessage("  after r-extend: %*s%s\t\t(score=%d)\n",a_l," ",dbg.substr(b_l,b_r-b_l+1).chars(),maxscore);
#ifdef DEBUG
  if (debug) GMessage("     extended to: %*s\n",a_r+1,dbg.substr(b_l,b_r-b_l+1).chars());
#endif
  if (maxscore>=a_min_score && mmovh3<2 && mmovh5<2) {
     if (a_l<a_ovh3) {
        //adapter closer to the left end (typical for 5' adapter)
        l5=a_r+1;
        l3=alen-1;
        }
      else {
        //adapter matching at the right end (typical for 3' adapter)
        l5=0;
        l3=a_l-1;
        }
     return true;
     }
 else { //keep this segment pair for later (gapped alignment)
   segs.addSegPair(new CSegPair(a_l+1, a_r+1, b_l+1, b_r+1, maxscore));
   //this will also update min & max coordinates in segs (segs.astart, .aend, .bstart, .bend)
   }
  //do not trim:
  l5=0;
  l3=alen-1;
  return false;
 }

/*
int getWordValue(const char* s, int wlen) {
 int r=0;
 while (wlen--) { r+=(((int)s[wlen])<<wlen) }
 return r;
 }
*/
int get3mer_value(const char* s) {
 return (s[0]<<16)+(s[1]<<8)+s[2];
 }
 
int w3_match(int qv, const char* str, int slen, int start_index=0) {
 if (start_index>=slen || start_index<0) return -1;
 for (int i=start_index;i<slen-3;i++) {
   int rv=get3mer_value(str+i);
   if (rv==qv) return i;
   }
 return -1;
 }

int w3_rmatch(int qv, const char* str, int slen, int end_index=-1) {
 if (end_index>=slen) return -1;
 if (end_index<0) end_index=slen-1;
 for (int i=end_index-2;i>=0;i--) {
   int rv=get3mer_value(str+i);
   if (rv==qv) return i;
   }
 return -1;
 }

int fast4match(int32 qv, const char* str, int slen, int start_index=0) {
 if (start_index>=slen || start_index<0) return -1;
 for (int i=start_index;i<slen-4;i++) {
   int32* rv=(int32*)(str+i);
   if (*rv==qv) return i;
   }
 return -1;
 }

int fast4rmatch(int32 qv, const char* str, int slen, int end_index=-1) {
 if (end_index>=slen) return -1;
 if (end_index<0) end_index=slen-1;
 for (int i=end_index-3;i>=0;i--) {
   int32* rv=(int32*)(str+i);
   if (*rv==qv) return i;
   }
 return -1;
 }

#ifdef DEBUG
void dbgPrintChain(CSegChain& chain, const char* aseq) {
  GStr s(aseq);
  for (int i=0;i<chain.Count();i++) {
   CSegPair& seg=*chain[i];
   GMessage("  dbg chain seg%d: %*s [%d-%d:%d-%d]\n",i,seg.a.start-1+seg.len(),
            s.substr(seg.b.start-1, seg.len()).chars(), seg.b.start,seg.b.end,seg.a.start,seg.a.end);
   }
}
#endif

bool trim_adapter3(GStr& seq, int&l5, int &l3) {
 int rlen=seq.length();
 l5=0;
 l3=rlen-1;
 //first try a full match, we might get lucky
 int fi=-1;
 if ((fi=seq.index(adapter3))>=0) {
   if (fi<rlen-fi-a3len) {//match is closer to the right end
      l5=fi+a3len;
      l3=rlen-1;
      }
    else {
      l5=0;
      l3=fi-1;
      }
   return true;
   }
 #ifdef DEBUG
 if (debug) GMessage(">TRIM3 >>   Read: %s\n",seq.chars());
 #endif

 //also, for fast detection of other adapter-only reads that start past 
 // the beginning of the adapter sequence, try to see if the first a3len-4 
 // bases of the read are a substring of the adapter
 if (rlen>a3len-3) {
   GStr rstart=seq.substr(1,a3len-4);
   if ((fi=adapter3.index(rstart))>=0) {
     l3=rlen-1;
     l5=a3len-4;
     while (fi+l5<a3len && l5<l3 && adapter3[fi+l5]==seq[l5]) l5++;
     return true;
     }
  }
 CSegChain a3segs; //no chains here, just an ordered collection of segment pairs
  //check the easy cases - 11 bases exact match at the end
 int fdlen=11;
  if (a3len<16) {
   fdlen=a3len>>1;
   }
 if (fdlen>4) {
     //check if we're lucky enough to have the last 11 bases of the read a part of the adapter
     GStr rstart=seq.substr(-fdlen-3,fdlen);
     if ((fi=adapter3.index(rstart))>=0) {
#ifdef DEBUG
       if (debug) GMessage("  W11match found: %*s\n", rlen-3, (adapter3.substr(fi,fdlen)).chars());
#endif
       if (extendMatch(seq.chars(), rlen, rlen-fdlen-3,
                     adapter3.chars(), a3len, fi,  fdlen, l5,l3, a3segs))
            return true;
       }
     //another easy case: first 11 characters of the adaptor found as a substring of the read
     GStr bstr=adapter3.substr(0, fdlen);
     if ((fi=seq.rindex(bstr))>=0) {
#ifdef DEBUG
       if (debug) GMessage("  A11match found: %*s\n", fi+fdlen, bstr.chars());
#endif
       if (extendMatch(seq.chars(), rlen, fi,
                     adapter3.chars(), a3len, 0,  fdlen, l5,l3, a3segs))
            return true;
       }
     } //tried to match 11 bases first 
     
 //no easy cases, so let's do the wmer hashing for the first 12 bases of the adaptor
 //-- only extend if the match is in the 3' (ending) region of the read
 int wordSize=3;
 int hlen=12;
 if (hlen>a3len-wordSize) hlen=a3len-wordSize;
 int imin=rlen>>1; //last half of the read, left boundary for the wmer match
 if (imin<a3len) { imin=GMIN(a3len, rlen-wordSize); }
 imin=rlen-imin;
 for (int iw=0;iw<hlen;iw++) {
   //int32* qv=(int32*)(adapter3.chars()+iw);
   int qv=get3mer_value(adapter3.chars()+iw);
   fi=-1;
   //while ((fi=fast4rmatch(*qv, seq.chars(), rlen, fi))>=0 && fi>=imin) {
   while ((fi=w3_rmatch(qv, seq.chars(), rlen, fi))>=0 && fi>=imin) {
     //GMessage(" ... fi=%d after w3_rmatch() (imin=%d)\n", fi, imin);

#ifdef DEBUG
     if (debug) GMessage("    Wmatch found: %*s\n", fi+wordSize, (adapter3.substr(iw,wordSize)).chars());
#endif
     if (extendMatch(seq.chars(), rlen, fi, adapter3.chars(), 
                   a3len, iw, wordSize, l5,l3, a3segs)) return true;
     fi--;
     if (fi<imin) break;
     }
   } //for each wmer in the first hlen bases of the adaptor
 /*
 //couldn't find a good trimming extension, hash 12 more bases of the adapter to collect more segment pairs there
 //but only do this if we already have segment pairs collected in the last 12 bases of the adapter
 if (a3segs.bstart>3 || a3segs.bend<(uint)(hlen-wordSize)) return false; 
 int hlen2=a3len-wordSize;
 //if (hlen2>a3len-4) hlen2=a3len-4;
 if (hlen2>hlen) {
#ifdef DEBUG
     if (debug && a3segs.Count()>0) {
        GMessage("  >>>>>2nd. hash: %s\n",seq.chars());
        }
#endif
     for (int iw=hlen;iw<hlen2;iw++) {
         //int* qv=(int32 *)(adapter3.chars()+iw);
         int qv=get3mer_value(adapter3.chars()+iw);
         fi=-1;
         //while ((fi=fast4rmatch(*qv, seq.chars(), rlen, fi))>=0 && fi>=imin) {
         while ((fi=w3_rmatch(qv, seq.chars(), rlen, fi))>=0 && fi>=imin) {
           extendMatch(seq.chars(), rlen, fi, adapter3.chars(), 
                         a3len, iw, wordSize, l5,l3, a3segs);
           fi--;
           if (fi<imin) break;
           }
         } //for each wmer between hlen2 and hlen bases of the adaptor
     }
 //lastly, analyze collected a3segs for a possible gapped alignment:
 GList<CSegChain> segchains(false,true,false);
#ifdef DEBUG
 if (debug && a3segs.Count()>0) {
   GMessage(">>>>>>>>>   Read: %s\n",seq.chars());
   }
#endif
 for (int i=0;i<a3segs.Count();i++) {
   if (a3segs[i]->chain==NULL) {
       if (a3segs[i]->b.start>3) continue; //don't start a hopeless chain
       CSegChain* newchain=new CSegChain();
       newchain->setFreeItem(false);
       newchain->addSegPair(a3segs[i]);
       a3segs[i]->chain=newchain;
       segchains.Add(newchain); //just to free them when done
       }
   for (int j=i+1;j<a3segs.Count();j++) {
      CSegChain* chain=a3segs[i]->chain;
      if (chain->extendChain(a3segs[j])) {
          a3segs[j]->chain=chain;
#ifdef DEBUG
          if (debug) dbgPrintChain(*chain, adapter3.chars());
#endif      
          //save time by checking here if the extended chain is already acceptable for trimming
          if (chain->aend>(uint)(rlen-4) && chain->bstart<4 && chain->score>a_min_chain_score) {
            l5=0;
            l3=chain->astart-2;
#ifdef DEBUG
          if (debug && a3segs.Count()>0) {
            GMessage(">>> >> trimmed-3: %*s\n",l3-l5+1,seq.substr(l5,l3-l5+1).chars());
            }
#endif
            return true;
            }
          } //chain can be extended
      }
   } //collect segment alignments into chains
 */  
 return false; //no adapter parts found 
 }

bool trim_adapter5(GStr& seq, int&l5, int &l3) {
 //if (debug) GMessage("trim_adapter5 on: %s\n", seq.chars());
 int rlen=seq.length();
 l5=0;
 l3=rlen-1;
 //try to see if adapter is fully included in the read
 int fi=-1;
 if ((fi=seq.index(adapter5))>=0) {
   if (fi<rlen-fi-a5len) {//match is closer to the right end
      l5=fi+a5len;
      l3=rlen-1;
      }
    else {
      l5=0;
      l3=fi-1; 
      }
   return true;
   }
#ifdef DEBUG
 if (debug) GMessage(">TRIM5 >>   Read: %s\n",seq.chars());
#endif

 CSegChain a5segs(true); //list of segment pairs to analyze later if no extendMatch succeeded

 //try the easy way out first - look for an exact match of 11 bases
 int fdlen=11;
  if (a5len<16) {
   fdlen=a5len>>1;
   }
 if (fdlen>4) {
     GStr rstart=seq.substr(1,fdlen); //skip the first base as it's sometimes bogus
     if ((fi=adapter5.index(rstart))>=0) {
#ifdef DEBUG
       if (debug) GMessage("  W11match found: %*s\n", 1+fdlen, (adapter3.substr(fi,fdlen)).chars());
#endif
       if (extendMatch(seq.chars(), rlen, 1,
                     adapter5.chars(), a5len, fi,  fdlen, l5,l3, a5segs, true))
           return true;
       }
     //another easy case: last 11 characters of the adaptor found as a substring of the read
     GStr bstr=adapter5.substr(-fdlen);
     if ((fi=seq.index(bstr))>=0) {
#ifdef DEBUG
       if (debug) GMessage("  A11match found: %*s\n", fi+fdlen, bstr.chars());
#endif
       if (extendMatch(seq.chars(), rlen, fi,
                     adapter5.chars(), a5len, a5len-fdlen,  fdlen, l5,l3,a5segs,true))
          return true;
       }
     } //tried to matching at most 11 bases first 
     
 //-- no easy cases, do the wmer hashing for the last 12 bases of the adaptor
 //-- only extend a wmer if it matches in the 5' (beginning) region of the read
 int wordSize=3; 
 int hlen=12;
 if (hlen>a5len-wordSize) hlen=a5len-wordSize;
 int imax=rlen>>1; //first half of the read, right boundary for the wmer match
 if (imax<a5len) { imax=GMIN(a5len, rlen-wordSize); }
 for (int iw=0;iw<=hlen;iw++) {
   int apstart=a5len-iw-wordSize;
   fi=0;
   //int* qv=(int32 *)(adapter5.chars()+apstart);
   int qv=get3mer_value(adapter5.chars()+apstart);
   //while ((fi=fast4match(*qv, seq.chars(), rlen, fi))>=0 && fi<=imax) {
   while ((fi=w3_match(qv, seq.chars(), rlen, fi))>=0 && fi<=imax) {
#ifdef DEBUG
     if (debug) GMessage("    Wmatch found: %*s\n", fi+wordSize, (adapter5.substr(apstart,wordSize)).chars());
#endif
     if (extendMatch(seq.chars(), rlen, fi, adapter5.chars(), 
                a5len, apstart, wordSize, l5,l3, a5segs, true)) return true;
     fi++;
     if (fi>imax) break;
     }
   } //for each wmer in the last hlen bases of the adaptor
 /* 

 //couldn't find a good trimming extension, hash 12 more bases of the adapter to collect more segment pairs there
 //but only do this if we already have segment pairs collected in the last 12 bases of the adapter
 if (a5segs.bend<(uint)(a5len-3) || a5segs.bstart>(uint)(a5len-hlen+4)) return false;
 int hlen2=a5len-wordSize;
 //if (hlen2>a5len-wordSize) hlen2=a5len-wordSize;
#ifdef DEBUG
      if (debug && a5segs.Count()>0) {
        GMessage("  >>>>>2nd. hash: %s\n",seq.chars());
        }
#endif
 if (hlen2>hlen) {
     for (int iw=hlen+1;iw<=hlen2;iw++) {
         int apstart=a5len-iw-wordSize;
         fi=0;
         //int* qv=(int32 *)(adapter5.chars()+apstart);
         int qv=get3mer_value(adapter5.chars()+apstart);
         //while ((fi=fast4match(*qv, seq.chars(), rlen, fi))>=0 && fi<=imax) {
         while ((fi=w3_match(qv, seq.chars(), rlen, fi))>=0 && fi<=imax) {
           extendMatch(seq.chars(), rlen, fi, adapter5.chars(), 
                      a5len, apstart, wordSize, l5,l3, a5segs, true);
           fi++;
           if (fi>imax) break;
           }
         } //for each wmer between hlen2 and hlen bases of the adaptor
     }
 if (a5segs.bend<(uint)(a5len-3) || a5segs.bstart>(uint)(a5len-hlen+4)) return false;
 // lastly, analyze collected a5segs for a possible gapped alignment:
 GList<CSegChain> segchains(false,true,false);
#ifdef DEBUG
 if (debug && a5segs.Count()>0) {
   GMessage(">>>>>>>>>   Read: %s\n",seq.chars());
   }
#endif
 for (int i=0;i<a5segs.Count();i++) {
   if (a5segs[i]->chain==NULL) {
       if (a5segs[i]->b.end<(int)(a5len-4)) continue; //don't start a hopeless chain
       CSegChain* newchain=new CSegChain(true);
       newchain->setFreeItem(false);
       newchain->addSegPair(a5segs[i]);
       a5segs[i]->chain=newchain;
       segchains.Add(newchain); //just to free them when done
       }
   for (int j=i+1;j<a5segs.Count();j++) {
      CSegChain* chain=a5segs[i]->chain;
      if (chain->extendChain(a5segs[j])) {
         a5segs[j]->chain=chain;
#ifdef DEBUG
         if (debug) dbgPrintChain(*chain, adapter5.chars());
#endif      
      //save time by checking here if the extended chain is already acceptable for trimming
         if (chain->bend>(uint)(a5len-3) && chain->astart<4 && chain->score>a_min_chain_score) {
            l5=chain->aend;
            l3=rlen-1;
            return true;
            }
         } //chain can be extended
      }
   } //collect segment alignments into chains
 */
 return false; //no adapter parts found 
 }

//convert qvs to/from phred64 from/to phread33 
void convertPhred(GStr& q) {
 for (int i=0;i<q.length();i++) q[i]+=qv_cvtadd;
}

void convertPhred(char* q, int len) {
 for (int i=0;i<len;i++) q[i]+=qv_cvtadd;
}

bool getFastxRec(GLineReader& fq, GStr& rseq, GStr& rqv, 
          GStr& rname, GStr& rinfo, GStr& infname) {
 rseq="";
 rqv="";
 rname="";
 rinfo="";
 if (fq.eof()) return false;
 char* l=fq.getLine();
 while (l!=NULL && (l[0]==0 || isspace(l[0]))) l=fq.getLine(); //ignore empty lines
 if (l==NULL) return false;
 /* if (rawFormat) {
      //TODO: implement raw qseq parsing here
      //if (raw type=N) then continue; //skip invalid/bad records
      } //raw qseq format
 else { // FASTQ or FASTA */
 isfasta=(l[0]=='>');
 if (!isfasta && l[0]!='@') GError("Error: fasta/fastq record marker not found(%s)\n%s\n", 
      infname.chars(), l);
 GStr s(l);
 rname=&(l[1]);
 for (int i=0;i<rname.length();i++)
    if (rname[i]<=' ') { 
       if (i<rname.length()-2) rinfo=rname.substr(i+1); 
       rname.cut(i); 
       break; 
       }
  //now get the sequence
 if ((l=fq.getLine())==NULL) 
      GError("Error: unexpected EOF after header for read %s (%s)\n",
                   rname.chars(), infname.chars());
 rseq=l; //this must be the DNA line
 while ((l=fq.getLine())!=NULL) {
      //seq can span multiple lines
      if (l[0]=='>' || l[0]=='+') {
           fq.pushBack(); 
           break; //
           }
      rseq+=l;
      } //check for multi-line seq 
 if (!isfasta) { //reading fastq quality values, which can also be multi-line
    if ((l=fq.getLine())==NULL)
        GError("Error: unexpected EOF after sequence for %s\n", rname.chars());
    if (l[0]!='+') GError("Error: fastq qv header marker not detected!\n");
    if ((l=fq.getLine())==NULL) 
        GError("Error: unexpected EOF after qv header for %s\n", rname.chars());
    rqv=l;
    //if (rqv.length()!=rseq.length()) 
    //  GError("Error: qv len != seq len for %s\n", rname.chars());
    while (rqv.length()<rseq.length() && ((l=fq.getLine())!=NULL)) {
      rqv+=l; //append to qv string
      }
    }// fastq
 // } //<-- FASTA or FASTQ
 rseq.upper(); //TODO: what if we care about masking?
 return true;
}

char process_read(GStr& rname, GStr& rseq, GStr& rqv, int &l5, int &l3) {
 //returns 0 if the read was untouched, 1 if it was just trimmed 
 // and a trash code if it was trashed
 l5=0;
 l3=rseq.length()-1;
 if (l3-l5+1<min_read_len) {
   return 's';
   }
GStr wseq(rseq.chars());
GStr wqv(rqv.chars());
int w5=l5;
int w3=l3;
//first do the q-based trimming
if (qvtrim_qmin!=0 && !wqv.is_empty() && qtrim(wqv, w5, w3)) { // qv-threshold trimming
   if (w3-w5+1<min_read_len) {
     return 'Q'; //invalid read
     }
    //-- keep only the w5..w3 range
   l5=w5;
   l3=w3;
   wseq=wseq.substr(w5, w3-w5+1);
   if (!wqv.is_empty())
      wqv=wqv.substr(w5, w3-w5+1);
   } //qv trimming
// N-trimming on the remaining read seq
if (ntrim(wseq, w5, w3)) {
   //GMessage("before: %s\n",wseq.chars());
   //GMessage("after : %s (%d)\n",wseq.substr(w5,w3-w5+1).chars(),w3-w5+1);
   l5+=w5;
   l3-=(wseq.length()-1-w3);
   if (w3-w5+1<min_read_len) {
     return 'N'; //to be trashed
     }
    //-- keep only the w5..w3 range
   wseq=wseq.substr(w5, w3-w5+1);
   if (!wqv.is_empty())
      wqv=wqv.substr(w5, w3-w5+1);
   w5=0;
   w3=wseq.length()-1;
   }
if (a3len>0) {
  if (trim_adapter3(wseq, w5, w3)) {
     int trimlen=wseq.length()-(w3-w5+1);
     num_trimmed3++;
     if (trimlen<min_trimmed3) 
         min_trimmed3=trimlen;
     l5+=w5;
     l3-=(wseq.length()-1-w3);
     if (w3-w5+1<min_read_len) {
         return '3';
         }
      //-- keep only the w5..w3 range
      wseq=wseq.substr(w5, w3-w5+1);
      if (!wqv.is_empty())
         wqv=wqv.substr(w5, w3-w5+1);
      }//some adapter was trimmed
   } //adapter trimming
if (a5len>0) {
  if (trim_adapter5(wseq, w5, w3)) {
     int trimlen=wseq.length()-(w3-w5+1);
     num_trimmed5++;
     if (trimlen<min_trimmed5)
         min_trimmed5=trimlen;
     l5+=w5;
     l3-=(wseq.length()-1-w3);
     if (w3-w5+1<min_read_len) {
         return '5';
         }
      //-- keep only the w5..w3 range
      wseq=wseq.substr(w5, w3-w5+1);
      if (!wqv.is_empty())
         wqv=wqv.substr(w5, w3-w5+1);
      }//some adapter was trimmed
   } //adapter trimming
if (doCollapse) {
   //keep read for later
   FqDupRec* dr=dhash.Find(wseq.chars());
   if (dr==NULL) { //new entry
          //if (prefix.is_empty()) 
             dhash.Add(wseq.chars(), 
                  new FqDupRec(&wqv, rname.chars()));
          //else dhash.Add(wseq.chars(), new FqDupRec(wqv.chars(),wqv.length()));
         }
      else
         dr->add(wqv);
   } //collapsing duplicates
 else { //not collapsing duplicates
   //apply the dust filter now
   if (doDust) {
     int dustbases=dust(wseq);
     if (dustbases>(wseq.length()>>1)) {
        return 'D';
        }
     }
   } //not collapsing duplicates
return (l5>0 || l3<rseq.length()-1) ? 1 : 0;
}

void printHeader(FILE* f_out, char recmarker, GStr& rname, GStr& rinfo) {
 //GMessage("printing Header..%c%s\n",recmarker, rname.chars());
 if (rinfo.is_empty()) fprintf(f_out, "%c%s\n",recmarker,rname.chars());
  else fprintf(f_out, "%c%s %s\n",recmarker, rname.chars(), rinfo.chars());
 }

void writeRead(FILE* f_out, GStr& rname, GStr& rinfo, GStr& rseq, GStr& rqv, int& outcounter) {
   outcounter++;
   bool asFasta=(rqv.is_empty() || fastaOutput);
   if (asFasta) {
    if (prefix.is_empty()) {
       printHeader(f_out, '>',rname,rinfo);
       fprintf(f_out, "%s\n", rseq.chars()); //plain one-line fasta for now
       }
      else {
       fprintf(f_out, ">%s%08d\n%s\n", prefix.chars(), outcounter, 
                          rseq.chars());
       }
     }
   else {  //fastq
    if (convert_phred) convertPhred(rqv);
    if (prefix.is_empty()) {
       printHeader(f_out, '@', rname, rinfo);
       fprintf(f_out, "%s\n+\n%s\n", rseq.chars(), rqv.chars());
       }
      else
       fprintf(f_out, "@%s_%08d\n%s\n+\n%s\n", prefix.chars(), outcounter, 
                          rseq.chars(),rqv.chars() );
     }
}

void trash_report(char trashcode, GStr& rname, FILE* freport) {
 if (freport==NULL || trashcode<=' ') return;
 if (trashcode=='3' || trashcode=='5') {
   fprintf(freport, "%s\tA%c\n",rname.chars(),trashcode);
   }
 else {
   fprintf(freport, "%s\t%c\n",rname.chars(),trashcode);
   }
 //tcounter++;
}

GStr getFext(GStr& s, int* xpos=NULL) {
 //s must be a filename without a path
 GStr r("");
 if (xpos!=NULL) *xpos=0;
 if (s.is_empty() || s=="-") return r;
 int p=s.rindex('.');
 int d=s.rindex('/');
 if (p<=0 || p>s.length()-2 || p<s.length()-7 || p<d) return r;
 r=s.substr(p+1);
 if (xpos!=NULL) *xpos=p+1;
 r.lower();
 return r;
 }

void baseFileName(GStr& fname) {
 //remove all known extensions, like .txt,fq,fastq,fasta,fa)(.gz .gzip .bz2 .bzip2) .
 int xpos=0;
 GStr fext=getFext(fname, &xpos);
 if (xpos<=1) return;
 bool extdel=false;
 GStr f2;
 if (fext=="z" || fext=="zip") {
   extdel=true;
   }
  else if (fext.length()>=2) {
   f2=fext.substr(0,2);
   extdel=(f2=="gz" || f2=="bz");
   }
 if (extdel) {
   fname.cut(xpos-1);
   fext=getFext(fname, &xpos);
   if (xpos<=1) return;
   }
 extdel=false;
 if (fext=="f" || fext=="fq" || fext=="txt" || fext=="seq" || fext=="sequence") {
   extdel=true;
   }
  else if (fext.length()>=2) {
   extdel=(fext.substr(0,2)=="fa");
   }
 if (extdel) fname.cut(xpos-1);
 GStr fncp(fname);
 fncp.lower();
 fncp.chomp("seq");
 fncp.chomp("sequence");
 fncp.trimR("_.");
 if (fncp.length()<fname.length()) fname.cut(fncp.length());
}

FILE* prepOutFile(GStr& infname, GStr& pocmd) {
  FILE* f_out=NULL;
  GStr fname(getFileName(infname.chars()));
  //eliminate known extensions
  baseFileName(fname);
  if (outsuffix.is_empty() || outsuffix=="-") { return stdout; }
    else if (pocmd.is_empty()) {
               GStr oname(fname);
               oname.append('.');
               oname.append(outsuffix);
               f_out=fopen(oname.chars(),"w");
               if (f_out==NULL) GError("Error: cannot create '%s'\n",oname.chars());
               }
            else {
              GStr oname(">");
              oname.append(fname);
              oname.append('.');
              oname.append(outsuffix);
              pocmd.append(oname);
              f_out=popen(pocmd.chars(), "w");
              if (f_out==NULL) GError("Error: cannot popen '%s'\n",pocmd.chars());
              }
 return f_out;
}

void guess_unzip(GStr& fname, GStr& picmd) {
 GStr fext=getFext(fname);
 if (fext=="gz" || fext=="gzip" || fext=="z") {
    picmd="gzip -cd ";
    }
   else if (fext=="bz2" || fext=="bzip2" || fext=="bz" || fext=="bzip") {
    picmd="bzip2 -cd ";
    }
}

void setupFiles(FILE*& f_in, FILE*& f_in2, FILE*& f_out, FILE*& f_out2, 
                       GStr& s, GStr& infname, GStr& infname2) { 
// uses outsuffix to generate output file names and open file handles as needed
 infname="";
 infname2="";
 f_in=NULL;
 f_in2=NULL;
 f_out=NULL;
 f_out2=NULL;
 //analyze outsuffix intent
 GStr pocmd;
 GStr ox=getFext(outsuffix);
 if (ox.length()>2) ox=ox.substr(0,2);
 if (ox=="gz") pocmd="gzip -9 -c ";
   else if (ox=="bz") pocmd="bzip2 -9 -c ";
 if (s=="-") {
    f_in=stdin;
    infname="stdin";
    f_out=prepOutFile(infname, pocmd);
    return;
    } // streaming from stdin
 s.startTokenize(",:");
 s.nextToken(infname);
 bool paired=s.nextToken(infname2);
 if (fileExists(infname.chars())==0) 
    GError("Error: cannot find file %s!\n",infname.chars());
 GStr fname(getFileName(infname.chars()));
 GStr picmd;
 guess_unzip(fname, picmd);
 if (picmd.is_empty()) {
   f_in=fopen(infname.chars(), "r");
   if (f_in==NULL) GError("Error opening file '%s'!\n",infname.chars());
   }
  else {
   picmd.append(infname);
   f_in=popen(picmd.chars(), "r");
   if (f_in==NULL) GError("Error at popen %s!\n", picmd.chars());
   }
 f_out=prepOutFile(infname, pocmd);
 if (!paired) return;
 if (doCollapse) GError("Error: sorry, -C option cannot be used with paired reads!\n");
 // ---- paired reads:-------------
 if (fileExists(infname2.chars())==0) 
     GError("Error: cannot find file %s!\n",infname2.chars());
 picmd="";
 GStr fname2(getFileName(infname2.chars()));
 guess_unzip(fname2, picmd);
 if (picmd.is_empty()) {
   f_in2=fopen(infname2.chars(), "r");
   if (f_in2==NULL) GError("Error opening file '%s'!\n",infname2.chars());
   }
  else {
   picmd.append(infname2);
   f_in2=popen(picmd.chars(), "r");
   if (f_in2==NULL) GError("Error at popen %s!\n", picmd.chars());
   }
 f_out2=prepOutFile(infname2, pocmd);
}
