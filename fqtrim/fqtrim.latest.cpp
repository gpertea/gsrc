#include "GArgs.h"
#include "GStr.h"
#include "GHash.hh"
#include "GList.hh"

#define USAGE "Usage:\n\
 fqtrim [-5 <5'adapter>] [-3 <3'adapter>] [-l <minlen>] [-C] [-D] [-Q] \\\n\
    [-n <rename_prefix>] [-o <trimmed.fq>] [-r <discarded.lst>] <input.fq>\n\
 \n\
 Trim low quality bases at 3' end, optional adapter, filter for low complexity and\n\
 optionally collapse duplicates in the input fastq data.\n\
\n\
Options:\n\
-n  rename all the reads using the <prefix> followed by a read counter;\n\
    if -C option was given, the suffix \"_x<N>\" is appended, with <N> being\n\
    the read duplication count\n\
-o  write the trimmed/filtered fastq into <trimmed.fq>(instead of stdout)\n\
-5  trim the given adapter or primer sequence at the 5' end of each read\n\
    (e.g. -5 CGACAGGTTCAGAGTTCTACAGTCCGACGATC)\n\
-3  trim the given adapter sequence at the 3' end of each read\n\
    (e.g. -3 TCGTATGCCGTCTTCTGCTTG)\n\
-l  minimum \"clean\" length at the high quality 5'end that a read must\n\
    have in order to pass the filter (default: 16)\n\
-r  write a simple \"trash report\" file listing the discarded reads with a\n\
    one letter code for the reason why they were trashed\n\
-C  collapse duplicate reads and add a prefix with their count to the read\n\
    name\n\
-D  apply a low-complexity (dust) filter and discard any read that has at\n\
    least half of it masked by this\n\
-Q  quality values in the input data are interpreted as phred64, convert\n\
    them to phred33\n\
 "
// example 3' adapter for miRNAs: TCGTATGCCGTCTTCTGCTTG

//For pair ends sequencing:
//3' : ACACTCTTTCCCTACACGACGCTCTTCCGATCT
//5' : GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG
FILE* f_out=NULL; //stdout if not provided
FILE* f_in=NULL; //input fastq (stdin if not provided)
FILE* freport=NULL;
bool debug=false;
bool doCollapse=false;
bool doDust=false;
bool trashReport=false;
bool rawFormat=false;
int min_read_len=16;
int dust_cutoff=16;
bool isfasta=false;
bool phred64=false;
GStr prefix;
GStr adapter3;
GStr adapter5;
int a3len=0;
int a5len=0;
// adaptor matching metrics -- for extendMatch() function
const int a_m_score=2; //match score
const int a_mis_score=-3; //mismatch
const int a_dropoff_score=7;
const int a_min_score=8; //an exact match of 4 bases at the proper ends WILL be trimmed
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
      return (a.start==d.a.start && a.end==d.a.end && b.start==d.b.start && b.end==d.b.end);
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
       f=fopen(s,"w");
       if (f==NULL) GError("Error creating file: %s\n", s.chars());
       }
     }
}

#define FWCLOSE(fh) if (fh!=NULL && fh!=stdout) fclose(fh)
#define FRCLOSE(fh) if (fh!=NULL && fh!=stdin) fclose(fh)

GHash<FqDupRec> dhash; //hash to keep track of duplicates

bool ntrim(GStr& rseq, int &l5, int &l3); //returns true if any trimming occured
 
int dust(GStr& seq);
bool trim_adapter3(GStr& seq, int &l5, int &l3); //returns true if any trimming occured
bool trim_adapter5(GStr& seq, int &l5, int &l3); //returns true if any trimming occured

void convertQ64(char* q, int len);
void convertQ64(GStr& q);

int main(int argc, char * const argv[]) {
  GArgs args(argc, argv, "YQRDCl:d:3:5:n:r:o:");
  int e;
  int icounter=0; //counter for input reads
  int outcounter=0; //counter for output reads
  if ((e=args.isError())>0) {
      GMessage("%s\nInvalid argument: %s\n", USAGE, argv[e]);
      exit(224);
      }
  debug=(args.getOpt('Y')!=NULL);
  phred64=(args.getOpt('Q')!=NULL);
  doCollapse=(args.getOpt('C')!=NULL);
  doDust=(args.getOpt('D')!=NULL);
  rawFormat=(args.getOpt('R')!=NULL);
  if (rawFormat) {
    GError("Sorry, raw qseq format parsing is not implemented yet!\n");
    }
  prefix=args.getOpt('n');
  GStr s=args.getOpt('l');
  if (!s.is_empty()) 
     min_read_len=s.asInt();
  s=args.getOpt('d');
  if (!s.is_empty()) {
     dust_cutoff=s.asInt();
     doDust=true;
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
  trashReport=  (args.getOpt('r')!=NULL);
  if (args.startNonOpt()==0) {
    GMessage(USAGE);
    exit(224);
    }

  openfw(f_out, args, 'o');
  if (f_out==NULL) f_out=stdout;
  if (trashReport)
    openfw(freport, args, 'r');
  char* infile=NULL;
  while ((infile=args.nextNonOpt())!=NULL) {
    GStr infname(infile);
    if (strcmp(infile,"-")==0) {
       f_in=stdin; infname="stdin"; }
     else {
        f_in=fopen(infile,"r");
        if (f_in==NULL)
            GError("Cannot open input file %s!\n",infile);
        }
     GLineReader fq(f_in);
     char* l=NULL;
     while ((l=fq.getLine())!=NULL) {
        GStr rname; //current read name
        GStr rseq;  //current read sequence
        GStr rqv;   //current read quality values
        GStr s;
        if (rawFormat) {
          //TODO: implement qseq parsing here
          //if (raw type=N) then continue; //skip invalid/bad records
          
          } //raw qseq format
         else { // FASTQ or FASTA
          isfasta=(l[0]=='>');
          if (!isfasta && l[0]!='@') GError("Error: fastq record marker not detected!\n");
          s=l;
          rname=&(l[1]);
          icounter++;
          //GMessage("readname=%s\n",rname.chars());
          for (int i=0;i<rname.length();i++)
            if (rname[i]<=' ') { rname.cut(i); break; }
          //now get the sequence
          if ((l=fq.getLine())==NULL) 
              GError("Error: unexpected EOF after header for %s\n",rname.chars());
          rseq=l; //this must be the DNA line
          if (isfasta) {
            while ((l=fq.getLine())!=NULL) {
              //fasta can have multiple sequence lines
              if (l[0]=='>') { 
                   fq.pushBack(); 
                   break; //
                   }
              rseq+=l;
              }
            } //multi-line fasta file
          if (!isfasta) {
            if ((l=fq.getLine())==NULL) 
                GError("Error: unexpected EOF after sequence for %s\n", rname.chars());
            if (l[0]!='+') GError("Error: fastq qv header marker not detected!\n");
            if ((l=fq.getLine())==NULL) 
                GError("Error: unexpected EOF after qv header for %s\n", rname.chars());
            rqv=l;
            if (rqv.length()!=rseq.length()) 
                GError("Error: qv len != seq len for %s\n", rname.chars());
            }
        } //<-- FASTA or FASTQ
        rseq.upper();
        int l5=0;
        int l3=rseq.length()-1;
        if (l3-l5+1<min_read_len) {
           if (trashReport) {
                  fprintf(freport, "%s\ts\t%s\n",rname.chars(), rseq.chars());
                  }
           continue;
           }
        if (ntrim(rseq, l5, l3)) {
           //GMessage("before: %s\n",rseq.chars());
           //GMessage("after : %s (%d)\n",rseq.substr(l5,l3-l5+1).chars(),l3-l5+1);
           if (l3-l5+1<min_read_len) {
             if (trashReport) {
                    fprintf(freport, "%s\tN\t%s\n", rname.chars(), rseq.chars());
                    }
             continue; //invalid read
             }
            //-- keep only the l5..l3 range
           rseq=rseq.substr(l5, l3-l5+1);

           if (!rqv.is_empty())
              rqv=rqv.substr(l5, l3-l5+1);
           }
           
        if (a3len>0) {
          if (trim_adapter3(rseq, l5, l3)) {
             if (l3-l5+1<min_read_len) {
                 if (trashReport) {
                     fprintf(freport, "%s\tA3\t%s\n",rname.chars(), rseq.chars());
                     }
                 continue;
                 }
              //-- keep only the l5..l3 range
              rseq=rseq.substr(l5, l3-l5+1);
              if (!rqv.is_empty())
                 rqv=rqv.substr(l5, l3-l5+1);
              }//some adapter was trimmed
           } //adapter trimming
        if (a5len>0) {
          if (trim_adapter5(rseq, l5, l3)) {
             if (l3-l5+1<min_read_len) {
                 if (trashReport) {
                     fprintf(freport, "%s\tA5\t%s\n",rname.chars(), rseq.chars());
                     }
                 continue;
                 }
              //-- keep only the l5..l3 range
              rseq=rseq.substr(l5, l3-l5+1);
              if (!rqv.is_empty())
                 rqv=rqv.substr(l5, l3-l5+1);
              }//some adapter was trimmed
           } //adapter trimming
        if (doCollapse) {
           //keep read for later
           FqDupRec* dr=dhash.Find(rseq.chars());
           if (dr==NULL) { //new entry
                  //if (prefix.is_empty()) 
                     dhash.Add(rseq.chars(), 
                          new FqDupRec(&rqv, rname.chars()));
                  //else dhash.Add(rseq.chars(), new FqDupRec(rqv.chars(),rqv.length()));
                 }
              else    
                 dr->add(rqv);
           } //collapsing duplicates
         else { //not collapsing duplicates
           //do the dust filter now
           if (doDust) {
             int dustbases=dust(rseq);
             if (dustbases>(rseq.length()>>1)) {
                if (trashReport) {
                  fprintf(freport, "%s\tD\t%s\n",rname.chars(),rseq.chars());
                  }
                continue;
                }
             }
           //print this record here  
           outcounter++;
           if (isfasta) {
            if (prefix.is_empty())
               fprintf(f_out, ">%s\n%s\n", rname.chars(), rseq.chars());
              else 
               fprintf(f_out, ">%s%08d\n%s\n", prefix.chars(), outcounter, 
                                  rseq.chars());
             }
           else {  //fastq
            if (phred64) convertQ64(rqv);
            if (prefix.is_empty())
               fprintf(f_out, "@%s\n%s\n+\n%s\n", rname.chars(), rseq.chars(),rqv.chars());
              else 
               fprintf(f_out, "@%s_%08d\n%s\n+\n%s\n", prefix.chars(), outcounter, 
                                  rseq.chars(),rqv.chars() );
             }
           } //not collapsing duplicates
        } //for each fastq record
   } //while each input file
 FRCLOSE(f_in);
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
              fprintf(freport, "%s_x%d\tD\t%s\n",qd->firstname, qd->count,seq);
              }
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
       if (phred64) convertQ64(qd->qv, qd->len);
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
   } //report collapsed dhash entries
 GMessage("Number of input reads: %9d\n", icounter);
 GMessage("       Output records: %9d\n", outcounter);
 if (trashReport) {
    FWCLOSE(freport);
    }
 
 FWCLOSE(f_out);
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
   NData():NCount(0),end5(0),end3(0),seqlen(0),perc_N(0),
  		 seq(NULL), valid(true) {   }
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


bool ntrim(GStr& rseq, int &l5, int &l3) {
 //count Ns in the sequence, trim N-rich ends
 feat.init(rseq);
 l5=feat.end5-1;
 l3=feat.end3-1;
 N_analyze(feat.end5-1, feat.end3-1, 0, feat.NCount-1); 
 if (l5==feat.end5-1 && l3==feat.end3-1) {
    if (feat.perc_N>6.2) {
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
 
 if (feat.perc_N>6.2) { //not valid if more than 1 N per 16 bases
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
 int imin=rlen<<1; //last half of the read, left boundary for the wmer match
 if (imin<a3len) { imin=GMIN(a3len, rlen-wordSize); }
 imin=rlen-imin;
 for (int iw=0;iw<hlen;iw++) {
   //int32* qv=(int32*)(adapter3.chars()+iw);
   int qv=get3mer_value(adapter3.chars()+iw);
   fi=-1;
   //while ((fi=fast4rmatch(*qv, seq.chars(), rlen, fi))>=0 && fi>=imin) {
   while ((fi=w3_rmatch(qv, seq.chars(), rlen, fi))>=0 && fi>=imin) {
#ifdef DEBUG
     if (debug) GMessage("    Wmatch found: %*s\n", fi+wordSize, (adapter3.substr(iw,wordSize)).chars());
#endif
     if (extendMatch(seq.chars(), rlen, fi, adapter3.chars(), 
                   a3len, iw, wordSize, l5,l3, a3segs)) return true;
     fi--;
     if (fi<imin) break;
     }
   } //for each wmer in the first hlen bases of the adaptor
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
 //int hlen=a5len;
 if (hlen>a5len-wordSize) hlen=a5len-wordSize;
 //int hlen=a5len-5;
 int imax=rlen<<1; //tfirst half of the read, right boundary for the wmer match
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

 //couldn't find a good trimming extension, hash 12 more bases of the adapter to collect more segment pairs there
 //but only do this if we already have segment pairs collected in the last 12 bases of the adapter
 if (a5segs.bend<(uint)(a5len-3) || a5segs.bstart>(uint)(a5len-hlen+4)) return false;
 //if (a5segs.bend>(uint)(a5len-3) && a5segs.bstart<(uint)(a5len-hlen+4)) {
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
 return false; //no adapter parts found 
 }

//conversion of phred64 to phread33
void convertQ64(GStr& q) {
 for (int i=0;i<q.length();i++) q[i]-=31;
}

void convertQ64(char* q, int len) {
 for (int i=0;i<len;i++) q[i]-=31;
}

