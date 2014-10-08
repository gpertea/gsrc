#include "GArgs.h"
#include "GStr.h"
#include "GHash.hh"
#include "GList.hh"

#define USAGE "Usage:\n\
 fqtrim [-t <3'adapter>] [-l <minlen>] [-n <rename_prefix>] [-C] [-D] [-Q] \\ \n\
    [-o <trimmed.fq>] <input.fq>\n\
 \n\
 Trim low quality bases at 3' end, optional adapter, filter for low complexity and\n\
 optionally collapse duplicates in the input fastq data.\n\
\n\
Options:\n\
-n  rename all the reads using the <prefix> followed by a read counter;\n\
    if -C option was given, the suffix \"_x<N>\" is appended, with <N> being\n\
    the read duplication count\n\
-o  write the trimmed/filtered fastq into <trimmed.fq>(instead of stdout)\n\
-t  check and trim the given adapter sequence at the 3' end or each read\n\
-l  minimum \"clean\" length at the high quality 5'end that a read must\n\
    have in order to pass the filter (default: 16)\n\
-C  collapse duplicate reads and add a prefix with their count to the read\n\
    name\n\
-D  apply a low-complexity (dust) filter and discard any read that has at\n\
    least half of it masked by this\n\
-Q  quality values in the input data are interpreted as phred64, convert \n\
    them to phred33\n\
 "
// example 3' adapter for miRNAs: TCGTATGCCGTCTTCTGCTTG
FILE* f_out=NULL; //stdout if not provided
FILE* f_in=NULL; //input fastq (stdin if not provided)
bool debug=false;
bool doCollapse=false;
bool doDust=false;
int min_read_len=16;
int dust_cutoff=16;
bool phred64=false;
GStr prefix;
char adapter3[256];
int a3len=0;
// element in hash
class FqDupRec {
 public:
   int count; //how many of these reads are the same
   int len; //length of qv
   char* qv;
   FqDupRec(const char* q=NULL, int l=0) {
     len=l;
     GMALLOC(qv, len+1);
     strncpy(qv, q, len);
     count=1;
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

void ntrim(GStr& rseq, GStr &rqv, int &l5, int &l3);
int dust(GStr& seq);
void trim_adapter3(GStr& seq, int &l5, int &l3);

int main(int argc, char * const argv[]) {
  GArgs args(argc, argv, "QDCl:d:t:n:o:");
  int e;
  if ((e=args.isError())>0)
      GError("%s\nInvalid argument: %s\n", USAGE, argv[e]);
  if (args.getOpt('h')!=NULL) GError("%s\n", USAGE);
  debug=(args.getOpt('D')!=NULL);
  phred64=(args.getOpt('Q')!=NULL);
  doCollapse=(args.getOpt('C')!=NULL);
  doDust=(args.getOpt('D')!=NULL);
  prefix=args.getOpt('n');
  GStr s=args.getOpt('l');
  if (!s.is_empty()) 
     min_read_len=s.asInt();
  s=args.getOpt('d');
  if (!s.is_empty()) {
     dust_cutoff=s.asInt();
     doDust=true;
     }
     
  if (args.getOpt('t')!=NULL) {
    a3len=strlen(args.getOpt('t'));
    if (a3len>255) a3len=255;
    strncpy(adapter3, args.getOpt('t'), a3len);
    }
  GStr rname; //current read name
  GStr rseq; //current read sequence
  GStr rqv;  //current read quality values
  args.startNonOpt();
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
        bool isfasta=(l[0]=='>');
        if (!isfasta && l[0]!='@') GError("Error: fastq record marker not detected!\n");
        GStr s(l);
        GStr rname(&(l[1]));
        GMessage("readname=%s\n",rname.chars());
        for (int i=0;i<rname.length();i++)
          if (rname[i]<=' ') { rname.cut(i); break; }
        //now get the sequence
        if ((l=fq.getLine())==NULL) 
            GError("Error: unexpected EOF after header for %s\n",rname.chars());
        rseq=l; //this must be the DNA line
        rseq.upper();
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
        int l5=0;
        int l3=rseq.length()-1;
        ntrim(rseq, rqv, l5, l3);
        if (l3-l5+1<min_read_len) continue; //invalid read
        rseq=rseq.substr(l5, l3-l5+1);
        if (!rqv.is_empty())
           rqv=rqv.substr(l5, l3-l5+1);
        if (a3len>0)
          trim_adapter3(rseq, l5, l3);
        if (l3-l5+1<min_read_len) continue;
        if (doCollapse) {
           //keep read for later
           FqDupRec* dr=dhash.Find(rseq.chars());
           if (dr==NULL) 
                 dhash.Add(rseq.chars(), new FqDupRec(rqv.chars(),rqv.length()));
              else    
                 dr->add(rqv);
           } //collapsing duplicates
         else { //not collapsing duplicates
           //do the dust filter now
           if (doDust) {
             int dustbases=dust(rseq);
             if (dustbases>(rseq.length()>>1)) continue;
             }
           //print this record here  
           } //not collapsing duplicates
        } //for each fastq record
   } //while each input file
 FRCLOSE(f_in);
 
 if (doCollapse) {
   //TODO: run dust() on each record here
   //      and print those that are not >50% dusted
   }
 FWCLOSE(f_out);
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
     seqlen=rseq.length();
     seq=rseq.chars();
     end5=1;
     end3=seqlen;
     for (int i=0;i<seqlen;i++)
        if (rseq[i]=='N') {// if (!ichrInStr(rseq[i], "ACGT")
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
 if (l3-l5<min_read_len || p5>p3 ) {
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


void ntrim(GStr& rseq, GStr &rqv, int &l5, int &l3) {
 //count Ns in the sequence
 feat.init(rseq);
 N_analyze(feat.end5-1, feat.end3-1, 0, feat.NCount-1); 
 l5=feat.end5-1;
 l3=feat.end3-1;
 if (l3-l5+1<min_read_len) {
   feat.valid=false;
   return;
   }
 feat.N_calc();
 if (feat.perc_N>6.2) { //not valid if more than 1 N per 16 bases
      feat.valid=false;
      l3=l5+1;
      return;
      }
 }
//--------------- dust functions ----------------
#define MAXREG    1001

static int dustword = 3; 
static int dustwindow = 32; 
static int dustwindow2 = 16; 
static int mv, iv, jv;

void set_dust_window(int value) {
   dustwindow = value;
   dustwindow2 = (dustwindow >> 1);
}

void set_dust_word(int value) {
   dustword = value;
}

static void dust_wo1(int len, const char* s, int ivv) {
      int i, ii, j, v, t, n, n1, sum;
      static int counts[32*32*32];
      static int iis[32*32*32];
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

int dust_wo(int len, const char* s, int* beg, int* end) {
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
            dust_wo1(len-i, s+i, i);
      }
      *beg = iv;
      *end = iv + jv;
      return mv;
}

void duststr(const char* seq, char* seqmsk, int seqlen, int level) { //, maskFunc maskfn) {
      int i, j, l, a, b, v;
      //DustRgn* first=NULL;
      //DustRgn* cur=NULL;
      //DustRgn* prev=NULL;
      a=0;b=0;
      for (i=0; i < seqlen; i += dustwindow2) {
            l = (seqlen > i+dustwindow) ? dustwindow : seqlen-i;
            v = dust_wo(l, seq+i, &a, &b);
            if (v > level) {
               //for (j = a; j <= b && j < dustwindow2; j++) {
               for (j = a; j <= b; j++) {
                        seqmsk[i+j]='N';
                        }
                }
      }
 //return first; 
}

int dust(GStr& rseq) {
 char* seq=Gstrdup(rseq.chars());
 duststr(rseq.chars(), seq, rseq.length(), dust_cutoff);
 //check the number of Ns:
 int ncount=0;
 for (int i=0;i<rseq.length();i++) {
   if (seq[i]=='N') ncount++;
   }
 fprintf(stderr,"orig: %s\n",rseq.chars());
 fprintf(stderr,"dust: %s\n",seq);
 GFREE(seq);
 return ncount;
 }

// ------------------ adapter matching

void trim_adapter3(GStr& seq, int&l5, int &l3) {
 
 }


